#!/usr/bin/env python3
"""
Incremental sync of ./dependencies to a Zenodo DRAFT deposition.

Usage:
    ZENODO_TOKEN=xxxx python3 zenodo_sync.py <deposition_id> [--prune] [--dry-run]

What it does
------------
Same hybrid packaging as zenodo_upload.py (big files individual, small files
bin-packed into <=PART_CAP zip parts, <=100 objects total), but it only uploads
what CHANGED since last time:

  * Each file's content MD5 is tracked in a local state file
    (dependencies_zenodo_staging/upload_state.json).
  * A big file (>= BIG_THRESHOLD) maps 1:1 to a Zenodo object -> if its MD5 is
    unchanged and it's present on the record, it's skipped; if changed/new, only
    that one object is re-uploaded.
  * Small files live inside zip "parts". Part membership is STABLE (stored in
    state): existing files keep their part; new small files fill a part that has
    room or start a new part. A part is rebuilt + re-uploaded ONLY if one of its
    members changed/was added/removed (detected via a per-part content signature).
  * manifest.json (key -> kind/target/members) is uploaded so the record
    self-documents how to rebuild ./dependencies.

Flags
-----
  --prune    : delete record objects (and state entries) whose source files no
               longer exist locally. Without it, removed files are left on the record.
  --dry-run  : print the planned actions (upload/rebuild/skip/delete) and exit.

First run after the initial bulk upload establishes the MD5 baseline: objects
already present on the record are recorded as in-sync without re-uploading.

NOTE: Zenodo published records are immutable. This targets a DRAFT (or a new
version's draft). For a published record, create a new version first, then sync
to that new version's deposition id.
"""
import os
import sys
import json
import time
import hashlib
import zipfile
import urllib.request
import urllib.parse
import urllib.error
import socket

socket.setdefaulttimeout(180)

BASE = "https://zenodo.org/api"
ROOT = os.path.dirname(os.path.abspath(__file__))
DEP_DIR = os.path.join(ROOT, "dependencies")
STAGE = os.path.join(ROOT, "dependencies_zenodo_staging")
STATE_PATH = os.path.join(STAGE, "upload_state.json")
BIG_THRESHOLD = 250_000_000
PART_CAP = 2_000_000_000
MAX_OBJECTS = 95


# --------------------------------------------------------------------------- http helpers
def _request(req, retries=6):
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(req) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as e:
            if e.code >= 500 and attempt < retries:
                time.sleep(5 * attempt); continue
            raise
        except urllib.error.URLError:
            if attempt < retries:
                time.sleep(5 * attempt); continue
            raise


def api_get(url, token):
    sep = "&" if "?" in url else "?"
    _, body = _request(urllib.request.Request(url + sep + "access_token=" + token))
    return json.loads(body)


def api_delete(url, token):
    sep = "&" if "?" in url else "?"
    status, _ = _request(urllib.request.Request(url + sep + "access_token=" + token,
                                                 method="DELETE"))
    return status


def put_file(bucket, key, path, token):
    size = os.path.getsize(path)
    url = f"{bucket}/{urllib.parse.quote(key)}?access_token={token}"
    with open(path, "rb") as fh:
        req = urllib.request.Request(url, data=fh, method="PUT")
        req.add_header("Content-Type", "application/octet-stream")
        req.add_header("Content-Length", str(size))
        for attempt in range(1, 11):
            try:
                with urllib.request.urlopen(req) as r:
                    return json.load(r)
            except Exception as e:
                if attempt == 10:
                    raise
                print(f"    retry {attempt} on {key}: {e}")
                time.sleep(5 * attempt)
                fh.seek(0)
                req.data = fh


# ------------------------------------------------------------------------------- helpers
def md5(path, _buf=1024 * 1024):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(_buf), b""):
            h.update(chunk)
    return h.hexdigest()


def key_for_path(rel):
    return rel.replace(os.sep, "__")


def collect_files():
    out = {}
    for r, _, names in os.walk(DEP_DIR):
        for n in names:
            if n == ".DS_Store":
                continue
            p = os.path.join(r, n)
            rel = os.path.relpath(p, DEP_DIR)
            try:
                out[rel] = os.path.getsize(p)
            except OSError:
                pass
    return out


def load_state():
    if os.path.exists(STATE_PATH):
        with open(STATE_PATH) as fh:
            return json.load(fh)
    return {"files": {}, "part_assignment": {}, "parts": {}}
    # files: relpath -> {size, mtime, md5}
    # part_assignment: relpath -> part_index
    # parts: key -> {"members": [relpath,...], "signature": "..."}


def save_state(state):
    os.makedirs(STAGE, exist_ok=True)
    tmp = STATE_PATH + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(state, fh, indent=1, sort_keys=True)
    os.rename(tmp, STATE_PATH)


def file_md5_cached(rel, size, state):
    """Return md5, using state cache when size+mtime match (avoids re-hashing)."""
    p = os.path.join(DEP_DIR, rel)
    mtime = os.path.getmtime(p)
    rec = state["files"].get(rel)
    if rec and rec.get("size") == size and abs(rec.get("mtime", -1) - mtime) < 1e-6:
        return rec["md5"], mtime
    return md5(p), mtime


def assign_parts(small, state):
    """Stable part assignment. small: dict relpath->size. Returns {part_idx:[relpaths]}."""
    assignment = dict(state.get("part_assignment", {}))
    # drop assignments for files that no longer exist or are no longer 'small'
    assignment = {r: i for r, i in assignment.items() if r in small}
    # current load per part
    load = {}
    for r, i in assignment.items():
        load[i] = load.get(i, 0) + small[r]
    # place new small files
    for r in sorted(small, key=lambda x: -small[x]):
        if r in assignment:
            continue
        placed = None
        for i in sorted(load):
            if load[i] + small[r] <= PART_CAP:
                placed = i
                break
        if placed is None:
            placed = (max(load) + 1) if load else 0
            load[placed] = 0
        assignment[r] = placed
        load[placed] = load.get(placed, 0) + small[r]
    # group
    parts = {}
    for r, i in assignment.items():
        parts.setdefault(i, []).append(r)
    state["part_assignment"] = assignment
    return parts


def build_part(key, members):
    os.makedirs(STAGE, exist_ok=True)
    zpath = os.path.join(STAGE, key)
    tmp = zpath + ".tmp"
    with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_STORED, allowZip64=True) as zf:
        for rel in sorted(members):
            zf.write(os.path.join(DEP_DIR, rel), rel)
    os.rename(tmp, zpath)
    return zpath


# ----------------------------------------------------------------------------------- main
def main():
    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        sys.exit("ERROR: set ZENODO_TOKEN env var")
    args = sys.argv[1:]
    prune = "--prune" in args
    dry = "--dry-run" in args
    ids = [a for a in args if not a.startswith("--")]
    if not ids:
        sys.exit("Usage: ZENODO_TOKEN=xxx python3 zenodo_sync.py <deposition_id> "
                 "[--prune] [--dry-run]")
    dep_id = ids[0]

    dep = api_get(f"{BASE}/deposit/depositions/{dep_id}", token)
    bucket = dep["links"].get("bucket")
    print(f"Deposition {dep_id}: '{dep.get('title')}'  state={dep.get('state')}  "
          f"submitted={dep.get('submitted')}")
    if dep.get("submitted") or not bucket:
        sys.exit("ERROR: record is published/immutable. Create a NEW VERSION first.")

    # current record contents: key -> md5 checksum
    record = {}
    for f in dep.get("files", []):
        ck = (f.get("checksum") or "").replace("md5:", "")
        record[f["filename"]] = ck

    state = load_state()
    files = collect_files()                       # relpath -> size
    big = {r: s for r, s in files.items() if s >= BIG_THRESHOLD}
    small = {r: s for r, s in files.items() if s < BIG_THRESHOLD}

    # ---- plan objects ----
    plan = []   # dicts: kind, key, members(optional), target
    # big files: 1:1
    for rel in sorted(big):
        plan.append({"kind": "file", "key": key_for_path(rel), "target": rel,
                     "rel": rel})
    # small files: stable parts
    parts = assign_parts(small, state)
    for idx in sorted(parts):
        key = f"dependencies_part_{idx:03d}.zip"
        plan.append({"kind": "zip", "key": key, "target": "",
                     "members": sorted(parts[idx])})

    n_obj = len(plan) + 1  # +1 for manifest.json
    if n_obj > MAX_OBJECTS:
        sys.exit(f"ABORT: {n_obj} objects exceeds safety limit {MAX_OBJECTS} "
                 f"(Zenodo cap 100). Raise BIG_THRESHOLD or PART_CAP.")

    # ---- decide actions ----
    to_upload, to_rebuild, skips = [], [], []
    new_files_state = {}

    for p in plan:
        if p["kind"] == "file":
            rel = p["rel"]
            digest, mtime = file_md5_cached(rel, files[rel], state)
            new_files_state[rel] = {"size": files[rel], "mtime": mtime, "md5": digest}
            # Skip if the record already holds this exact content (Zenodo stores md5).
            if digest and record.get(p["key"]) == digest:
                skips.append(p["key"])
            else:
                to_upload.append((p, digest))
        else:  # zip part
            sig_parts = []
            for rel in p["members"]:
                digest, mtime = file_md5_cached(rel, files[rel], state)
                new_files_state[rel] = {"size": files[rel], "mtime": mtime,
                                        "md5": digest}
                sig_parts.append(f"{rel}:{digest}")
            signature = hashlib.md5("\n".join(sorted(sig_parts)).encode()).hexdigest()
            prev = state.get("parts", {}).get(p["key"], {})
            present = p["key"] in record
            # Skip if unchanged, OR first-sync baseline (part already on record and we
            # have no prior signature -> trust it matches the current member set).
            if present and prev.get("signature") in (signature, None):
                skips.append(p["key"])
                state.setdefault("parts", {})[p["key"]] = {
                    "members": p["members"], "signature": signature}
            else:
                to_rebuild.append((p, signature))

    # ---- deletions (objects/files no longer in plan) ----
    planned_keys = {p["key"] for p in plan} | {"manifest.json",
                                               "manifest.csv", "unpack_dependencies.sh"}
    orphan_keys = [k for k in record if k not in planned_keys]

    print(f"\nPlan: {len(plan)} objects (+manifest). "
          f"upload={len(to_upload)} files, rebuild={len(to_rebuild)} parts, "
          f"skip={len(skips)}, orphan={len(orphan_keys)}")
    for p, _ in to_upload:
        print(f"  UPLOAD  {p['key']}")
    for p, _ in to_rebuild:
        print(f"  REBUILD {p['key']}  ({len(p['members'])} files)")
    if orphan_keys:
        print(f"  ORPHANS {'(will delete with --prune)' if prune else '(kept)'}: "
              f"{', '.join(orphan_keys)}")

    if dry:
        print("\n--dry-run: no changes made.")
        return

    # ---- execute ----
    for p, digest in to_upload:
        print(f"UPLOAD {p['key']} ...")
        put_file(bucket, p["key"], os.path.join(DEP_DIR, p["rel"]), token)

    for p, signature in to_rebuild:
        print(f"REBUILD {p['key']} ({len(p['members'])} files) ...")
        z = build_part(p["key"], p["members"])
        try:
            put_file(bucket, p["key"], z, token)
        finally:
            try:
                os.remove(z)
            except OSError:
                pass
        state.setdefault("parts", {})[p["key"]] = {"members": p["members"],
                                                    "signature": signature}

    if prune:
        for k in orphan_keys:
            print(f"DELETE {k} ...")
            # find file id
            fid = next((f["id"] for f in dep.get("files", []) if f["filename"] == k),
                       None)
            if fid:
                api_delete(f"{BASE}/deposit/depositions/{dep_id}/files/{fid}", token)
        state["part_assignment"] = {r: i for r, i in state["part_assignment"].items()
                                    if r in small}

    state["files"] = new_files_state
    save_state(state)

    # ---- manifest.json (key -> kind/target/members) ----
    manifest = {"objects": [{"key": p["key"], "kind": p["kind"],
                             "target": p["target"],
                             "members": p.get("members", [])} for p in plan]}
    mpath = os.path.join(STAGE, "manifest.json")
    with open(mpath, "w") as fh:
        json.dump(manifest, fh, indent=1)
    put_file(bucket, "manifest.json", mpath, token)

    print(f"\nSync complete. {len(to_upload)} files + {len(to_rebuild)} parts "
          f"uploaded, {len(skips)} unchanged. State: {STATE_PATH}")
    print("Review the draft on zenodo.org and publish when ready.")


if __name__ == "__main__":
    main()
