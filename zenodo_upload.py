#!/usr/bin/env python3
"""
Package ./dependencies and upload it to a Zenodo DRAFT deposition.

Usage:
    ZENODO_TOKEN=xxxx python3 zenodo_upload.py <deposition_id> [--purge]

Constraints handled
--------------------
* Zenodo allows a MAXIMUM OF 100 FILES per record  -> we must produce <=100 objects.
* Zenodo storage is FLAT (filenames cannot contain "/").

Scheme
------
* Files >= BIG_THRESHOLD are uploaded INDIVIDUALLY (per-file resume), with key =
  relative path, "/" -> "__"  (e.g. shared/RV_data.rds -> shared__RV_data.rds).
* All smaller files are BIN-PACKED into zip "parts" of <= PART_CAP bytes
  (dependencies_part_000.zip, ...). Each file is stored in the zip under its path
  relative to ./dependencies, so `unzip -d dependencies` reconstructs the tree.
  A crash loses at most one part (<= PART_CAP).
* manifest.csv + unpack_dependencies.sh are uploaded too, so the full tree can be
  rebuilt from the downloaded record.

Resume: existing bucket contents are listed; objects already present (files: key+size
match; zip parts: key present) are skipped. Re-run the same command after any crash.
--purge deletes ALL files currently on the draft first (use when changing the scheme).
"""
import os
import sys
import time
import json
import csv
import zipfile
import socket
import urllib.request
import urllib.parse
import urllib.error

# Per-socket-operation inactivity timeout: a stalled/half-open TCP connection
# raises socket.timeout (caught by retry logic) instead of hanging forever.
# Safe for slow-but-steady transfers since data keeps flowing within this window.
socket.setdefaulttimeout(180)

# Optional: force IPv4 (ZENODO_FORCE_IPV4=1). IPv6 paths via a VPN often cause
# TLS handshake timeouts / mid-transfer drops on large uploads; restricting DNS
# resolution to A records routes everything over IPv4.
if os.environ.get("ZENODO_FORCE_IPV4") == "1":
    _orig_getaddrinfo = socket.getaddrinfo
    def _getaddrinfo_v4(host, port, family=0, *args, **kwargs):
        res = _orig_getaddrinfo(host, port, socket.AF_INET, *args, **kwargs)
        return res or _orig_getaddrinfo(host, port, family, *args, **kwargs)
    socket.getaddrinfo = _getaddrinfo_v4
    print("ZENODO_FORCE_IPV4=1 -> using IPv4 only")

BASE = "https://zenodo.org/api"
ROOT = os.path.dirname(os.path.abspath(__file__))
DEP_DIR = os.path.join(ROOT, "dependencies")
STAGE = os.path.join(ROOT, "dependencies_zenodo_staging")
BIG_THRESHOLD = 250_000_000   # files >= this go up individually
PART_CAP = 2_000_000_000      # max bytes per zip part (max loss on a crash)
MAX_OBJECTS = 95              # safety: refuse to start if plan exceeds this


# --------------------------------------------------------------------------- http helpers
def _request(req, retries=5):
    """Run a urllib request, retrying transient 5xx / network errors."""
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(req) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as e:
            if e.code >= 500 and attempt < retries:
                print(f"  transient HTTP {e.code} on {req.get_method()} "
                      f"(attempt {attempt}); retrying...")
                time.sleep(5 * attempt)
                continue
            raise
        except urllib.error.URLError as e:
            if attempt < retries:
                print(f"  network error {e} (attempt {attempt}); retrying...")
                time.sleep(5 * attempt)
                continue
            raise


def api_get(url, token):
    sep = "&" if "?" in url else "?"
    req = urllib.request.Request(url + sep + "access_token=" + token)
    _, body = _request(req)
    return json.loads(body)


def api_delete(url, token):
    sep = "&" if "?" in url else "?"
    req = urllib.request.Request(url + sep + "access_token=" + token, method="DELETE")
    status, _ = _request(req)
    return status


def put_file(bucket, key, path, token, size):
    url = f"{bucket}/{urllib.parse.quote(key)}?access_token={token}"
    with open(path, "rb") as fh:
        req = urllib.request.Request(url, data=fh, method="PUT")
        req.add_header("Content-Type", "application/octet-stream")
        req.add_header("Content-Length", str(size))
        with urllib.request.urlopen(req) as r:
            return json.load(r)


# ------------------------------------------------------------------------------- planning
def key_for_path(rel):
    return rel.replace(os.sep, "__")


def collect_files():
    out = []
    for r, _, names in os.walk(DEP_DIR):
        for n in names:
            if n == ".DS_Store":
                continue
            p = os.path.join(r, n)
            try:
                out.append((p, os.path.relpath(p, DEP_DIR), os.path.getsize(p)))
            except OSError:
                pass
    return out


def build_plan():
    """Return list of plan dicts. Bin-packs small files into <=PART_CAP zip parts."""
    files = collect_files()
    big = [(p, rel, s) for (p, rel, s) in files if s >= BIG_THRESHOLD]
    small = sorted((f for f in files if f[2] < BIG_THRESHOLD),
                   key=lambda x: -x[2])  # first-fit-decreasing

    bins = []  # each: [total, [(abspath, relpath), ...]]
    for p, rel, s in small:
        placed = False
        for b in bins:
            if b[0] + s <= PART_CAP:
                b[0] += s
                b[1].append((p, rel))
                placed = True
                break
        if not placed:
            bins.append([s, [(p, rel)]])

    plan = []
    for p, rel, s in sorted(big, key=lambda x: x[2]):
        plan.append({"kind": "file", "key": key_for_path(rel), "target": rel,
                     "src": p, "size": s})
    for i, b in enumerate(bins):
        plan.append({"kind": "zip", "key": f"dependencies_part_{i:03d}.zip",
                     "target": "", "members": b[1], "size": b[0]})
    return plan


def build_part(key, members):
    """Create a zip part in STAGE; arcnames are relative to ./dependencies."""
    os.makedirs(STAGE, exist_ok=True)
    zpath = os.path.join(STAGE, key)
    tmp = zpath + ".tmp"
    with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_STORED, allowZip64=True) as zf:
        for abspath, rel in members:
            zf.write(abspath, rel)
    os.rename(tmp, zpath)
    return zpath


def write_meta(plan):
    os.makedirs(STAGE, exist_ok=True)
    man = os.path.join(STAGE, "manifest.csv")
    with open(man, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["zenodo_key", "kind", "target_path"])
        for p in plan:
            w.writerow([p["key"], p["kind"], p["target"]])
    script = os.path.join(STAGE, "unpack_dependencies.sh")
    with open(script, "w") as fh:
        fh.write(r"""#!/usr/bin/env bash
# Reconstruct ./dependencies from a downloaded Zenodo record.
# Usage: ./unpack_dependencies.sh <downloaded_dir> [output_dir]   (default ./dependencies)
set -euo pipefail
DL="${1:?Usage: unpack_dependencies.sh <downloaded_dir> [output_dir]}"
OUT="${2:-./dependencies}"
mkdir -p "$OUT"
tail -n +2 "$DL/manifest.csv" | while IFS=, read -r key kind target; do
  [ -z "${key:-}" ] && continue
  case "$kind" in
    file)
      mkdir -p "$OUT/$(dirname "$target")"
      cp "$DL/$key" "$OUT/$target" ;;
    zip)
      unzip -o -q "$DL/$key" -d "$OUT" ;;
  esac
done
echo "Reconstructed tree under: $OUT"
""")
    os.chmod(script, 0o755)
    plan.append({"kind": "file", "key": "manifest.csv", "target": "_meta_manifest.csv",
                 "src": man, "size": os.path.getsize(man)})
    plan.append({"kind": "file", "key": "unpack_dependencies.sh",
                 "target": "_meta_unpack.sh", "src": script,
                 "size": os.path.getsize(script)})
    return plan


# ----------------------------------------------------------------------------------- main
def main():
    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        sys.exit("ERROR: set ZENODO_TOKEN env var")
    args = sys.argv[1:]
    purge = "--purge" in args
    ids = [a for a in args if not a.startswith("--")]
    if not ids:
        sys.exit("Usage: ZENODO_TOKEN=xxx python3 zenodo_upload.py <deposition_id> [--purge]")
    dep_id = ids[0]

    dep = api_get(f"{BASE}/deposit/depositions/{dep_id}", token)
    bucket = dep["links"].get("bucket")
    print(f"Deposition {dep_id}: '{dep.get('title')}'  state={dep.get('state')}  "
          f"submitted={dep.get('submitted')}")
    if dep.get("submitted") or not bucket:
        sys.exit("ERROR: record is published/immutable. Create a NEW VERSION first.")

    # NOTE: the bucket-listing GET (api/files/<id>) is unreliable (intermittent 500s),
    # so we use the deposition files API for listing/deleting; PUT still uses the bucket.
    def list_dep_files():
        d = api_get(f"{BASE}/deposit/depositions/{dep_id}", token)
        return d.get("files", [])

    if purge:
        cur = list_dep_files()
        print(f"Purging {len(cur)} existing files from draft...")
        for f in cur:
            api_delete(f"{BASE}/deposit/depositions/{dep_id}/files/{f['id']}", token)
        print("Purge complete.")

    print("Building plan...")
    plan = build_plan()
    plan = write_meta(plan)
    nz = sum(1 for p in plan if p["kind"] == "zip")
    total = len(plan)
    total_bytes = sum(p["size"] for p in plan)
    print(f"Plan: {total} objects ({total-nz} individual files + {nz} zip parts), "
          f"{total_bytes/1e9:.1f} GB")
    if total > MAX_OBJECTS:
        sys.exit(f"ABORT: {total} objects exceeds safety limit {MAX_OBJECTS} "
                 f"(Zenodo hard cap is 100). Raise thresholds.")

    existing = {}
    try:
        for f in list_dep_files():
            existing[f["filename"]] = f.get("filesize")
    except Exception as e:
        print(f"  (could not list existing files: {e}; assuming none)")

    plan.sort(key=lambda p: p["size"])  # small first
    done = 0
    for i, p in enumerate(plan, 1):
        key = p["key"]
        if (p["kind"] == "file" and existing.get(key) == p["size"]) or \
           (p["kind"] == "zip" and key in existing):
            done += p["size"]
            print(f"[{i}/{total}] SKIP {key}")
            continue
        if p["kind"] == "zip":
            src = build_part(key, p["members"])
            size = os.path.getsize(src)
        else:
            src, size = p["src"], p["size"]
        for attempt in range(1, 11):  # up to 10 tries per object (flaky-link resilience)
            try:
                t0 = time.time()
                put_file(bucket, key, src, token, size)
                dt = time.time() - t0
                done += size
                print(f"[{i}/{total}] OK {key}  {size/1e6:.1f} MB  "
                      f"{size/1e6/max(dt,0.01):.1f} MB/s  ({done/total_bytes*100:.0f}%)")
                break
            except urllib.error.HTTPError as e:
                body = e.read().decode("utf-8", "replace")[:300]
                print(f"[{i}/{total}] HTTP {e.code} {key} (try {attempt}): {body}")
                if e.code in (400, 401, 403):
                    sys.exit("Fatal HTTP error -- stopping.")
                time.sleep(5 * attempt)
            except Exception as e:
                print(f"[{i}/{total}] ERR {key} (try {attempt}): {e}")
                time.sleep(5 * attempt)
        else:
            sys.exit(f"Failed {key} after retries. Re-run to resume.")
        if p["kind"] == "zip":
            try:
                os.remove(src)
            except OSError:
                pass

    print(f"\nAll {total} objects uploaded. Review the draft on zenodo.org and "
          f"publish manually when ready.")


if __name__ == "__main__":
    main()
