##############################################################################
# Generate_Supplementary_Data.R
#
# Builds the supplementary data files (Data S1-S22) for the RV atlas
# paper, following the v60 manuscript numbering (see the Data S legends in
# RV_snRNASeq_v60).  v60 has 22 cited Data S files; the section order below
# is physical (legacy order) but each section now writes its v60 number.
# Data S14 is CellChat spatial niche communication (per-niche ligand-receptor
# probabilities across NF/pRV/RVF); its section reads the shipped table at
# dependencies/shared/Xenium/xenium_cellchat_niche_communication.csv.
#
# REMOVED FROM v60 (sections commented out below; not cited as a Data S):
#   - ChEA TF enrichment              (was Data S21; xlsx deleted)
#   - Bulk RNA-seq per-sample QC       (was Data S24, extra)
#   - Differential cell-type abundance (was Data S25, extra)
#   - Module-eigengene projection      (was Data S26, extra)
#   - CollecTRI/decoupleR TF activity  (was Data S27, extra)
#   - Cardiac fibrosis quantification  (was Data S28, extra)
# The 5 extra xlsx remain on disk (Data_S24..S28) but are no longer regenerated.
#
# TWO MODES
#   (default) REBUILD = FALSE
#     Load cached CSVs / RDS under ./dependencies/shared/ and emit xlsx.
#
#   REBUILD = TRUE
#     Re-derive every CSV de novo from raw inputs (Seurat RDS objects, kallisto
#     tximport, hdWGCNA pipeline) and THEN emit xlsx.  Per-section toggles
#     (REBUILD_BULK_DEGS, REBUILD_SN_PSEUDOBULK, etc.) let you rebuild a
#     single supplement without rebuilding the rest.
#
# Cached CSV outputs are written back to ./dependencies/shared/ so future runs
# are fast.  xlsx outputs land in ./output/Supp_Data/.
#
# Inputs (for rebuild)
#   ./dependencies/shared/BulkRNA/{counts,metadata}.csv               (bulk kallisto)
#   ./dependencies/shared/RV_bulkRNASeq_seurat.rds                    (hdWGCNA host)
#   ./dependencies/shared/snRV_ref.rds                                (snRNA Seurat)
#   ./dependencies/shared/{cm,fb,ec,myeloid_new,mural}_subclust*.rds  (sn subclust)
#   ./dependencies/shared/Xenium/                                     (Xenium subclust)
#   ./dependencies/shared/PAB_data_clean.rds                          (murine PAB)
#   ./dependencies/shared/all_peds_data.rds                           (pediatric)
#   ./dependencies/shared/GSE240921_processed-data-human.csv          (PAH)
#
# Helpers are recycled from
#   new_scripts/Figure_1.R        (bulk DESeq2 + hdWGCNA + trend genes)
#   new_scripts/Figure_3.R        (.run_pseudobulk_deseq2)
#   new_scripts/Figure_6.R        (ChEA / enrichR)
#   additional_scripts/snRNAReanalysis.r (run_pseudobulk_deseq2, 'RNA' assay)
#   additional_scripts/AnalysisPAH.R      (GSE240921 PAH DESeq2)
##############################################################################

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

# ==========================================================================
# REBUILD FLAGS
# ==========================================================================
REBUILD                   <- FALSE   # master switch
REBUILD_BULK_DEGS         <- REBUILD # S1  — NF/pRV/RVF DESeq tables
REBUILD_BULK_ENRICH       <- REBUILD # S2
REBUILD_BULK_WGCNA        <- REBUILD # S4  — bulk_heart_modules
REBUILD_BULK_TRENDS       <- REBUILD # S3  — NF_2_pRV_*_2_RVF_*
REBUILD_PAH               <- REBUILD # S5
REBUILD_SN_ANNOT          <- REBUILD # S6
REBUILD_SN_SUBCLUST_MARKS <- REBUILD # S7  — FindAllMarkers per subclust obj
REBUILD_SN_PSEUDOBULK_L1  <- REBUILD # S8  — sn lineage pseudobulk
REBUILD_XEN_PSEUDOBULK_L1 <- REBUILD # S11
REBUILD_XEN_SUBCLUST_MARKS<- REBUILD # S12
REBUILD_XEN_NICHE         <- REBUILD # S13 — composition tables
REBUILD_EC_HDWGCNA        <- REBUILD # S18
REBUILD_PAB               <- REBUILD # S19
REBUILD_CX3CR1            <- TRUE    # S20 — FORCE: add LV_PAB_vs_LV_Sham contrast
REBUILD_PEDS              <- REBUILD # S21 / S22
REBUILD_CHEA              <- REBUILD # (removed in v60; section commented out)
REBUILD_SN_PSEUDOBULK_L2  <- REBUILD # S9

# Parallelism ------------------------------------------------------------- #
PARALLEL      <- FALSE
N_CORES       <- 4

# ==========================================================================
# PATHS
# ==========================================================================
.shared   <- '/Users/ikuz/Documents/RV_Atlas/dependencies/shared'
.bulkdir  <- .shared
.outdir   <- '/Users/ikuz/Documents/RV_Atlas/output/Supp_Data'
.cachedir <- '/Users/ikuz/Documents/RV_Atlas/output'
.xcache   <- file.path(.shared, 'Xenium')
dir.create(.outdir,   showWarnings = FALSE, recursive = TRUE)
dir.create(.cachedir, showWarnings = FALSE, recursive = TRUE)

# ==========================================================================
# HELPERS
# ==========================================================================
.write_xlsx <- function(sheets, out_name) {
  path <- file.path(.outdir, out_name)
  sheets <- Filter(function(x) !is.null(x) && (is.data.frame(x) || is.list(x)) && NROW(x) > 0, sheets)
  if (length(sheets) == 0) { message('SKIP empty: ', out_name); return(invisible(NULL)) }
  openxlsx::write.xlsx(sheets, file = path, overwrite = TRUE)
  message('WROTE ', out_name, '  (', length(sheets), ' sheet(s))')
}
.read_csv_safe <- function(p) {
  p <- path.expand(p)
  if (!file.exists(p)) { message('MISSING: ', p); return(NULL) }
  tryCatch(readr::read_csv(p, show_col_types = FALSE),
           error = function(e) { message('READ FAIL: ', p, '  (', conditionMessage(e), ')'); NULL })
}
.section <- function(n, title) message(sprintf('\n============ Data S%02d - %s ============', n, title))
.read_panel <- function() {
  p <- .read_csv_safe(file.path(.shared, 'XeniumPanel.csv'))
  if (is.null(p)) p <- .read_csv_safe(file.path(.shared, 'panel_genes.csv'))
  p
}
# Canonical subname mapping (manuscript naming, space style).
#  - FB: static Fb#->descriptive map per Supplementary_Figure_3.R (the
#    `.cm_subname_map`/`Fb` block); fb_subclust.rds still uses old Fb#.
#  - CM: per-cell barcode lookup from cm_subclust_new_new.rds (10 labels:
#    CM Baseline / BetaMHC / HAND2 / HMGCS2 / HTR4 / KCNJ3 / MYH6 / NPP /
#    RORA / XIRP), precomputed once into dependencies/shared/
#    cm_subclust_labels.csv to avoid re-loading the 1.3 GB rds.
.FB_SUBNAME_RMAP <- c(
  Fb1 = 'Fb Resident',     Fb2 = 'Fb Adventitial',
  Fb3 = 'Fb Elastogenic',  Fb4 = 'Fb Interstitial',
  Fb5 = 'Fb Stressed',     Fb6 = 'Fb Pro-fibrotic',
  Fb7 = 'Fb Anti-fibrotic')
.CM_LABELS_CSV <- file.path(.shared, 'cm_subclust_labels.csv')
.CM_OLD2NEW_CSV <- file.path(.shared, 'cm_oldnum_to_new.csv')
.CM_ORPHAN_LABEL <- 'CM unassigned'
.cm_lab_env <- new.env(parent = emptyenv())
.load_cm_labels <- function() {
  if (!is.null(.cm_lab_env$lab)) return(invisible(.cm_lab_env$lab))
  v <- character(0); o2n <- character(0)
  if (file.exists(.CM_LABELS_CSV)) {
    d <- utils::read.csv(.CM_LABELS_CSV, stringsAsFactors = FALSE)
    v <- as.character(d$Subnames); names(v) <- as.character(d$cell_id)
  } else message('  cm_subclust_labels.csv missing — CM barcode lookup disabled')
  if (file.exists(.CM_OLD2NEW_CSV)) {
    d2 <- utils::read.csv(.CM_OLD2NEW_CSV, stringsAsFactors = FALSE)
    o2n <- as.character(d2$new); names(o2n) <- as.character(d2$old)
  } else message('  cm_oldnum_to_new.csv missing — CM cluster-level fallback disabled')
  .cm_lab_env$lab <- v; .cm_lab_env$o2n <- o2n
  invisible(v)
}
## Relabel rules for snRNA subname strings (manuscript naming, space style):
##   Fb1..Fb7  -> .FB_SUBNAME_RMAP (static, from SuppFig3)
##   Cm#  + cell_id, in cm_subclust cache         -> new label from cache
##   Cm#  + cell_id, NOT in cache (orphan)        -> CM_ORPHAN_LABEL
##   Cm#  without cell_id (cluster-level table)   -> mode of new labels
##                                                   per old Cm# (cm_oldnum_to_new.csv)
.relabel_subnames <- function(raw, cell_id = NULL) {
  raw <- as.character(raw); new <- raw
  fb <- raw %in% names(.FB_SUBNAME_RMAP)
  if (any(fb)) new[fb] <- .FB_SUBNAME_RMAP[raw[fb]]
  cm <- grepl('^Cm[0-9]+$', raw)
  if (any(cm)) {
    .load_cm_labels()
    if (!is.null(cell_id)) {
      bc  <- as.character(cell_id)[cm]
      hit <- bc %in% names(.cm_lab_env$lab)
      new[cm][hit]  <- .cm_lab_env$lab[bc[hit]]
      new[cm][!hit] <- .CM_ORPHAN_LABEL
    } else if (length(.cm_lab_env$o2n)) {
      ok <- raw[cm] %in% names(.cm_lab_env$o2n)
      new[cm][ok]  <- .cm_lab_env$o2n[raw[cm][ok]]
      new[cm][!ok] <- .CM_ORPHAN_LABEL
    }
  }
  new
}
.rename_subnames <- function(df, src_col = 'Subnames_manual') {
  if (is.null(df) || !src_col %in% names(df)) return(df)
  df[[src_col]] <- .relabel_subnames(df[[src_col]],
    if ('cell_id' %in% names(df)) df$cell_id else NULL)
  df
}
.apply_rmap_inplace <- function(df, src_col) {
  if (is.null(df) || !src_col %in% names(df)) return(df)
  df[[src_col]] <- .relabel_subnames(df[[src_col]],
    if ('cell_id' %in% names(df)) df$cell_id else NULL)
  df
}
.stub <- function(n, title, note) {
  .section(n, title); message('TODO: ', note)
  .write_xlsx(list(README = data.frame(note = note, title = title)),
              sprintf('Data_S%02d_TODO.xlsx', n))
}
.require <- function(pkgs) {
  ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (!all(ok)) { message('MISSING packages: ', paste(pkgs[!ok], collapse = ', ')); return(FALSE) }
  for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
  TRUE
}

# --------------------------------------------------------------------------
# BUILDERS
# --------------------------------------------------------------------------

# ---- Bulk DEGs (S1) ------------------------------------------------------ #
# Recycles Figure_1.R lines ~47-286 (tximport -> SVA -> DESeq2 -> lfcShrink).
.build_bulk_degs <- function(force = FALSE) {
  out_nf_prv <- file.path(.bulkdir, 'NF_vs_pRV_deseq.csv')
  out_nf_rvf <- file.path(.bulkdir, 'NF_vs_RVF_deseq.csv')
  out_prv_rvf<- file.path(.bulkdir, 'pRV_vs_RVF_deseq.csv')
  if (!force && all(file.exists(c(out_nf_prv, out_nf_rvf, out_prv_rvf)))) return(invisible())
  message('BUILD: bulk DEGs (S1)')
  if (!.require(c('DESeq2', 'tximport', 'biomaRt', 'sva'))) return(invisible())

  bulk <- read.csv(file.path(.shared, 'BulkRNA/counts.csv'))
  meta <- read.csv(file.path(.shared, 'BulkRNA/metadata.csv'))
  toDel <- seq(1, nrow(meta), 2); meta <- meta[-toDel, ]

  .ctx <- file.path(.cachedir, 'fig1_tx2gene_cache.rds')
  if (file.exists(.ctx)) {
    tx2gene <- readRDS(.ctx)
  } else {
    mart <- biomaRt::useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl',
                             host = 'useast.ensembl.org')
    res  <- biomaRt::getBM(
      attributes = c('ensembl_transcript_id_version', 'ensembl_gene_id',
                     'external_transcript_name', 'external_gene_name'),
      filters = 'ensembl_transcript_id_version', values = bulk[, 1], mart = mart)
    tx2gene <- res[, c(1, 4)]; tx2gene <- tx2gene[tx2gene$external_gene_name != '', ]
    saveRDS(tx2gene, .ctx)
  }
  .cxi <- file.path(.cachedir, 'fig1_txi_cache.rds')
  if (file.exists(.cxi)) {
    txi.kallisto <- readRDS(.cxi)
  } else {
    message('tximport pass requires original kallisto h5s; using counts matrix')
    mat <- as.matrix(bulk[, -1]); rownames(mat) <- bulk[, 1]
    mat <- stats::aggregate(mat, by = list(gene = tx2gene$external_gene_name[match(rownames(mat), tx2gene$ensembl_transcript_id_version)]), FUN = sum)
    rn <- mat$gene; mat <- as.matrix(mat[, -1]); rownames(mat) <- rn
    txi.kallisto <- list(counts = round(mat), abundance = mat, length = mat + 1,
                         countsFromAbundance = 'no')
  }
  bulk.meta <- meta
  rownames(bulk.meta) <- colnames(txi.kallisto$counts)[seq_len(nrow(bulk.meta))]

  # SVA
  .csv_cache <- file.path(.cachedir, 'fig1_svseq_cache.rds')
  if (file.exists(.csv_cache)) {
    svseq <- readRDS(.csv_cache)
  } else {
    dat <- round(txi.kallisto$counts); dat <- dat[rowMeans(dat) > 1, ]
    mod  <- stats::model.matrix(~ category, data = bulk.meta)
    mod0 <- stats::model.matrix(~ 1, data = bulk.meta)
    svseq <- sva::svaseq(dat, mod, mod0, n.sv = 21)
    saveRDS(svseq, .csv_cache)
  }
  for (i in seq_len(ncol(svseq$sv))) bulk.meta[[paste0('SV', i)]] <- svseq$sv[, i]

  ddsSE <- DESeq2::DESeqDataSetFromTximport(
    txi.kallisto, bulk.meta,
    design = as.formula(paste0('~ category + ', paste0('SV', 1:21, collapse = ' + '))))
  .ddc <- file.path(.cachedir, 'fig1_deseq_fitted_cache.rds')
  if (file.exists(.ddc)) { ddsSE <- readRDS(.ddc) } else { ddsSE <- DESeq2::DESeq(ddsSE); saveRDS(ddsSE, .ddc) }

  sh <- function(a, b) as.data.frame(DESeq2::lfcShrink(ddsSE, contrast = c('category', a, b), type = 'ashr')) %>%
    tibble::rownames_to_column('gene')
  write.csv(sh('NF',  'pRV'),  out_nf_prv,  row.names = FALSE)
  write.csv(sh('NF',  'RVF'),  out_nf_rvf,  row.names = FALSE)
  write.csv(sh('pRV', 'RVF'),  out_prv_rvf, row.names = FALSE)
  message('  wrote ', out_nf_prv, ', ', out_nf_rvf, ', ', out_prv_rvf)
}

# ---- Bulk WGCNA modules (S3) --------------------------------------------- #
# Recycles Figure_1.R lines ~330-395.
.build_bulk_wgcna <- function(force = FALSE) {
  out <- file.path(.shared, 'bulk_heart_modules.csv')
  if (!force && file.exists(out)) return(invisible())
  message('BUILD: bulk WGCNA (S3)')
  if (!.require(c('Seurat', 'hdWGCNA'))) return(invisible())
  seurat_obj <- readRDS(file.path(.shared, 'RV_bulkRNASeq_seurat.rds'))
  seurat_obj <- hdWGCNA::SetupForWGCNA(seurat_obj, gene_select = 'fraction',
                                       fraction = 0.05, wgcna_name = 'bulk_RV')
  seurat_obj <- hdWGCNA::MetacellsByGroups(seurat_obj, group.by = 'category',
                                           k = 20, min_cells = 10, ident.group = 'category')
  seurat_obj <- hdWGCNA::NormalizeMetacells(seurat_obj)
  seurat_obj <- hdWGCNA::SetDatExpr(seurat_obj, assay = 'RNA', slot = 'data')
  seurat_obj <- hdWGCNA::TestSoftPowers(seurat_obj, networkType = 'signed')
  seurat_obj <- hdWGCNA::ConstructNetwork(seurat_obj, tom_name = 'bulk_RV',
                                          tom_dir = file.path(.cachedir, 'hdWGCNA_TOM'),
                                          overwrite_tom = TRUE)
  seurat_obj <- hdWGCNA::ModuleEigengenes(seurat_obj, group.by.vars = NULL)
  modules <- hdWGCNA::GetModules(seurat_obj)
  write.csv(modules, row.names = FALSE, quote = FALSE, file = out)
  message('  wrote ', out)
}

# ---- Bulk trend genes (S4) ---------------------------------------------- #
# Recycles Figure_1.R lines ~855-860.  Needs the three bulk DEG tables.
# Thresholds match Supplementary_Figure_1.R (padj < 0.1, sign of log2FoldChange only,
# no magnitude cutoff) so the S4 trend CSVs reproduce the shipped/figure gene sets.
.build_bulk_trends <- function(force = FALSE, padj_cut = 0.1, lfc = 0) {
  targets <- c('NF_2_pRV_up_2_RVF_up.csv', 'NF_2_pRV_down_2_RVF_down.csv',
               'NF_2_pRV_up_2_RVF_down.csv', 'NF_2_pRV_down_2_RVF_up.csv',
               'NF_2_pRV_flat_2_RVF_up.csv', 'NF_2_pRV_flat_2_RVF_down.csv')
  all_exist <- all(file.exists(file.path(.bulkdir, targets)))
  if (!force && all_exist) return(invisible())
  message('BUILD: bulk stepwise trends (S4)')
  .build_bulk_degs(force = FALSE)
  d1 <- .read_csv_safe(file.path(.bulkdir, 'NF_vs_pRV_deseq.csv'))
  d2 <- .read_csv_safe(file.path(.bulkdir, 'NF_vs_RVF_deseq.csv'))
  d3 <- .read_csv_safe(file.path(.bulkdir, 'pRV_vs_RVF_deseq.csv'))
  gcol <- intersect(c('gene', 'Gene', 'symbol'), names(d1))[1]
  if (is.na(gcol)) gcol <- names(d1)[1]   # readr names the unnamed gene column '...1'
  .sig <- function(d, dir) {
    d %>% filter(!is.na(padj), padj < padj_cut, {if (dir == 'up') log2FoldChange > lfc else log2FoldChange < -lfc}) %>%
      pull(!!sym(gcol))
  }
  nf.prv.up  <- .sig(d1, 'up');  nf.prv.dn  <- .sig(d1, 'down')
  nf.rvf.up  <- .sig(d2, 'up');  nf.rvf.dn  <- .sig(d2, 'down')
  prv.rvf.up <- .sig(d3, 'up');  prv.rvf.dn <- .sig(d3, 'down')

  sets <- list(
    NF_2_pRV_up_2_RVF_up     = intersect(intersect(nf.prv.dn, nf.rvf.dn), prv.rvf.dn),
    NF_2_pRV_down_2_RVF_down = intersect(intersect(nf.prv.up, nf.rvf.up), prv.rvf.up),
    NF_2_pRV_up_2_RVF_down   = setdiff(intersect(nf.prv.dn, prv.rvf.up), nf.rvf.dn),
    NF_2_pRV_down_2_RVF_up   = setdiff(intersect(nf.prv.up, prv.rvf.dn), nf.rvf.up),
    NF_2_pRV_flat_2_RVF_up   = setdiff(intersect(nf.rvf.up, prv.rvf.up), nf.prv.up),
    NF_2_pRV_flat_2_RVF_down = setdiff(intersect(nf.rvf.dn, prv.rvf.dn), nf.prv.dn)
  )
  for (nm in names(sets)) {
    write.csv(data.frame(gene = sets[[nm]]), file.path(.bulkdir, paste0(nm, '.csv')), row.names = FALSE)
  }
  message('  wrote ', length(sets), ' trend CSVs')
}

# ---- PAH DEGs (S5) ------------------------------------------------------- #
# Recycles AnalysisPAH.R.
.build_pah <- function(force = FALSE) {
  out1 <- file.path(.shared, 'pah.nf.vs.prv.csv')
  out2 <- file.path(.shared, 'pah.nf.vs.rvf.csv')
  out3 <- file.path(.shared, 'pah.prv.vs.rvf.csv')
  if (!force && all(file.exists(c(out1, out2, out3)))) return(invisible())
  message('BUILD: PAH DEGs (S5)')
  if (!.require(c('DESeq2', 'sva', 'stringr'))) return(invisible())
  pah_raw    <- read.csv(file.path(.shared, 'GSE240921_processed-data-human.csv'))
  gene.annot <- read.csv(file.path(.shared, 'Human.GRCh38.p13.annot.tsv'), sep = '\t')
  pah_raw$names <- gene.annot$Symbol[match(pah_raw$id, gene.annot$EnsemblGeneID)]
  pah_raw$names[is.na(pah_raw$names)] <- pah_raw$id[is.na(pah_raw$names)]
  subj.names <- colnames(pah_raw)[2:41]
  subj.group <- stringr::str_split_i(subj.names, '_', 1)
  subj.group[subj.group == 'RV.Normal']  <- 'Control'
  subj.group[subj.group == 'RV.Compen']  <- 'Compensated'
  subj.group[subj.group == 'RV.Failing'] <- 'Decompensated'
  pah_proc <- stats::aggregate(pah_raw[, 2:41], by = list(gene = pah_raw$names), FUN = sum)
  rn <- pah_proc$gene; pah_proc <- as.matrix(pah_proc[, -1]); rownames(pah_proc) <- rn
  coldata <- data.frame(group = factor(subj.group),
                        row.names = subj.names, stringsAsFactors = FALSE)
  dat  <- round(pah_proc); dat <- dat[rowMeans(dat) > 1, ]
  mod  <- stats::model.matrix(~ group, data = coldata)
  mod0 <- stats::model.matrix(~ 1,     data = coldata)
  sv   <- sva::svaseq(dat, mod, mod0, n.sv = 5)
  for (i in 1:5) coldata[[paste0('SV', i)]] <- sv$sv[, i]
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(pah_proc), colData = coldata,
                                        design = ~ group + SV1 + SV2 + SV3 + SV4 + SV5)
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- dds[rowSums(DESeq2::counts(dds, normalized = TRUE) >= 5) >= 3, ]
  dds <- DESeq2::DESeq(dds)
  .ss <- function(a, b) as.data.frame(DESeq2::lfcShrink(dds, contrast = c('group', a, b), type = 'ashr')) %>%
    tibble::rownames_to_column('gene')
  write.csv(.ss('Control',     'Compensated'),   out1, row.names = FALSE)
  write.csv(.ss('Control',     'Decompensated'), out2, row.names = FALSE)
  write.csv(.ss('Compensated', 'Decompensated'), out3, row.names = FALSE)
  message('  wrote 3 PAH CSVs')
}

# ---- PAB pseudobulk (S21) ----------------------------------------------- #
# Per-mouse, per-lineage aggregation -> DESeq2 with `group` (Nor/Mod/Sev) and
# all 3 pairwise contrasts. Writes a long CSV to
# dependencies/shared/pab_pseudobulk_lineage_deseq2.csv that S21 picks up.
.build_pab_pseudobulk <- function(force = FALSE) {
  out <- file.path(.shared, 'pab_pseudobulk_lineage_deseq2.csv')
  if (!force && file.exists(out)) return(invisible())
  rds <- file.path(.shared, 'PAB_data_clean.rds')
  if (!file.exists(rds)) { message('  PAB_data_clean.rds not found'); return(invisible()) }
  if (!.require(c('Seurat', 'DESeq2', 'Matrix'))) return(invisible())
  message('BUILD: PAB pseudobulk (lineage, group=Nor/Mod/Sev)')
  obj <- readRDS(rds)
  md  <- obj@meta.data
  ct  <- md$Names; grp <- md$orig.ident; pat <- md$patient
  keep <- !is.na(ct) & !is.na(grp) & !is.na(pat) & nchar(as.character(ct)) > 0
  ct <- as.character(ct[keep]); grp <- as.character(grp[keep]); pat <- as.character(pat[keep])
  cells <- colnames(obj)[keep]
  cnts <- Seurat::GetAssayData(obj, assay = 'RNA', layer = 'counts')[, cells, drop = FALSE]
  rm(obj); gc()
  contrasts <- list(
    Mod_vs_Nor = c('Mod', 'Nor'),
    Sev_vs_Nor = c('Sev', 'Nor'),
    Sev_vs_Mod = c('Sev', 'Mod'))
  rows <- list()
  for (this_ct in sort(unique(ct))) {
    sel <- ct == this_ct
    if (sum(sel) < 30) next
    key <- paste(pat[sel], grp[sel], sep = '__')
    fk  <- factor(key)
    pb  <- as.matrix(cnts[, sel, drop = FALSE] %*%
                     Matrix::sparse.model.matrix(~ 0 + fk))
    colnames(pb) <- sub('^fk', '', colnames(pb))
    meta <- data.frame(
      sample = colnames(pb),
      patient = sub('__.*$', '', colnames(pb)),
      group   = factor(sub('^.*__', '', colnames(pb)),
                       levels = c('Nor', 'Mod', 'Sev')),
      stringsAsFactors = FALSE)
    rownames(meta) <- meta$sample
    pb <- pb[, meta$sample, drop = FALSE]
    if (length(levels(droplevels(meta$group))) < 2) next
    dds <- tryCatch(DESeq2::DESeqDataSetFromMatrix(countData = round(pb), colData = meta,
                                                   design = ~ group),
                    error = function(e) NULL)
    if (is.null(dds)) next
    dds <- tryCatch(DESeq2::DESeq(dds, quiet = TRUE), error = function(e) NULL)
    if (is.null(dds)) next
    for (cn in names(contrasts)) {
      lv <- contrasts[[cn]]
      if (!all(lv %in% levels(droplevels(meta$group)))) next
      res <- tryCatch(as.data.frame(DESeq2::results(dds, contrast = c('group', lv[1], lv[2]))),
                      error = function(e) NULL)
      if (is.null(res) || nrow(res) == 0) next
      res$gene <- rownames(res)
      res$cell_type <- this_ct
      res$comparison <- cn
      rows[[paste(this_ct, cn, sep = '__')]] <-
        res[, c('gene','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj',
                'cell_type','comparison')]
    }
    message('  ', this_ct, ': ', sum(sel), ' cells -> n_samples=', ncol(pb))
  }
  rm(cnts); gc()
  if (length(rows) == 0) { message('  no contrasts computed'); return(invisible()) }
  df <- do.call(rbind, rows); rownames(df) <- NULL
  write.csv(df, out, row.names = FALSE)
  message('  wrote ', out, ' (', nrow(df), ' rows)')
}

# ---- snRNA pseudobulk (S8 / S22) ---------------------------------------- #
# Recycles additional_scripts/snRNAReanalysis.r run_pseudobulk_deseq2.
.run_pseudobulk_deseq <- function(obj, ident_col, assay_name = 'RNA') {
  if (!.require(c('Seurat', 'DESeq2'))) return(NULL)
  obj$pseudobulk_id <- paste0(obj$patient, '..', obj@meta.data[[ident_col]])
  pseudo_meta <- unique(obj@meta.data[, c('patient', 'group', ident_col, 'pseudobulk_id')])
  pseudo_meta$group <- factor(pseudo_meta$group, levels = c('NF', 'pRV', 'RVF'))
  rownames(pseudo_meta) <- pseudo_meta$pseudobulk_id
  pseudo_counts <- Seurat::AggregateExpression(
    obj, assays = assay_name, group.by = 'pseudobulk_id',
    slot = 'counts', return.seurat = FALSE)[[assay_name]]
  colnames(pseudo_counts) <- gsub('-', '_', sub('^g', '', colnames(pseudo_counts)))
  pseudo_meta$pseudobulk_id <- gsub('-', '_', sub('^g', '', pseudo_meta$pseudobulk_id))
  rownames(pseudo_meta)     <- pseudo_meta$pseudobulk_id
  out <- list()
  for (ct in unique(pseudo_meta[[ident_col]])) {
    ct_s <- pseudo_meta$pseudobulk_id[pseudo_meta[[ident_col]] == ct]
    ct_s <- ct_s[ct_s %in% colnames(pseudo_counts)]
    if (length(ct_s) < 2) next
    cc <- pseudo_counts[, ct_s, drop = FALSE]; cm <- pseudo_meta[ct_s, , drop = FALSE]
    if (length(unique(cm$group)) < 2) next
    cc <- cc[rowSums(cc >= 5) >= 2, , drop = FALSE]; if (nrow(cc) < 10) next
    dds <- tryCatch(DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(countData = round(cc), colData = cm, design = ~ group)),
                    error = function(e) NULL)
    if (is.null(dds)) next
    for (comp in list(c('RVF', 'NF'), c('pRV', 'NF'), c('RVF', 'pRV'))) {
      r <- tryCatch(as.data.frame(DESeq2::results(dds, contrast = c('group', comp[1], comp[2]))) %>%
                      tibble::rownames_to_column('gene') %>%
                      dplyr::mutate(cell_type = ct, comparison = paste0(comp[1], '_vs_', comp[2])) %>%
                      dplyr::arrange(padj),
                    error = function(e) NULL)
      if (!is.null(r)) out[[paste0(ct, '_', comp[1], '_vs_', comp[2])]] <- r
    }
  }
  dplyr::bind_rows(out)
}
.build_sn_pseudobulk <- function(level = c('lineage', 'sublineage'), force = FALSE) {
  level <- match.arg(level)
  label_col <- if (level == 'lineage') 'Names' else 'Subnames_manual'
  out <- file.path(.shared, sprintf('sn_pseudobulk_%s_deseq2.csv', level))
  if (!force && file.exists(out)) return(invisible())
  message('BUILD: snRNA pseudobulk (', level, ')')
  rds <- file.path(.shared, 'snRV_ref.rds')
  if (!file.exists(rds)) { message('  snRV_ref.rds not found'); return(invisible()) }
  obj <- readRDS(rds)
  df  <- .run_pseudobulk_deseq(obj, ident_col = label_col, assay_name = 'RNA')
  rm(obj); gc()
  if (!is.null(df)) { write.csv(df, out, row.names = FALSE); message('  wrote ', out) }
}

# ---- snRNA subclustering markers (S7) ----------------------------------- #
.build_sn_subclust_marks <- function(force = FALSE) {
  targets <- list(
    cm      = list(rds = file.path(.shared, 'cm_subclust_new_new.rds'),
                   csv = file.path(.shared, 'cm_subclust_marks.csv'),
                   col = 'Subnames',
                   findargs = list(only.pos = FALSE, min.pct = 0.01,
                                   logfc.threshold = 0.25)),
    fb      = list(rds = file.path(.shared, 'fb_subclust.rds'),
                   csv = file.path(.shared, 'fb_subclust_marks.csv'),
                   col = 'fb_subtype'),
    ec      = list(rds = file.path(.shared, 'ec_subclust.rds'),
                   csv = file.path(.shared, 'ec_subclust_marks.csv'),
                   col = 'ec_subtype'),
    myeloid = list(rds = file.path(.shared, 'myeloid_subclust_new.rds'),
                   csv = file.path(.shared, 'myeloid_subclust_marks.csv'),
                   col = 'myeloid_subtype'),
    mural   = list(rds = file.path(.shared, 'pc_sm_subclust.rds'),
                   csv = file.path(.shared, 'mural_subclust_marks.csv'),
                   col = 'mural_subtype')
  )
  if (!.require(c('Seurat'))) return(invisible())
  for (nm in names(targets)) {
    t <- targets[[nm]]
    if (!force && file.exists(t$csv)) next
    if (!file.exists(t$rds)) { message('  skip ', nm, ' (', t$rds, ' missing)'); next }
    message('BUILD: ', nm, ' FindAllMarkers')
    obj <- readRDS(t$rds)
    # Choose identity:
    #   1. If t$col is a meta.data column and looks like cluster labels, use it
    #   2. Else use the object's existing active.ident (these subclust rds files
    #      were saved with cluster labels as the active ident — e.g. Fb1..Fb7).
    # Reject choices with > 30 unique levels (likely patient IDs).
    set_from_meta <- FALSE
    if (!is.null(t$col) && t$col %in% names(obj@meta.data)) {
      vals <- as.character(obj@meta.data[[t$col]])
      if (length(unique(vals)) <= 30) {
        Seurat::Idents(obj) <- vals; set_from_meta <- TRUE
      }
    }
    if (!set_from_meta) {
      cur <- as.character(Seurat::Idents(obj))
      if (length(unique(cur)) > 30) {
        message('  refusing ', nm, ': active.ident has ',
                length(unique(cur)), ' levels (looks like patient/cell IDs)')
        rm(obj); gc(); next
      }
      message('  using active.ident (', length(unique(cur)), ' levels) for ', nm)
    }
    default_args <- list(only.pos = TRUE, min.pct = 0.25,
                         logfc.threshold = 0.25, verbose = FALSE)
    final_args <- modifyList(default_args,
                             if (is.null(t$findargs)) list() else t$findargs)
    final_args$object <- obj
    mk <- do.call(Seurat::FindAllMarkers, final_args)
    write.csv(mk, t$csv, row.names = FALSE); rm(obj); gc()
    message('  wrote ', t$csv)
  }
}

# ---- Xenium pseudobulk (S9) --------------------------------------------- #
# Recycles Figure_3.R lines ~1773-1820.
.build_xenium_pseudobulk <- function(level = c('lineage', 'sublineage'), force = FALSE) {
  level <- match.arg(level)
  out   <- file.path(.shared, sprintf('xenium_pseudobulk_%s_deseq2.csv', level))
  if (!force && file.exists(out)) return(invisible())
  message('BUILD: Xenium pseudobulk (', level, ')')
  rds <- file.path(.xcache, 'xenium_obj_subclustered.rds')
  if (!file.exists(rds)) rds <- file.path(.xcache, 'xenium_corrected.rds')
  if (!file.exists(rds)) { message('  xenium corrected/subclustered rds not found'); return(invisible()) }
  obj <- readRDS(rds)
  # Pseudobulk uses expression counts only; drop stale spatial FOVs so subset() does
  # not trip the FOV-class validation error this object throws under current Seurat.
  obj@images <- list()
  # Drop subtype-Unassigned cells for BOTH levels, matching the Data_S11 generator.
  # (This was previously gated on level=='sublineage'; omitting it at lineage level
  # left ~54% Unassigned cells in each lineage's pseudobulk, so the lineage CSV did
  # not reproduce -- log2FC Spearman ~0.71 vs the stored table. With the filter it
  # reproduces at ~0.999. See analysis/figure_reproducibility/reports/DEG_ACCURACY_AND_METHODS.md.)
  obj <- subset(obj, subset = cell_types_subclustering != 'Unassigned')
  # The Xenium assay now also carries SpaGE-imputed genes (~17k features); restrict to
  # the measured probe panel (477 genes) so DE is not run on imputed features. The panel
  # list is a fixed design fact; if the object is ever restored to the 477-gene assay
  # this intersect is a harmless no-op.
  .panel_f <- file.path(.xcache, 'xenium_panel_genes.txt')
  if (file.exists(.panel_f)) {
    panel  <- readLines(.panel_f)
    keep_g <- intersect(rownames(obj), panel)
    if (length(keep_g) >= 100) obj <- subset(obj, features = keep_g)
  } else message('  NOTE: xenium_panel_genes.txt missing -> DE may include imputed genes')
  id_col <- if (level == 'lineage') 'cell_type_rctd_doublet' else 'cell_types_subclustering'
  if (!id_col %in% names(obj@meta.data)) id_col <- names(obj@meta.data)[grepl('cell_type', names(obj@meta.data))][1]
  df <- .run_pseudobulk_deseq(obj, ident_col = id_col, assay_name = 'Xenium')
  rm(obj); gc()
  if (!is.null(df)) { write.csv(df, out, row.names = FALSE); message('  wrote ', out) }
}

# ---- Xenium subclustering markers (S10) --------------------------------- #
.build_xenium_subclust_marks <- function(force = FALSE) {
  out <- file.path(.xcache, 'xenium_corrected_findallmarkers.csv')
  if (!force && file.exists(out) && file.size(out) > 100) return(invisible())
  message('BUILD: Xenium FindAllMarkers (S10)')
  rds_marks <- file.path(.xcache, 'xenium_corrected_findallmarkers.rds')
  if (file.exists(rds_marks)) {
    mk <- readRDS(rds_marks)
    if (is.data.frame(mk)) {
      write.csv(mk, out, row.names = FALSE)
      message('  converted ', rds_marks, ' -> ', out)
      return(invisible())
    }
  }
  rds <- file.path(.xcache, 'xenium_obj_subclustered.rds')
  if (!file.exists(rds)) { message('  xenium subclust rds not found'); return(invisible()) }
  if (!.require(c('Seurat'))) return(invisible())
  obj <- readRDS(rds)
  .sub <- subset(obj, subset = cell_types_subclustering != 'Unassigned')
  Seurat::Idents(.sub) <- .sub$cell_types_subclustering
  mk <- Seurat::FindAllMarkers(.sub, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = 0.25, verbose = FALSE)
  write.csv(mk, out, row.names = FALSE); rm(obj, .sub); gc()
  message('  wrote ', out)
}

# ==========================================================================
# MAIN PIPELINE — Data S01..S22 (v60); removed sections commented out
# ==========================================================================

# ---------- S01 : Bulk RNA-seq DEGs ---------------------------------------- #
try({
  .section(01, 'Bulk RNA-seq DEGs')
  if (REBUILD_BULK_DEGS) .build_bulk_degs(force = TRUE)
  sheets <- list(
    NF_vs_pRV  = .read_csv_safe(file.path(.bulkdir, 'NF_vs_pRV_deseq.csv')),
    NF_vs_RVF  = .read_csv_safe(file.path(.bulkdir, 'NF_vs_RVF_deseq.csv')),
    pRV_vs_RVF = .read_csv_safe(file.path(.bulkdir, 'pRV_vs_RVF_deseq.csv'))
  )
  ## Ensure unnamed first column (gene IDs from readr's `...1` artifact)
  ## carries an explicit header.
  sheets <- lapply(sheets, function(d) {
    if (is.null(d) || ncol(d) == 0) return(d)
    if (names(d)[1] %in% c('..1', '...1', '', 'X', 'X.1', 'V1'))
      names(d)[1] <- 'gene'
    d
  })
  .write_xlsx(sheets, 'Data_S01_bulk_DEGs.xlsx')
})

# ---------- S02 : Bulk pathway enrichment ---------------------------------- #
try({
  .section(02, 'Bulk pathway enrichment')
  cache_path <- file.path(.cachedir, 'S02_bulk_enrichment_cache.rds')
  if (!REBUILD_BULK_ENRICH && file.exists(cache_path)) {
    sheets <- readRDS(cache_path)
  } else {
    if (!.require('enrichR')) { .stub(02, 'Bulk pathway enrichment', 'install enrichR then rerun'); sheets <- NULL }
    else {
      # v57 promises: Reactome, GO Biological Processes, GTEx Aging Signatures
      dbs <- c('GO_Biological_Process_2023', 'Reactome_2022',
               'GTEx_Aging_Signatures_2021')
      run_enrich <- function(deg, lfc_col='log2FoldChange', padj_col='padj', gene_col='gene', n=300) {
        if (is.null(deg)) return(NULL)
        if (!gene_col %in% names(deg)) gene_col <- names(deg)[1]
        sig <- deg %>% filter(!is.na(.data[[padj_col]]), .data[[padj_col]] < 0.05,
                              abs(.data[[lfc_col]]) > 0.5) %>%
          arrange(.data[[padj_col]]) %>% head(n)
        if (nrow(sig) == 0) return(NULL)
        up   <- sig %>% filter(.data[[lfc_col]] > 0) %>% pull(.data[[gene_col]])
        down <- sig %>% filter(.data[[lfc_col]] < 0) %>% pull(.data[[gene_col]])
        out <- list()
        if (length(up)   > 0) { e <- tryCatch(enrichR::enrichr(up,   dbs), error = function(e) NULL)
          e <- Filter(function(x) is.data.frame(x) && nrow(x) > 0 && is.character(x$Term), e)
          if (length(e) > 0) out$UP   <- bind_rows(e, .id = 'db') }
        if (length(down) > 0) { e <- tryCatch(enrichR::enrichr(down, dbs), error = function(e) NULL)
          e <- Filter(function(x) is.data.frame(x) && nrow(x) > 0 && is.character(x$Term), e)
          if (length(e) > 0) out$DOWN <- bind_rows(e, .id = 'db') }
        out
      }
      r1 <- run_enrich(.read_csv_safe(file.path(.bulkdir, 'NF_vs_pRV_deseq.csv')))
      r2 <- run_enrich(.read_csv_safe(file.path(.bulkdir, 'pRV_vs_RVF_deseq.csv')))
      r3 <- run_enrich(.read_csv_safe(file.path(.bulkdir, 'NF_vs_RVF_deseq.csv')))
      sheets <- list(NF_vs_pRV_UP   = r1$UP, NF_vs_pRV_DOWN   = r1$DOWN,
                     pRV_vs_RVF_UP  = r2$UP, pRV_vs_RVF_DOWN  = r2$DOWN,
                     NF_vs_RVF_UP   = r3$UP, NF_vs_RVF_DOWN   = r3$DOWN)
      saveRDS(sheets, cache_path)
    }
  }
  if (!is.null(sheets)) {
    ## Drop enrichR legacy "Old P-value" columns across all sheets.
    .drop_old_p <- c('Old.P.value', 'Old.Adjusted.P.value',
                      'Old P-value', 'Old Adjusted P-value')
    sheets <- lapply(sheets, function(d) {
      if (is.null(d)) return(d)
      d[, setdiff(names(d), .drop_old_p), drop = FALSE]
    })
    .write_xlsx(sheets, 'Data_S02_bulk_enrichment.xlsx')
  }
})

# ---------- S04 : Bulk WGCNA modules --------------------------------------- #
try({
  .section(04, 'Bulk WGCNA modules')
  if (REBUILD_BULK_WGCNA) .build_bulk_wgcna(force = TRUE)
  mods <- .read_csv_safe(file.path(.shared, 'bulk_heart_modules.csv'))
  if (is.null(mods)) mods <- .read_csv_safe(file.path(.bulkdir, 'bulk_heart_modules.csv'))
  if (!is.null(mods)) {
    color_col <- intersect(c('color', 'Color', 'module', 'Module'), names(mods))
    if (length(color_col) > 0 && .require('WGCNA')) {
      cc          <- color_col[1]
      colors_vec  <- as.character(mods[[cc]])
      std_colors  <- WGCNA::standardColors()
      mods$module_num <- ifelse(tolower(colors_vec) %in% c('grey', 'gray'),
                                0L, match(colors_vec, std_colors))
      mods$module     <- paste0('M', mods$module_num)
      sheets <- list(all_modules = mods)
      for (n in sort(unique(mods$module_num))) {
        key <- paste0('M', n)
        sheets[[key]] <- mods[mods$module_num == n, , drop = FALSE]
      }
    } else {
      sheets <- list(all_modules = mods)
    }
    .write_xlsx(sheets, 'Data_S04_bulk_WGCNA_modules.xlsx')
  }
})

# ---------- S03 : Bulk stepwise trends ------------------------------------- #
try({
  .section(03, 'Bulk stepwise trend genes')
  if (REBUILD_BULK_TRENDS) .build_bulk_trends(force = TRUE)
  trend_files <- c(down_down = 'NF_2_pRV_down_2_RVF_down.csv',
                   down_up   = 'NF_2_pRV_down_2_RVF_up.csv',
                   flat_down = 'NF_2_pRV_flat_2_RVF_down.csv',
                   flat_up   = 'NF_2_pRV_flat_2_RVF_up.csv',
                   up_down   = 'NF_2_pRV_up_2_RVF_down.csv',
                   up_up     = 'NF_2_pRV_up_2_RVF_up.csv')
  sheets <- lapply(trend_files, function(f) {
    df <- .read_csv_safe(file.path(.bulkdir, f))
    if (is.null(df)) return(NULL)
    gene_col <- if ('x' %in% names(df)) 'x' else names(df)[ncol(df)]
    data.frame(gene = df[[gene_col]])
  })
  names(sheets) <- names(trend_files)
  .write_xlsx(sheets, 'Data_S03_bulk_stepwise_trends.xlsx')
})

# ---------- S05 : PAH comparison (concordant / discordant + enrichment) ---- #
# Manuscript claim: concordant DEGs (shared mitochondrial / IFN suppression)
# and discordant DEGs (CD163, mitochondrial translation) between our RV cohort
# and the GSE240921 PAH cohort, with pathway enrichment of each set.
try({
  .section(05, 'PAH comparison: concordance + discordance')
  .build_pah(force = REBUILD_PAH)
  cache_path <- file.path(.cachedir, 'S05_PAH_compare_cache.rds')
  if (!REBUILD_PAH && file.exists(cache_path)) {
    sheets <- readRDS(cache_path)
  } else {
    pairings <- list(
      NF_vs_pRV  = list(ours = 'NF_vs_pRV_deseq.csv',  pah = 'pah.nf.vs.prv.csv'),
      NF_vs_RVF  = list(ours = 'NF_vs_RVF_deseq.csv',  pah = 'pah.nf.vs.rvf.csv'),
      pRV_vs_RVF = list(ours = 'pRV_vs_RVF_deseq.csv', pah = 'pah.prv.vs.rvf.csv')
    )
    .pick <- function(df) {
      gc  <- intersect(c('gene','Gene','symbol'), names(df))[1]
      if (is.na(gc)) gc <- names(df)[1]
      df <- df[, c(gc, 'log2FoldChange','padj')]
      names(df)[1] <- 'gene'
      df[!is.na(df$padj), , drop = FALSE]
    }
    # Criterion follows additional_scripts/AnalysisPAH.R:413-423:
    #   padj < 0.1  AND  |log2FoldChange| > 0.1  in BOTH cohorts.
    .merge_pair <- function(ours_f, pah_f, padj_thresh = 0.1, lfc_thresh = 0.1) {
      a <- .read_csv_safe(file.path(.bulkdir, ours_f))
      b <- .read_csv_safe(file.path(.shared,  pah_f))
      if (is.null(a) || is.null(b)) return(NULL)
      m <- merge(.pick(a), .pick(b), by = 'gene', suffixes = c('_RV', '_PAH'))
      m$sig_RV  <- !is.na(m$padj_RV)  & m$padj_RV  < padj_thresh &
                   abs(m$log2FoldChange_RV)  > lfc_thresh
      m$sig_PAH <- !is.na(m$padj_PAH) & m$padj_PAH < padj_thresh &
                   abs(m$log2FoldChange_PAH) > lfc_thresh
      m$concordant <- m$sig_RV & m$sig_PAH &
                      sign(m$log2FoldChange_RV) == sign(m$log2FoldChange_PAH)
      m$discordant <- m$sig_RV & m$sig_PAH &
                      sign(m$log2FoldChange_RV) != sign(m$log2FoldChange_PAH)
      m
    }
    .enrich_set <- function(genes,
                            dbs = c('GO_Biological_Process_2023','Reactome_2022',
                                    'KEGG_2021_Human','MSigDB_Hallmark_2020')) {
      genes <- unique(genes[!is.na(genes) & nzchar(genes)])
      if (length(genes) < 5 || !.require('enrichR')) return(NULL)
      e <- tryCatch(enrichR::enrichr(genes, dbs), error = function(e) NULL)
      if (is.null(e)) return(NULL)
      e <- Filter(function(x) is.data.frame(x) && nrow(x) > 0 && is.character(x$Term), e)
      if (length(e) > 0) bind_rows(e, .id = 'db') else NULL
    }
    sheets <- list(
      GSE240921_raw  = .read_csv_safe(file.path(.shared, 'GSE240921_processed-data-human.csv')),
      NF_vs_pRV_PAH  = .read_csv_safe(file.path(.shared, 'pah.nf.vs.prv.csv')),
      NF_vs_RVF_PAH  = .read_csv_safe(file.path(.shared, 'pah.nf.vs.rvf.csv')),
      pRV_vs_RVF_PAH = .read_csv_safe(file.path(.shared, 'pah.prv.vs.rvf.csv'))
    )
    for (nm in names(pairings)) {
      m <- .merge_pair(pairings[[nm]]$ours, pairings[[nm]]$pah)
      if (is.null(m)) next
      conc_up <- m$gene[m$concordant & m$log2FoldChange_RV > 0]
      conc_dn <- m$gene[m$concordant & m$log2FoldChange_RV < 0]
      disc_up <- m$gene[m$discordant & m$log2FoldChange_RV > 0]  # up in RV, down in PAH
      disc_dn <- m$gene[m$discordant & m$log2FoldChange_RV < 0]  # down in RV, up in PAH
      sheets[[paste0(nm, '_concordant')]]      <- m[m$concordant, , drop = FALSE]
      sheets[[paste0(nm, '_discordant')]]      <- m[m$discordant, , drop = FALSE]
      sheets[[paste0(nm, '_enrich_conc_UP')]]  <- .enrich_set(conc_up)
      sheets[[paste0(nm, '_enrich_conc_DN')]]  <- .enrich_set(conc_dn)
      sheets[[paste0(nm, '_enrich_disc_RVup')]]<- .enrich_set(disc_up)
      sheets[[paste0(nm, '_enrich_disc_RVdn')]]<- .enrich_set(disc_dn)
    }
    saveRDS(sheets, cache_path)
  }
  ## Drop enrichR legacy 'Old.P.value' / 'Old.Adjusted.P.value' columns.
  sheets <- lapply(sheets, function(d) {
    if (is.null(d)) return(d)
    d[, !names(d) %in% c('Old.P.value', 'Old.Adjusted.P.value'), drop = FALSE]
  })
  ## --- S1I program membership: genes in the four coloured programs of Fig S1I
  ##     (PAH vs RV gene-set scatter). log2FC/padj from the standard DEG tables
  ##     (NF vs RVF; positive = higher in NF), consistent with the other S5 sheets.
  sheets[['S1I_program_genes']] <- local({
    .rv  <- tryCatch(read.csv(file.path(.bulkdir, 'NF_vs_RVF_deseq.csv'), row.names = 1), error = function(e) NULL)
    .pah <- tryCatch(read.csv(file.path(.shared,  'pah.nf.vs.rvf.csv'),   row.names = 1), error = function(e) NULL)
    .univ <- intersect(rownames(.rv), rownames(.pah))   # genes shared by both cohorts = those shown in S1I
    .progs <- list(
      'ETC / OXPHOS'   = list(re = '^(NDUF|SDH[ABCD]|SDHAF|UQCR|COX(4|5|6|7|8|10|11|15|17)|ATP5|CYC1|CYCS|ATPAF|ETFA|ETFB|ETFDH)',
                              x  = c('MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6','MT-CO1','MT-CO2','MT-CO3','MT-CYB','MT-ATP6','MT-ATP8')),
      'Interferon'     = list(re = NULL, x = c('ISG15','IFI6','IFI27','IFI35','IFI44','IFI44L','IFIH1','IFIT1','IFIT2','IFIT3','IFIT5','IFITM1','IFITM2','IFITM3','IRF1','IRF7','IRF9','MX1','MX2','OAS1','OAS2','OAS3','OASL','RSAD2','STAT1','STAT2','USP18','DDX58','DDX60','DDX60L','HERC5','HERC6','XAF1','BST2','GBP1','GBP2','GBP4','GBP5','SAMD9','SAMD9L','EIF2AK2','PARP9','PARP12','PARP14','DTX3L','CMPK2','EPSTI1','LY6E','PLSCR1','TRIM22','SP100','SP110','UBE2L6','LGALS3BP','IFI16','ISG20','HELZ2')),
      'Cardiac stress' = list(re = NULL, x = c('NPPA','NPPB','MYH7','ACTA1','ACTA2','ANKRD1','XIRP2','MYL7','CMYA5')),
      'Complement'     = list(re = NULL, x = c('C1QA','C1QB','C1QC','C1R','C1S','C2','C3','C4A','C4B','C5','C6','C7','C8A','C8B','C8G','C9','CFB','CFH','CFD','CFI','CFP','CR1','SERPING1','MASP1','MASP2','MBL2','C3AR1','C5AR1'))
    )
    .val <- function(tab, g, col) if (!is.null(tab) && g %in% rownames(tab)) tab[g, col] else NA_real_
    do.call(rbind, lapply(names(.progs), function(nm) {
      s   <- .progs[[nm]]
      mem <- s$x
      if (!is.null(s$re)) mem <- c(mem, grep(s$re, .univ, value = TRUE))
      mem <- sort(unique(mem[mem %in% .univ]))
      if (length(mem) == 0) return(NULL)
      data.frame(program = nm, gene = mem,
                 log2FC_NF_vs_RVF_RV  = vapply(mem, .val, numeric(1), tab = .rv,  col = 'log2FoldChange'),
                 padj_RV              = vapply(mem, .val, numeric(1), tab = .rv,  col = 'padj'),
                 log2FC_NF_vs_RVF_PAH = vapply(mem, .val, numeric(1), tab = .pah, col = 'log2FoldChange'),
                 padj_PAH             = vapply(mem, .val, numeric(1), tab = .pah, col = 'padj'),
                 row.names = NULL, stringsAsFactors = FALSE)
    }))
  })
  .write_xlsx(sheets, 'Data_S05_PAH_comparison.xlsx')
})

# ---------- S06 : snRNA annotations --------------------------------------- #
try({
  .section(06, 'snRNA cell type annotations')
  rds <- file.path(.shared, 'snRV_ref.rds')
  if (!file.exists(rds)) rds <- file.path(.shared, 'Post_R3_FINAL_with_counts.rds')
  cache_path <- file.path(.cachedir, 'S06_sn_annotations.csv')
  if (!REBUILD_SN_ANNOT && file.exists(cache_path)) {
    annot <- .rename_subnames(.read_csv_safe(cache_path))
    ## Drop the legacy duplicate column if the cache was written by an
    ## older script version that added Subnames_renamed alongside Subnames_manual.
    if ('Subnames_renamed' %in% names(annot))
      annot$Subnames_renamed <- NULL
    .write_xlsx(list(annotations = annot), 'Data_S06_snRNA_annotations.xlsx')
  } else if (file.exists(rds) && .require('Seurat')) {
    obj  <- readRDS(rds)
    meta <- obj@meta.data %>% tibble::rownames_to_column('cell_id')
    keep <- intersect(c('cell_id', 'patient', 'group', 'orig.ident', 'seurat_clusters',
                        'Names', 'Subnames_manual', 'lineage', 'subtype', 'subname',
                        'cm_subtype', 'fb_subtype', 'ec_subtype', 'myeloid_subtype',
                        'mural_subtype', 'neuron_subtype', 'epi_subtype', 'adipo_subtype',
                        # QC metrics (v57 S6 promise)
                        'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.exon',
                        'entropy', 'scrublet_score'),
                      names(meta))
    annot <- meta[, keep, drop = FALSE]
    # UMAP coordinates (v57 S6 promise)
    umap_red <- tryCatch(Seurat::Embeddings(obj, 'umap'), error = function(e) NULL)
    if (!is.null(umap_red) && nrow(umap_red) == nrow(annot)) {
      annot$UMAP_1 <- umap_red[, 1]
      annot$UMAP_2 <- umap_red[, 2]
    }
    annot <- .rename_subnames(annot)
    write.csv(annot, cache_path, row.names = FALSE)
    .write_xlsx(list(annotations = annot), 'Data_S06_snRNA_annotations.xlsx')
    rm(obj); gc()
  } else {
    .stub(06, 'snRNA annotations', paste('Seurat RDS not found at', rds))
  }
})

# ---------- S07 : snRNA subclustering markers ----------------------------- #
try({
  .section(07, 'snRNA subclustering markers')
  .build_sn_subclust_marks(force = REBUILD_SN_SUBCLUST_MARKS)
  sheets <- list(
    cm      = .read_csv_safe(file.path(.shared, 'cm_subclust_marks.csv')),
    fb      = .read_csv_safe(file.path(.shared, 'fb_subclust_marks.csv')),
    ec      = .read_csv_safe(file.path(.shared, 'ec_subclust_marks.csv')),
    myeloid = .read_csv_safe(file.path(.shared, 'myeloid_subclust_marks.csv')),
    mural   = .read_csv_safe(file.path(.shared, 'mural_subclust_marks.csv'))
  )
  ## Rename cluster ident -> manuscript subname (cluster-level: FB static
  ## map + CM old-Cm# fallback).
  sheets <- lapply(sheets, function(d) {
    if (!is.null(d) && 'cluster' %in% names(d))
      .apply_rmap_inplace(as.data.frame(d), 'cluster') else d })
  ## EC numeric cluster -> subname via majority join (cm_subclust_new_new
  ## style) cached at dependencies/shared/ec_cluster_to_subname.csv.
  .ec_map_csv <- file.path(.shared, 'ec_cluster_to_subname.csv')
  if (!is.null(sheets$ec) && 'cluster' %in% names(sheets$ec) &&
      file.exists(.ec_map_csv)) {
    em <- utils::read.csv(.ec_map_csv, stringsAsFactors = FALSE)
    nm <- as.character(em$subname); names(nm) <- as.character(em$cluster)
    if (anyDuplicated(nm)) {
      d <- duplicated(nm) | duplicated(nm, fromLast = TRUE)
      seq_within <- ave(seq_along(nm[d]), nm[d], FUN = seq_along)
      nm[d] <- paste0(nm[d], '_', seq_within)
    }
    raw <- as.character(sheets$ec$cluster)
    sheets$ec$cluster <- ifelse(raw %in% names(nm), nm[raw], raw)
  }
  ## Mural numeric cluster -> PC/SM lineage label via cached marker-derived
  ## map (dependencies/shared/mural_cluster_to_subname.csv): cluster 0 -> SM_1,
  ## 1 -> PC_1, 2 -> PC_2, 3 -> SM_2 (PC: ABCC9/RGS5/THBS4 markers; SM:
  ## MYH11/TPM2/KCNMA1 markers). Suffixes preserve cluster identity.
  .mu_map_csv <- file.path(.shared, 'mural_cluster_to_subname.csv')
  if (!is.null(sheets$mural) && 'cluster' %in% names(sheets$mural) &&
      file.exists(.mu_map_csv)) {
    mm <- utils::read.csv(.mu_map_csv, stringsAsFactors = FALSE)
    nm <- as.character(mm$subname); names(nm) <- as.character(mm$cluster)
    raw <- as.character(sheets$mural$cluster)
    sheets$mural$cluster <- ifelse(raw %in% names(nm), nm[raw], raw)
  }
  .write_xlsx(sheets, 'Data_S07_snRNA_subclustering_markers.xlsx')
})

# ---------- S08 : snRNA lineage pseudobulk ------------------------------- #
try({
  .section(08, 'snRNA lineage pseudobulk DEGs')
  if (REBUILD_SN_PSEUDOBULK_L1) .build_sn_pseudobulk('lineage', force = TRUE)
  df <- .read_csv_safe(file.path(.shared, 'sn_pseudobulk_lineage_deseq2.csv'))
  if (!is.null(df)) {
    ct_col <- intersect(c('cell_type', 'lineage'), names(df))[1]
    cp_col <- intersect(c('comparison', 'contrast'), names(df))[1]
    sheets <- list(all = df)
    if (!is.na(ct_col) && !is.na(cp_col)) {
      for (ct in sort(unique(df[[ct_col]]))) for (cp in sort(unique(df[[cp_col]]))) {
        sub <- df[df[[ct_col]] == ct & df[[cp_col]] == cp, , drop = FALSE]
        if (nrow(sub) == 0) next
        nm <- substr(gsub('[^A-Za-z0-9_]+', '_', paste0(ct, '__', cp)), 1, 31)
        sheets[[nm]] <- sub
      }
    }
    .write_xlsx(sheets, 'Data_S08_snRNA_lineage_pseudobulk.xlsx')
  }
})

# ---------- S11 : Xenium lineage pseudobulk + sn concordance ------------ #
# v57 promise: "differential expression results within each spatially resolved
# cell type, including cross-platform concordance with snRNA-seq".  We attach
# snRNA log2FC and padj for each (gene, cell_type, comparison) tuple where
# matched at the lineage level.
try({
  .section(11, 'Xenium lineage pseudobulk DEGs')
  if (REBUILD_XEN_PSEUDOBULK_L1) .build_xenium_pseudobulk('lineage', force = TRUE)
  df <- .read_csv_safe(file.path(.xcache, 'xenium_pseudobulk_lineage_deseq2.csv'))
  if (!is.null(df)) {
    ct_col <- intersect(c('cell_type', 'lineage', 'subtype'), names(df))[1]
    cp_col <- intersect(c('comparison', 'contrast'), names(df))[1]
    # Cross-platform concordance: join with snRNA lineage pseudobulk and
    # compute per-(cell-type x contrast) sign-concordance + Spearman of LFC.
    concordance_summary <- NULL
    sn <- .read_csv_safe(file.path(.shared, 'sn_pseudobulk_lineage_deseq2.csv'))
    if (!is.null(sn) && !is.na(ct_col) && !is.na(cp_col)) {
      sn_ct <- intersect(c('cell_type', 'lineage', 'subtype'), names(sn))[1]
      sn_cp <- intersect(c('comparison', 'contrast'), names(sn))[1]
      if (!is.na(sn_ct) && !is.na(sn_cp) && 'gene' %in% names(sn)) {
        sn_pick <- sn[, c('gene', sn_ct, sn_cp, 'log2FoldChange', 'padj')]
        names(sn_pick) <- c('gene', '_ct', '_cp',
                            'snRNA_log2FoldChange', 'snRNA_padj')
        df$`_ct` <- df[[ct_col]]; df$`_cp` <- df[[cp_col]]
        df <- merge(df, sn_pick, by = c('gene', '_ct', '_cp'), all.x = TRUE)
        df$`_ct` <- NULL; df$`_cp` <- NULL
        df$concordant_sign <- !is.na(df$snRNA_log2FoldChange) &
                              !is.na(df$log2FoldChange) &
                              sign(df$snRNA_log2FoldChange) == sign(df$log2FoldChange)
        ## Per (cell-type x contrast) summary
        concordance_summary <- do.call(rbind, by(
          df, list(df[[ct_col]], df[[cp_col]]),
          function(d) {
            v <- !is.na(d$snRNA_log2FoldChange) & !is.na(d$log2FoldChange)
            if (sum(v) < 5) return(NULL)
            sig <- v & !is.na(d$padj) & !is.na(d$snRNA_padj) &
                   d$padj < 0.05 & d$snRNA_padj < 0.05
            sp <- tryCatch(stats::cor(d$log2FoldChange[v],
                                      d$snRNA_log2FoldChange[v],
                                      method = 'spearman'),
                           error = function(e) NA_real_)
            ## Sig-only Spearman + sign-concordance (FDR<0.05 in BOTH platforms).
            ## With ~470 panel genes the all-gene Spearman is dragged toward 0
            ## by non-DE noise; the both-sig subset is the load-bearing metric.
            sp_sig <- if (sum(sig) >= 5)
              tryCatch(stats::cor(d$log2FoldChange[sig],
                                  d$snRNA_log2FoldChange[sig],
                                  method = 'spearman'),
                       error = function(e) NA_real_) else NA_real_
            pct_sig_conc <- if (sum(sig) > 0)
              round(100 * mean(d$concordant_sign[sig], na.rm = TRUE), 1) else NA_real_
            data.frame(
              cell_type   = d[[ct_col]][1],
              contrast    = d[[cp_col]][1],
              n_overlap   = sum(v),
              n_both_sig  = sum(sig),
              spearman_lfc       = round(sp, 3),
              spearman_lfc_sig   = round(sp_sig, 3),
              pct_sign_concordant     = round(100 * mean(
                d$concordant_sign[v], na.rm = TRUE), 1),
              pct_sign_concordant_sig = pct_sig_conc,
              stringsAsFactors = FALSE)
          }))
        if (!is.null(concordance_summary))
          concordance_summary <- concordance_summary[order(
            concordance_summary$cell_type, concordance_summary$contrast), ]
      }
    }
    sheets <- list(all = df)
    if (!is.null(concordance_summary))
      sheets$cross_platform_concordance <- concordance_summary
    if (!is.na(ct_col) && !is.na(cp_col)) {
      for (ct in sort(unique(df[[ct_col]]))) for (cp in sort(unique(df[[cp_col]]))) {
        sub <- df[df[[ct_col]] == ct & df[[cp_col]] == cp, , drop = FALSE]
        if (nrow(sub) == 0) next
        nm <- substr(gsub('[^A-Za-z0-9_]+', '_', paste0(ct, '__', cp)), 1, 31)
        sheets[[nm]] <- sub
      }
    }
    .write_xlsx(sheets, 'Data_S11_xenium_lineage_pseudobulk.xlsx')
  }
})

# ---------- S12 : Xenium subclustering (assignments + markers) --------- #
# v57 promise: "subtype assignments and marker genes for spatially resolved
# cell populations within each of the 12 major cell types after Proseg
# resegmentation, RCTD annotation, and SPLIT correction."
try({
  .section(12, 'Xenium subclustering assignments + markers')
  .build_xenium_subclust_marks(force = REBUILD_XEN_SUBCLUST_MARKS)
  sheets <- list()
  # 1. Subtype assignments (per-cell barcode -> subtype) from Xenium_metadata.csv.
  # CSV was written with row.names but its header lacks a leading slot, so
  # readr shifts every column by 1. Base read.csv auto-detects this; restore
  # the barcode as an explicit `cell_id` column.
  meta_csv <- file.path(.shared, 'Xenium_metadata.csv')
  if (file.exists(meta_csv)) {
    md <- tryCatch(
      utils::read.csv(meta_csv, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL)
    if (!is.null(md)) {
      if (!is.null(rownames(md)) && !all(rownames(md) == as.character(seq_len(nrow(md)))))
        md <- tibble::rownames_to_column(md, 'cell_id')
      ## Drop redundant/internal Xenium QC + cluster cols (user request).
      drop_cols <- c('cluster_kmeans_10', 'orig.ident', 'cluster_split_10',
                     'cell_type_seurat', 'roi_Sample', 'Xenium_snn_res.1',
                     'seurat_clusters', 'subnames', 'names',
                     'cell_types_cm_collapsed', 'nCount_niche_broad',
                     'nFeature_niche_broad', 'kmeans_12', 'kmeans_13',
                     'kmeans_14', 'kmeans_15', 'kmeans_16')
      md <- md[, !(names(md) %in% drop_cols), drop = FALSE]
      sheets[['subtype_assignments']] <- md
    }
  }
  # 2. FindAllMarkers output (gene-level marker tables per subtype)
  marks_csv <- file.path(.xcache, 'xenium_corrected_findallmarkers.csv')
  if (file.exists(marks_csv) && file.size(marks_csv) > 100) {
    d <- .read_csv_safe(marks_csv)
    if (!is.null(d) && 'gene' %in% names(d) && names(d)[1] != 'gene')
      d <- d[, c('gene', setdiff(names(d), 'gene')), drop = FALSE]
    sheets[['markers']] <- d
  }
  # (Reference snRNA CM/EC marker sheets removed -- they live in Data S7
  # (snRNA_subclustering_markers); PAB markers belong in S25 (PAB_murine).
  # S14 is Xenium-only by intent.)
  ## Manuscript-naming rename on cluster ident / cell_type columns.
  sheets <- lapply(sheets, function(d) {
    if (is.null(d)) return(d)
    d <- as.data.frame(d)
    for (cc in intersect(c('cell_type_rctd_doublet','cell_type','cluster','subtype'),
                          names(d)))
      d <- .apply_rmap_inplace(d, cc)
    d
  })
  if (length(sheets) == 0) {
    .stub(12, 'Xenium subclustering',
          'no marker / metadata CSVs found; set REBUILD_XEN_SUBCLUST_MARKS=TRUE')
  } else {
    .write_xlsx(sheets, 'Data_S12_xenium_subclustering.xlsx')
  }
})

# ---------- S13 : Xenium niche composition ----------------------------- #
# v57 promise: cell type proportions within each spatial niche (k=100 NN,
# CLR-transformed, collapsed CM subtypes), niche cluster assignments, and
# niche abundance across NF / pRV / RVF disease states.
try({
  .section(13, 'Xenium niche composition')
  ## Use base read.csv (auto-detects the leading row.names slot the CSV omits;
  ## readr loses the last column on this file).
  md <- tryCatch(
    utils::read.csv(file.path(.shared, 'Xenium_metadata.csv'),
                    stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL)
  if (is.null(md)) {
    .stub(13, 'Xenium niche composition', 'Xenium_metadata.csv not found')
  } else {
    ## Xenium_metadata.csv columns (verified May 20 v58):
    ##   `patient`    -> patient ID (1343/1467/...)
    ##   `group`      -> disease group (NF/pRV/RVF) <-- use this for stratification
    ##   `roi_Sample` -> tissue-section spatial position (Top/Middle/Bottom)
    ##   `cell_type_rctd_doublet` -> coarse lineage
    ##   `cell_types_subclustering` / `cell_types_manual` -> fine subtype
    ## The `niche_manual` column is serialized as `<kmeans_id>,"<label>"`,
    ## so the same label appears N times (one per kmeans cluster). Extract
    ## just the trailing label.
    if ('niche_manual' %in% names(md))
      md$niche_manual <- sub('^[^,]*,"?([^"]+?)"?$', '\\1', as.character(md$niche_manual))
    ## Filter out 'Unassigned' cells so proportions are computed on a
    ## meaningful denominator (per user request).
    keep_rows <- !is.na(md$niche_manual) & !is.na(md$cell_types_manual) &
                 md$cell_types_manual != 'Unassigned' &
                 md$niche_manual      != 'Unassigned'
    md <- md[keep_rows, , drop = FALSE]
    sheets <- list()
    if (all(c('niche_manual', 'cell_types_manual') %in% names(md))) {
      tab <- table(md$niche_manual, md$cell_types_manual)
      prop_per_niche <- as.data.frame.matrix(prop.table(tab, margin = 1))
      prop_per_niche$niche <- rownames(prop_per_niche)
      prop_per_niche <- prop_per_niche[, c('niche',
                          setdiff(names(prop_per_niche), 'niche'))]
      sheets$cell_type_proportions_per_niche <- prop_per_niche
    }
    if (all(c('niche_manual', 'group') %in% names(md))) {
      tab2 <- table(md$niche_manual, md$group)
      ab <- as.data.frame.matrix(prop.table(tab2, margin = 2))   # within each disease group
      # Order columns NF / pRV / RVF if present
      group_order <- intersect(c('NF', 'pRV', 'RVF'), names(ab))
      if (length(group_order)) ab <- ab[, group_order, drop = FALSE]
      ab$niche <- rownames(ab)
      ab <- ab[, c('niche', setdiff(names(ab), 'niche'))]
      sheets$niche_abundance_by_group <- ab
    }
    .write_xlsx(sheets, 'Data_S13_xenium_niche.xlsx')
  }
})

# ---------- S15 : respirometry ---------------------------------------- #
try({
  .section(15, 'Mitochondrial respirometry')
  mito_path <- file.path(.shared, 'mito_data.xlsx')
  if (file.exists(mito_path) && .require('openxlsx')) {
    sheet_names <- openxlsx::getSheetNames(mito_path)
    .clean_mito_cols <- function(nm) {
      nm <- gsub('\\.(?=\\()', ' ', nm, perl = TRUE) # ".(" -> " ("
      nm <- gsub('^X([0-9]+)$', 'sample_\\1', nm)    # X1/X2 placeholders
      nm
    }
    sheets <- setNames(
      lapply(sheet_names, function(s) {
        d <- openxlsx::read.xlsx(mito_path, sheet = s)
        if (!is.null(d)) names(d) <- .clean_mito_cols(names(d))
        d
      }),
      substr(sheet_names, 1, 31)
    )
    .write_xlsx(sheets, 'Data_S15_mito_respirometry.xlsx')
  } else {
    .stub(15, 'Mitochondrial respirometry',
          'mito_data.xlsx not found at dependencies/shared/')
  }
})

# ---------- S17 : Myeloid transcriptional program gene lists ----------- #
# Sourced from Figure_6.R's curated 4-program tibble (GR_targets, MHCII,
# Inflammasome, TypeI_IFN) written to output/Figure_6/.  If Figure_6
# hasn't run yet, fall back to a stub.
try({
  .section(17, 'Myeloid transcriptional program gene lists')
  myeloid_csv <- Find(file.exists, c(
    '/Users/ikuz/Documents/RV_Atlas/output/Figure_6/SuppTable_Myeloid_Programs.csv',
    '/Users/ikuz/Documents/RV_Atlas/output/SuppTable_Myeloid_Programs.csv'))
  if (!is.null(myeloid_csv) && file.exists(myeloid_csv)) {
    df <- .read_csv_safe(myeloid_csv)
    ## Override the per-program `source` strings with literature-free
    ## descriptions (avoids any specific-paper attribution claim that
    ## would need per-gene PDF verification).
    if (!is.null(df) && all(c('program','source') %in% names(df))) {
      .src_override <- c(
        GR_targets = paste0(
          'Canonical glucocorticoid-receptor (GR) negative-feedback ',
          'transactivation reporters (FKBP5, TSC22D3/GILZ, ZBTB16, KLF9, ',
          'ANGPTL4); each gene independently confirmed downregulated ',
          'NF->disease by per-patient pseudobulk in pooled myeloid.'),
        MHCII = paste0(
          'MHC class II / antigen-presentation machinery driven by the ',
          'master transactivator CIITA; confirmed upregulated NF->disease ',
          'in pooled myeloid.'),
        Inflammasome = paste0(
          'NLRP3 inflammasome machinery (sensors NLRP3/NLRP1/AIM2/NAIP, ',
          'adaptor PYCARD, caspases CASP1/CASP4, effector GSDMD, output ',
          'cytokines IL1B/IL18); pertinent negative - unchanged NF->disease.'),
        TypeI_IFN = paste0(
          'Type-I interferon / antiviral interferon-stimulated gene (ISG) ',
          'signature; pertinent negative - unchanged NF->disease in pooled ',
          'myeloid.'))
      for (p in names(.src_override))
        df$source[df$program == p] <- .src_override[[p]]
    }
    sheets <- list(programs = df)
    if (!is.null(df) && 'program' %in% names(df)) {
      for (p in sort(unique(df$program))) {
        sheets[[substr(p, 1, 31)]] <- df[df$program == p, , drop = FALSE]
      }
    }
    .write_xlsx(sheets, 'Data_S17_myeloid_programs.xlsx')
  } else {
    .stub(17, 'Myeloid transcriptional program gene lists',
          'run Figure_6.R first to generate output/Figure_6/SuppTable_Myeloid_Programs.csv')
  }
})

# ---------- S18 : EC hdWGCNA ------------------------------------------- #
try({
  .section(18, 'EC hdWGCNA modules')
  rds <- file.path(.shared, 'EC_hdWGCNA_by_celltype.rds')
  if (REBUILD_EC_HDWGCNA) message('  (rebuild EC hdWGCNA is upstream in Figure_5/7.R; keeping cached rds)')
  if (file.exists(rds)) {
    ec <- readRDS(rds); sheets <- list()
    .extract_hdwgcna <- function(x, prefix = '') {
      if (!.require('hdWGCNA')) return(list())
      out <- list()
      mods <- tryCatch(as.data.frame(hdWGCNA::GetModules(x)),
                       error = function(e) NULL)
      ## Build color -> ecM# map (size descending, excluding grey).
      cmap <- character(0)
      if (!is.null(mods) && 'module' %in% names(mods)) {
        cnt <- sort(table(mods$module[mods$module != 'grey']),
                    decreasing = TRUE)
        cmap <- setNames(paste0('ec', prefix, 'M', seq_along(cnt)), names(cnt))
        cmap <- c(cmap, grey = 'grey')
        mods$module <- cmap[as.character(mods$module)]
        # rename kME_<color> -> kME_<ecM#>
        kidx <- grep('^kME_', names(mods))
        if (length(kidx)) {
          nn <- names(mods)[kidx]
          for (col in names(cmap))
            nn <- sub(paste0('^kME_', col, '$'),
                      paste0('kME_', cmap[[col]]), nn)
          names(mods)[kidx] <- nn
        }
      }
      if (!is.null(mods) && nrow(mods) > 0)
        out[[paste0(prefix, 'modules')]] <- mods
      mes  <- tryCatch(as.data.frame(hdWGCNA::GetMEs(x, harmonized = TRUE)),
                       error = function(e) NULL)
      if (!is.null(mes) && nrow(mes) > 0) {
        if (length(cmap)) {
          # GetMEs columns are color names (or "<color>") -> rename
          nm <- names(mes)
          for (col in names(cmap)) nm <- sub(paste0('^', col, '$'), cmap[[col]], nm)
          names(mes) <- nm
        }
        mes <- cbind(cell_id = rownames(mes), mes); rownames(mes) <- NULL
        out[[paste0(prefix, 'module_eigengenes')]] <- mes
      }
      hubs <- tryCatch(as.data.frame(hdWGCNA::GetHubGenes(x, n_hubs = 25)),
                       error = function(e) NULL)
      if (!is.null(hubs) && nrow(hubs) > 0) {
        if (length(cmap) && 'module' %in% names(hubs))
          hubs$module <- cmap[as.character(hubs$module)]
        out[[paste0(prefix, 'hub_genes')]] <- hubs
      }
      out
    }
    if (is.data.frame(ec)) {
      sheets$modules <- ec
    } else if (inherits(ec, 'Seurat')) {
      sheets <- modifyList(sheets, .extract_hdwgcna(ec))
    } else if (is.list(ec)) {
      for (nm in names(ec)) {
        x <- ec[[nm]]
        if (is.data.frame(x)) { sheets[[substr(nm, 1, 31)]] <- x }
        else if (inherits(x, 'Seurat')) {
          add <- .extract_hdwgcna(x, prefix = paste0(nm, '_'))
          names(add) <- vapply(names(add), function(k) substr(k, 1, 31), character(1))
          sheets <- modifyList(sheets, add)
        }
      }
    }
    if (length(sheets) == 0)
      sheets$README <- data.frame(note = 'EC_hdWGCNA object contained no extractable tables')
    .write_xlsx(sheets, 'Data_S18_EC_hdWGCNA.xlsx')
  } else {
    .stub(18, 'EC hdWGCNA', 'EC_hdWGCNA_by_celltype.rds not found')
  }
})

# ---------- S14 : CellChat spatial niche communication ---------------- #
# Per-niche ligand-receptor communication (CellChat) across NF / pRV / RVF,
# from the spatial CellChat pipeline (replicability source:
# additional_scripts/spatial_cellchat_niche_communication.R). One row per
# ligand-receptor interaction x spatial niche: per-condition mean CellChat
# communication probability + patient counts, unpaired Wilcoxon for NF->pRV
# and pRV->RVF with rank-biserial r effect sizes, confidence tiers, and the
# directional pattern. The README sheet documents every column.
try({
  .section(14, 'CellChat spatial niche communication')
  cc <- tryCatch(
    utils::read.csv(file.path(.xcache, 'xenium_cellchat_niche_communication.csv'),
                    stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL)
  if (is.null(cc)) {
    .stub(14, 'CellChat spatial niche communication',
          'place xenium_cellchat_niche_communication.csv in dependencies/shared/Xenium/')
  } else {
    readme <- data.frame(
      column = c('niche', 'interaction_name', 'pathway_name', 'ligand', 'receptor',
                 'source', 'target', 'mean_prob_NF', 'mean_prob_pRV', 'mean_prob_RVF',
                 'n_patients_NF', 'n_patients_pRV', 'n_patients_RVF',
                 'wilcox_p_NF_vs_pRV', 'effect_size_r_NF_vs_pRV', 'tier_NF_vs_pRV', 'pattern_NF_vs_pRV',
                 'wilcox_p_pRV_vs_RVF', 'effect_size_r_pRV_vs_RVF', 'tier_pRV_vs_RVF', 'pattern_pRV_vs_RVF',
                 'has_absent_niche'),
      description = c(
        'Spatial tissue niche the interaction was analyzed in (e.g. Myocardium, Fibrotic remodeling stroma)',
        'Ligand-receptor pair (e.g. SEMA3C_PLXND1)',
        'Signaling pathway the LR pair belongs to (e.g. SEMA3, COLLAGEN)',
        'Ligand gene',
        'Receptor gene',
        'Sending cell-type group (e.g. Neural, EC, CM)',
        'Receiving cell-type group',
        'Mean CellChat communication probability in NF, averaged across patients',
        'Mean CellChat communication probability in pRV, averaged across patients',
        'Mean CellChat communication probability in RVF, averaged across patients',
        'Number of patients (of 3) in NF with detectable signal for this interaction',
        'Number of patients (of 3) in pRV with detectable signal for this interaction',
        'Number of patients (of 3) in RVF with detectable signal for this interaction',
        'Wilcoxon rank-sum p-value, NF vs pRV',
        'Rank-biserial r effect size, NF vs pRV (-1..1; positive = higher in pRV)',
        'Confidence tier, NF vs pRV: 1_top (|r|=1, p<0.1), 2_large_effect (|r|>=0.67, p<0.2), or 3_exploratory',
        'Directional pattern, NF vs pRV: gained_in_pRV, lost_in_pRV, or quantitative_change',
        'Wilcoxon rank-sum p-value, pRV vs RVF',
        'Rank-biserial r effect size, pRV vs RVF (-1..1; positive = higher in RVF)',
        'Confidence tier, pRV vs RVF (same scheme as tier_NF_vs_pRV)',
        'Directional pattern, pRV vs RVF: gained_in_RVF, lost_in_RVF, or quantitative_change',
        'TRUE if any patient lacked this niche entirely (vs niche present but interaction undetected); aids zero interpretation'),
      stringsAsFactors = FALSE)
    .write_xlsx(list(README = readme, niche_communication = cc),
                'Data_S14_CellChat_niche_communication.xlsx')
  }
})

# ---------- S19 : Murine PAB ------------------------------------------- #
try({
  .section(19, 'Murine PAB dataset')
  rds <- file.path(.shared, 'PAB_data_clean.rds')
  sheets <- list()
  if (file.exists(rds) && .require('Seurat')) {
    pab <- readRDS(rds)
    if (inherits(pab, 'Seurat')) {
      md <- pab@meta.data %>% tibble::rownames_to_column('cell_id')
      ## Drop QC/cluster cols + rename patient -> mouse (user request).
      ## Also: the existing `group` column has NAs; `orig.ident` is the
      ## complete Nor/Mod/Sev label -> drop `group`, rename orig.ident -> group.
      drop_cols <- c('scrublet_call', 'RNA_snn_res.2', 'seurat_clusters',
                     'SCT_snn_res.0.8', 'Names', 'SCT_snn_res.1', 'group')
      md <- md[, !(names(md) %in% drop_cols), drop = FALSE]
      if ('orig.ident' %in% names(md)) names(md)[names(md) == 'orig.ident'] <- 'group'
      if ('patient' %in% names(md)) names(md)[names(md) == 'patient'] <- 'mouse'
      sheets$PAB_annotations <- md
      if (REBUILD_PAB) {
        # rebuild PAB marker CSVs from Seurat object
        Seurat::Idents(pab) <- md[[intersect(c('cell_type','Names','lineage'), names(md))[1]]]
        mk <- Seurat::FindAllMarkers(pab, only.pos = TRUE, min.pct = 0.25,
                                     logfc.threshold = 0.25, verbose = FALSE)
        write.csv(mk, file.path(.bulkdir, 'PAB_marks.csv'), row.names = FALSE)
      }
    }
    rm(pab); gc()
  }
  .clean_pab_marker <- function(d) {
    if (is.null(d)) return(d)
    cn <- names(d)
    drop_idx <- which(cn == '' | cn == 'X' | grepl('^\\.\\.\\.?\\d?$', cn))
    if (length(drop_idx)) d <- d[, -drop_idx, drop = FALSE]
    if ('gene' %in% names(d) && names(d)[1] != 'gene')
      d <- d[, c('gene', setdiff(names(d), 'gene')), drop = FALSE]
    d
  }
  pab_ec   <- .clean_pab_marker(.read_csv_safe(file.path(.bulkdir, 'PAB_EC.csv')))
  pab_mark <- .clean_pab_marker(.read_csv_safe(file.path(.bulkdir, 'PAB_marks.csv')))
  if (!is.null(pab_ec))   sheets$PAB_EC      <- pab_ec
  if (!is.null(pab_mark)) sheets$PAB_markers <- pab_mark
  ## Per-mouse pseudobulk DESeq2 contrasts (Mod_vs_Nor / Sev_vs_Nor / Sev_vs_Mod)
  ## across all 11 lineages in `Names`. Built once per PAB_data_clean.rds.
  if (REBUILD_PAB) .build_pab_pseudobulk(force = TRUE)
  if (!file.exists(file.path(.shared, 'pab_pseudobulk_lineage_deseq2.csv')))
    .build_pab_pseudobulk(force = FALSE)
  pab_deg <- .read_csv_safe(file.path(.shared, 'pab_pseudobulk_lineage_deseq2.csv'))
  if (!is.null(pab_deg)) {
    sheets$PAB_pseudobulk_all <- pab_deg
    if (all(c('cell_type', 'comparison') %in% names(pab_deg))) {
      for (ct in sort(unique(pab_deg$cell_type))) for (cp in sort(unique(pab_deg$comparison))) {
        sub <- pab_deg[pab_deg$cell_type == ct & pab_deg$comparison == cp, , drop = FALSE]
        if (nrow(sub) == 0) next
        nm <- substr(gsub('[^A-Za-z0-9_]+', '_', paste0(ct, '__', cp)), 1, 31)
        sheets[[nm]] <- sub
      }
    }
  }
  .write_xlsx(sheets, 'Data_S19_PAB_murine.xlsx')
})

# ---------- S20 : CX3CR1 ----------------------------------------------- #
# Derived in-script from kallisto quants at dependencies/shared/Mouse_PAB_Myeloid/
# (same raw inputs as Supplementary_Figure_5.R S5I/J).  Three contrasts on the
# combined Origin x Type group factor.
try({
  .section(20, 'CX3CR1-tdTomato lineage tracing bulk DEGs')
  .cx3dir  <- file.path(.shared, 'Mouse_PAB_Myeloid')
  .cx3meta <- file.path(.cx3dir, 'meta.csv')
  .cx3t2g  <- file.path(.cx3dir, 't2g.txt')
  cx3_outs <- list(
    RV_PAB_vs_RV_Sham  = file.path(.shared, 'CX3CR1_RV_PAB_vs_RV_Sham.csv'),
    LV_PAB_vs_LV_Sham  = file.path(.shared, 'CX3CR1_LV_PAB_vs_LV_Sham.csv'),
    RV_PAB_vs_LV_PAB   = file.path(.shared, 'CX3CR1_RV_PAB_vs_LV_PAB.csv'),
    RV_Sham_vs_LV_Sham = file.path(.shared, 'CX3CR1_RV_Sham_vs_LV_Sham.csv')
  )
  rebuild_needed <- REBUILD_CX3CR1 || !all(file.exists(unlist(cx3_outs)))
  if (rebuild_needed) {
    if (!file.exists(.cx3meta) || !file.exists(.cx3t2g)) {
      .stub(20, 'CX3CR1 lineage tracing bulk',
            'Mouse_PAB_Myeloid/{meta.csv,t2g.txt} not found')
    } else if (!.require(c('DESeq2', 'tximport'))) {
      message('  DESeq2/tximport unavailable; skipping')
    } else {
      samples <- read.csv(.cx3meta, header = TRUE, check.names = TRUE)
      samples$group <- factor(paste(samples$Origin, samples$Type, sep = '_'),
                              levels = c('LV_Sham','LV_PAB','RV_Sham','RV_PAB'))
      files <- file.path(.cx3dir, 'nascent', samples$ID, 'abundance.h5')
      names(files) <- samples$ID
      ok <- file.exists(files)
      if (!all(ok)) {
        message('  missing kallisto h5s for: ', paste(samples$ID[!ok], collapse = ', '))
        files   <- files[ok]
        samples <- samples[ok, , drop = FALSE]
      }
      t2g <- read.table(.cx3t2g, fill = TRUE)
      tx2gene <- data.frame(TXNAME = t2g$V1, GENEID = t2g$V3)
      txi <- tximport::tximport(files, type = 'kallisto', tx2gene = tx2gene)
      rownames(samples) <- samples$ID
      dds <- DESeq2::DESeqDataSetFromTximport(txi, samples, design = ~ group)
      dds <- DESeq2::DESeq(dds)
      sh <- function(a, b) as.data.frame(
        DESeq2::lfcShrink(dds, contrast = c('group', a, b), type = 'ashr')
      ) %>% tibble::rownames_to_column('gene')
      write.csv(sh('RV_PAB',  'RV_Sham'),  cx3_outs$RV_PAB_vs_RV_Sham,  row.names = FALSE)
      write.csv(sh('LV_PAB',  'LV_Sham'),  cx3_outs$LV_PAB_vs_LV_Sham,  row.names = FALSE)
      write.csv(sh('RV_PAB',  'LV_PAB' ),  cx3_outs$RV_PAB_vs_LV_PAB,   row.names = FALSE)
      write.csv(sh('RV_Sham', 'LV_Sham'),  cx3_outs$RV_Sham_vs_LV_Sham, row.names = FALSE)
      message('  wrote ', paste(unlist(cx3_outs), collapse = ', '))
    }
  }
  ## Excel auto-corrupts Sept1-Sept12 -> "1".."12" when the kallisto-output
  ## CSVs are opened by Excel/Numbers. Restore the canonical gene symbols.
  .fix_sept_ids <- function(d) {
    if (is.null(d) || !'gene' %in% names(d)) return(d)
    g <- as.character(d$gene)
    bad <- g %in% as.character(1:12)
    if (any(bad)) d$gene[bad] <- paste0('Sept', g[bad])
    d
  }
  sheets <- Filter(Negate(is.null), lapply(cx3_outs, .read_csv_safe))
  sheets <- lapply(sheets, .fix_sept_ids)
  if (length(sheets) > 0)
    .write_xlsx(sheets, 'Data_S20_CX3CR1_tracing_bulk.xlsx')
  else
    .stub(20, 'CX3CR1 lineage tracing bulk', 'no contrasts computed')
})

# ---------- S21 : Pediatric annotations -------------------------------- #
try({
  .section(21, 'Pediatric annotations')
  rds <- file.path(.shared, 'all_peds_data.rds')
  if (file.exists(rds) && .require('Seurat')) {
    peds <- readRDS(rds)
    if (inherits(peds, 'Seurat')) {
      md <- peds@meta.data %>% tibble::rownames_to_column('cell_id')
      ## Drop QC/cluster cols that aren't useful to readers (user request).
      drop_cols <- c('orig.ident', 'bh_pval', 'scrublet_cluster_score',
                     'SCT_snn_res.0.1', 'SCT_snn_res.0.2', 'SCT_snn_res.0.3',
                     'SCT_snn_res.0.4', 'SCT_snn_res.0.5',
                     'seurat_clusters', 'id')
      md <- md[, !(names(md) %in% drop_cols), drop = FALSE]
      .write_xlsx(list(peds_annotations = md), 'Data_S21_peds_annotations.xlsx')
    }
    rm(peds); gc()
  } else {
    .stub(21, 'Pediatric annotations', 'all_peds_data.rds not found')
  }
})

# ---------- S22 : Pediatric DEGs --------------------------------------- #
try({
  .section(22, 'Pediatric DEGs')
  rds        <- file.path(.shared, 'all_peds_data.rds')
  cache_path <- file.path(.cachedir, 'S19_peds_DEGs.rds')
  .gene_first <- function(d) {
    if (is.null(d) || !'gene' %in% names(d)) return(d)
    d[, c('gene', setdiff(names(d), 'gene')), drop = FALSE]
  }
  if (!REBUILD_PEDS && file.exists(cache_path)) {
    sheets <- lapply(readRDS(cache_path), .gene_first)
    .write_xlsx(sheets, 'Data_S22_peds_DEGs.xlsx')
  } else if (!file.exists(rds)) {
    .stub(22, 'Pediatric DEGs', 'all_peds_data.rds not found')
  } else if (!.require('Seurat')) {
    .stub(22, 'Pediatric DEGs', 'Seurat not installed')
  } else {
    peds <- readRDS(rds); md <- peds@meta.data
    group_col <- intersect(c('group', 'condition', 'disease'), names(md))[1]
    ct_col    <- intersect(c('cell_type', 'lineage', 'subtype', 'celltype',
                             'cell.type', 'Names', 'subname'), names(md))[1]
    sheets <- list()
    if (!is.na(group_col) && !is.na(ct_col) && length(unique(md[[group_col]])) >= 2) {
      Seurat::Idents(peds) <- md[[ct_col]]
      for (ct in unique(md[[ct_col]])) {
        sub <- tryCatch(subset(peds, idents = ct), error = function(e) NULL); if (is.null(sub)) next
        Seurat::Idents(sub) <- sub@meta.data[[group_col]]
        if (length(unique(Seurat::Idents(sub))) < 2) next
        mk <- tryCatch(Seurat::FindAllMarkers(sub, only.pos = FALSE, min.pct = 0.1,
                                              logfc.threshold = 0.25, verbose = FALSE),
                       error = function(e) NULL)
        if (!is.null(mk) && nrow(mk) > 0)
          sheets[[substr(gsub('[^A-Za-z0-9_]+', '_', ct), 1, 31)]] <- mk
      }
    }
    if (length(sheets) == 0)
      sheets$README <- data.frame(note = 'no usable group/celltype columns in peds metadata')
    saveRDS(sheets, cache_path)
    sheets <- lapply(sheets, .gene_first)
    .write_xlsx(sheets, 'Data_S22_peds_DEGs.xlsx')
    rm(peds); gc()
  }
})

# ---------- S10 : Xenium gene panel + imputation provenance ----------- #
# Merged file (panel sheet + provenance sheet + summary + historical) per v57
# S20 promise.  The standalone S19 panel xlsx was removed in this version.
try({
  .section(10, 'Xenium gene panel + provenance')
  ## Two sheets: summary (panel + imputed gene counts) + provenance
  ## (all genes used in Xenium analyses, labelled panel_measured vs
  ## SpaGE_imputed). Imputed-gene catalog is read from the imputed
  ## Xenium object's gene-id cache (gene_ids_Xenium_resegmented_imputed
  ## _final_1343.csv); negative-control probes are filtered out.
  panel_csv <- file.path(.shared, 'panel_genes.csv')
  imp_csv   <- file.path(.shared,
                  'gene_ids_Xenium_resegmented_imputed_final_1343.csv')
  if (!file.exists(panel_csv) || !file.exists(imp_csv)) {
    .stub(10, 'Xenium gene provenance',
          'panel_genes.csv or gene_ids_*_imputed_*.csv missing')
  } else {
    p <- suppressMessages(readr::read_csv(panel_csv, show_col_types = FALSE))
    panel_genes <- unique(trimws(as.character(
      p[[intersect(c('x','gene','symbol'), names(p))[1]]])))
    panel_genes <- panel_genes[nzchar(panel_genes)]

    gids   <- as.character(utils::read.csv(imp_csv, stringsAsFactors = FALSE)[[1]])
    all_g  <- unique(trimws(gids))
    real_g <- all_g[!grepl('^(NEGCONTROL|UNASSIGNED|BLANK)', all_g,
                           ignore.case = TRUE)]
    imp_off <- setdiff(real_g, panel_genes)

    provenance <- data.frame(
      gene = c(panel_genes, imp_off),
      on_panel = c(rep(TRUE,  length(panel_genes)),
                   rep(FALSE, length(imp_off))),
      provenance = c(rep('panel_measured', length(panel_genes)),
                     rep('SpaGE_imputed',   length(imp_off))),
      stringsAsFactors = FALSE)
    provenance <- provenance[order(provenance$gene), ]
    rownames(provenance) <- NULL
    summary_tbl <- data.frame(
      metric = c('panel_gene_count', 'imputed_gene_count'),
      value  = c(length(panel_genes), length(imp_off)),
      stringsAsFactors = FALSE)
    .write_xlsx(list(summary = summary_tbl, provenance = provenance),
                'Data_S10_xenium_gene_provenance.xlsx')
  }
})

# [REMOVED in v60 — not cited as a Data S; code retained, commented]
# # ---------- S21 : ChEA TF enrichment ---------------------------------- #
# try({
#   .section(21, 'ChEA TF enrichment')
#   ## Drop enrichR legacy 'Old.P.value' / 'Old.Adjusted.P.value' columns that
#   ## carry no usable information (pre-FDR-correction values from the legacy API).
#   .drop_old_pvals <- function(d) {
#     if (is.null(d)) return(d)
#     d[, !names(d) %in% c('Old.P.value', 'Old.Adjusted.P.value'), drop = FALSE]
#   }
#   cache_path <- file.path(.cachedir, 'S21_ChEA_cache.rds')
#   if (!REBUILD_CHEA && file.exists(cache_path)) {
#     sheets <- lapply(readRDS(cache_path), .drop_old_pvals)
#     .write_xlsx(sheets, 'Data_S21_ChEA_TF_enrichment.xlsx')
#   } else if (!.require('enrichR')) {
#     .stub(21, 'ChEA TF enrichment', 'install enrichR')
#   } else {
#     dbs <- c('ChEA_2022', 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP_X',
#              'TRRUST_Transcription_Factors_2019')
#     deg_files <- list(NF_vs_pRV  = file.path(.bulkdir, 'NF_vs_pRV_deseq.csv'),
#                       pRV_vs_RVF = file.path(.bulkdir, 'pRV_vs_RVF_deseq.csv'),
#                       NF_vs_RVF  = file.path(.bulkdir, 'NF_vs_RVF_deseq.csv'))
#     sheets <- list()
#     for (nm in names(deg_files)) {
#       d <- .read_csv_safe(deg_files[[nm]]); if (is.null(d)) next
#       gcol <- intersect(c('gene', 'Gene', 'symbol'), names(d))[1]
#       if (is.na(gcol)) gcol <- names(d)[1]  # CSV with unnamed rownames column
#       sig <- d %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5)
#       for (dir in c('UP', 'DOWN')) {
#         gs <- if (dir == 'UP') sig %>% filter(log2FoldChange > 0) %>% pull(!!sym(gcol))
#               else             sig %>% filter(log2FoldChange < 0) %>% pull(!!sym(gcol))
#         if (length(gs) == 0) next
#         e <- tryCatch(enrichR::enrichr(gs, dbs), error = function(e) NULL)
#         e <- Filter(function(x) is.data.frame(x) && nrow(x) > 0 && is.character(x$Term), e)
#         if (length(e) > 0)
#           sheets[[substr(paste0(nm, '_', dir), 1, 31)]] <- bind_rows(e, .id = 'db')
#       }
#     }
#     saveRDS(sheets, cache_path)
#     sheets <- lapply(sheets, .drop_old_pvals)
#     .write_xlsx(sheets, 'Data_S21_ChEA_TF_enrichment.xlsx')
#   }
# })

# ---------- S09 : snRNA sublineage pseudobulk ------------------------- #
try({
  .section(09, 'snRNA sublineage pseudobulk DEGs')
  if (REBUILD_SN_PSEUDOBULK_L2) .build_sn_pseudobulk('sublineage', force = TRUE)
  df <- .read_csv_safe(file.path(.shared, 'sn_pseudobulk_sublineage_deseq2.csv'))
  if (!is.null(df)) {
    ct_col <- intersect(c('cell_type', 'subtype', 'subname'), names(df))[1]
    if (!is.na(ct_col)) df <- .apply_rmap_inplace(df, ct_col)
    ## Excel's single-sheet limit is 1,048,576 rows. The 34-subtype x
    ## 3-contrast pseudobulk exceeds this (~1.14 M rows), which causes
    ## Excel's "needs to be repaired" error. Split by contrast (each
    ## contrast ~380 k rows, well under the limit). subtype column is
    ## preserved within each sheet for downstream slicing.
    split_col <- intersect(c('contrast', 'comparison'), names(df))[1]
    sheets <- if (!is.na(split_col) && nrow(df) > 1048575)
      split(df, df[[split_col]]) else list(all = df)
    .write_xlsx(sheets, 'Data_S09_snRNA_sublineage_pseudobulk.xlsx')
  }
})

# ---------- S16 : WGA cardiomyocyte cross-sectional area --------------- #
# Per-cell Minimum Feret diameter (proxy for CSA) measured by wheat-germ
# agglutinin staining of human RV sections (NF n=4, pRV n=4, RVF n=4;
# 17,864 CMs total).  Source: dependencies/shared/wga_minferet_human_rv.csv
try({
  .section(16, 'WGA cardiomyocyte cross-sectional area')
  wga_csv <- file.path(.shared, 'wga_minferet_human_rv.csv')
  if (file.exists(wga_csv)) {
    df <- .read_csv_safe(wga_csv)
    if (!is.null(df) && nrow(df) > 0) {
      idx_col <- which(names(df) %in% c('', '...1'))
      if (length(idx_col) > 0) df <- df[, -idx_col, drop = FALSE]
      sheets <- list(per_cell = df)
      if ('Group' %in% names(df) && 'AreaShape_MinFeretDiameter' %in% names(df)) {
        med_tbl <- aggregate(AreaShape_MinFeretDiameter ~ HHTB_ID + Group, data = df,
                             FUN = median, na.rm = TRUE)
        names(med_tbl)[3] <- 'median_MinFeret'
        sheets$per_patient_median <- med_tbl
      }
      .write_xlsx(sheets, 'Data_S16_WGA_CM_cross_sectional_area.xlsx')
    } else {
      .stub(16, 'WGA cross-sectional area', 'wga_minferet_human_rv.csv empty')
    }
  } else {
    .stub(16, 'WGA cross-sectional area',
          'wga_minferet_human_rv.csv not found in dependencies/shared/')
  }
})

# [REMOVED in v60 — not cited as a Data S; code retained, commented]
# # ---------- S27 : CollecTRI / decoupleR TF activity ------------------- #
# # Per-cell ULM transcription-factor activity (CollecTRI net, decoupleR
# # run_ulm) from additional_scripts/mito_TF.R, cached at
# # dependencies/shared/TF_activity.rds (738 TFs x ~21.5k CM nuclei).
# # Per-cell is ~16M rows (won't fit a sheet) so we emit a per-TF summary
# # and a compact TF x patient mean-activity matrix; join patient->disease
# # group via the cohort table (Data S5 / Supp Table S3).
# try({
#   .section(27, 'CollecTRI/decoupleR TF activity (snRNA CM)')
#   tf_rds <- file.path(.shared, 'TF_activity.rds')
#   if (file.exists(tf_rds)) {
#     acts <- readRDS(tf_rds)
#     acts <- acts[acts$statistic == 'ulm' & !is.na(acts$score), , drop = FALSE]
#     acts$patient <- sub('_.*$', '', as.character(acts$condition))
#     summ <- data.frame(
#       TF           = names(tapply(acts$score, acts$source, length)),
#       n_cells      = as.integer(tapply(acts$score, acts$source, length)),
#       mean_score   = as.numeric(tapply(acts$score, acts$source, mean)),
#       sd_score     = as.numeric(tapply(acts$score, acts$source, stats::sd)),
#       frac_sig_p05 = as.numeric(tapply(acts$p_value < 0.05, acts$source,
#                                        mean, na.rm = TRUE)),
#       stringsAsFactors = FALSE)
#     summ <- summ[order(-abs(summ$mean_score)), ]
#     m <- tapply(acts$score, list(acts$source, acts$patient), mean)
#     mat <- data.frame(TF = rownames(m), round(as.data.frame(m), 4),
#                        check.names = FALSE, stringsAsFactors = FALSE)
#     rownames(mat) <- NULL
#     .write_xlsx(list(per_TF_summary = summ, TF_by_patient_mean = mat),
#                 'Data_S27_TF_activity_CollecTRI.xlsx')
#   } else {
#     .stub(27, 'CollecTRI/decoupleR TF activity',
#           'TF_activity.rds not found (run additional_scripts/mito_TF.R)')
#   }
# })

# [REMOVED in v60 — not cited as a Data S; code retained, commented]
# # ---------- S25 : Differential cell-type abundance -------------------- #
# # Per-patient cell-type proportions + Kruskal-Wallis across disease
# # groups, for the human snRNA atlas (lineage + subtype) and Xenium.
# try({
#   .section(25, 'Differential cell-type abundance (snRNA & Xenium)')
#   out <- list()
#   .abund <- function(md, idc, pc, gcc, tag) {
#     sub <- md[!is.na(md[[idc]]) & md[[idc]] != 'Unassigned', , drop = FALSE]
#     tab <- as.data.frame(table(patient = as.character(sub[[pc]]),
#                                cell_type = as.character(sub[[idc]])),
#                          stringsAsFactors = FALSE)
#     tot <- tapply(tab$Freq, tab$patient, sum)
#     tab$prop <- tab$Freq / tot[tab$patient]
#     g <- unique(data.frame(patient = as.character(sub[[pc]]),
#                            group = as.character(sub[[gcc]]),
#                            stringsAsFactors = FALSE))
#     tab <- merge(tab, g, by = 'patient', all.x = TRUE)
#     tab$cell_type <- .relabel_subnames(tab$cell_type)   # cluster-level rename
#     out[[paste0(tag, '_proportions')]] <<- tab
#     st <- do.call(rbind, by(tab, tab$cell_type, function(d) data.frame(
#       cell_type = d$cell_type[1], n_patients = nrow(d),
#       kruskal_p = tryCatch(stats::kruskal.test(
#         d$prop, factor(d$group))$p.value, error = function(e) NA_real_),
#       stringsAsFactors = FALSE)))
#     st$kruskal_fdr <- stats::p.adjust(st$kruskal_p, 'BH')
#     out[[paste0(tag, '_KW_stats')]] <<- st[order(st$kruskal_p), ]
#   }
#   rds <- file.path(.shared, 'snRV_ref.rds')
#   if (file.exists(rds)) {
#     obj <- readRDS(rds); md <- obj@meta.data; rm(obj); invisible(gc(FALSE))
#     pc  <- intersect(c('patient', 'Patient', 'orig.ident'), names(md))[1]
#     gcc <- intersect(c('group', 'category', 'disease'), names(md))[1]
#     ## snRNA uses `Names` (lineage) and `Subnames_manual` (subtype) -- not the
#     ## Xenium column names. Iterate over whichever snRNA cell-type cols exist.
#     for (lv in list(c('lineage', 'Names'),
#                     c('subtype', 'Subnames_manual'))) {
#       if (lv[2] %in% names(md) && !is.na(pc) && !is.na(gcc))
#         .abund(md, lv[2], pc, gcc, paste0('snRNA_', lv[1]))
#     }
#   } else {
#     .stub(25, 'Differential abundance', 'snRV_ref.rds not found')
#   }
#   xm <- tryCatch(
#     utils::read.csv(file.path(.shared, 'Xenium_metadata.csv'),
#                     stringsAsFactors = FALSE, check.names = FALSE),
#     error = function(e) NULL)
#   if (!is.null(xm)) {
#     ## Xenium_metadata.csv columns (verified May 20 v58):
#     ##   `patient`                  -> patient ID (1343/1467/...)
#     ##   `group`                    -> disease group (NF/pRV/RVF)
#     ##   `roi_Sample`               -> tissue position (Top/Middle/Bottom)
#     ##   `cell_type_rctd_doublet`   -> coarse lineage
#     ##   `cell_types_subclustering` -> fine subtype
#     xm <- as.data.frame(xm)
#     if (all(c('patient','group','cell_type_rctd_doublet') %in% names(xm)))
#       .abund(xm, 'cell_type_rctd_doublet',   'patient', 'group', 'Xenium_lineage')
#     if (all(c('patient','group','cell_types_subclustering') %in% names(xm)))
#       .abund(xm, 'cell_types_subclustering', 'patient', 'group', 'Xenium_subtype')
#   }
#   if (length(out)) .write_xlsx(out, 'Data_S25_celltype_abundance.xlsx')
# })

# [REMOVED in v60 — not cited as a Data S; code retained, commented]
# # ---------- S28 : Cardiac fibrosis quantification --------------------- #
# # Picrosirius-red / Masson-trichrome % fibrosis. Sources (group-blocked
# # wide sheets) tidied to long: murine PAB (PAB_Sirius.xlsx: Sham /
# # PAB.Mod / PAB.Sev) and human RV (Adult_Trichrome_Sirius.xlsx: NF /
# # pRV / RVF, with patient IDs). '*' marks flagged measurements.
# try({
#   .section(28, 'Cardiac fibrosis quantification (Sirius / Trichrome)')
#   .ffill <- function(nm) {                       # group = last non-filler header
#     g <- nm
#     for (i in seq_along(g))
#       if (i > 1 && (g[i] == '' || grepl('^X[0-9]+$|^\\.\\.\\.', g[i])))
#         g[i] <- g[i - 1]
#     g
#   }
#   .num <- function(x) suppressWarnings(as.numeric(
#     gsub('[^0-9.eE+-]', '', trimws(as.character(x)))))
#   sheets <- list()
#   p1 <- file.path(.shared, 'PAB_Sirius.xlsx')
#   if (file.exists(p1)) {
#     d <- openxlsx::read.xlsx(p1, 1); g <- .ffill(names(d))
#     v <- unlist(d[1, ], use.names = FALSE)
#     long <- data.frame(Group = g, fibrosis_pct = .num(v),
#                        stringsAsFactors = FALSE)
#     sheets$PAB_Sirius <- long[!is.na(long$fibrosis_pct), ]
#   }
#   p2 <- file.path(.shared, 'Adult_Trichrome_Sirius.xlsx')
#   if (file.exists(p2)) {
#     d <- openxlsx::read.xlsx(p2, 1); g <- .ffill(names(d))
#     pid <- unlist(d[1, ], use.names = FALSE)
#     val <- unlist(d[2, ], use.names = FALSE)
#     pid_clean <- as.character(pid)
#     pid_clean[is.na(pid_clean) | !nzchar(trimws(pid_clean))] <- 'Out of cohort'
#     long <- data.frame(Group = g, PatientID = pid_clean,
#                        fibrosis_pct = .num(val),
#                        stringsAsFactors = FALSE)
#     sheets$Adult_Trichrome_Sirius <- long[!is.na(long$fibrosis_pct), ]
#   }
#   if (length(sheets))
#     .write_xlsx(sheets, 'Data_S28_fibrosis_quantification.xlsx')
#   else
#     .stub(28, 'Fibrosis quantification',
#           'PAB_Sirius.xlsx / Adult_Trichrome_Sirius.xlsx not found')
# })

# [REMOVED in v60 — not cited as a Data S; code retained, commented]
# # ---------- S26 : Module-eigengene matrix (RV WGCNA modules) --------- #
# # Per-cell bulk->snRNA WGCNA module eigengenes (committed dependency
# # scWGCNA_bulk2sn_MEs.rds: 61,398 RV nuclei x 29 colour-named modules),
# # joined to the RV snRNA metadata (snRV_ref.rds, barcodes match) and
# # summarised per disease group / per patient / per cell type. NOTE: the
# # committed MEs dependency is RV-ONLY; the pediatric & PAB cross-
# # projection (Fig 8C/E) is produced by the OOM-prone ProjectModules
# # path and is NOT a committed input, so it cannot be emitted deps-only
# # here -- it needs the same kind of one-time precompute as S27 (promote
# # the F8-projected MEs to a dependency). Module colour<->M# mapping is
# # in Data S4 (bulk WGCNA module memberships).
# try({
#   .section(26, 'Module-eigengene matrix (RV WGCNA modules)')
#   me_rds <- file.path(.shared, 'scWGCNA_bulk2sn_MEs.rds')
#   sn_rds <- file.path(.shared, 'snRV_ref.rds')
#   if (file.exists(me_rds) && file.exists(sn_rds)) {
#     MEs <- as.data.frame(readRDS(me_rds), check.names = FALSE)
#     MEs <- MEs[, colnames(MEs) != 'grey', drop = FALSE]
#     obj <- readRDS(sn_rds); md <- obj@meta.data; rm(obj); invisible(gc(FALSE))
#     grc <- intersect(c('group', 'category', 'condition'), names(md))[1]
#     smc <- intersect(c('patient', 'orig.ident', 'sample'), names(md))[1]
#     ctc <- names(md)[grepl('cell.?type|subname|cluster|annot',
#                            names(md), ignore.case = TRUE) &
#                      !grepl('orig', names(md), ignore.case = TRUE)][1]
#     bc  <- intersect(rownames(MEs), rownames(md))
#     if (length(bc) > 0 && !is.na(grc)) {
#       M    <- MEs[bc, , drop = FALSE]; mods <- colnames(M)
#       base <- data.frame(
#         group     = as.character(md[bc, grc]),
#         patient   = if (!is.na(smc)) as.character(md[bc, smc]) else NA,
#         cell_type = if (!is.na(ctc)) as.character(md[bc, ctc]) else NA,
#         stringsAsFactors = FALSE)
#       base$cell_type <- .relabel_subnames(base$cell_type, bc)
#       long <- reshape(cbind(base, M), varying = mods, v.names = 'ME',
#                       timevar = 'module', times = mods, direction = 'long')
#       sh <- list(
#         module_by_group   = aggregate(ME ~ group + module, long, mean),
#         module_by_patient = aggregate(ME ~ patient + group + module,
#                                       long, mean))
#       if (!is.na(ctc))
#         sh$module_by_group_celltype <-
#           aggregate(ME ~ group + cell_type + module, long, mean)
#       .write_xlsx(sh, 'Data_S26_module_eigengene_projection.xlsx')
#     } else {
#       .stub(26, 'Module-eigengene matrix',
#             'MEs/snRV_ref barcodes or group column not resolvable')
#     }
#   } else {
#     .stub(26, 'Module-eigengene matrix',
#           'scWGCNA_bulk2sn_MEs.rds or snRV_ref.rds not found')
#   }
# })

# [REMOVED in v60 — not cited as a Data S; code retained, commented]
# # ---------- S24 : Bulk RNA-seq per-sample QC --------------------------- #
# # Per-patient bulk RNA-seq data-quality overview (addresses reviewer
# # comment on mapped reads / gene detection / RIN).  Parses kallisto
# # pseudoalignment stats from dependencies/shared/BulkRNA/30-*/P*_quant/
# # run_info.json (n_processed, n_pseudoaligned, p_pseudoaligned, n_unique,
# # p_unique), library size + genes detected from BulkRNA/counts.csv, and
# # joins disease group / sex / age / RIN from BulkRNA/metadata.csv.
# try({
#   .section(24, 'Bulk RNA-seq per-sample QC')
#   bulk_dir <- file.path(.shared, 'BulkRNA')
#   ji <- list.files(bulk_dir, 'run_info\\.json$',
#                    recursive = TRUE, full.names = TRUE)
#   if (!length(ji) || !file.exists(file.path(bulk_dir, 'counts.csv'))) {
#     .stub(24, 'Bulk RNA-seq QC',
#           'BulkRNA/run_info.json files or counts.csv not found')
#   } else {
#     .req_json <- .require('jsonlite')
#     .getn <- function(txt, k) {
#       m <- regmatches(txt, regexpr(paste0('"', k,
#                        '"\\s*:\\s*[-0-9.eE+]+'), txt))
#       if (length(m)) as.numeric(sub('.*:\\s*', '', m)) else NA_real_
#     }
#     qc <- do.call(rbind, lapply(ji, function(p) {
#       pid <- regmatches(basename(dirname(p)),
#                         regexpr('P[0-9]+', basename(dirname(p))))[1]
#       if (is.na(pid)) return(NULL)
#       if (.req_json) {
#         j <- tryCatch(jsonlite::fromJSON(p), error = function(e) NULL)
#       } else {
#         txt <- paste(readLines(p, warn = FALSE), collapse = ' ')
#         j <- list(n_processed = .getn(txt, 'n_processed'),
#                   n_pseudoaligned = .getn(txt, 'n_pseudoaligned'),
#                   n_unique = .getn(txt, 'n_unique'),
#                   p_pseudoaligned = .getn(txt, 'p_pseudoaligned'),
#                   p_unique = .getn(txt, 'p_unique'),
#                   n_targets = .getn(txt, 'n_targets'))
#       }
#       if (is.null(j)) return(NULL)
#       data.frame(patient = pid,
#                  n_processed     = as.numeric(j$n_processed),
#                  n_pseudoaligned = as.numeric(j$n_pseudoaligned),
#                  p_pseudoaligned = as.numeric(j$p_pseudoaligned),
#                  n_unique        = as.numeric(j$n_unique),
#                  p_unique        = as.numeric(j$p_unique),
#                  n_targets       = as.numeric(j$n_targets),
#                  stringsAsFactors = FALSE)
#     }))
#     co <- read.csv(file.path(bulk_dir, 'counts.csv'),
#                    check.names = FALSE, row.names = 1)
#     lib <- data.frame(
#       patient          = colnames(co),
#       library_size     = colSums(co),
#       n_genes_detected = colSums(co > 0),
#       median_count_per_detected_gene = apply(co, 2, function(x) {
#         x <- x[x > 0]; if (length(x)) stats::median(x) else NA_real_ }),
#       stringsAsFactors = FALSE)
#     md <- tryCatch(read.csv(file.path(bulk_dir, 'metadata.csv'),
#                             check.names = FALSE, stringsAsFactors = FALSE),
#                    error = function(e) NULL)
#     if (!is.null(md)) {
#       pid_col <- intersect(c('sample_name', 'patient'), names(md))[1]
#       keep    <- c(pid_col,
#                    intersect(c('var_disease', 'var_sex',
#                                 'var_age',     'HHT_age', 'HHT_eth',
#                                 'RIN_score'),
#                               names(md)))
#       md <- md[, keep, drop = FALSE]; names(md)[1] <- 'patient'
#       out <- merge(merge(qc, lib, by = 'patient', all = TRUE),
#                     md, by = 'patient', all.x = TRUE)
#     } else out <- merge(qc, lib, by = 'patient', all = TRUE)
#     out <- out[order(out$patient), , drop = FALSE]; rownames(out) <- NULL
#     .write_xlsx(list(per_sample_QC = out),
#                 'Data_S24_bulk_RNAseq_QC.xlsx')
#   }
# })

message('\nDone. Output directory: ', .outdir)
