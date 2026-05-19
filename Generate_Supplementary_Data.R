##############################################################################
# Generate_Supplementary_Data.R
#
# Builds all supplementary data files (Data S1-S22) for v53 of the RV atlas
# paper, using the new numbering scheme established during the v52 -> v53
# reorganization.
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
REBUILD_BULK_WGCNA        <- REBUILD # S3  — bulk_heart_modules
REBUILD_BULK_TRENDS       <- REBUILD # S4  — NF_2_pRV_*_2_RVF_*
REBUILD_PAH               <- REBUILD # S5
REBUILD_SN_ANNOT          <- REBUILD # S6
REBUILD_SN_SUBCLUST_MARKS <- REBUILD # S7  — FindAllMarkers per subclust obj
REBUILD_SN_PSEUDOBULK_L1  <- REBUILD # S8  — sn lineage pseudobulk
REBUILD_XEN_PSEUDOBULK_L1 <- REBUILD # S9
REBUILD_XEN_SUBCLUST_MARKS<- REBUILD # S10
REBUILD_XEN_NICHE         <- REBUILD # S11 — composition tables
REBUILD_EC_HDWGCNA        <- REBUILD # S13
REBUILD_PAB               <- REBUILD # S15
REBUILD_CX3CR1            <- REBUILD # S16 — CX3CR1-tdTomato bulk DEGs
REBUILD_PEDS              <- REBUILD # S17 / S18
REBUILD_CHEA              <- REBUILD # S21
REBUILD_SN_PSEUDOBULK_L2  <- REBUILD # S22

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
# Map snRNA subcluster labels (Cm#, Fb#) to v57 manuscript names.  CM mapping is
# derived from top-marker LFC in CM_marks.csv.  FB mapping is partial: Fb1 is the
# most-abundant baseline (Resident), Fb7 is too rare for label transfer
# (Anti-fibrotic per manuscript); Fb2-Fb6 await marker confirmation.
.SUBNAME_RMAP <- c(
  Cm2 = 'CM MYH6',   Cm4 = 'CM HMGCS2', Cm5 = 'CM HTR4',
  Cm6 = 'CM HAND2',  Cm8 = 'CM NPP',    Cm9 = 'CM KCNJ3',
  Cm10 = 'CM BetaMHC',
  Fb1 = 'Fb Resident', Fb7 = 'Fb Anti-fibrotic'
)
.rename_subnames <- function(df, src_col = 'Subnames_manual') {
  if (is.null(df) || !src_col %in% names(df)) return(df)
  raw <- as.character(df[[src_col]])
  df$Subnames_renamed <- ifelse(raw %in% names(.SUBNAME_RMAP),
                                .SUBNAME_RMAP[raw], raw)
  df
}
.apply_rmap_inplace <- function(df, src_col) {
  if (is.null(df) || !src_col %in% names(df)) return(df)
  raw <- as.character(df[[src_col]])
  df[[src_col]] <- ifelse(raw %in% names(.SUBNAME_RMAP),
                          .SUBNAME_RMAP[raw], raw)
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
.build_bulk_trends <- function(force = FALSE, padj = 0.05, lfc = 0.5) {
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
  .sig <- function(d, dir) {
    d %>% filter(!is.na(padj), padj < padj, {if (dir == 'up') log2FoldChange > lfc else log2FoldChange < -lfc}) %>%
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
  if (level == 'sublineage')
    obj <- subset(obj, subset = cell_types_subclustering != 'Unassigned')
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
# MAIN PIPELINE — Data S01..S22
# ==========================================================================

# ---------- S1 : Bulk RNA-seq DEGs ---------------------------------------- #
try({
  .section(1, 'Bulk RNA-seq DEGs')
  if (REBUILD_BULK_DEGS) .build_bulk_degs(force = TRUE)
  sheets <- list(
    NF_vs_pRV  = .read_csv_safe(file.path(.bulkdir, 'NF_vs_pRV_deseq.csv')),
    NF_vs_RVF  = .read_csv_safe(file.path(.bulkdir, 'NF_vs_RVF_deseq.csv')),
    pRV_vs_RVF = .read_csv_safe(file.path(.bulkdir, 'pRV_vs_RVF_deseq.csv'))
  )
  .write_xlsx(sheets, 'Data_S01_bulk_DEGs.xlsx')
})

# ---------- S2 : Bulk pathway enrichment ---------------------------------- #
try({
  .section(2, 'Bulk pathway enrichment')
  cache_path <- file.path(.cachedir, 'S02_bulk_enrichment_cache.rds')
  if (!REBUILD_BULK_ENRICH && file.exists(cache_path)) {
    sheets <- readRDS(cache_path)
  } else {
    if (!.require('enrichR')) { .stub(2, 'Bulk pathway enrichment', 'install enrichR then rerun'); sheets <- NULL }
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
  if (!is.null(sheets)) .write_xlsx(sheets, 'Data_S02_bulk_enrichment.xlsx')
})

# ---------- S3 : Bulk WGCNA modules --------------------------------------- #
try({
  .section(3, 'Bulk WGCNA modules')
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
    .write_xlsx(sheets, 'Data_S03_bulk_WGCNA_modules.xlsx')
  }
})

# ---------- S4 : Bulk stepwise trends ------------------------------------- #
try({
  .section(4, 'Bulk stepwise trend genes')
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
  .write_xlsx(sheets, 'Data_S04_bulk_stepwise_trends.xlsx')
})

# ---------- S5 : PAH comparison (concordant / discordant + enrichment) ---- #
# Manuscript claim: concordant DEGs (shared mitochondrial / IFN suppression)
# and discordant DEGs (CD163, mitochondrial translation) between our RV cohort
# and the GSE240921 PAH cohort, with pathway enrichment of each set.
try({
  .section(5, 'PAH comparison: concordance + discordance')
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
  .write_xlsx(sheets, 'Data_S05_PAH_comparison.xlsx')
})

# ---------- S6 : snRNA annotations --------------------------------------- #
try({
  .section(6, 'snRNA cell type annotations')
  rds <- file.path(.shared, 'snRV_ref.rds')
  if (!file.exists(rds)) rds <- file.path(.shared, 'Post_R3_FINAL_with_counts.rds')
  cache_path <- file.path(.cachedir, 'S06_sn_annotations.csv')
  if (!REBUILD_SN_ANNOT && file.exists(cache_path)) {
    annot <- .rename_subnames(.read_csv_safe(cache_path))
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
    .stub(6, 'snRNA annotations', paste('Seurat RDS not found at', rds))
  }
})

# ---------- S7 : snRNA subclustering markers ----------------------------- #
try({
  .section(7, 'snRNA subclustering markers')
  .build_sn_subclust_marks(force = REBUILD_SN_SUBCLUST_MARKS)
  sheets <- list(
    cm      = .read_csv_safe(file.path(.shared, 'cm_subclust_marks.csv')),
    fb      = .read_csv_safe(file.path(.shared, 'fb_subclust_marks.csv')),
    ec      = .read_csv_safe(file.path(.shared, 'ec_subclust_marks.csv')),
    myeloid = .read_csv_safe(file.path(.shared, 'myeloid_subclust_marks.csv')),
    mural   = .read_csv_safe(file.path(.shared, 'mural_subclust_marks.csv'))
  )
  .write_xlsx(sheets, 'Data_S07_snRNA_subclustering_markers.xlsx')
})

# ---------- S8 : snRNA lineage pseudobulk ------------------------------- #
try({
  .section(8, 'snRNA lineage pseudobulk DEGs')
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

# ---------- S9 : Xenium lineage pseudobulk + sn concordance ------------ #
# v57 promise: "differential expression results within each spatially resolved
# cell type, including cross-platform concordance with snRNA-seq".  We attach
# snRNA log2FC and padj for each (gene, cell_type, comparison) tuple where
# matched at the lineage level.
try({
  .section(9, 'Xenium lineage pseudobulk DEGs')
  if (REBUILD_XEN_PSEUDOBULK_L1) .build_xenium_pseudobulk('lineage', force = TRUE)
  df <- .read_csv_safe(file.path(.xcache, 'xenium_pseudobulk_lineage_deseq2.csv'))
  if (!is.null(df)) {
    ct_col <- intersect(c('cell_type', 'lineage'), names(df))[1]
    cp_col <- intersect(c('comparison', 'contrast'), names(df))[1]
    # Cross-platform concordance: join with snRNA lineage pseudobulk
    sn <- .read_csv_safe(file.path(.shared, 'sn_pseudobulk_lineage_deseq2.csv'))
    if (!is.null(sn) && !is.na(ct_col) && !is.na(cp_col)) {
      sn_ct <- intersect(c('cell_type', 'lineage'), names(sn))[1]
      sn_cp <- intersect(c('comparison', 'contrast'), names(sn))[1]
      if (!is.na(sn_ct) && !is.na(sn_cp) && 'gene' %in% names(sn)) {
        sn_pick <- sn[, c('gene', sn_ct, sn_cp, 'log2FoldChange', 'padj')]
        names(sn_pick) <- c('gene', '_ct', '_cp', 'snRNA_log2FoldChange', 'snRNA_padj')
        df$`_ct` <- df[[ct_col]]; df$`_cp` <- df[[cp_col]]
        df <- merge(df, sn_pick, by = c('gene', '_ct', '_cp'), all.x = TRUE)
        df$`_ct` <- NULL; df$`_cp` <- NULL
        df$concordant_sign <- !is.na(df$snRNA_log2FoldChange) & !is.na(df$log2FoldChange) &
                              sign(df$snRNA_log2FoldChange) == sign(df$log2FoldChange)
      }
    }
    sheets <- list(all = df)
    if (!is.na(ct_col) && !is.na(cp_col)) {
      for (ct in sort(unique(df[[ct_col]]))) for (cp in sort(unique(df[[cp_col]]))) {
        sub <- df[df[[ct_col]] == ct & df[[cp_col]] == cp, , drop = FALSE]
        if (nrow(sub) == 0) next
        nm <- substr(gsub('[^A-Za-z0-9_]+', '_', paste0(ct, '__', cp)), 1, 31)
        sheets[[nm]] <- sub
      }
    }
    .write_xlsx(sheets, 'Data_S09_xenium_lineage_pseudobulk.xlsx')
  }
})

# ---------- S10 : Xenium subclustering (assignments + markers) --------- #
# v57 promise: "subtype assignments and marker genes for spatially resolved
# cell populations within each of the 12 major cell types after Proseg
# resegmentation, RCTD annotation, and SPLIT correction."
try({
  .section(10, 'Xenium subclustering assignments + markers')
  .build_xenium_subclust_marks(force = REBUILD_XEN_SUBCLUST_MARKS)
  sheets <- list()
  # 1. Subtype assignments (per-cell barcode -> subtype) from Xenium_metadata.csv
  meta_csv <- file.path(.shared, 'Xenium_metadata.csv')
  if (file.exists(meta_csv)) {
    md <- .read_csv_safe(meta_csv)
    if (!is.null(md)) sheets[['subtype_assignments']] <- md
  }
  # 2. FindAllMarkers output (gene-level marker tables per subtype)
  marks_csv <- file.path(.xcache, 'xenium_corrected_findallmarkers.csv')
  if (file.exists(marks_csv) && file.size(marks_csv) > 100) {
    sheets[['markers']] <- .read_csv_safe(marks_csv)
  }
  # 3. Reference marker tables (snRNA lineage-level for cross-platform mapping)
  for (extra in c('CM_marks.csv', 'EC_marks.csv', 'PAB_marks.csv')) {
    p <- file.path(.shared, extra)
    if (file.exists(p)) sheets[[substr(sub('\\.csv$', '', extra), 1, 31)]] <- .read_csv_safe(p)
  }
  if (length(sheets) == 0) {
    .stub(10, 'Xenium subclustering',
          'no marker / metadata CSVs found; set REBUILD_XEN_SUBCLUST_MARKS=TRUE')
  } else {
    .write_xlsx(sheets, 'Data_S10_xenium_subclustering.xlsx')
  }
})

# ---------- S11 : Xenium niche composition ----------------------------- #
# v57 promise: cell type proportions within each spatial niche (k=100 NN,
# CLR-transformed, collapsed CM subtypes), niche cluster assignments, and
# niche abundance across NF / pRV / RVF disease states.
try({
  .section(11, 'Xenium niche composition')
  md <- .read_csv_safe(file.path(.shared, 'Xenium_metadata.csv'))
  if (is.null(md)) {
    .stub(11, 'Xenium niche composition', 'Xenium_metadata.csv not found')
  } else {
    # Rename the unnamed first column (cell barcode) -> 'barcode'
    nm1 <- names(md)[1]
    if (is.na(nm1) || nm1 == '' || nm1 == '...1') names(md)[1] <- 'barcode'
    keep <- intersect(c('barcode', 'x_centroid', 'y_centroid',
                        'cell_type_rctd_doublet', 'group', 'patient',
                        'cell_types_manual', 'cell_types_subclustering',
                        'niche_manual'), names(md))
    assignments <- md[, keep, drop = FALSE]
    sheets <- list(per_cell_assignments = assignments)

    if (all(c('niche_manual', 'cell_types_manual') %in% names(md))) {
      tab <- table(md$niche_manual, md$cell_types_manual)
      prop_per_niche <- as.data.frame.matrix(prop.table(tab, margin = 1))
      prop_per_niche$niche <- rownames(prop_per_niche)
      prop_per_niche <- prop_per_niche[, c('niche', setdiff(names(prop_per_niche), 'niche'))]
      sheets$cell_type_proportions_per_niche <- prop_per_niche
    }

    if (all(c('niche_manual', 'group') %in% names(md))) {
      tab2 <- table(md$niche_manual, md$group)
      ab <- as.data.frame.matrix(prop.table(tab2, margin = 2))  # proportion within each group
      ab$niche <- rownames(ab)
      ab <- ab[, c('niche', setdiff(names(ab), 'niche'))]
      sheets$niche_abundance_by_group <- ab
    }

    .write_xlsx(sheets, 'Data_S11_xenium_niche.xlsx')
  }
})

# ---------- S12 : respirometry ---------------------------------------- #
try({
  .section(12, 'Mitochondrial respirometry')
  mito_path <- file.path(.shared, 'mito_data.xlsx')
  if (file.exists(mito_path) && .require('openxlsx')) {
    sheet_names <- openxlsx::getSheetNames(mito_path)
    sheets <- setNames(
      lapply(sheet_names, function(s) openxlsx::read.xlsx(mito_path, sheet = s)),
      substr(sheet_names, 1, 31)
    )
    .write_xlsx(sheets, 'Data_S12_mito_respirometry.xlsx')
  } else {
    .stub(12, 'Mitochondrial respirometry',
          'mito_data.xlsx not found at dependencies/shared/')
  }
})

# ---------- S13 : Myeloid transcriptional program gene lists ----------- #
# Sourced from Figure_6.R's curated 4-program tibble (GR_homeostatic,
# HIF_vascular, NFkB_MHCII, IFNg_AP) written to output/Figure_6/.  If Figure_6
# hasn't run yet, fall back to a stub.
try({
  .section(13, 'Myeloid transcriptional program gene lists')
  myeloid_csv <- '/Users/ikuz/Documents/RV_Atlas/output/Figure_6/SuppTable_Myeloid_Programs.csv'
  if (file.exists(myeloid_csv)) {
    df <- .read_csv_safe(myeloid_csv)
    sheets <- list(programs = df)
    if (!is.null(df) && 'program' %in% names(df)) {
      for (p in sort(unique(df$program))) {
        sheets[[substr(p, 1, 31)]] <- df[df$program == p, , drop = FALSE]
      }
    }
    .write_xlsx(sheets, 'Data_S13_myeloid_programs.xlsx')
  } else {
    .stub(13, 'Myeloid transcriptional program gene lists',
          'run Figure_6.R first to generate output/Figure_6/SuppTable_Myeloid_Programs.csv')
  }
})

# ---------- S14 : EC hdWGCNA ------------------------------------------- #
try({
  .section(14, 'EC hdWGCNA modules')
  rds <- file.path(.shared, 'EC_hdWGCNA_by_celltype.rds')
  if (REBUILD_EC_HDWGCNA) message('  (rebuild EC hdWGCNA is upstream in Figure_5/7.R; keeping cached rds)')
  if (file.exists(rds)) {
    ec <- readRDS(rds); sheets <- list()
    if (is.data.frame(ec)) {
      sheets$modules <- ec
    } else if (is.list(ec)) {
      for (nm in names(ec)) {
        x <- ec[[nm]]
        if (is.data.frame(x)) { sheets[[substr(nm, 1, 31)]] <- x }
        else if (.require('hdWGCNA') && inherits(x, 'Seurat')) {
          mods <- tryCatch(hdWGCNA::GetModules(x), error = function(e) NULL)
          if (!is.null(mods)) sheets[[substr(paste0(nm, '_mods'), 1, 31)]] <- as.data.frame(mods)
        }
      }
    }
    if (length(sheets) == 0)
      sheets$README <- data.frame(note = 'EC_hdWGCNA object contained no extractable tables')
    .write_xlsx(sheets, 'Data_S14_EC_hdWGCNA.xlsx')
  } else {
    .stub(14, 'EC hdWGCNA', 'EC_hdWGCNA_by_celltype.rds not found')
  }
})

# ---------- S15 : CellNEST --------------------------------------------- #
.stub(15, 'CellNEST LR interactions',
      'export NEST output to dependencies/shared/CellNEST_LR.csv and rerun')

# ---------- S16 : Murine PAB ------------------------------------------- #
try({
  .section(16, 'Murine PAB dataset')
  rds <- file.path(.shared, 'PAB_data_clean.rds')
  sheets <- list()
  if (file.exists(rds) && .require('Seurat')) {
    pab <- readRDS(rds)
    if (inherits(pab, 'Seurat')) {
      md <- pab@meta.data %>% tibble::rownames_to_column('cell_id')
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
  pab_ec   <- .read_csv_safe(file.path(.bulkdir, 'PAB_EC.csv'))
  pab_mark <- .read_csv_safe(file.path(.bulkdir, 'PAB_marks.csv'))
  if (!is.null(pab_ec))   sheets$PAB_EC      <- pab_ec
  if (!is.null(pab_mark)) sheets$PAB_markers <- pab_mark
  .write_xlsx(sheets, 'Data_S16_PAB_murine.xlsx')
})

# ---------- S17 : CX3CR1 ----------------------------------------------- #
# Derived in-script from kallisto quants at dependencies/shared/Mouse_PAB_Myeloid/
# (same raw inputs as Supplementary_Figure_5.R S5I/J).  Three contrasts on the
# combined Origin x Type group factor.
try({
  .section(17, 'CX3CR1-tdTomato lineage tracing bulk DEGs')
  .cx3dir  <- file.path(.shared, 'Mouse_PAB_Myeloid')
  .cx3meta <- file.path(.cx3dir, 'meta.csv')
  .cx3t2g  <- file.path(.cx3dir, 't2g.txt')
  cx3_outs <- list(
    RV_PAB_vs_RV_Sham  = file.path(.shared, 'CX3CR1_RV_PAB_vs_RV_Sham.csv'),
    RV_PAB_vs_LV_PAB   = file.path(.shared, 'CX3CR1_RV_PAB_vs_LV_PAB.csv'),
    RV_Sham_vs_LV_Sham = file.path(.shared, 'CX3CR1_RV_Sham_vs_LV_Sham.csv')
  )
  rebuild_needed <- REBUILD_CX3CR1 || !all(file.exists(unlist(cx3_outs)))
  if (rebuild_needed) {
    if (!file.exists(.cx3meta) || !file.exists(.cx3t2g)) {
      .stub(17, 'CX3CR1 lineage tracing bulk',
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
      write.csv(sh('RV_PAB',  'LV_PAB' ),  cx3_outs$RV_PAB_vs_LV_PAB,   row.names = FALSE)
      write.csv(sh('RV_Sham', 'LV_Sham'),  cx3_outs$RV_Sham_vs_LV_Sham, row.names = FALSE)
      message('  wrote ', paste(unlist(cx3_outs), collapse = ', '))
    }
  }
  sheets <- Filter(Negate(is.null), lapply(cx3_outs, .read_csv_safe))
  if (length(sheets) > 0)
    .write_xlsx(sheets, 'Data_S17_CX3CR1_tracing_bulk.xlsx')
  else
    .stub(17, 'CX3CR1 lineage tracing bulk', 'no contrasts computed')
})

# ---------- S18 : Pediatric annotations -------------------------------- #
try({
  .section(18, 'Pediatric annotations')
  rds <- file.path(.shared, 'all_peds_data.rds')
  if (file.exists(rds) && .require('Seurat')) {
    peds <- readRDS(rds)
    if (inherits(peds, 'Seurat')) {
      md <- peds@meta.data %>% tibble::rownames_to_column('cell_id')
      .write_xlsx(list(peds_annotations = md), 'Data_S18_peds_annotations.xlsx')
    }
    rm(peds); gc()
  } else {
    .stub(18, 'Pediatric annotations', 'all_peds_data.rds not found')
  }
})

# ---------- S19 : Pediatric DEGs --------------------------------------- #
try({
  .section(19, 'Pediatric DEGs')
  rds        <- file.path(.shared, 'all_peds_data.rds')
  cache_path <- file.path(.cachedir, 'S19_peds_DEGs.rds')
  if (!REBUILD_PEDS && file.exists(cache_path)) {
    .write_xlsx(readRDS(cache_path), 'Data_S19_peds_DEGs.xlsx')
  } else if (!file.exists(rds)) {
    .stub(19, 'Pediatric DEGs', 'all_peds_data.rds not found')
  } else if (!.require('Seurat')) {
    .stub(19, 'Pediatric DEGs', 'Seurat not installed')
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
    .write_xlsx(sheets, 'Data_S19_peds_DEGs.xlsx')
    rm(peds); gc()
  }
})

# ---------- S20 : Xenium gene panel + imputation provenance ----------- #
# Merged file (panel sheet + provenance sheet + summary + historical) per v57
# S20 promise.  The standalone S19 panel xlsx was removed in this version.
try({
  .section(20, 'Xenium gene panel + provenance')
  panel <- .read_panel()
  xen_l <- .read_csv_safe(file.path(.xcache, 'xenium_pseudobulk_lineage_deseq2.csv'))
  xen_s <- .read_csv_safe(file.path(.xcache, 'xenium_pseudobulk_sublineage_deseq2.csv'))

  if (is.null(panel) || (is.null(xen_l) && is.null(xen_s))) {
    .stub(20, 'Xenium gene provenance',
          'XeniumPanel.csv or xenium pseudobulk CSVs missing')
  } else {
    panel_genes <- unique(trimws(panel[[1]]))
    panel_genes <- panel_genes[nzchar(panel_genes)]

    .xen_genes <- function(df) if (is.null(df)) character() else
      unique(trimws(df[[intersect(c('gene','Gene','symbol'), names(df))[1]]]))

    g_l <- .xen_genes(xen_l)
    g_s <- .xen_genes(xen_s)
    all_xen <- sort(unique(c(g_l, g_s)))

    provenance <- tibble::tibble(
      gene = all_xen,
      on_panel         = gene %in% panel_genes,
      provenance       = ifelse(on_panel, 'panel_measured', 'SpaGE_imputed'),
      in_lineage_DE    = gene %in% g_l,
      in_sublineage_DE = gene %in% g_s
    )

    # Off-panel genes referenced in v51/v52 Xenium analyses (historical transparency)
    historical_imputed <- tibble::tibble(
      gene = c('POSTN','COL1A1','FKBP5','TGFBR3','MLIP','VEGFA','MDM2'),
      note = c('v51/v52 stromal Xenium only; snRNA+bulk independent',
               'v51/v52 stromal Xenium only; snRNA+bulk independent',
               'CM module gene; off-panel, snRNA+bulk only',
               'stromal DEG; off-panel, snRNA+bulk only',
               'CM module gene; off-panel, snRNA+bulk only',
               'pericyte Phase 2 hit (v51/v52); imputed, not retained in v53',
               'pericyte Phase 2 hit (v51/v52); imputed, not retained in v53')
    )
    historical_imputed$still_used_in_v53 <- historical_imputed$gene %in% all_xen

    panel_sheet <- tibble::tibble(gene = panel_genes, in_panel = TRUE)

    summary_tbl <- tibble::tibble(
      metric = c('panel_gene_count',
                 'genes_in_xenium_lineage_DE',
                 'genes_in_xenium_sublineage_DE',
                 'total_unique_xenium_DE_genes',
                 'xenium_DE_genes_on_panel',
                 'xenium_DE_genes_off_panel'),
      value  = c(length(panel_genes),
                 length(g_l),
                 length(g_s),
                 length(all_xen),
                 sum(provenance$on_panel),
                 sum(!provenance$on_panel))
    )

    sheets <- list(
      summary            = summary_tbl,
      provenance         = provenance,
      panel_genes        = panel_sheet,
      historical_imputed = historical_imputed
    )
    .write_xlsx(sheets, 'Data_S20_xenium_gene_provenance.xlsx')
  }
})

# ---------- S21 : ChEA TF enrichment ---------------------------------- #
try({
  .section(21, 'ChEA TF enrichment')
  cache_path <- file.path(.cachedir, 'S21_ChEA_cache.rds')
  if (!REBUILD_CHEA && file.exists(cache_path)) {
    .write_xlsx(readRDS(cache_path), 'Data_S21_ChEA_TF_enrichment.xlsx')
  } else if (!.require('enrichR')) {
    .stub(21, 'ChEA TF enrichment', 'install enrichR')
  } else {
    dbs <- c('ChEA_2022', 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP_X',
             'TRRUST_Transcription_Factors_2019')
    deg_files <- list(NF_vs_pRV  = file.path(.bulkdir, 'NF_vs_pRV_deseq.csv'),
                      pRV_vs_RVF = file.path(.bulkdir, 'pRV_vs_RVF_deseq.csv'),
                      NF_vs_RVF  = file.path(.bulkdir, 'NF_vs_RVF_deseq.csv'))
    sheets <- list()
    for (nm in names(deg_files)) {
      d <- .read_csv_safe(deg_files[[nm]]); if (is.null(d)) next
      gcol <- intersect(c('gene', 'Gene', 'symbol'), names(d))[1]
      if (is.na(gcol)) gcol <- names(d)[1]  # CSV with unnamed rownames column
      sig <- d %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5)
      for (dir in c('UP', 'DOWN')) {
        gs <- if (dir == 'UP') sig %>% filter(log2FoldChange > 0) %>% pull(!!sym(gcol))
              else             sig %>% filter(log2FoldChange < 0) %>% pull(!!sym(gcol))
        if (length(gs) == 0) next
        e <- tryCatch(enrichR::enrichr(gs, dbs), error = function(e) NULL)
        e <- Filter(function(x) is.data.frame(x) && nrow(x) > 0 && is.character(x$Term), e)
        if (length(e) > 0)
          sheets[[substr(paste0(nm, '_', dir), 1, 31)]] <- bind_rows(e, .id = 'db')
      }
    }
    saveRDS(sheets, cache_path)
    .write_xlsx(sheets, 'Data_S21_ChEA_TF_enrichment.xlsx')
  }
})

# ---------- S22 : snRNA sublineage pseudobulk ------------------------- #
try({
  .section(22, 'snRNA sublineage pseudobulk DEGs')
  if (REBUILD_SN_PSEUDOBULK_L2) .build_sn_pseudobulk('sublineage', force = TRUE)
  df <- .read_csv_safe(file.path(.shared, 'sn_pseudobulk_sublineage_deseq2.csv'))
  if (!is.null(df)) {
    ct_col <- intersect(c('cell_type', 'subtype', 'subname'), names(df))[1]
    if (!is.na(ct_col)) df <- .apply_rmap_inplace(df, ct_col)
    .write_xlsx(list(all = df), 'Data_S22_snRNA_sublineage_pseudobulk.xlsx')
  }
})

# ---------- S23 : WGA cardiomyocyte cross-sectional area --------------- #
# Per-cell Minimum Feret diameter (proxy for CSA) measured by wheat-germ
# agglutinin staining of human RV sections (NF n=4, pRV n=4, RVF n=4;
# 17,864 CMs total).  Source: dependencies/shared/wga_minferet_human_rv.csv
try({
  .section(23, 'WGA cardiomyocyte cross-sectional area')
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
      .write_xlsx(sheets, 'Data_S23_WGA_CM_cross_sectional_area.xlsx')
    } else {
      .stub(23, 'WGA cross-sectional area', 'wga_minferet_human_rv.csv empty')
    }
  } else {
    .stub(23, 'WGA cross-sectional area',
          'wga_minferet_human_rv.csv not found in dependencies/shared/')
  }
})

message('\nDone. Output directory: ', .outdir)
