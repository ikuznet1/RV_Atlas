##############################################
#### snRNAReanalysis.r
#### Pseudobulk DESeq2, volcano plots, enrichR
#### on the snRNA-seq RV reference object
####
#### Level 1: broad cell types (Names)
#### Level 2: sub cell types   (Subnames_manual)
##############################################

library(Seurat)
library(DESeq2)
library(EnhancedVolcano)
library(enrichR)
library(tidyverse)
library(tibble)
library(dplyr)

dir.create('./output/snRNA', showWarnings = FALSE, recursive = TRUE)


##############################################
#### Load data
##############################################

sn <- readRDS('~/Documents/XeniumWorkflow/snRV_ref.rds')
stopifnot(all(c('Names','Subnames_manual','patient','group') %in% colnames(sn@meta.data)))

# Seurat v5: join layers if split
if (inherits(sn[['RNA']], 'Assay5')) sn <- JoinLayers(sn)

sn$group   <- factor(sn$group,   levels = c('NF','pRV','RVF'))
sn$patient <- as.character(sn$patient)
DefaultAssay(sn) <- 'RNA'

cat('Cells per group:\n');  print(table(sn$group))
cat('Cells per patient:\n'); print(table(sn$patient, sn$group))
cat('Names labels:\n');      print(sort(table(sn$Names)))
cat('Subnames_manual labels:\n');   print(sort(table(sn$Subnames_manual)))


##############################################
#### Shared helper functions
##############################################

cap_pval <- function(df) {
  mutate(df, padj = ifelse(is.na(padj), NA_real_, pmax(padj, 1e-5)))
}

sig_labs <- function(df, n = 30) {
  df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5) %>%
    arrange(padj) %>%
    head(n) %>%
    pull(gene)
}

make_volcano_colors_sn <- function(df) {
  is_sig <- !is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 0.5
  cols <- ifelse(is_sig & df$log2FoldChange  > 0, 'red3',
           ifelse(is_sig & df$log2FoldChange < 0, 'royalblue3', 'grey80'))
  names(cols) <- df$gene
  cols
}

volcano_args <- list(
  x              = 'log2FoldChange',
  y              = 'padj',
  pCutoff        = 0.05,
  FCcutoff       = 0.5,
  pointSize      = 2,
  labSize        = 3.5,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  xlim           = c(-5, 5),
  ylim           = c(0, 5),
  legendPosition = 'none'
)

make_volcano_pdf_sn <- function(dat, filepath, title_str) {
  if (is.null(dat) || nrow(dat) == 0) {
    cat('  Skipping', title_str, '- no data\n'); return(invisible(NULL))
  }
  pdf(filepath, width = 10, height = 8)
  tryCatch(
    print(do.call(EnhancedVolcano,
                  c(list(toptable  = dat,
                         lab       = dat$gene,
                         title     = title_str,
                         selectLab = sig_labs(dat),
                         colCustom = make_volcano_colors_sn(dat)),
                    volcano_args))),
    error = function(e) cat('  Volcano error for', title_str, ':', conditionMessage(e), '\n')
  )
  dev.off()
}

# enrichR
enrichR_dbs <- c('GO_Biological_Process_2023', 'Reactome_2022')
setEnrichrSite('Enrichr')

combine_enrich <- function(enrich_res, direction, comparison) {
  if (is.null(enrich_res)) return(NULL)
  lapply(names(enrich_res), function(db) {
    df <- enrich_res[[db]]
    if (nrow(df) == 0) return(NULL)
    df %>% mutate(Term       = as.character(Term),
                  database   = db,
                  direction  = direction,
                  comparison = comparison)
  }) %>% bind_rows()
}

plot_enrich <- function(enrich_df, title_str) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
  enrich_df %>%
    arrange(Adjusted.P.value) %>%
    head(10) %>%
    mutate(Term = factor(Term, levels = rev(Term))) %>%
    ggplot(aes(x = -log10(Adjusted.P.value), y = Term)) +
    geom_col(fill = 'steelblue') +
    theme_bw() +
    labs(x = '-log10(adj. P)', y = NULL, title = title_str)
}

run_enrich <- function(dat) {
  if (is.null(dat) || nrow(dat) == 0) return(list(up = NULL, down = NULL))
  sig  <- dat %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5)
  up   <- sig %>% filter(log2FoldChange > 0) %>% pull(gene)
  down <- sig %>% filter(log2FoldChange < 0) %>% pull(gene)
  list(
    up   = if (length(up)   > 0) tryCatch(enrichr(up,   enrichR_dbs), error = function(e) NULL) else NULL,
    down = if (length(down) > 0) tryCatch(enrichr(down, enrichR_dbs), error = function(e) NULL) else NULL
  )
}

# Shared DESeq2 runner — called for both levels
run_pseudobulk_deseq2 <- function(obj, label_col, assay_name = 'RNA') {
  obj@meta.data[[label_col]] <- as.character(obj@meta.data[[label_col]])
  obj$pseudobulk_id <- paste0(obj$patient, '..', obj@meta.data[[label_col]])

  pseudo_meta <- unique(obj@meta.data[, c('patient','group', label_col, 'pseudobulk_id')])
  pseudo_meta$group <- factor(pseudo_meta$group, levels = c('NF','pRV','RVF'))
  rownames(pseudo_meta) <- pseudo_meta$pseudobulk_id

  pseudo_counts <- AggregateExpression(
    obj, assays = assay_name, group.by = 'pseudobulk_id',
    slot = 'counts', return.seurat = FALSE
  )[[assay_name]]
  colnames(pseudo_counts) <- gsub('-', '_', sub('^g', '', colnames(pseudo_counts)))
  pseudo_counts <- round(as.matrix(pseudo_counts))
  storage.mode(pseudo_counts) <- 'integer'

  cell_types <- unique(pseudo_meta[[label_col]])
  results_list <- list()

  for (ct in cell_types) {
    ct_samples <- pseudo_meta$pseudobulk_id[pseudo_meta[[label_col]] == ct]
    ct_samples <- ct_samples[ct_samples %in% colnames(pseudo_counts)]
    if (length(ct_samples) < 2) next

    ct_counts <- pseudo_counts[, ct_samples, drop = FALSE]
    ct_meta   <- pseudo_meta[ct_samples, , drop = FALSE]
    if (length(unique(ct_meta$group)) < 2) next

    keep <- rowSums(ct_counts >= 5) >= 2
    ct_counts <- ct_counts[keep, , drop = FALSE]
    if (nrow(ct_counts) < 10) next

    dds <- DESeqDataSetFromMatrix(countData = ct_counts, colData = ct_meta, design = ~ group)
    dds <- tryCatch(DESeq(dds), error = function(e) { cat('DESeq2 failed for', ct, '\n'); NULL })
    if (is.null(dds)) next

    for (comp in list(c('RVF','NF'), c('pRV','NF'), c('RVF','pRV'))) {
      comp_label <- paste0(comp[1], '_vs_', comp[2])
      res <- tryCatch(
        results(dds, contrast = c('group', comp[1], comp[2])) %>%
          as.data.frame() %>%
          tibble::rownames_to_column('gene') %>%
          mutate(cell_type = ct, comparison = comp_label) %>%
          arrange(padj),
        error = function(e) NULL
      )
      if (!is.null(res)) results_list[[paste0(ct, '_', comp_label)]] <- res
    }
    cat('Done:', ct, '\n')
  }
  bind_rows(results_list)
}

# Shared volcano + enrichR loop
run_volcano_enrich <- function(all_deseq, label_col, prefix, enrich_list) {
  cell_types <- unique(all_deseq$cell_type)
  for (ct in cell_types) {
    ct_slug <- gsub('[^A-Za-z0-9]', '_', ct)
    dat_rvf     <- all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_NF')  %>% cap_pval()
    dat_prv     <- all_deseq %>% filter(cell_type == ct, comparison == 'pRV_vs_NF')  %>% cap_pval()
    dat_rvf_prv <- all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_pRV') %>% cap_pval()

    make_volcano_pdf_sn(dat_rvf,
      paste0('./output/snRNA/', prefix, ct_slug, '_volcano_RVF_vs_NF.pdf'),
      paste0(ct, ': RVF vs NF'))
    make_volcano_pdf_sn(dat_prv,
      paste0('./output/snRNA/', prefix, ct_slug, '_volcano_pRV_vs_NF.pdf'),
      paste0(ct, ': pRV vs NF'))
    make_volcano_pdf_sn(dat_rvf_prv,
      paste0('./output/snRNA/', prefix, ct_slug, '_volcano_RVF_vs_pRV.pdf'),
      paste0(ct, ': RVF vs pRV'))

    # enrichR
    er_rvf     <- run_enrich(dat_rvf)
    er_prv     <- run_enrich(dat_prv)
    er_rvf_prv <- run_enrich(dat_rvf_prv)

    ct_enrich <- bind_rows(
      combine_enrich(er_rvf$up,       'up',   'RVF_vs_NF'),
      combine_enrich(er_rvf$down,     'down', 'RVF_vs_NF'),
      combine_enrich(er_prv$up,       'up',   'pRV_vs_NF'),
      combine_enrich(er_prv$down,     'down', 'pRV_vs_NF'),
      combine_enrich(er_rvf_prv$up,   'up',   'RVF_vs_pRV'),
      combine_enrich(er_rvf_prv$down, 'down', 'RVF_vs_pRV')
    )

    if (!is.null(ct_enrich) && nrow(ct_enrich) > 0) {
      ct_enrich$cell_type <- ct
      enrich_list[[ct]]   <- ct_enrich

      # Per cell type enrichR plots
      pdf(paste0('./output/snRNA/', prefix, ct_slug, '_enrichr.pdf'), width = 10, height = 6)
      for (comp_label in c('RVF_vs_NF','pRV_vs_NF','RVF_vs_pRV')) {
        for (dir_label in c('up','down')) {
          er_sub <- ct_enrich %>%
            filter(comparison == comp_label, direction == dir_label)
          for (db in enrichR_dbs) {
            p <- plot_enrich(er_sub %>% filter(database == db),
                             paste(ct, comp_label, toupper(dir_label), db))
            if (!is.null(p)) print(p)
          }
        }
      }
      dev.off()
    }
  }
  enrich_list
}


##############################################
#### Level 1: DESeq2 by broad cell type (Names)
##############################################

cat('\n--- Level 1: pseudobulk DESeq2 by Names ---\n')
all_deseq_l1 <- run_pseudobulk_deseq2(sn, label_col = 'Names', assay_name = 'RNA')
write.csv(all_deseq_l1, './output/snRNA/names_pseudobulk_deseq2.csv', row.names = FALSE)
cat('Level 1 DESeq2 done:', nrow(all_deseq_l1), 'rows\n')

cat('\n--- Level 1: Volcano plots + enrichR by Names ---\n')
enrich_l1 <- list()
enrich_l1 <- run_volcano_enrich(all_deseq_l1, label_col = 'Names',
                                 prefix = '', enrich_list = enrich_l1)

if (length(enrich_l1) > 0) {
  all_enrich_l1 <- bind_rows(enrich_l1)
  write.csv(all_enrich_l1, './output/snRNA/names_all_enrichr.csv', row.names = FALSE)
}


##############################################
#### Level 2: DESeq2 by sub cell type (Subnames_manual)
##############################################

cat('\n--- Level 2: pseudobulk DESeq2 by Subnames_manual ---\n')
all_deseq_l2 <- run_pseudobulk_deseq2(sn, label_col = 'Subnames_manual', assay_name = 'RNA')
write.csv(all_deseq_l2, './output/snRNA/subnames_pseudobulk_deseq2.csv', row.names = FALSE)
cat('Level 2 DESeq2 done:', nrow(all_deseq_l2), 'rows\n')

cat('\n--- Level 2: Volcano plots + enrichR by Subnames_manual ---\n')
enrich_l2 <- list()
enrich_l2 <- run_volcano_enrich(all_deseq_l2, label_col = 'Subnames_manual',
                                 prefix = 'sub_', enrich_list = enrich_l2)

if (length(enrich_l2) > 0) {
  all_enrich_l2 <- bind_rows(enrich_l2)
  write.csv(all_enrich_l2, './output/snRNA/subnames_all_enrichr.csv', row.names = FALSE)
}

##############################################
#### Marker genes by Subnames_manual
##############################################

cat('\n--- Marker genes by Subnames_manual (one-vs-rest) ---\n')

Idents(sn) <- 'Subnames_manual'

# Run FindAllMarkers: Wilcoxon one-vs-rest, only positive markers, min.pct 0.25
all_markers <- FindAllMarkers(
  sn,
  assay        = 'RNA',
  only.pos     = TRUE,
  min.pct      = 0.25,
  logfc.threshold = 0.25,
  test.use     = 'wilcox',
  return.thresh = 0.05
)

all_markers <- all_markers %>%
  arrange(cluster, p_val_adj, desc(avg_log2FC))

write.csv(all_markers, './output/snRNA/subnames_markers_all.csv', row.names = FALSE)
cat('Total markers found:', nrow(all_markers), '\n')

# Top 20 markers per cluster for quick reference
top20 <- all_markers %>%
  group_by(cluster) %>%
  slice_head(n = 20) %>%
  ungroup()

write.csv(top20, './output/snRNA/subnames_markers_top20.csv', row.names = FALSE)

# Summary table: top 5 genes per cluster (for annotation)
top5_summary <- all_markers %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  summarise(
    n_markers   = n(),
    top5_genes  = paste(gene, collapse = ', '),
    top5_log2FC = paste(round(avg_log2FC, 2), collapse = ', '),
    top5_padj   = paste(format(p_val_adj, digits = 2), collapse = ', '),
    .groups     = 'drop'
  )

write.csv(top5_summary, './output/snRNA/subnames_markers_top5_summary.csv', row.names = FALSE)
cat('Top 5 summary:\n')
print(as.data.frame(top5_summary), right = FALSE)

# Heatmap of top 5 markers per cluster
top5_genes <- all_markers %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  pull(gene) %>%
  unique()

# Downsample for heatmap (max 100 cells per cluster)
set.seed(42)
cells_ds <- unlist(lapply(levels(Idents(sn)), function(cl) {
  cl_cells <- WhichCells(sn, idents = cl)
  if (length(cl_cells) > 100) sample(cl_cells, 100) else cl_cells
}))

pdf('./output/snRNA/subnames_markers_heatmap.pdf', width = 20, height = 16)
tryCatch({
  p <- DoHeatmap(
    subset(sn, cells = cells_ds),
    features = top5_genes,
    group.by = 'Subnames_manual',
    size     = 3,
    angle    = 90
  ) + theme(axis.text.y = element_text(size = 6))
  print(p)
}, error = function(e) cat('Heatmap error:', conditionMessage(e), '\n'))
dev.off()

# Dot plot of top 3 markers per cluster
top3_genes <- all_markers %>%
  group_by(cluster) %>%
  slice_head(n = 3) %>%
  pull(gene) %>%
  unique()

pdf('./output/snRNA/subnames_markers_dotplot.pdf', width = 18, height = 10)
tryCatch({
  p <- DotPlot(sn, features = top3_genes, group.by = 'Subnames_manual') +
    RotatedAxis() +
    theme(axis.text.x = element_text(size = 7))
  print(p)
}, error = function(e) cat('DotPlot error:', conditionMessage(e), '\n'))
dev.off()

cat('\nMarker gene analysis complete.\n')
cat('Outputs in ./output/snRNA/\n')

cat('\nsnRNAReanalysis complete.\n')
cat('Outputs in ./output/snRNA/\n')
