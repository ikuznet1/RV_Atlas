###############################################################################
## Test script: re-run F5 Panel C M2 enrichR using the v52 CM object
## (cm_new_subclust.rds) to see if it replicates the published v57 image.
##
## Reads cm_new_subclust.rds directly from atlas backup (no copy required).
## Writes test PDFs to ./output/_test_F5C/ so nothing in the normal outputs is
## clobbered.
###############################################################################
suppressPackageStartupMessages({
  library(Seurat)
  library(WGCNA)
  library(enrichR)
  library(ggplot2)
  library(forcats)
  library(viridis)
  library(scales)
})

source('./helper_scripts/_shared_helpers.R')

OUT <- './output/_test_F5C'
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

CM_PATH <- '~/Downloads/hdWGCNA_TOM/RV_Atlas_Backup/dependencies/shared/cm_new_subclust.rds'
message('Loading cm_new_subclust.rds (~1.6 GB)...')
seurat_ref <- readRDS(path.expand(CM_PATH))
message('Cells: ', ncol(seurat_ref), '   Genes: ', nrow(seurat_ref))
message('group levels: ', paste(unique(seurat_ref$group), collapse = ', '))

## Build module → integer mapping (WGCNA colour names → M-numbers)
mapping <- WGCNA::labels2colors(1:100)
bulk_modules <- read.csv('./dependencies/shared/bulk_heart_modules.csv')
bulk_modules$module <- match(bulk_modules$module, mapping)

Idents(seurat_ref) <- 'group'

dbs <- c('ChEA_2022', 'WikiPathway_2023_Human', 'Reactome_2016',
         'GO_Biological_Process_2023')

parse_ratio <- function(ratio) {
  ratio <- sub("^\\s*", "", as.character(ratio))
  ratio <- sub("\\s*$", "", ratio)
  numerator   <- as.numeric(sub("/\\d+$", "", ratio))
  denominator <- as.numeric(sub("^\\d+/", "", ratio))
  numerator / denominator
}
wrapText <- function(x, len) {
  sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"),
         USE.NAMES = FALSE)
}

run_contrast <- function(ident_a, ident_b, tag) {
  message(sprintf('\n=== %s (%s vs %s) ===', tag, ident_a, ident_b))

  ## FindMarkers per-module
  combined_set <- data.frame()
  mods_idx <- c(2, 12, 28, 10, 25, 26)
  for (i in mods_idx) {
    key_genes <- subset(bulk_modules, module %in% c(i))$gene_name
    key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
    if (length(key_genes) < 3) next
    gene_set <- FindMarkers(seurat_ref, ident.1 = ident_a, ident.2 = ident_b,
                             features = key_genes)
    gene_set <- subset(gene_set, p_val_adj < 0.05)
    if (nrow(gene_set) == 0) next
    gene_set$module <- paste0('M', i)
    gene_set$color  <- mapping[i]
    combined_set <- if (nrow(combined_set) == 0) gene_set
                    else rbind(combined_set, gene_set)
  }

  M2_genes_up <- rownames(subset(combined_set, module == 'M2' & avg_log2FC > 0))
  message('  M2_genes_up: n=', length(M2_genes_up), ';  first 10: ',
          paste(head(M2_genes_up, 10), collapse = ', '))

  ## enrichR
  enriched <- enrichr(M2_genes_up, dbs)
  go <- enriched[[4]]   # GO_Biological_Process_2023
  go_sig <- subset(go, Adjusted.P.value < 0.05)
  message('  GO_BP rows total: ', nrow(go), '   padj<0.05: ', nrow(go_sig))

  ## Top 5 by Combined Score (after filter), reversed for plot order
  go_sig <- go_sig[order(go_sig$Combined.Score, decreasing = TRUE), ]
  cat('\n  Top 5 (padj<0.05, sorted by Combined.Score):\n')
  print(head(go_sig[, c('Term','Adjusted.P.value','Combined.Score','Overlap')], 5))

  ## Plot exactly like F5
  if (nrow(go_sig) > 0) {
    top5 <- head(go_sig, 5)[rev(seq_len(min(5, nrow(go_sig)))), ]
    p <- ggplot(top5,
                aes(x = Combined.Score,
                    y = fct_inorder(Term),
                    color = as.numeric(Adjusted.P.value),
                    size  = parse_ratio(Overlap))) +
      geom_point() + xlab('Combined Score') + ylab('Term') +
      labs(color = 'P value', size = 'Overlap') +
      theme_classic() + ggtitle(paste0('M2 Up - ', tag)) +
      scale_y_discrete(labels = fct_inorder(
        wrapText(sapply(strsplit(top5$Term, ' \\(GO'), `[`, 1), 35))) +
      theme(axis.text = element_text(colour = 'black')) +
      scale_color_stepsn(colors = rev(magma(256)))
    out_pdf <- file.path(OUT, paste0('TEST_M2_enrichr_up_', gsub(' ', '_', tag), '.pdf'))
    pdf(out_pdf, width = 5, height = 2.5); print(p); dev.off()
    message('  wrote ', out_pdf)
  } else {
    message('  no significant terms — no PDF written')
  }

  invisible(list(M2_genes_up = M2_genes_up, top5 = head(go_sig, 5)))
}

res_NF  <- run_contrast('RVF', 'NF',  'RVF vs NF')
res_pRV <- run_contrast('RVF', 'pRV', 'RVF vs pRV')

cat('\n\n=== SUMMARY ===\n')
cat('Test PDFs in ', OUT, '\n')
cat('RVF vs NF  top 5 terms:\n'); print(res_NF$top5$Term)
cat('RVF vs pRV top 5 terms:\n'); print(res_pRV$top5$Term)
