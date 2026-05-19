##############################################################################
# Generate_Supplementary_Tables.R
#
# Builds Tables S1-S5 (clinical / experimental metadata) per v57 dictionary,
# emitted as a SINGLE workbook (per Nature CVR convention: one combined
# Supplementary_Tables.xlsx with one sheet per Table + a leading README).
#
#   S1  Adult RV cohort demographics + hemodynamics + platform allocation
#   S2  Hemodynamic / echo parameters used for disease-state stratification
#   S3  snRNA-seq cohort demographics (n=11, balanced)
#   S4  Murine PAB snRNA-seq cohort echo + group (n=10)
#   S5  Per-group n for adult / PAB / peds respirometry
#
# Sources
#   dependencies/shared/BulkRNA/metadata.csv     (adult bulk: var_disease, HHT_*,
#                                                  var_RAP/PCWP/RAP_PCWP/PAS/D/M,
#                                                  var_CI, var_PVRI, var_wt/ht,
#                                                  var_BSA/BMI, var_thyroid,
#                                                  var_pacer, prot_sn, var_sc_fin,
#                                                  RIN_score, var_batch)
#   dependencies/shared/PAB_data_clean.rds       (PAB snRNA mouse IDs + group)
#   dependencies/shared/snPAB_Echo.xlsx          (PAB echo measurements)
#   dependencies/shared/mito_data.xlsx           (respirometry; per-group n)
#
# Output
#   output/Supp_Tables/Supplementary_Tables.xlsx
##############################################################################

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# Render NA values as literal "NA" in xlsx cells (not as blank).
options(openxlsx.na.string = 'NA')

.shared <- '/Users/ikuz/Documents/RV_Atlas/dependencies/shared'
.outdir <- '/Users/ikuz/Documents/RV_Atlas/output/Supp_Tables'
dir.create(.outdir, showWarnings = FALSE, recursive = TRUE)

.section <- function(n, title) message(sprintf('\n============ Table S%d - %s ============', n, title))

# Curated echo / metric leaf labels (abbreviations + units; shared by the
# murine PAB tables S4 & S5). Keys cover both the snRNA-cache column names
# and the wk2_/wk8_-stripped PAB+respir column names.
.ECHO_ABBR <- c(
  Body_Weight_g              = 'BW (g)',
  Surgical_Body_Weight_g_    = 'BW (g)',
  RV_Area_Diastole_mm2       = 'RVAD (mm²)',
  RV_Area_Diastole_mm_2_     = 'RVAD (mm²)',
  RV_Area_Systole_mm2        = 'RVAS (mm²)',
  RV_Area_Systole_mm_2_      = 'RVAS (mm²)',
  FAC                        = 'FAC (%)',
  RV_Free_Wall_mm            = 'RVFW (mm)',
  RV_Free_Wall_Thickness_mm_ = 'RVFW (mm)',
  TAPSE_mm                   = 'TAPSE (mm)',
  TAPSE_mm_                  = 'TAPSE (mm)',
  TDI_S_Wave                 = 'TDI S′',
  TDI_S_Wave_mm_s_           = 'TDI S′',
  Pk_Gradient_mmHg           = 'Pk Grad (mmHg)',
  Pk_Velocity_mm_sec         = 'Pk Vel (mm/s)',
  Tricuspid_Regurgitation    = 'TR',
  TV_E_wave_mm_s_            = 'TV E (mm/s)',
  TDI_e_mm_s_                = 'TDI e′ (mm/s)',
  PA_VTI                     = 'PA VTI',
  Group                      = 'Group',
  Cohort                     = 'Cohort',
  Echo_Date                  = 'Echo Date',
  Echo_Weight_g_             = 'Echo Weight (g)',
  Sx_Date                    = 'Surgery Date'
)
# Registry of sheets that need a 2-row (grouped) header. Keyed by sheet
# name; value = list(leaf = <display leaf labels>, group = <top-row span
# label per column, '' = no group). Read by the xlsx writer and .render_pdf.
.header2 <- list()
# Registry of sheets rendered as several vertically-stacked sub-tables in
# one sheet/section. Keyed by sheet name; value = list of parts, each
# list(title=, cols=<original col names>, labels=<display headers>). The
# parent sheet's data.frame (post .round_sf3) is the data source.
.stacked <- list()

# Canonical snRNA-seq cohort: 11 patient IDs (from snRV_ref.rds$patient unique).
SN_PATIENT_IDS <- c('1343','1392','1467','1561','1567','1618',
                    '1632','1681','1691','1692','1697')

# Canonical Xenium cohort: 9 patient IDs (from Xenium_metadata.csv$patient unique).
XENIUM_PATIENT_IDS <- c('1343','1467','1561','1567','1618',
                        '1632','1691','1692','1697')

# Collect all tables in a single named list; we'll write one xlsx at the end.
all_sheets <- list()

# ===========================================================================
# Adult bulk metadata loader + tidy
# ===========================================================================
.load_adult_meta <- function() {
  raw <- read.csv(file.path(.shared, 'BulkRNA/metadata.csv'),
                  check.names = TRUE, fileEncoding = 'UTF-8-BOM',
                  stringsAsFactors = FALSE)
  raw <- raw[seq(2, nrow(raw), 2), , drop = FALSE]   # drop empty padding rows
  rownames(raw) <- NULL
  raw
}

.tidy_adult <- function(df) {
  rmap <- c(
    title              = 'patient_id',
    maineffect         = 'group',       # disease state: NF / pRV / RVF
    var_disease        = 'etiology',    # NF / DCM / ICM
    HHT_sex            = 'sex',
    HHT_age            = 'age',
    HHT_eth            = 'ethnicity',
    var_RAP            = 'RAP_mmHg',
    var_PAS            = 'PA_systolic_mmHg',
    var_PAD            = 'PA_diastolic_mmHg',
    var_PAM            = 'PA_mean_mmHg',
    var_PCWP           = 'PCWP_mmHg',
    var_CI             = 'cardiac_index',
    var_RAP_PCWP       = 'RAP_PCWP_ratio',
    var_PVRI           = 'PVRI',
    var_wt             = 'weight_kg',
    var_ht             = 'height_cm',
    var_BSA            = 'BSA',
    var_BMI            = 'BMI',
    var_pacer          = 'pacemaker',
    RIN_score          = 'RIN'
  )
  keep_raw <- intersect(names(rmap), names(df))
  out <- df[, keep_raw, drop = FALSE]
  names(out) <- rmap[keep_raw]
  out$in_snRNA_cohort  <- as.character(out$patient_id) %in% SN_PATIENT_IDS
  out$in_Xenium_cohort <- as.character(out$patient_id) %in% XENIUM_PATIENT_IDS

  # Column ordering
  ord <- c('patient_id','group','etiology',
           'sex','age','ethnicity',
           'RAP_mmHg','RAP_PCWP_ratio',
           'PCWP_mmHg','PA_systolic_mmHg','PA_diastolic_mmHg','PA_mean_mmHg',
           'cardiac_index','PVRI',
           'weight_kg','height_cm','BSA','BMI',
           'pacemaker',
           'in_snRNA_cohort','in_Xenium_cohort','RIN')
  out <- out[, intersect(ord, names(out)), drop = FALSE]
  # Round RAP:PCWP ratio to 3 significant figures
  if ('RAP_PCWP_ratio' %in% names(out))
    out$RAP_PCWP_ratio <- signif(suppressWarnings(as.numeric(out$RAP_PCWP_ratio)), 3)

  # Order: NF -> pRV -> RVF, then by patient_id
  if ('group' %in% names(out)) {
    out$group <- factor(out$group,
                        levels = c('NF','pRV','RVF',
                                   setdiff(unique(out$group), c('NF','pRV','RVF'))))
    out <- out[order(out$group, out$patient_id), , drop = FALSE]
    out$group <- as.character(out$group)
    rownames(out) <- NULL
  }
  # Empty strings -> NA
  for (c in names(out)) {
    if (is.character(out[[c]])) out[[c]][!nzchar(out[[c]])] <- NA
  }
  # Drop columns that are entirely NA AND not in the v57 must-keep list
  must_keep <- c('RAP_mmHg','RAP_PCWP_ratio',
                 'underlying_diagnosis',
                 'in_snRNA_cohort','in_Xenium_cohort')
  all_na <- vapply(out, function(x) all(is.na(x)) || all(!nzchar(as.character(x)) | is.na(x)),
                   logical(1))
  drop_cols <- names(out)[all_na & !names(out) %in% must_keep]
  if (length(drop_cols) > 0) {
    message('  dropping all-NA columns: ', paste(drop_cols, collapse = ', '))
    out <- out[, setdiff(names(out), drop_cols), drop = FALSE]
  }
  out
}

# ===========================================================================
# Table S1 — Adult RV cohort
# ===========================================================================
.section(1, 'Adult RV cohort demographics + platform allocation')
adult <- .load_adult_meta()
s1_full <- .tidy_adult(adult)
# Split into S1 (demographics + platform allocation) and S2 (hemodynamics).
.HEMO_COLS <- c('RAP_mmHg','RAP_PCWP_ratio','PCWP_mmHg',
                'PA_systolic_mmHg','PA_diastolic_mmHg','PA_mean_mmHg',
                'cardiac_index','PVRI')
s1 <- s1_full[, !names(s1_full) %in% .HEMO_COLS, drop = FALSE]
all_sheets[['Table_S1']] <- s1

.section(2, 'Adult RV cohort hemodynamics')
s2 <- s1_full[, c('patient_id', 'group',
                  intersect(.HEMO_COLS, names(s1_full))),
              drop = FALSE]
all_sheets[['Table_S2']] <- s2

# ===========================================================================
# Table S3 — snRNA-seq cohort (n=11) demographics + QC
# Demographics joined from S1; QC summarized per-patient from snRV_ref.rds
# (mean % mito, mean entropy, mean % intronic, median nCount / nFeature).
# ===========================================================================
.section(3, 'snRNA-seq cohort (n=11) demographics + QC')
sn_rds <- file.path(.shared, 'snRV_ref.rds')
if (file.exists(sn_rds)) {
  sn <- readRDS(sn_rds)
  sm <- sn@meta.data
  rm(sn); invisible(gc(verbose = FALSE))
  if (!'patient' %in% names(sm)) {
    message('  snRV_ref.rds metadata lacks patient column; skipping Table S3')
  } else {
    sm$patient <- as.character(sm$patient)
    # Normalise percentage columns to 0-100 scale before aggregating; many
    # snRNA-seq objects store percent.* values as 0-1 fractions, which
    # produced spurious >99% intronic values before this fix.
    .as_pct <- function(v) {
      v <- suppressWarnings(as.numeric(v))
      if (suppressWarnings(max(v, na.rm = TRUE)) <= 1) v <- v * 100
      v
    }
    has_exon <- 'percent.exon'     %in% names(sm)
    has_intr <- 'percent.intronic' %in% names(sm)
    if ('percent.mt' %in% names(sm)) sm$percent_mt_pct <- .as_pct(sm$percent.mt)
    if (has_intr) {
      sm$percent_intronic_pct <- .as_pct(sm$percent.intronic)
    } else if (has_exon) {
      exon_pct <- .as_pct(sm$percent.exon)
      sm$percent_intronic_pct <- 100 - exon_pct
    }
    qc <- sm %>%
      dplyr::group_by(patient_id = patient) %>%
      dplyr::summarise(
        n_nuclei              = dplyr::n(),
        median_nCount_RNA     = stats::median(nCount_RNA,   na.rm = TRUE),
        median_nFeature_RNA   = stats::median(nFeature_RNA, na.rm = TRUE),
        mean_percent_mt       = if ('percent_mt_pct' %in% names(sm))
                                  mean(percent_mt_pct, na.rm = TRUE) else NA_real_,
        mean_entropy          = if ('entropy' %in% names(sm))
                                  mean(entropy, na.rm = TRUE) else NA_real_,
        mean_percent_intronic = if ('percent_intronic_pct' %in% names(sm))
                                  mean(percent_intronic_pct, na.rm = TRUE) else NA_real_,
        .groups = 'drop'
      )
    demo <- s1[s1$in_snRNA_cohort == TRUE,
               intersect(c('patient_id','group','sex','age',
                           'ethnicity','underlying_diagnosis',
                           'RAP_mmHg','RAP_PCWP_ratio'),
                         names(s1)),
               drop = FALSE]
    demo$patient_id <- as.character(demo$patient_id)
    s3 <- merge(demo, qc, by = 'patient_id', all = TRUE)
    s3 <- s3[order(s3$group, s3$patient_id), , drop = FALSE]
    rownames(s3) <- NULL
    all_sheets[['Table_S3']] <- s3
    message('  snRNA cohort rows: ', nrow(s3),
            if (nrow(s3) != 11) '  (NOTE: v57 expects n=11)' else '')
  }
  rm(sm); invisible(gc(verbose = FALSE))
} else {
  message('  snRV_ref.rds not found; skipping Table S3')
}

# ===========================================================================
# Table S4 — Murine PAB snRNA-seq cohort
# ===========================================================================
.section(4, 'Murine PAB snRNA-seq cohort (n=10)')
## NOTE: the snRNA-seq mice and the echo/respirometry mice are DISTINCT
## animals that merely reuse ID integers (per the wet-lab). They are kept
## as two SEPARATE tables: Table S4 = the n=10 sequenced cohort (this
## block); Table S5 = the 2-wk + 2-mo echo / respirometry cohort (below).
## Reads the small precomputed cohort summary (deps-compliant; avoids
## loading the 55k-cell PAB_data_clean.rds Seurat object at table-gen
## time — that readRDS peaks memory and OOMs on constrained machines).
## Regenerate the cache from PAB_data_clean.rds + 'PAB snRNAseq
## samples.xlsx' if the upstream snRNA object changes.
sn_cache <- file.path(.shared, 'PAB_snRNA_cohort_summary.csv')
if (file.exists(sn_cache)) {
  s4 <- read.csv(sn_cache, stringsAsFactors = FALSE, check.names = FALSE)
  s4 <- s4[order(factor(s4$group, c('Sham','Moderate','Severe')), s4$ID), ,
           drop = FALSE]
  s4$HR <- NULL                                    # drop HR per spec
  .g4 <- table(s4$group)                           # tally before rename
  ## Final abbreviated headers (units; ² superscript, mm/s parenthesised).
  s4_lab <- c(ID = 'ID', group = 'Group', n_nuclei = 'Nuclei', .ECHO_ABBR)
  names(s4) <- ifelse(names(s4) %in% names(s4_lab),
                      s4_lab[names(s4)], names(s4))
  rownames(s4) <- NULL
  all_sheets[['Table_S4']] <- s4
  message('  Table S4 (snRNA cohort): ', nrow(s4), ' mice | cols: ',
          ncol(s4), ' | ',
          paste(names(.g4), .g4, sep = '=', collapse = ', '))
} else {
  message('  PAB_snRNA_cohort_summary.csv not found; skipping Table S4')
}

# ===========================================================================
# Table S5 — Murine PAB 2-wk + 2-mo echo / respirometry cohort:
#   the full echo/respirometry animal set (DISTINCT mice from Table S4),
#   2-week + 8-week (2-month) echocardiography. Detailed respirometry flux
#   columns and admin date/weight columns are dropped per spec; sorted by
#   PAB_vs_Sham.
# ===========================================================================
.section(5, 'Murine PAB 2-wk + 2-mo echo / respirometry cohort')
pab_resp_xlsx <- file.path(.shared, 'PAB + respirometry data.xlsx')
if (file.exists(pab_resp_xlsx)) {
  raw <- as.data.frame(suppressMessages(readxl::read_excel(
           pab_resp_xlsx, sheet = 1, col_names = FALSE, .name_repair = 'minimal')))
  hdr <- as.character(unlist(raw[2, ]))
  dat <- raw[-(1:2), , drop = FALSE]; rownames(dat) <- NULL
  nm <- hdr
  nm[1:21]  <- ifelse(seq_len(21) %in% 1:6, hdr[1:21], paste0('wk2_', hdr[1:21]))
  nm[22:40] <- paste0('wk8_', hdr[22:40])
  nm[41:49] <- paste0('resp_', ifelse(is.na(hdr[41:49]) | hdr[41:49] == '',
                                      paste0('c', 41:49), hdr[41:49]))
  nm <- make.unique(gsub('[^A-Za-z0-9]+', '_', trimws(nm)))
  colnames(dat) <- nm
  exp_col <- nm[1]; id_col <- nm[2]
  dat[[id_col]] <- suppressWarnings(as.integer(dat[[id_col]]))
  dat <- dat[!is.na(dat[[id_col]]), , drop = FALSE]
  grp_cands <- grep('PAB_vs_SHAM|^wk8_PAB|^resp_PAB', colnames(dat),
                    value = TRUE, ignore.case = TRUE)
  pab_vs <- if (length(grp_cands))
    apply(dat[, grp_cands, drop = FALSE], 1,
          function(v){ v <- v[!is.na(v) & v != '']; if (length(v)) v[1] else NA })
    else NA
  exp_set <- as.character(dat[[exp_col]])           # sort key only (dropped)
  s5 <- data.frame(
    ID          = as.integer(dat[[id_col]]),
    PAB_vs_Sham = pab_vs,
    dat[, setdiff(colnames(dat), c(exp_col, id_col)), drop = FALSE],
    check.names = FALSE, stringsAsFactors = FALSE)
  ## Drop per spec: Exp Set (sort-only), HR, detailed respirometry flux
  ## cols, admin date/weight cols, pure-duplicate key cols.
  drop5 <- c('wk2_HR', 'wk8_HR',
             'wk2_TV_E_wave_mm_s_', 'wk8_TV_E_wave_mm_s_',
             'wk2_TDI_e_mm_s_', 'wk8_TDI_e_mm_s_',
             'wk2_PA_VTI', 'wk8_PA_VTI',
             'Echo_Date', 'Echo_Weight_g_', 'Sx_Date',
             'resp_Batch_ID', 'resp_PM', 'resp_ADP', 'resp_Cyt_C',
             'resp_OC', 'resp_Succ', 'resp_Cyt_C_', 'resp_RCR',
             'resp_PAB_vs_SHAM', 'wk8_ID', 'wk8_PAB_vs_SHAM')
  s5 <- s5[, setdiff(colnames(s5), drop5), drop = FALSE]
  ord <- order(s5$PAB_vs_Sham, exp_set, s5$ID)     # exp_set: sort key only
  s5 <- s5[ord, , drop = FALSE]
  rownames(s5) <- NULL
  ## "Dead 12/4/25" (and any "Dead <date>") -> "Died".
  for (cn in names(s5)) if (is.character(s5[[cn]]))
    s5[[cn]] <- sub('^\\s*Dead\\b.*$', 'Died', s5[[cn]], ignore.case = TRUE)
  ## readxl returns the metric columns as character (header sat in row 2);
  ## coerce numeric-looking cols to numeric so .round_sf3 caps echo
  ## parameters (RVAD/RVAS/FAC/RVFW/TAPSE/...) at 3 significant figures.
  for (cn in names(s5)) {
    x <- s5[[cn]]
    if (!is.character(x)) next
    real <- !is.na(x) & nzchar(trimws(x)) &
            !(tolower(trimws(x)) %in% c('na', 'n/a'))
    v <- suppressWarnings(as.numeric(x))
    if (any(real) && all(!is.na(v[real]))) s5[[cn]] <- v
  }
  ## Split into two vertically-stacked sub-tables (one Table S5 heading):
  ## Week 2 then Week 8. Each repeats ID / Category / BW (g) and that
  ## week's echo params; the Week 8 BW is the wk8 surgical body weight.
  echo_base <- c('RV_Area_Diastole_mm_2_', 'RV_Area_Systole_mm_2_', 'FAC',
                 'RV_Free_Wall_Thickness_mm_', 'TAPSE_mm_',
                 'TDI_S_Wave_mm_s_', 'Pk_Gradient_mmHg',
                 'Pk_Velocity_mm_sec', 'Tricuspid_Regurgitation')
  echo_lab <- unname(.ECHO_ABBR[echo_base])
  lead_lab <- c('ID', 'Category', 'BW (g)')
  ## Force the known NUMERIC echo metrics to numeric and cap at 3 sig
  ## figs here: the generic all-or-nothing coercion above leaves a column
  ## as character if ANY cell is non-numeric (e.g. a stray '-'/'ND'),
  ## which then bypasses .round_sf3 and prints full precision. TR is
  ## text (grades) so it is excluded.
  num_base <- setdiff(echo_base, 'Tricuspid_Regurgitation')
  num_cols <- c('Surgical_Body_Weight_g_', 'wk8_Surgical_Body_Weight_g_',
                paste0('wk2_', num_base), paste0('wk8_', num_base))
  for (cn in intersect(num_cols, names(s5))) {
    s5[[cn]] <- signif(suppressWarnings(as.numeric(
      gsub('[^0-9.eE+-]', '', trimws(as.character(s5[[cn]]))))), 3)
  }
  ## Cohort group: a mouse is "Week 8" if it has any 8-week echo datum,
  ## else "Week 2 only" (respirometry-only mice with no 8-week echo).
  .wk8_metric <- intersect(paste0('wk8_', num_base), names(s5))
  s5$wk_group <- ifelse(
    rowSums(!is.na(s5[, .wk8_metric, drop = FALSE])) > 0,
    'Week 8', 'Week 2 only')
  ## nlead = leading non-echo cols; rows whose echo columns are ALL NA
  ## (respirometry-only mice with no 8-week echo) are dropped per part.
  ## Week 2 carries the Group column + ID; Week 8 drops ID.
  .stacked[['Table_S5']] <- list(
    list(title = 'Week 2', nlead = 4L,
         cols = c('ID', 'PAB_vs_Sham', 'wk_group',
                  'Surgical_Body_Weight_g_', paste0('wk2_', echo_base)),
         labels = c('ID', 'Category', 'Group', 'BW (g)', echo_lab)),
    list(title = 'Week 8', nlead = 3L,
         cols = c('ID', 'PAB_vs_Sham', 'wk8_Surgical_Body_Weight_g_',
                  paste0('wk8_', echo_base)),
         labels = c('ID', 'Category', 'BW (g)', echo_lab)))
  all_sheets[['Table_S5']] <- s5
  message('  Table S5 (echo/respir cohort): ', nrow(s5), ' mice | cols: ',
          ncol(s5), ' | ',
          paste(names(table(s5$PAB_vs_Sham)), table(s5$PAB_vs_Sham),
                sep = '=', collapse = ', '))
} else {
  message('  PAB + respirometry data.xlsx not found; skipping Table S5')
}

# ===========================================================================
# Table S6 — Xenium spatial transcriptomics cohort (n=9):
#   per-patient demographics (joined from S1), tissue section details, QC.
# ===========================================================================
.section(6, 'Xenium spatial transcriptomics cohort (n=9)')
xen_meta_path <- file.path(.shared, 'Xenium_metadata.csv')
if (file.exists(xen_meta_path)) {
  xm <- read.csv(xen_meta_path, stringsAsFactors = FALSE)
  # Per-patient QC + tissue section summary
  xm$patient <- as.character(xm$patient)
  qc <- xm %>%
    dplyr::group_by(patient_id = patient) %>%
    dplyr::summarise(
      n_cells           = dplyr::n(),
      median_nCount     = stats::median(nCount_Xenium,   na.rm = TRUE),
      median_nFeature   = stats::median(nFeature_Xenium, na.rm = TRUE),
      median_cell_area  = stats::median(cell_area,       na.rm = TRUE),
      x_min             = min(x_centroid, na.rm = TRUE),
      x_max             = max(x_centroid, na.rm = TRUE),
      y_min             = min(y_centroid, na.rm = TRUE),
      y_max             = max(y_centroid, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::mutate(
      # x/y_centroid are in micrometers (10x Xenium standard output);
      # convert bounding-box area to mm^2.
      section_area_mm2 = ((x_max - x_min) * (y_max - y_min)) / 1e6,
      cells_per_mm2    = n_cells / section_area_mm2
    ) %>%
    dplyr::select(-x_min, -x_max, -y_min, -y_max)

  # Demographics from S1 (cohort already filtered by in_Xenium_cohort)
  demo <- s1[s1$in_Xenium_cohort == TRUE,
             intersect(c('patient_id','group','sex','age',
                         'ethnicity','underlying_diagnosis',
                         'RAP_mmHg','RAP_PCWP_ratio'),
                       names(s1)),
             drop = FALSE]
  demo$patient_id <- as.character(demo$patient_id)
  s6 <- merge(demo, qc, by = 'patient_id', all = TRUE)
  s6 <- s6[order(s6$group, s6$patient_id), , drop = FALSE]
  rownames(s6) <- NULL
  all_sheets[['Table_S6']] <- s6
  message('  Xenium cohort rows: ', nrow(s6),
          if (nrow(s6) != 9) '  (NOTE: v57 expects n=9)' else '')
} else {
  message('  Xenium_metadata.csv not found; skipping Table S6')
}

# ===========================================================================
# README sheet (leading) — one-line description per Table
# ===========================================================================
index_df <- data.frame(
  `Table #` = c('Table S1', 'Table S2', 'Table S3', 'Table S4',
                'Table S5', 'Table S6'),
  Title = c(
    'Adult RV cohort: demographics & platform allocation (n=142).',
    'Adult RV cohort: hemodynamics (n=142).',
    'snRNA-seq cohort: demographics & QC (n=11).',
    'Murine PAB snRNA-seq cohort: echo (n=10).',
    'Murine PAB 2-wk + 2-mo echo / respirometry cohort.',
    'Xenium spatial cohort: demographics & QC (n=9).'
  ),
  check.names = FALSE
)
all_sheets <- c(list(INDEX = index_df), all_sheets)

# Columns that should be rendered as whole integers (no sig-fig cap)
.INT_COLS <- c(
  'patient_id', 'mouse_id', 'sample_num',
  'n_nuclei', 'n_cells',
  'median_nCount_RNA', 'median_nFeature_RNA',
  'median_nCount', 'median_nFeature',
  'median_cell_area', 'cells_per_mm2',
  'Pk_Velocity_mm/sec'
)

# Round numeric columns: integers stay whole, others capped at 3 sig figs.
.round_sf3 <- function(df) {
  for (c in names(df)) {
    if (!is.numeric(df[[c]])) next
    if (c %in% .INT_COLS) {
      df[[c]] <- as.integer(round(df[[c]]))
    } else if (!is.integer(df[[c]])) {
      df[[c]] <- signif(df[[c]], 3)
    }
  }
  df
}
all_sheets <- lapply(all_sheets, .round_sf3)

out_path <- file.path(.outdir, 'Supplementary_Tables.xlsx')
{
  wb  <- openxlsx::createWorkbook()
  hs1 <- openxlsx::createStyle(textDecoration = 'bold', halign = 'center',
                               valign = 'center', border = 'bottom')
  for (nm in names(all_sheets)) {
    df <- all_sheets[[nm]]
    openxlsx::addWorksheet(wb, nm)
    h2 <- .header2[[nm]]
    st <- .stacked[[nm]]
    if (!is.null(st)) {
      ## Vertically-stacked sub-tables: bold title row, header row, data,
      ## blank-row separator; repeat for each part.
      cur <- 1L
      for (p in st) {
        sub <- df[, p$cols, drop = FALSE]; names(sub) <- p$labels
        ec <- seq_len(ncol(sub))[-seq_len(p$nlead)]   # echo columns
        sub <- sub[rowSums(!is.na(sub[, ec, drop = FALSE])) > 0, ,
                   drop = FALSE]
        openxlsx::writeData(wb, nm, p$title, startRow = cur, startCol = 1)
        openxlsx::addStyle(wb, nm,
          openxlsx::createStyle(textDecoration = 'bold'),
          rows = cur, cols = 1)
        openxlsx::writeData(wb, nm, sub, startRow = cur + 1L,
                            colNames = TRUE)
        openxlsx::addStyle(wb, nm, hs1, rows = cur + 1L,
                           cols = seq_along(p$cols), gridExpand = TRUE)
        cur <- cur + 1L + 1L + nrow(sub) + 2L
      }
    } else if (!is.null(h2) && length(h2$leaf) == ncol(df)) {
      ## Row 1 = merged group spans; row 2 = leaf headers; data from row 3.
      dd <- df; names(dd) <- h2$leaf
      openxlsx::writeData(wb, nm, dd, startRow = 2, colNames = TRUE)
      g <- h2$group; r <- rle(g); end <- cumsum(r$lengths)
      start <- end - r$lengths + 1
      for (k in seq_along(r$values)) {
        if (nzchar(r$values[k])) {
          openxlsx::mergeCells(wb, nm, cols = start[k]:end[k], rows = 1)
          openxlsx::writeData(wb, nm, r$values[k], startRow = 1,
                              startCol = start[k])
        }
      }
      openxlsx::addStyle(wb, nm, hs1, rows = 1,
                         cols = seq_len(ncol(df)), gridExpand = TRUE)
      openxlsx::freezePane(wb, nm, firstActiveRow = 3)
    } else {
      openxlsx::writeData(wb, nm, df, colNames = TRUE)
      openxlsx::freezePane(wb, nm, firstActiveRow = 2)
    }
  }
  openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)
}
message('\nWROTE ', out_path, '  (', length(all_sheets), ' sheet(s))')

# ===========================================================================
# Formatted PDF rendering (for submission)
# ===========================================================================
# Pretty-print column names: replace _ with spaces, parenthesize known unit
# suffixes, and use curated labels for common fields.
.PRETTY_COL_MAP <- c(
  sample_num             = 'Sample #',
  patient_id             = 'Patient ID',
  group                  = 'Group',
  etiology               = 'Etiology',
  sex                    = 'Sex',
  age                    = 'Age (yrs)',
  ethnicity              = 'Ethnicity',
  underlying_diagnosis   = 'Underlying Diagnosis',
  high_impact_variants_WES = 'High-impact Variants (WES)',
  RAP_mmHg               = 'RAP (mmHg)',
  RAP_PCWP_ratio         = 'RAP:PCWP',
  TAPSE_cm               = 'TAPSE (cm)',
  TR_grade               = 'TR Grade',
  PCWP_mmHg              = 'PCWP (mmHg)',
  PA_systolic_mmHg       = 'PA Systolic (mmHg)',
  PA_diastolic_mmHg      = 'PA Diastolic (mmHg)',
  PA_mean_mmHg           = 'PA Mean (mmHg)',
  cardiac_index          = 'Cardiac Index (L/min/m²)',
  PVRI                   = 'PVRI',
  weight_kg              = 'Weight (kg)',
  height_cm              = 'Height (cm)',
  BSA                    = 'BSA (m²)',
  BMI                    = 'BMI (kg/m²)',
  thyroid_status         = 'Thyroid Status',
  pacemaker              = 'Pacemaker',
  exogenous_glucocorticoid = 'Exogenous Glucocorticoid',
  in_snRNA_cohort        = 'In snRNA Cohort',
  in_Xenium_cohort       = 'In Xenium Cohort',
  RIN                    = 'RIN',
  mouse_id               = 'Mouse ID',
  Cage                   = 'Cage',
  PAB_vs_SHAM            = 'PAB vs Sham',
  Sx_Date                = 'Surgery Date',
  Date_of_Echo           = 'Echo Date',
  Echo_Body_Weight       = 'Echo Body Weight (g)',
  Surgical_Body_Weight   = 'Surgical Body Weight (g)',
  HR                     = 'HR (bpm)',
  RV_Area_Diastole       = 'RVAD (mm²)',
  RV_Area_Systole        = 'RVAS (mm²)',
  FAC                    = 'FAC (%)',
  RV_Free_Wall_Thickness_mm = 'RVFWT (mm)',
  TAPSE_mm               = 'TAPSE (mm)',
  TDI_S_Wave             = 'TDI S\' (mm/s)',
  Pk_Gradient_mmHg       = 'Pk Gradient (mmHg)',
  `Pk_Velocity_mm/sec`   = 'Pk Velocity (mm/s)',
  Tricuspid_Regurgitation = 'Tricuspid Regurgitation',
  group_xenium           = 'Group (Xenium)',
  n_cells                = 'Cells',
  median_nCount          = 'Median Count',
  median_nFeature        = 'Median Features',
  median_cell_area       = 'Median Cell Area (µm²)',
  section_area_mm2       = 'Section Area (mm²)',
  cells_per_mm2          = 'Cell Density (cells/mm²)',
  n_nuclei               = 'Nuclei',
  median_nCount_RNA      = 'Median Count',
  median_nFeature_RNA    = 'Median Features',
  mean_percent_mt        = '% Mitochondrial',
  mean_entropy           = 'Entropy',
  mean_percent_intronic  = '% Intronic'
)
.pretty_cols <- function(df) {
  new <- vapply(names(df), function(n) {
    if (n %in% names(.PRETTY_COL_MAP)) .PRETTY_COL_MAP[[n]]
    else gsub('_', ' ', n)
  }, character(1))
  names(df) <- new
  df
}

.render_pdf <- function(sheets, pdf_path) {
  if (!requireNamespace('rmarkdown', quietly = TRUE) ||
      !requireNamespace('kableExtra', quietly = TRUE) ||
      !requireNamespace('knitr', quietly = TRUE)) {
    message('  rmarkdown / kableExtra / knitr not installed — skipping PDF')
    return(invisible(FALSE))
  }
  if (requireNamespace('tinytex', quietly = TRUE) && !tinytex::is_tinytex()) {
    message('  tinytex not detected — attempting xelatex via system')
  }
  pretty_sheets <- lapply(sheets, .pretty_cols)
  rds_path <- file.path(tempdir(), 'supp_tables_sheets.rds')
  saveRDS(list(pretty = pretty_sheets, raw = sheets), rds_path)

  rmd_lines <- c(
    '---',
    'title: "Supplementary Tables"',
    'output:',
    '  pdf_document:',
    '    latex_engine: xelatex',
    'geometry: paperwidth=8.5in, paperheight=11in, margin=1in',
    'mainfont: Arial',
    'header-includes:',
    '  - \\usepackage{longtable}',
    '  - \\usepackage{booktabs}',
    '  - \\usepackage{caption}',
    '  - \\usepackage{array}',
    '  - \\setlength{\\tabcolsep}{2pt}',
    '  - \\renewcommand{\\arraystretch}{0.85}',
    '---',
    '',
    paste0('```{r setup, include=FALSE}'),
    'library(knitr); library(kableExtra)',
    paste0('.obj <- readRDS("', rds_path,
           '"); sheets <- .obj$pretty; raw <- .obj$raw'),
    'knitr::opts_chunk$set(echo = FALSE, results = "asis",',
    '                      warning = FALSE, message = FALSE)',
    '```',
    ''
  )
  # Table captions (updated from v57 to reflect actual column composition)
  v57_captions <- list(
    Table_S1 = 'Table S1. Patient demographics for the adult RV cohort (n=142), including disease state classification (NF, pRV, RVF), etiology (DCM/ICM/NF), age, sex, ethnicity, anthropometrics, pacemaker status, RNA integrity (RIN), and platform allocation flags (in_snRNA_cohort, in_Xenium_cohort).',
    Table_S2 = 'Table S2. Hemodynamic parameters for the adult RV cohort (n=142) used for disease-state stratification: right atrial pressure (RAP), pulmonary capillary wedge pressure (PCWP), RA:PCWP ratio (primary stratification variable), pulmonary arterial pressures, cardiac index, and pulmonary vascular resistance index (PVRI).',
    Table_S3 = 'Table S3. snRNA-seq cohort (n=11): per-patient demographics joined from Table S1, plus per-nucleus quality control metrics — nuclei recovered, median UMI count and feature count, mean mitochondrial percent, mean transcriptional entropy, and mean intronic read percent.',
    Table_S4 = 'Table S4. Murine PAB snRNA-seq cohort (n=10: 3 Sham, 3 Moderate-RVF, 4 Severe-RVF), with nuclei recovered and per-mouse echocardiographic parameters. Column abbreviations: BW, body weight (g); RVAD, RV area at diastole (mm²); RVAS, RV area at systole (mm²); FAC, fractional area change (%); RVFW, RV free-wall thickness (mm); TAPSE, tricuspid annular plane systolic excursion (mm); TDI S′, tissue-Doppler systolic wave; Pk Grad, peak gradient (mmHg); Pk Vel, peak velocity (mm/s); TR, tricuspid regurgitation. These are DISTINCT animals from the echo/respirometry cohort (Table S5) despite shared ID integers.',
    Table_S5 = 'Table S5. Murine PAB 2-week + 2-month (8-week) echocardiography / respirometry cohort — the full echo-phenotyped animal set (distinct mice from the snRNA-seq cohort in Table S4), categorised as PAB vs Sham (Category) and sorted by Category. Presented as two vertically-stacked sub-tables: Week 2 and Week 8. ID, Category and BW (g) are repeated in each; the Week 8 BW (g) is the 8-week surgical body weight. The Week 2 Group column flags whether a mouse is "Week 2 only" (respirometry-only, no 8-week echo) or "Week 8" (also echoed at 8 weeks); the Week 8 sub-table lists exactly the latter. Echo abbreviations match Table S4 (BW, body weight (g); RVAD, RV area at diastole (mm²); RVAS, RV area at systole (mm²); FAC, fractional area change (%); RVFW, RV free-wall thickness (mm); TAPSE, tricuspid annular plane systolic excursion (mm); TDI S′; Pk Grad, peak gradient (mmHg); Pk Vel, peak velocity (mm/s); TR, tricuspid regurgitation); echo values are shown to 3 significant figures. Week 8 mouse ID 17 died during the echocardiogram, so only the parameters acquired before death are reported (subsequent metrics — TAPSE, TDI S′, Pk Grad, Pk Vel, TR — are shown as NA).',
    Table_S6 = 'Table S6. Patient demographics for the Xenium spatial transcriptomics cohort (n=9), with tissue section details and per-section quality control metrics (cells captured, median UMI/features per cell, median cell area, section area, cell density).'
  )

  for (nm in names(pretty_sheets)) {
    display_nm <- gsub('_', ' ', nm)
    fs <- if (nm == 'INDEX') 11 else 6
    h2 <- .header2[[nm]]
    st <- .stacked[[nm]]
    if (!is.null(st)) {
      body <- character(0)
      for (p in st) {
        cols_lit <- paste0('c(', paste(sprintf('"%s"', p$cols),
                                       collapse = ', '), ')')
        lab_lit  <- paste0('c(', paste(sprintf('"%s"', p$labels),
                                       collapse = ', '), ')')
        body <- c(body,
          paste0('\\textbf{', p$title, '}'), '',
          '```{r}',
          paste0('p <- raw[["', nm, '"]][, ', cols_lit, ', drop = FALSE]'),
          paste0('colnames(p) <- ', lab_lit),
          paste0('p <- p[rowSums(!is.na(p[, -(1:', p$nlead,
                 '), drop = FALSE])) > 0, , drop = FALSE]'),
          'rownames(p) <- NULL',
          'kbl(p, format = "latex", booktabs = TRUE, longtable = TRUE,',
          '    row.names = FALSE, na.string = "NA", linesep = "") %>%',
          '  kable_styling(latex_options = c("repeat_header","HOLD_position"),',
          paste0('                font_size = ', fs, ')'),
          '```', '')
      }
      rmd_lines <- c(rmd_lines, paste0('## ', display_nm), '', body)
    } else {
      if (!is.null(h2) && length(h2$leaf) == ncol(pretty_sheets[[nm]])) {
        leaf_lit <- paste0('c(', paste(sprintf('"%s"', h2$leaf),
                                       collapse = ', '), ')')
        r <- rle(h2$group)
        span_lit <- paste0('c(', paste(sprintf('"%s" = %d',
                              ifelse(nzchar(r$values), r$values, ' '),
                              r$lengths), collapse = ', '), ')')
        kbl_block <- c(
          paste0('df <- sheets[["', nm, '"]]'),
          paste0('colnames(df) <- ', leaf_lit),
          'kbl(df, format = "latex", booktabs = TRUE, longtable = TRUE,',
          '    na.string = "NA", linesep = "") %>%',
          paste0('  add_header_above(', span_lit, ') %>%'),
          '  kable_styling(latex_options = c("repeat_header","HOLD_position"),',
          paste0('                font_size = ', fs, ')'))
      } else {
        kbl_block <- c(
          paste0('df <- sheets[["', nm, '"]]'),
          'kbl(df, format = "latex", booktabs = TRUE, longtable = TRUE,',
          '    na.string = "NA", linesep = "") %>%',
          '  kable_styling(latex_options = c("repeat_header","HOLD_position"),',
          paste0('                font_size = ', fs, ')'))
      }
      rmd_lines <- c(rmd_lines,
        paste0('## ', display_nm), '',
        '```{r}', kbl_block, '```', ''
      )
    }
    if (!is.null(v57_captions[[nm]])) {
      rmd_lines <- c(rmd_lines,
        paste0('\\small\\textbf{', strsplit(v57_captions[[nm]], '\\. ')[[1]][1],
               '.}\\ ',
               sub('^[^.]+\\. ', '', v57_captions[[nm]])),
        ''
      )
    }
    rmd_lines <- c(rmd_lines, '\\newpage', '')
  }
  rmd_path <- file.path(tempdir(), '_supp_tables.Rmd')
  writeLines(rmd_lines, rmd_path)
  ok <- tryCatch({
    rmarkdown::render(rmd_path, output_file = pdf_path,
                      quiet = TRUE, clean = TRUE)
    TRUE
  }, error = function(e) {
    message('  PDF render failed: ', conditionMessage(e))
    FALSE
  })
  if (ok) message('  WROTE ', pdf_path)
  ok
}

pdf_path <- file.path(.outdir, 'Supplementary_Tables.pdf')
invisible(.render_pdf(all_sheets, pdf_path))
message('Done. Output directory: ', .outdir)
