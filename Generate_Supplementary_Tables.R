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
pab_rds      <- file.path(.shared, 'PAB_data_clean.rds')
sn_echo_xlsx <- file.path(.shared, 'PAB snRNAseq samples.xlsx')
if (file.exists(pab_rds)) {
  pab <- readRDS(pab_rds); md <- pab@meta.data
  pc  <- intersect(c('patient','Patient'),  names(md))[1]
  gcl <- intersect(c('orig.ident','group'), names(md))[1]
  agg <- aggregate(rep(1L, nrow(md)),
                   by = list(mouse_id = as.character(md[[pc]]),
                             grp_raw  = as.character(md[[gcl]])), FUN = sum)
  names(agg)[3] <- 'n_nuclei'
  gmap <- c(Nor = 'Sham', Mod = 'Moderate', Sev = 'Severe')
  s4 <- data.frame(
    ID       = suppressWarnings(as.integer(agg$mouse_id)),
    group    = ifelse(agg$grp_raw %in% names(gmap),
                      gmap[agg$grp_raw], agg$grp_raw),
    n_nuclei = as.integer(agg$n_nuclei),
    stringsAsFactors = FALSE)
  if (file.exists(sn_echo_xlsx)) {
    e <- openxlsx::read.xlsx(sn_echo_xlsx, sheet = 1)
    e$.mid <- suppressWarnings(as.integer(e$ID))
    pick <- c(Body_Weight_g           = 'Surgical.Body.Weight',
              HR                      = 'HR',
              RV_Area_Diastole_mm2    = 'RV.Area.Diastole',
              RV_Area_Systole_mm2     = 'RV.Area.Systole',
              FAC                     = 'FAC',
              RV_Free_Wall_mm         = 'RV.Free.Wall.Thickness.(mm)',
              TAPSE_mm                = 'TAPSE.(mm)',
              TDI_S_Wave              = 'TDI.S.Wave',
              Pk_Gradient_mmHg        = 'Pk.Gradient.mmHg',
              Pk_Velocity_mm_sec      = 'Pk.Velocity.mm/sec',
              Tricuspid_Regurgitation = 'Tricuspid.Regurgitation')
    for (k in names(pick))
      if (pick[[k]] %in% names(e))
        s4[[k]] <- e[[pick[[k]]]][match(s4$ID, e$.mid)]
  }
  s4 <- s4[order(factor(s4$group, c('Sham','Moderate','Severe')), s4$ID), ,
           drop = FALSE]
  rownames(s4) <- NULL
  all_sheets[['Table_S4']] <- s4
  message('  Table S4 (snRNA cohort): ', nrow(s4), ' mice | cols: ',
          ncol(s4), ' | ',
          paste(names(table(s4$group)), table(s4$group),
                sep = '=', collapse = ', '))
  rm(pab, md); invisible(gc(verbose = FALSE))
} else {
  message('  PAB_data_clean.rds not found; skipping Table S4')
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
  s5 <- data.frame(
    Exp_set     = as.character(dat[[exp_col]]),
    ID          = as.integer(dat[[id_col]]),
    PAB_vs_Sham = pab_vs,
    dat[, setdiff(colnames(dat), c(exp_col, id_col)), drop = FALSE],
    check.names = FALSE, stringsAsFactors = FALSE)
  ## Drop per spec: detailed respirometry flux cols + admin date/weight
  ## cols + pure-duplicate key cols (ID / group are already up front).
  drop5 <- c('Echo_Date', 'Echo_Weight_g_', 'Sx_Date',
             'resp_Batch_ID', 'resp_PM', 'resp_ADP', 'resp_Cyt_C',
             'resp_OC', 'resp_Succ', 'resp_Cyt_C_', 'resp_RCR',
             'resp_PAB_vs_SHAM', 'wk8_ID', 'wk8_PAB_vs_SHAM')
  s5 <- s5[, setdiff(colnames(s5), drop5), drop = FALSE]
  s5 <- s5[order(s5$PAB_vs_Sham, s5$Exp_set, s5$ID), , drop = FALSE]
  rownames(s5) <- NULL
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
  s5 <- merge(demo, qc, by = 'patient_id', all = TRUE)
  s5 <- s5[order(s5$group, s5$patient_id), , drop = FALSE]
  rownames(s5) <- NULL
  all_sheets[['Table_S5']] <- s5
  message('  Xenium cohort rows: ', nrow(s5),
          if (nrow(s5) != 9) '  (NOTE: v57 expects n=9)' else '')
} else {
  message('  Xenium_metadata.csv not found; skipping Table S5')
}

# ===========================================================================
# README sheet (leading) — one-line description per Table
# ===========================================================================
index_df <- data.frame(
  `Table #` = c('Table S1', 'Table S2', 'Table S3', 'Table S4', 'Table S5'),
  Title = c(
    'Adult RV cohort: demographics & platform allocation (n=142).',
    'Adult RV cohort: hemodynamics (n=142).',
    'snRNA-seq cohort: demographics & QC (n=11).',
    'Murine PAB snRNA-seq cohort: echo (n=10).',
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
openxlsx::write.xlsx(all_sheets, file = out_path, overwrite = TRUE)
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
  saveRDS(pretty_sheets, rds_path)

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
    paste0('sheets <- readRDS("', rds_path, '")'),
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
    Table_S4 = 'Table S4. Murine PAB snRNA-seq cohort (n=10: 3 Sham, 3 Moderate-RVF, 4 Severe-RVF), with per-mouse echocardiographic parameters — body weight at echo, heart rate (HR), RV areas at diastole (RVAD) and systole (RVAS), fractional area change (FAC), RV free wall thickness (RVFWT), tricuspid annular plane systolic excursion (TAPSE), tissue Doppler S′ (TDI S′), peak gradient and peak velocity across the PA band.',
    Table_S5 = 'Table S5. Patient demographics for the Xenium spatial transcriptomics cohort (n=9), with tissue section details and per-section quality control metrics (cells captured, median UMI/features per cell, median cell area, section area, cell density).'
  )

  for (nm in names(pretty_sheets)) {
    display_nm <- gsub('_', ' ', nm)
    fs <- if (nm == 'INDEX') 11 else 6
    rmd_lines <- c(rmd_lines,
      paste0('## ', display_nm), '',
      '```{r}',
      paste0('df <- sheets[["', nm, '"]]'),
      'kbl(df, format = "latex", booktabs = TRUE, longtable = TRUE,',
      '    na.string = "NA", linesep = "") %>%',
      paste0('  kable_styling(latex_options = c("repeat_header","HOLD_position"),'),
      paste0('                font_size = ', fs, ')'),
      '```',
      ''
    )
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
