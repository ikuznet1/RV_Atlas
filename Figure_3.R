###############################################################################
## Figure 3 (v53 final) — Xenium spatial transcriptomic atlas of the RV
##
## Panels (final, derived from new_scripts/Figure_3.png):
##   (A) Xenium UMAP, 625,035 cells, 12 lineage clusters
##       (CM, FB, EC, Adipo, Endo, Epi, LEC, Myeloid, Neuron, NKT, PC, SM)
##   (B) Marker dotplot — lineage × canonical gene panel
##   (C) Cross-platform concordance scatter:
##       Xenium log₂FC vs snRNA-seq log₂FC, coloured by lineage,
##       symbols mark concordant vs discordant calls
##   (D) Spatial Xenium maps (NF / pRV / RVF) coloured by lineage
##   (E) Spatial Xenium maps (NF / pRV / RVF) coloured by niche
##   (F) Stacked bar of niche composition by lineage (fraction-of-niche)
##   (G) Niche frequency boxplots by NF / pRV / RVF with non-parametric tests
##
## Outputs (./output/Figure_3/v52_figures/):
##   Figure_3A.pdf, Figure_3B.pdf, Figure_3C.pdf, Figure_3D.pdf,
##   Figure_3E.pdf, Figure_3F.pdf, Figure_3G.pdf
##   Supplementary: Figure_3_supp_panel_A_patient.pdf,
##                  Figure_3_supp_panel_F_heatmap.pdf,
##                  Figure_3_celltype_*.pdf,
##                  Figure_3_niche_clust_freq_stats.pdf
###############################################################################

source('./helper_scripts/_shared_helpers.R')

## Per-figure output directory (introduced for consistent output paths)
V52_FIG_DIR <- './output/Figure_3'
dir.create(V52_FIG_DIR, showWarnings = FALSE, recursive = TRUE)



## Per-figure Xenium WRITE directory (CSVs / RDS caches go to ./output, not dependencies)
.x_out <- file.path(V52_FIG_DIR, 'Xenium')
dir.create(.x_out, showWarnings = FALSE, recursive = TRUE)
## Suppress R's default Rplots.pdf in cwd when Rscript hits a plot call
## that's outside an explicit pdf() ... dev.off() envelope.
pdf(NULL)

## Capture a full traceback when F3 errors (the prior queue run died with
## a cryptic "Error: unexpected ')' in: ... abels)" with no call stack).
options(error = function() {
  message('\n=== F3 error traceback (innermost first) ===')
  tb <- sys.calls()
  for (i in rev(seq_along(tb))) {
    message(sprintf('  [%d] %s', i, paste(deparse(tb[[i]]), collapse = ' ')))
  }
  message('=== end traceback ===\n')
  quit(save = 'no', status = 1)
})
dir.create(file.path(V52_FIG_DIR, 'Xenium'), showWarnings = FALSE, recursive = TRUE)
## Composite figure dimensions (inches) — used by theme_v52() scaling
COMP_W <- 14
COMP_H <- 22

## Publication scales (geom widths, point sizes, text sizes) scaled to COMP_W.
PS <- pub_scales(COMP_W)

## ── Panel A cartoon (Endothelium schematic at top-left of Figure_3) ───────
## To regenerate cropped cartoon from the published Figure_3.pdf, run in R:
##   library(magick)
##   img <- image_read_pdf('~/Downloads/hdWGCNA_TOM/Manuscripts/Figure_3.pdf',
##                         pages = 1, density = 300)
##   ## Approx. geometry for the endothelium cartoon at top-left (17 cm ≈ 2008 px):
##   cart <- image_crop(img, '320x120+40+20')
##   image_write(cart, './new_scripts/assets/Figure_3_panel_A_endothelium_cartoon.png',
##               format = 'png', density = 300)
## Then compose with `insert_asset('Figure_3_panel_A_endothelium_cartoon.png')` in the layout.

########################################################
########################################################
########################################################
#### LIBRARIES AND HELPER FUNCTIONS
########################################################
########################################################
########################################################



library(Seurat)
library(arrow)
library(dplyr)
library(sf)
library(harmony)



#' @importFrom magrittr %>% %<>%
NULL

# TODO:
  # - edit args

#' Intermediate solution to \code{subset()}:
#' subset FOVs/centroids if selected cells are NOT found in each FOV
#' NOTE: some code parts and args are taken from SeuratObject

#' Function params/args:
#' @param object An S4 object or A \code{FOV} object
#' @param subset Logical expression indicating features/variables to keep
#' @param cells A vector of cells to keep; if \code{NULL}, defaults to all cells
#' @param idents A vector of identity classes to keep
#' @param Update.slots If to update slots of an object
#' @param Update.object If to update final object, default to TRUE.
#' @param ... Arguments passed to \code{subset()} and other methods


subset_opt <- function(
    object = NULL, 
    subset = NULL, 
    cells = NULL, 
    idents = NULL,
    features = NULL,
    Update.slots = TRUE,
    Update.object = TRUE,
    ...)
{
  
  if (Update.slots) { 
    message("Updating object slots..")
    object %<>% UpdateSlots()
  }
  
  message("Cloing object..")
  obj_subset <- object
  
  # sanity check - use only cell ids (no indices)
  if (all(is.integer(cells))) { 
    cells <- Cells(obj_subset)[cells]
  }
  
  if (!missing(subset) || !is.null(idents)) {
    message("Extracting cells matched to `subset` and/or `idents`")
  }
  
  if (class(obj_subset) == "FOV") {
    message("object class is `FOV` ")
    cells <- Cells(obj_subset)
  } else if (!class(obj_subset) == "FOV" && !missing(subset)) {
    subset <- enquo(arg = subset)
    # cells to keep in the object
    cells <-
      WhichCells(object = obj_subset, 
                 cells = cells,
                 idents = idents,
                 expression = subset,
                 return.null = TRUE, ...)
  } else if (!class(obj_subset) == "FOV" && !is.null(idents)) {
    cells <-
      WhichCells(object = obj_subset, 
                 cells = cells,
                 idents = idents,
                 return.null = TRUE, ...)
  } else if (is.null(cells)) {
    cells <- Cells(obj_subset)
  }
  
  # added support for object class `FOV`
  if (class(obj_subset) == "FOV") {
    message("Matching cells for object class `FOV`..")
    cells_check <- any(obj_subset %>% Cells %in% cells)
  } else { 
    # check if cells are present in all FOV
    message("Matching cells in FOVs..")
    cells_check <-
      lapply(Images(obj_subset) %>% seq, 
             function(i) { 
               any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells) 
             }) %>% unlist
  }
  
  if (all(cells_check)) { 
    message("Cell subsets are found in all FOVs!", "\n",
            "Subsetting object..")
    obj_subset %<>% base::subset(cells = cells, 
                                 idents = idents,
                                 features = features,
                                 ...)
    # subset FOVs
    message("Subsetting FOVs..")
    fovs <- 
      lapply(Images(obj_subset) %>% seq, function(i) {
          base::subset(x = obj_subset[[Images(obj_subset)[i]]],
                       cells = cells, 
                       idents = idents, 
                       features = features, 
                       ...)
      })
    # replace subsetted FOVs
    for (i in fovs %>% seq) { obj_subset[[Images(object)[i]]] <- fovs[[i]] }
    
  } else { 
    # if cells are present only in one or several FOVs:
    # subset FOVs
    fovs <- 
      lapply(Images(obj_subset) %>% seq, function(i) {
        if (any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells)) {
          message("Cell subsets are found only in FOV: ", "\n", Images(obj_subset)[i])
          message("Subsetting Centroids..")
          base::subset(x = obj_subset[[Images(obj_subset)[i]]],
                       cells = cells, 
                       idents = idents, 
                       features = features, 
                       ...)
        }
      })
    # remove FOVs with no matching cells
    message("Removing FOVs where cells are NOT found: ", "\n", 
            paste0(Images(object)[which(!cells_check == TRUE)], "\n"))
    # replace subsetted FOVs
    for (i in fovs %>% seq) { obj_subset[[Images(object)[i]]] <- fovs[[i]] } 
    
    # subset final object
    message("..subset final object")
    obj_subset %<>% 
      base::subset(cells = cells,
                   idents = idents,
                   features = features, 
                   ...)
  }
  
  if (Update.object && !class(obj_subset) == "FOV") { 
    message("Updating object..")
    obj_subset %<>% UpdateSeuratObject() }
  
  message("Object is ready!")
  return(obj_subset)
  
}


#Redefine ReadXenium()
ReadXenium <- function (data.dir, outs = c("matrix", "microns"), type = "centroids", 
          mols.qv.threshold = 20) 
{
  type <- match.arg(arg = type, choices = c("centroids", "segmentations"), 
                    several.ok = TRUE)
  outs <- match.arg(arg = outs, choices = c("matrix", "microns"), 
                    several.ok = TRUE)
  outs <- c(outs, type)
  has_dt <- requireNamespace("data.table", quietly = TRUE) && 
    requireNamespace("R.utils", quietly = TRUE)
  data <- sapply(outs, function(otype) {
    switch(EXPR = otype, matrix = {
      matrix <- suppressWarnings(Read10X(data.dir = file.path(data.dir, 
                                                              "cell_feature_matrix/")))
      matrix
    }, centroids = {
      if (has_dt) {
        cell_info <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                               "cells.csv.gz")))
      } else {
        cell_info <- read.csv(file.path(data.dir, "cells.csv.gz"))
      }
      cell_centroid_df <- data.frame(x = cell_info$x_centroid, 
                                     y = cell_info$y_centroid, cell = cell_info$cell_id, 
                                     stringsAsFactors = FALSE)
      cell_centroid_df
    }, segmentations = {
      if (has_dt) {
        cell_boundaries_df <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                                        "cell_boundaries.csv.gz")))
      } else {
        cell_boundaries_df <- read.csv(file.path(data.dir, 
                                                 "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
      }
      names(cell_boundaries_df) <- c("cell", "x", "y")
      cell_boundaries_df
    }, microns = {
      
      transcripts <- arrow::read_parquet(file.path(data.dir, "transcripts.parquet"))
      transcripts <- subset(transcripts, qv >= mols.qv.threshold)
      
      df <- data.frame(x = transcripts$x_location, y = transcripts$y_location, 
                       gene = transcripts$feature_name, stringsAsFactors = FALSE)
      df
    }, stop("Unknown Xenium input type: ", otype))
  }, USE.NAMES = TRUE)
  return(data)
}


#' Custom Crop function to remove any tissue area
#' @param x A data.frame with \code{x}, \code{y} and eg \code{cell} (cell barcode or ID) variable
#' @param object A \code{Seurat} object. NOTE, currently works on object with single FOV only. 
#' @param col_id A \code{character} vector specifing which variable has the cell ids
#' @param xy_pts A data.frame of xy point coordinates, must have 3 or more xy points, 
#'  it can also be a list of data frames. This will generate a convex hull for cropping
#' @param c_hull_include Everything under (\code{TRUE}) convex hull polygon is included or cropped out (\code{FALSE})
#' @param crop_molecules When \code{Seurat} is present, if to crop molecule cooridnates
#'  NOTE, this can take time especially when there are many molecules.
#' @param BPPARAM description
#' @return Returns cropped data.frame with same variables as input.
#'  Or cropped \code{Seurat} object.
#' @import sf
#' @import BiocParallel
#' @import SeuratObject
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr filter pull transmute
#' @examples
#' # NOTE, `sf` package must be installed!
#' # provide xy coords to make convex hull polygon
#' xy_pts <- 
#' data.frame("x" = c(5400, 11000, 5400, 6100, 7500, 11000),
#' "y" = c(0, 0, 1500, 2000, 4100, 4100))
#' # run crop function on Seurat object (image-based spatial object, like Xenium, Vizgen, etc..)
#' obj_crop <- 
#' Crop_custom(object = seurat_object, col_id = c("cell"), xy_pts = xy_pts,
#' c_hull_include = TRUE,
#' crop_molecules = TRUE,
#' BPPARAM = BiocParallel::MulticoreParam(5, tasks = 10L,force.GC = FALSE, progressbar = TRUE))
#' # run crop function on any dataframe that contains xy coords
#' df_crop <- 
#' Crop_custom(x = df_mols, col_id = c("molecule"), xy_pts = xy_pts,
#' c_hull_include = TRUE,
#' BPPARAM = BiocParallel::MulticoreParam(5, tasks = 10L,force.GC = FALSE, progressbar = TRUE))
#'
Crop_custom <- function(x = NULL, 
           object = NULL,
           col_id = c("cell"),
           xy_pts = NULL,
           c_hull_include = TRUE,
           crop_molecules = TRUE,
           BPPARAM = BiocParallel::SerialParam()) {

    # check packages
    pkgs <- c("data.table", "sf", "BiocParallel", 
              "tidyverse", "magrittr")
    lapply(pkgs %>% seq(), function(i)
    { !requireNamespace(pkgs[i], quietly = TRUE) } ) %>% 
      unlist() %>% 
      { if (c(which(.) > 0) %>% any()) 
      { c("Please install ->", "\n",
          paste0("'", pkgs[which(.)], "'", collapse = ", "), " for this function") %>% 
          stop(., call. = FALSE) } }
    
    # check inputs
    if (is.null(xy_pts)) {
      stop(">>> Please provide xy point coordinates in `xy_pts`")
    }
    
    if (!is.null(object)) {
      if (is(object, "Seurat")) {
        message(">>> Using `Seurat` object")
        x <- NULL
        df_xy <- object[[Images(object)[1]]] %>% GetTissueCoordinates()
      }
    } else if (!is.null(x)) {
      # check x
      is_x <-
      inherits(x = x, 
               what = c("data.frame", "data.table", "tibble", "matrix"))
      #grep("data.frame|data.table|tibble|matrix", 
      #       class(x)) %>% any()
      if (is_x) { 
        df_xy <- x
        object <- NULL
      }
    } else if (is.null(object) && is.null(x)) {
      stop(">>> Please provide either `object` or `x` data.frame")
    }
    
    # make sf data.frame
    sf_df <- st_as_sf(df_xy, coords = c("x", "y"))
    # make convex hull
    if (is(xy_pts, "list")) {
      xy_pts <- data.table::rbindlist(xy_pts)
    }
    c_hull <- 
      st_as_sf(xy_pts, coords = c("x", "y")) %>% 
      st_combine() %>% st_convex_hull()
    crop_df <-
      st_intersection(sf_df, c_hull)
    
    # subset data.frame given the cell ids
    # use col_id instead of $cell
    df_xy %<>% {
      if (c_hull_include) { 
        filter(., !!as.symbol(col_id) %in% pull(crop_df, col_id))
      } else {
        filter(., !(!!as.symbol(col_id)) %in% pull(crop_df, col_id))
      }
    }
    
    # TODO make more cleaner code for cropping df ----
    # output data.frame ----
    if (is.null(object) && !is.null(x)) {
      message(">>> Cropping `data.frame`")
      if (c_hull_include) {
        # faster with `st_join`
        mols <-
          st_join(x = sf_df, 
                  join = st_within, 
                  left = FALSE,
                  y = st_sf(geometry = c_hull))
      } else {
        #mols <- st_difference(sf_df, c_hull)
        mols <-
          st_join(x = sf_df, 
                  join = st_disjoint,
                  left = FALSE,
                  y = st_sf(geometry = c_hull))
      }
      
      genes <- mols %>% pull(col_id) %>% unique()
      mols <-
        bplapply(genes %>% seq(), function(i) {
          mols %>%
            
            # TODO: use col_id instead of molecule, eg !!as.symbol(col_id)
            filter(molecule == genes[i]) %>%
            st_geometry() %>%
            st_coordinates() %>%
            as.data.frame() %>%
            transmute(x = X, y = Y, 
                      molecule = genes[i])
        }, BPPARAM = BPPARAM) %>%
        data.table::rbindlist()
      message(">>> Return output: `data.frame`")
      if (crop_molecules) {
        return(mols)
      } else { return(df_xy) }
    }
    
  
    # cropping Seurat object ----
    # output Seurat
    if (!is.null(object)) {
      message(">>> Return output - `Seurat` object")
      # subset object, centroids and segmentations as well
      object %<>%
        subset(x = ., 
               cells = intersect(x = colnames(x = .),
                                 y =  pull(df_xy, col_id)))
      
      # crop molecules ----
      fov <- object[[Images(object)[1]]]
      if (crop_molecules && grep("molecules", names(fov)) != 0) {
        message(">>> Cropping molecule coordinates - might take time!")
        # using previously made convex hull
        sf_df_mols <- st_as_sf(fov[["molecule"]] %>% GetTissueCoordinates(), 
                               coords = c("x", "y"))
        if (c_hull_include) {
          #mols <- st_intersection(sf_df_mols, c_hull)
          # faster with `st_join`
          mols <-
            st_join(x = sf_df_mols, 
                    join = st_within, 
                    left = FALSE,
                    y = st_sf(geometry = c_hull))
        } else {
          #mols <- st_difference(sf_df_mols, c_hull)
          mols <-
            st_join(x = sf_df_mols, 
                    join = st_disjoint,
                    left = FALSE,
                    y = st_sf(geometry = c_hull))
        }
        genes <- mols$molecule %>% unique()
        mols <-
          bplapply(genes %>% seq(), function(i) {
            mols %>%
              filter(molecule == genes[i]) %>%
              st_geometry() %>%
              st_coordinates() %>%
              as.data.frame() %>%
              transmute(x = X, y = Y, 
                        gene = genes[i])
          }, BPPARAM = BPPARAM)
        # Create Molecule FOV only
        mols %<>% 
          data.table::rbindlist() %>%
          CreateMolecules()
        # replace and add to FOV of the object
        object[[Images(object)[1]]][["molecule"]] <- mols
        
        # TODO: make sure that cropped mols are added to object? ----
        #..and no mols are present in cropped out regions, since few mols were still present.
        # test with GetTissueCoordinates() and plot them!
        
        
        
      }
      validObject(object = object)
      return(object)
    }
  }



########################################################
########################################################
########################################################
#### SPATIAL TILING HELPER
########################################################
########################################################
########################################################

## offset_coords — tile multiple Xenium tissue slices within each disease-state facet.
## Each patient's tissue is normalised to origin, then arranged:
##   2 patients → side by side (gap_um apart)
##   3 patients → two on bottom row, one centred above (row_gap_um above)
##   4+ patients → 2-column grid
offset_coords <- function(meta, gap_um = 500, row_gap_um = 10000) {
  for (g in unique(meta$group)) {
    idx_g <- which(meta$group == g)
    patients <- unique(meta$patient[idx_g])
    for (p in patients) {
      idx_p <- which(meta$group == g & meta$patient == p)
      meta$x[idx_p] <- meta$x[idx_p] - min(meta$x[idx_p])
      meta$y[idx_p] <- meta$y[idx_p] - min(meta$y[idx_p])
    }
    if (length(patients) <= 1) next
    bboxes <- lapply(patients, function(p) {
      idx <- which(meta$group == g & meta$patient == p)
      list(w = diff(range(meta$x[idx])), h = diff(range(meta$y[idx])))
    })
    names(bboxes) <- patients
    areas <- sapply(bboxes, function(b) b$w * b$h)
    patients <- patients[order(-areas)]
    bboxes <- bboxes[as.character(patients)]
    if (length(patients) == 2) {
      idx2 <- which(meta$group == g & meta$patient == patients[2])
      meta$x[idx2] <- meta$x[idx2] + bboxes[[1]]$w + gap_um
    } else if (length(patients) == 3) {
      bot1 <- patients[1]; bot2 <- patients[2]; top1 <- patients[3]
      idx_bot2 <- which(meta$group == g & meta$patient == bot2)
      idx_top  <- which(meta$group == g & meta$patient == top1)
      bot_total_w <- bboxes[[1]]$w + gap_um + bboxes[[2]]$w
      meta$x[idx_bot2] <- meta$x[idx_bot2] + bboxes[[1]]$w + gap_um
      meta$x[idx_top]  <- meta$x[idx_top] + (bot_total_w - bboxes[[3]]$w) / 2
      meta$y[idx_top]  <- meta$y[idx_top] + row_gap_um
    } else {
      ncols <- 2
      nrows <- ceiling(length(patients) / ncols)
      col_widths    <- rep(0, ncols)
      row_heights_v <- rep(0, nrows)
      for (i in seq_along(patients)) {
        r  <- ceiling(i / ncols)
        cc <- ((i - 1) %% ncols) + 1
        col_widths[cc]    <- max(col_widths[cc],    bboxes[[i]]$w)
        row_heights_v[r]  <- max(row_heights_v[r],  bboxes[[i]]$h)
      }
      x_starts <- c(0, cumsum(col_widths[-ncols]    + gap_um))
      y_starts <- c(0, cumsum(row_heights_v[-nrows] + gap_um))
      for (i in seq_along(patients)) {
        p    <- patients[i]
        r    <- ceiling(i / ncols)
        cc   <- ((i - 1) %% ncols) + 1
        idx_p <- which(meta$group == g & meta$patient == p)
        meta$x[idx_p] <- meta$x[idx_p] + x_starts[cc]
        meta$y[idx_p] <- meta$y[idx_p] + y_starts[r]
      }
    }
  }
  meta
}

########################################################
########################################################
########################################################
#### DATA LOADING
#### Source: XeniumReanalysis.r lines 196-617 (additional_scripts)
#### Cached outputs: ~/Documents/XeniumWorkflow/functions/
########################################################

options(future.globals.maxSize = 48012896000)

.wf <- './dependencies/shared/Xenium/'
.xenium_corrected_path <- file.path(.wf, 'Xenium_resegmented_corrected.rds')
.xenium_imputed_path   <- file.path(.wf, 'Xenium_resegmented_imputed_final.rds')

##############################################
#### Corrected (non-imputed)
##############################################

if (file.exists(.xenium_corrected_path)) {
  message('Loading cached Xenium corrected object from ', .xenium_corrected_path)
  xenium.obj <- readRDS(.xenium_corrected_path)
  xenium.obj$names <- xenium.obj$cell_type_rctd_doublet  # alias for downstream code
} else {
  # Per-sample corrected .rds files produced by spatialdata2seurat.py
  xenium.obj.1 <- readRDS(file.path(.wf, 'proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb_corrected.rds'))
  xenium.obj.1$patient <- '1467'; xenium.obj.1$group <- 'RVF'
  idx <- rownames(xenium.obj.1)[1:477]; xenium.obj.1 <- subset(xenium.obj.1, features=idx)

  xenium.obj.2 <- readRDS(file.path(.wf, 'proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb_corrected.rds'))
  xenium.obj.2$patient <- '1632'; xenium.obj.2$group <- 'RVF'
  idx <- rownames(xenium.obj.2)[1:477]; xenium.obj.2 <- subset(xenium.obj.2, features=idx)

  xenium.obj.3 <- readRDS(file.path(.wf, 'proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb_corrected.rds'))
  xenium.obj.3$patient <- '1561'; xenium.obj.3$group <- 'NF'
  idx <- rownames(xenium.obj.3)[1:477]; xenium.obj.3 <- subset(xenium.obj.3, features=idx)

  xenium.obj.4 <- readRDS(file.path(.wf, 'proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb_corrected.rds'))
  xenium.obj.4$patient <- '1343'; xenium.obj.4$group <- 'RVF'
  idx <- rownames(xenium.obj.4)[1:477]; xenium.obj.4 <- subset(xenium.obj.4, features=idx)

  xenium.obj.5 <- readRDS(file.path(.wf, 'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_top_corrected.rds'))
  xenium.obj.5$patient <- '1697'; xenium.obj.5$group <- 'NF'
  idx <- rownames(xenium.obj.5)[1:477]; xenium.obj.5 <- subset(xenium.obj.5, features=idx)

  xenium.obj.6 <- readRDS(file.path(.wf, 'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_middle_corrected.rds'))
  xenium.obj.6$patient <- '1691'; xenium.obj.6$group <- 'NF'
  idx <- rownames(xenium.obj.6)[1:477]; xenium.obj.6 <- subset(xenium.obj.6, features=idx)

  xenium.obj.7 <- readRDS(file.path(.wf, 'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_bottom_corrected.rds'))
  xenium.obj.7$patient <- '1618'; xenium.obj.7$group <- 'pRV'
  idx <- rownames(xenium.obj.7)[1:477]; xenium.obj.7 <- subset(xenium.obj.7, features=idx)

  xenium.obj.8 <- readRDS(file.path(.wf, 'proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_top_corrected.rds'))
  xenium.obj.8$patient <- '1567'; xenium.obj.8$group <- 'pRV'
  idx <- rownames(xenium.obj.8)[1:477]; xenium.obj.8 <- subset(xenium.obj.8, features=idx)

  xenium.obj.9 <- readRDS(file.path(.wf, 'proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_bottom_corrected.rds'))
  xenium.obj.9$patient <- '1692'; xenium.obj.9$group <- 'pRV'
  idx <- rownames(xenium.obj.9)[1:477]; xenium.obj.9 <- subset(xenium.obj.9, features=idx)

  xenium.obj <- merge(xenium.obj.1, xenium.obj.2)
  xenium.obj <- merge(xenium.obj, xenium.obj.3)
  xenium.obj <- merge(xenium.obj, xenium.obj.4)
  xenium.obj <- merge(xenium.obj, xenium.obj.5)
  xenium.obj <- merge(xenium.obj, xenium.obj.6)
  xenium.obj <- merge(xenium.obj, xenium.obj.7)
  xenium.obj <- merge(xenium.obj, xenium.obj.8)
  xenium.obj <- merge(xenium.obj, xenium.obj.9)

  xenium.obj <- NormalizeData(xenium.obj)
  xenium.obj <- FindVariableFeatures(xenium.obj)
  xenium.obj <- ScaleData(xenium.obj, split.by = 'patient')
  xenium.obj <- RunPCA(xenium.obj, npcs = 30)
  xenium.obj <- JoinLayers(xenium.obj)
  xenium.obj <- RunHarmony(xenium.obj, group.by.vars = 'patient', reduction.save = "harmony",
                            kmeans_init_nstart = 20, kmeans_init_iter_max = 100)
  xenium.obj <- FindNeighbors(xenium.obj, dims = 1:30, reduction = "harmony")
  xenium.obj <- FindClusters(xenium.obj, resolution = 1)
  xenium.obj <- RunUMAP(xenium.obj, reduction = "pca", dims = 1:30, reduction.name = 'umap_orig')
  xenium.obj <- RunUMAP(xenium.obj, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  saveRDS(xenium.obj, .xenium_corrected_path)

  rm(xenium.obj.1, xenium.obj.2, xenium.obj.3, xenium.obj.4, xenium.obj.5,
     xenium.obj.6, xenium.obj.7, xenium.obj.8, xenium.obj.9); gc()
  xenium.obj$names <- xenium.obj$cell_type_rctd_doublet
}

##############################################
#### Imputed — XeniumReanalysis.r lines 534-617
##############################################

if (file.exists(.xenium_imputed_path)) {
  message('Loading cached Xenium imputed object from ', .xenium_imputed_path)
  xenium.imp <- readRDS(.xenium_imputed_path)
} else {
  library(BPCells)

  imputed_files <- file.path(.wf, c(
    'proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb_imputed.rds',
    'proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb_imputed.rds',
    'proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb_imputed.rds',
    'proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb_imputed.rds',
    'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_top_imputed.rds',
    'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_middle_imputed.rds',
    'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_bottom_imputed.rds',
    'proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_top_imputed.rds',
    'proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_bottom_imputed.rds'
  ))
  patients <- c('1467','1632','1561','1343','1697','1691','1618','1567','1692')
  groups   <- c('RVF','RVF','NF','RVF','NF','NF','pRV','pRV','pRV')

  bp_dir <- file.path(.wf, 'bpcells_imputed')
  dir.create(bp_dir, showWarnings = FALSE, recursive = TRUE)

  # Pass 1: discover shared genes
  message("=== Pass 1: discovering shared gene set ===")
  gene_sets <- list()
  for (i in seq_along(imputed_files)) {
    message("  Reading gene names from sample ", i, "/", length(imputed_files))
    obj <- readRDS(imputed_files[i])
    gene_sets[[i]] <- rownames(obj)
    rm(obj); gc(verbose = FALSE)
  }
  shared_genes <- Reduce(intersect, gene_sets)
  message("Shared genes: ", length(shared_genes))
  rm(gene_sets); gc(verbose = FALSE)

  # Pass 2: convert to on-disk BPCells format
  message("=== Pass 2: converting to on-disk BPCells format ===")
  obj_list <- vector("list", length(imputed_files))
  for (i in seq_along(imputed_files)) {
    message("Processing sample ", i, "/", length(imputed_files), ": ", patients[i])
    obj <- readRDS(imputed_files[i])
    obj$patient <- patients[i]
    obj$group   <- groups[i]
    obj <- subset(obj, features = shared_genes)
    assay_name <- DefaultAssay(obj)
    counts_path <- file.path(bp_dir, paste0("sample_", i, "_counts"))
    mat <- GetAssayData(obj, assay = assay_name, layer = "counts")
    BPCells::write_matrix_dir(mat, dir = counts_path, overwrite = TRUE)
    obj[[assay_name]]$counts <- BPCells::open_matrix_dir(counts_path)
    obj_list[[i]] <- obj
    rm(obj); gc(verbose = FALSE)
  }

  message("=== Merging ===")
  xenium.imp <- merge(obj_list[[1]], obj_list[2:length(obj_list)])
  rm(obj_list); gc(verbose = FALSE)

  xenium.imp <- NormalizeData(xenium.imp)
  xenium.imp <- FindVariableFeatures(xenium.imp)
  xenium.imp <- ScaleData(xenium.imp, split.by = 'patient')
  xenium.imp <- RunPCA(xenium.imp, npcs = 30)
  xenium.imp <- JoinLayers(xenium.imp)
  xenium.imp <- RunHarmony(xenium.imp, group.by.vars = 'patient', reduction.save = "harmony")
  xenium.imp <- FindNeighbors(xenium.imp, dims = 1:30, reduction = "harmony")
  xenium.imp <- FindClusters(xenium.imp, resolution = 1)
  xenium.imp <- RunUMAP(xenium.imp, reduction = "pca", dims = 1:30, reduction.name = 'umap_orig')
  xenium.imp <- RunUMAP(xenium.imp, reduction = "harmony", dims = 1:30, reduction.name = 'umap')

  saveRDS(xenium.imp, file.path(.wf, 'Xenium_resegmented_imputed.rds'))

  # Filter: remove low-count cells and unassigned types
  assay_name <- DefaultAssay(xenium.imp)
  xenium.imp <- xenium.imp[, colSums(GetAssayData(xenium.imp, assay = assay_name, layer = "counts")[rownames(xenium.obj), ]) > 10]
  xenium.imp <- xenium.imp[, !(xenium.imp$cell_type_rctd_doublet %in% c('None', 'Unknown'))]
  saveRDS(xenium.imp, .xenium_imputed_path)
  message("=== Done — saved Xenium_resegmented_imputed_final.rds ===")
}



VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

xenium.obj$subnames <- as.character(xenium.obj$names)
xenium.obj$markernames <- as.character(xenium.obj$names)


#ImageDimPlot(xenium.obj, fov = "fov.2",split.by='ident',cols = "polychrome", axes = TRUE) # FOV name may differ in new object

##############################################
##############################################
#### FINE-GRAINED CELL TYPE ASSIGNMENT
##############################################
##############################################

# ── Subclustering cache root ──────────────────────────────────────────────────
.xcache <- './dependencies/shared/Xenium'
dir.create(.xcache, showWarnings = FALSE, recursive = TRUE)

# Initialize subnames to broad cell type label
xenium.obj$subnames <- xenium.obj$cell_type_rctd_doublet

##############################################
#### CM subclustering
##############################################
.cache_cm <- file.path(.xcache, 'cm_clean_clean.rds')
if (file.exists(.cache_cm)) {
  message('Loading cached CM subclustering...')
  cm.obj.clean.clean <- readRDS(.cache_cm)
} else {
  ref    <- readRDS('./dependencies/shared/snRV_ref.rds')
  cm.ref <- subset(ref, subset = Subnames_manual %in%
    c('Cm1','Cm2','Cm3','Cm4','Cm5','Cm6','Cm7','Cm8','Cm9','Cm10'))

  cm.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet == 'CM')
  cm.obj <- JoinLayers(cm.obj)
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(cm.obj@meta.data), value = TRUE)
  if (length(umap_cols)) cm.obj@meta.data[umap_cols] <- NULL
  cm.obj[["pca"]] <- NULL; cm.obj[["harmony"]] <- NULL
  cm.obj[["umap"]] <- NULL; cm.obj[["umap_orig"]] <- NULL

  cm.obj <- FindVariableFeatures(cm.obj)
  cm.obj <- ScaleData(cm.obj)
  cm.obj <- RunPCA(cm.obj, npcs = 30)
  cm.obj <- RunHarmony(cm.obj, "patient")
  cm.obj <- RunUMAP(cm.obj, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  cm.obj <- RunUMAP(cm.obj, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  cm.obj <- FindNeighbors(cm.obj, reduction = "harmony", dims = 1:30)
  cm.obj <- FindClusters(cm.obj, resolution = 0.5)

  cm_labels <- c(
    '0'='CM_Metabolic','1'='CM_Stressed','2'='Fibroblast','3'='CM_Ventricular',
    '4'='CM_Contractile','5'='Endothelial','6'='CM_Ventricular_2','7'='Macrophage',
    '8'='Fibroblast_2','9'='Fibroblast_3','10'='Endothelial_2','11'='Pericyte',
    '12'='Fibroblast_4','13'='Mast_Cell','14'='CM_Ventricular_4','15'='T_NK_Cell',
    '16'='Dendritic','17'='Adipocyte','18'='Macrophage_2','19'='Fibroblast_5',
    '20'='Lymphatic_EC','21'='Schwann','22'='Proliferating','23'='CM_Ventricular_5')
  cm.obj$cm_subtype <- unname(cm_labels[as.character(cm.obj$seurat_clusters)])

  cm.obj.clean <- subset(cm.obj, subset = cm_subtype %in% c(
    'CM_Metabolic','CM_Stressed','CM_Ventricular','CM_Contractile',
    'CM_Ventricular_2','CM_Ventricular_4','CM_Ventricular_5'))
  cm.obj.clean[["pca"]] <- NULL; cm.obj.clean[["harmony"]] <- NULL
  cm.obj.clean[["umap"]] <- NULL; cm.obj.clean[["umap_orig"]] <- NULL
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(cm.obj.clean@meta.data), value = TRUE)
  if (length(umap_cols)) cm.obj.clean@meta.data[umap_cols] <- NULL

  cm.obj.clean <- JoinLayers(cm.obj.clean)
  cm.obj.clean <- FindVariableFeatures(cm.obj.clean)
  cm.obj.clean <- ScaleData(cm.obj.clean)
  cm.obj.clean <- RunPCA(cm.obj.clean, npcs = 30)
  cm.obj.clean <- RunHarmony(cm.obj.clean, "patient")
  cm.obj.clean <- RunUMAP(cm.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  cm.obj.clean <- RunUMAP(cm.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  cm.obj.clean <- FindNeighbors(cm.obj.clean, reduction = "harmony", dims = 1:30)
  cm.obj.clean <- FindClusters(cm.obj.clean, resolution = 0.5)

  cm_clean_labels <- c(
    '0'='CM_Ventricular','1'='CM_Ventricular_2','2'='CM_RV',
    '3'='CM_RV_2','4'='CM_Ventricular_3','5'='CM_Ventricular_4','6'='CM_RV_3',
    '7'='CM_Ventricular_5','8'='CM_Contractile','9'='CM_Ventricular_6',
    '10'='CM_Ventricular_7','11'='CM_Stressed')
  cm.obj.clean$cm_subtype <- unname(cm_clean_labels[as.character(cm.obj.clean$seurat_clusters)])

  cm.obj.clean.clean <- cm.obj.clean  # all clusters are true CMs

  cm.ref <- JoinLayers(cm.ref)
  cm.anchors <- FindTransferAnchors(
    reference = cm.ref, query = cm.obj.clean.clean,
    dims = 1:30, reference.reduction = 'pca')
  cm.predictions <- TransferData(
    anchorset = cm.anchors, refdata = cm.ref$Subnames_manual, dims = 1:30)
  cm.obj.clean.clean <- AddMetaData(cm.obj.clean.clean, metadata = cm.predictions)

  saveRDS(cm.obj.clean.clean, .cache_cm)
  rm(cm.obj, cm.obj.clean, cm.ref, ref); gc()
}
idx <- match(colnames(cm.obj.clean.clean), colnames(xenium.obj))
.keep <- !is.na(idx)
xenium.obj$subnames[idx[.keep]] <- cm.obj.clean.clean$cm_subtype[.keep]
rm(cm.obj.clean.clean); gc()

##############################################
#### FB subclustering
##############################################
.cache_fb <- file.path(.xcache, 'fb_clean_clean.rds')
if (file.exists(.cache_fb)) {
  message('Loading cached FB subclustering...')
  fb.obj.clean.clean <- readRDS(.cache_fb)
} else {
  # FB subclustering — full pipeline (from XeniumReanalysis.r lines 3410-3792)
  ref    <- readRDS('./dependencies/shared/snRV_ref.rds')
  fb.ref <- subset(ref, subset = Names %in% c('FB'))

  fb.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet %in% c('FB'))
  fb.obj <- JoinLayers(fb.obj)

  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(fb.obj@meta.data), value = TRUE)
  if (length(umap_cols)) fb.obj@meta.data[umap_cols] <- NULL
  fb.obj[["pca"]] <- NULL; fb.obj[["harmony"]] <- NULL
  fb.obj[["umap"]] <- NULL; fb.obj[["umap_orig"]] <- NULL

  fb.obj <- FindVariableFeatures(fb.obj)
  fb.obj <- ScaleData(fb.obj)
  fb.obj <- RunPCA(fb.obj, npcs = 30)
  fb.obj <- RunHarmony(fb.obj, group.by.vars = "patient",
                       theta = 8, lambda = 0.5, sigma = 0.1,
                       max_iter = 30, epsilon.cluster = -Inf,
                       epsilon.harmony = -Inf, plot_convergence = FALSE)

  library(uwot)
  .build_umap_reduc <- function(obj, reduction, name, key_prefix, dims = NULL) {
    emb <- Embeddings(obj, reduction)
    if (!is.null(dims)) emb <- emb[, dims, drop = FALSE]
    um <- uwot::umap(emb, n_neighbors = 30, min_dist = 0.3, metric = "cosine")
    rownames(um) <- rownames(emb)
    colnames(um) <- paste0(key_prefix, 1:2)
    obj[[name]] <- CreateDimReducObject(embeddings = um, key = key_prefix,
                                        assay = DefaultAssay(obj))
    obj
  }
  fb.obj <- .build_umap_reduc(fb.obj, "harmony", "umap",      "umap_",     dims = 1:30)
  fb.obj <- .build_umap_reduc(fb.obj, "pca",     "umap_orig", "umaporig_", dims = 1:30)
  fb.obj <- FindNeighbors(fb.obj, reduction = "harmony", dims = 1:30)
  fb.obj <- FindClusters(fb.obj, resolution = 0.5)

  fb_labels <- c(
    '0'  = 'FB_PCOLCE2',    '1'  = 'FB_Epicardial',  '2'  = 'FB_Homeostatic',
    '3'  = 'drop_CM',       '4'  = 'FB_Activated',   '5'  = 'FB_Resting',
    '6'  = 'drop_LowConf',  '7'  = 'drop_Ambig',     '8'  = 'drop_Macrophage',
    '9'  = 'drop_Mast',     '10' = 'drop_Pericyte',  '11' = 'drop_Mural',
    '12' = 'drop_Ambig',    '13' = 'drop_EC_Art',    '14' = 'drop_cDC_Mono',
    '15' = 'drop_Adipocyte','16' = 'drop_Tcell',     '17' = 'drop_CM2',
    '18' = 'drop_CM_stressed','19'= 'drop_Schwann',  '20' = 'drop_LEC',
    '21' = 'drop_Mural2',   '22' = 'drop_Myeloid_Prolif'
  )
  fb.obj$fb_subtype <- unname(fb_labels[as.character(fb.obj$seurat_clusters)])

  fb.obj.clean <- subset(fb.obj, subset = fb_subtype %in% c(
    'FB_PCOLCE2', 'FB_Epicardial', 'FB_Homeostatic', 'FB_Activated', 'FB_Resting'))
  fb.obj.clean[["pca"]] <- NULL; fb.obj.clean[["harmony"]] <- NULL
  fb.obj.clean[["umap"]] <- NULL; fb.obj.clean[["umap_orig"]] <- NULL
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(fb.obj.clean@meta.data), value = TRUE)
  if (length(umap_cols)) fb.obj.clean@meta.data[umap_cols] <- NULL

  fb.obj.clean <- JoinLayers(fb.obj.clean)
  fb.obj.clean <- FindVariableFeatures(fb.obj.clean)
  fb.obj.clean <- ScaleData(fb.obj.clean)
  fb.obj.clean <- RunPCA(fb.obj.clean, npcs = 30)
  fb.obj.clean <- RunHarmony(fb.obj.clean, "patient")
  fb.obj.clean <- RunUMAP(fb.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  fb.obj.clean <- RunUMAP(fb.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  fb.obj.clean <- FindNeighbors(fb.obj.clean, reduction = "harmony", dims = 1:30)
  fb.obj.clean <- FindClusters(fb.obj.clean, resolution = 0.5)

  fb_clean_labels <- c(
    '0'  = 'FB_PCOLCE2', '1'  = 'FB_Resting',  '2'  = 'FB_NOX4',
    '3'  = 'FB_Activated','4'  = 'FB_Homeostatic','5' = 'drop_EC',
    '6'  = 'FB_Resting', '7'  = 'drop_Ambig',  '8'  = 'drop_Mural',
    '9'  = 'FB_Stress',  '10' = 'FB_NTN4'
  )
  fb.obj.clean$fb_subtype <- unname(fb_clean_labels[as.character(fb.obj.clean$seurat_clusters)])

  fb.obj.clean.clean <- subset(fb.obj.clean, subset = fb_subtype %in% c(
    'FB_PCOLCE2', 'FB_Resting', 'FB_NOX4', 'FB_Activated',
    'FB_Homeostatic', 'FB_Stress', 'FB_NTN4'))
  fb.obj.clean.clean[["pca"]] <- NULL; fb.obj.clean.clean[["harmony"]] <- NULL
  fb.obj.clean.clean[["umap"]] <- NULL; fb.obj.clean.clean[["umap_orig"]] <- NULL
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(fb.obj.clean.clean@meta.data), value = TRUE)
  if (length(umap_cols)) fb.obj.clean.clean@meta.data[umap_cols] <- NULL

  fb.obj.clean.clean <- JoinLayers(fb.obj.clean.clean)
  fb.obj.clean.clean <- FindVariableFeatures(fb.obj.clean.clean)
  fb.obj.clean.clean <- ScaleData(fb.obj.clean.clean)
  fb.obj.clean.clean <- RunPCA(fb.obj.clean.clean, npcs = 30)
  fb.obj.clean.clean <- RunHarmony(fb.obj.clean.clean, "patient")
  fb.obj.clean.clean <- RunUMAP(fb.obj.clean.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  fb.obj.clean.clean <- RunUMAP(fb.obj.clean.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

  fb.ref <- JoinLayers(fb.ref)
  fb.anchors <- FindTransferAnchors(reference = fb.ref, query = fb.obj.clean.clean,
                                    dims = 1:30, reference.reduction = 'pca')
  fb.predictions <- TransferData(anchorset = fb.anchors, refdata = fb.ref$Subnames, dims = 1:30)
  fb.obj.clean.clean <- AddMetaData(fb.obj.clean.clean, metadata = fb.predictions)

  saveRDS(fb.obj.clean.clean, .cache_fb)
  rm(fb.obj, fb.obj.clean, fb.ref, ref); gc()
}
idx <- match(colnames(fb.obj.clean.clean), colnames(xenium.obj))
.keep <- !is.na(idx)
xenium.obj$subnames[idx[.keep]] <- fb.obj.clean.clean$fb_subtype[.keep]
rm(fb.obj.clean.clean); gc()

##############################################
#### EC subclustering (EC + Endo + LEC)
##############################################
.cache_ec <- file.path(.xcache, 'ec_clean_clean.rds')
if (file.exists(.cache_ec)) {
  message('Loading cached EC subclustering...')
  ec.obj.clean.clean <- readRDS(.cache_ec)
} else {
  # EC subclustering — full pipeline (from XeniumReanalysis.r lines 2941-3182)
  ref    <- readRDS('./dependencies/shared/snRV_ref.rds')
  ec.ref <- subset(ref, subset = Subnames %in% c(
    'EC_Arterial','EC_Capillary','EC_Endocardial','EC_Lymph','EC_Venous'))

  ec.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet %in% c('EC', 'Endo', 'LEC'))
  ec.obj <- JoinLayers(ec.obj)
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(ec.obj@meta.data), value = TRUE)
  if (length(umap_cols)) ec.obj@meta.data[umap_cols] <- NULL
  ec.obj[["pca"]] <- NULL; ec.obj[["harmony"]] <- NULL
  ec.obj[["umap"]] <- NULL; ec.obj[["umap_orig"]] <- NULL

  ec.obj <- FindVariableFeatures(ec.obj)
  ec.obj <- ScaleData(ec.obj)
  ec.obj <- RunPCA(ec.obj, npcs = 30)
  ec.obj <- RunHarmony(ec.obj, "patient")
  ec.obj <- RunUMAP(ec.obj, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  ec.obj <- RunUMAP(ec.obj, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  ec.obj <- FindNeighbors(ec.obj, reduction = "harmony", dims = 1:30)
  ec.obj <- FindClusters(ec.obj, resolution = 0.5)

  ec_labels <- c(
    '0'  = 'Capillary_EC', '1'  = 'Cardiomyocyte','2'  = 'Arterial_EC',
    '3'  = 'Fibroblast',   '4'  = 'Capillary_EC', '5'  = 'Neutrophil',
    '6'  = 'Myeloid',      '7'  = 'Venous_EC',    '8'  = 'Pericyte',
    '9'  = 'Endocardial',  '10' = 'Pericyte',     '11' = 'Endocardial',
    '12' = 'Unknown',      '13' = 'LEC',           '14' = 'TNK_Cell',
    '15' = 'Pericyte',     '16' = 'Unknown',       '17' = 'Venous_EC',
    '18' = 'Venous_EC'
  )
  ec.obj$ec_subtype <- unname(ec_labels[as.character(ec.obj$seurat_clusters)])

  ec.obj.clean.clean <- subset(ec.obj, subset = ec_subtype %in% c(
    'Capillary_EC', 'Arterial_EC', 'Venous_EC', 'Endocardial', 'LEC'))
  ec.obj.clean.clean[["pca"]] <- NULL; ec.obj.clean.clean[["harmony"]] <- NULL
  ec.obj.clean.clean[["umap"]] <- NULL; ec.obj.clean.clean[["umap_orig"]] <- NULL
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(ec.obj.clean.clean@meta.data), value = TRUE)
  if (length(umap_cols)) ec.obj.clean.clean@meta.data[umap_cols] <- NULL

  ec.obj.clean.clean <- JoinLayers(ec.obj.clean.clean)
  ec.obj.clean.clean <- FindVariableFeatures(ec.obj.clean.clean)
  ec.obj.clean.clean <- ScaleData(ec.obj.clean.clean)
  ec.obj.clean.clean <- RunPCA(ec.obj.clean.clean, npcs = 30)
  ec.obj.clean.clean <- RunHarmony(ec.obj.clean.clean, "patient")
  ec.obj.clean.clean <- RunUMAP(ec.obj.clean.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  ec.obj.clean.clean <- RunUMAP(ec.obj.clean.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  ec.obj.clean.clean <- FindNeighbors(ec.obj.clean.clean, reduction = "harmony", dims = 1:30)
  ec.obj.clean.clean <- FindClusters(ec.obj.clean.clean, resolution = 0.5)

  ec_clean_labels <- c(
    '0'  = 'Capillary_EC',    '1'  = 'Capillary_EC',    '2'  = 'Arterial_EC',
    '3'  = 'Arterial_EC',     '4'  = 'Capillary_EC',    '5'  = 'Venous_EC',
    '6'  = 'Endocardial',     '7'  = 'Endocardial',     '8'  = 'LEC',
    '9'  = 'Endocardial',     '10' = 'Cardiomyocyte',   '11' = 'Schwann_Cell',
    '12' = 'Proliferating_EC','13' = 'Activated_EC'
  )
  ec.obj.clean.clean$ec_subtype <- unname(ec_clean_labels[as.character(ec.obj.clean.clean$seurat_clusters)])

  # Remove contamination (CM, Schwann) and proliferating cluster 12
  ec.obj.clean.clean <- subset(ec.obj.clean.clean,
    subset = ec_subtype %in% c('Capillary_EC', 'Arterial_EC', 'Venous_EC',
                                'Endocardial', 'LEC', 'Activated_EC'))

  ec.ref <- JoinLayers(ec.ref)
  ec.anchors <- FindTransferAnchors(reference = ec.ref, query = ec.obj.clean.clean,
                                    dims = 1:30, reference.reduction = 'pca')
  ec.predictions <- TransferData(anchorset = ec.anchors, refdata = ec.ref$Subnames, dims = 1:30)
  ec.obj.clean.clean <- AddMetaData(ec.obj.clean.clean, metadata = ec.predictions)

  saveRDS(ec.obj.clean.clean, .cache_ec)
  rm(ec.obj, ec.ref, ref); gc()
}
idx <- match(colnames(ec.obj.clean.clean), colnames(xenium.obj))
.keep <- !is.na(idx)
xenium.obj$subnames[idx[.keep]] <- ec.obj.clean.clean$ec_subtype[.keep]
rm(ec.obj.clean.clean); gc()

##############################################
#### Myeloid subclustering
##############################################
.cache_myeloid <- file.path(.xcache, 'myeloid_clean_clean.rds')
if (file.exists(.cache_myeloid)) {
  message('Loading cached Myeloid subclustering...')
  myeloid.obj.clean.clean <- readRDS(.cache_myeloid)
} else {
  # Myeloid subclustering — full 3-round pipeline (from XeniumReanalysis.r lines 1102-1548)
  # NOTE: also produces tnk.obj.clean used by the NKT section below
  ref         <- readRDS('./dependencies/shared/snRV_ref.rds')
  myeloid.ref <- subset(ref, subset = Names %in% c('Myeloid'))

  myeloid.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet == "Myeloid")
  myeloid.obj <- JoinLayers(myeloid.obj)
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(myeloid.obj@meta.data), value = TRUE)
  if (length(umap_cols)) myeloid.obj@meta.data[umap_cols] <- NULL
  myeloid.obj[["pca"]] <- NULL; myeloid.obj[["harmony"]] <- NULL
  myeloid.obj[["umap"]] <- NULL; myeloid.obj[["umap_orig"]] <- NULL

  myeloid.obj <- FindVariableFeatures(myeloid.obj)
  myeloid.obj <- ScaleData(myeloid.obj)
  myeloid.obj <- RunPCA(myeloid.obj, npcs = 30)
  myeloid.obj <- RunHarmony(myeloid.obj, "patient")
  myeloid.obj <- RunUMAP(myeloid.obj, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  myeloid.obj <- RunUMAP(myeloid.obj, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  myeloid.obj <- FindNeighbors(myeloid.obj, reduction = "harmony", dims = 1:30)
  myeloid.obj <- FindClusters(myeloid.obj, resolution = 0.5)

  myeloid_labels <- c(
    '0'  = 'Macrophage_Resident',    '1'  = 'Macrophage_C1q',
    '2'  = 'Fibroblast',             '3'  = 'Cardiomyocyte',
    '4'  = 'Unknown',                '5'  = 'Monocyte',
    '6'  = 'Fibroblast_Doublet',     '7'  = 'Neutrophil',
    '8'  = 'Pericyte',               '9'  = 'Pericyte',
    '10' = 'Fibroblast',             '11' = 'Cardiomyocyte',
    '12' = 'Adipocyte',              '13' = 'TNK_Cell',
    '14' = 'Basophil',               '15' = 'Macrophage_Inflammatory',
    '16' = 'Myeloid_Proliferating',  '17' = 'Unknown',
    '18' = 'Schwann_Cell',           '19' = 'LymphEndothelial',
    '20' = 'Neutrophil',             '21' = 'Unknown'
  )
  myeloid.obj$myeloid_subtype <- unname(myeloid_labels[as.character(myeloid.obj$seurat_clusters)])

  myeloid.obj.clean <- subset(myeloid.obj, subset = myeloid_subtype %in% c(
    'Macrophage_Resident','Macrophage_C1q','Monocyte','Macrophage_Inflammatory',
    'Myeloid_Proliferating','Neutrophil','Basophil'))
  tnk.obj.clean <- subset(myeloid.obj, subset = myeloid_subtype %in% c('TNK_Cell'))

  myeloid.obj.clean[["pca"]] <- NULL; myeloid.obj.clean[["harmony"]] <- NULL
  myeloid.obj.clean[["umap"]] <- NULL
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(myeloid.obj.clean@meta.data), value = TRUE)
  if (length(umap_cols)) myeloid.obj.clean@meta.data[umap_cols] <- NULL

  myeloid.obj.clean <- JoinLayers(myeloid.obj.clean)
  myeloid.obj.clean <- FindVariableFeatures(myeloid.obj.clean)
  myeloid.obj.clean <- ScaleData(myeloid.obj.clean, split.by = 'patient')
  myeloid.obj.clean <- RunPCA(myeloid.obj.clean, npcs = 30)
  myeloid.obj.clean <- RunHarmony(myeloid.obj.clean, "patient")
  myeloid.obj.clean <- RunUMAP(myeloid.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  myeloid.obj.clean <- RunUMAP(myeloid.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  myeloid.obj.clean <- FindNeighbors(myeloid.obj.clean, reduction = "harmony", dims = 1:30)
  myeloid.obj.clean <- FindClusters(myeloid.obj.clean, resolution = 0.5)

  myeloid_clean_labels <- c(
    '0'  = 'Macrophage_Resident', '1'  = 'Macrophage_C1q',
    '2'  = 'Monocyte',            '3'  = 'Dendritic_Cell',
    '4'  = 'Doublet',             '5'  = 'Doublet',
    '6'  = 'Macrophage_Resident', '7'  = 'Macrophage_Resident',
    '8'  = 'Dendritic_Cell',      '9'  = 'Myeloid_Proliferating',
    '10' = 'Macrophage_Inflammatory','11'= 'Mast_Cell',
    '12' = 'Macrophage_TREM2',    '13' = 'CD8_T'
  )
  myeloid.obj.clean$myeloid_subtype <- unname(
    myeloid_clean_labels[as.character(myeloid.obj.clean$seurat_clusters)])

  myeloid.obj.clean.clean <- subset(myeloid.obj.clean,
    subset = myeloid_subtype %in% c('Macrophage_Resident', 'Macrophage_C1q', 'Monocyte',
                                     'Dendritic_Cell', 'Myeloid_Proliferating',
                                     'Macrophage_Inflammatory', 'Mast_Cell', 'Macrophage_TREM2'))
  myeloid.obj.clean.clean[["pca"]] <- NULL; myeloid.obj.clean.clean[["harmony"]] <- NULL
  myeloid.obj.clean.clean[["umap"]] <- NULL; myeloid.obj.clean.clean[["umap_orig"]] <- NULL
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(myeloid.obj.clean.clean@meta.data), value = TRUE)
  if (length(umap_cols)) myeloid.obj.clean.clean@meta.data[umap_cols] <- NULL

  myeloid.obj.clean.clean <- JoinLayers(myeloid.obj.clean.clean)
  myeloid.obj.clean.clean <- FindVariableFeatures(myeloid.obj.clean.clean)
  myeloid.obj.clean.clean <- ScaleData(myeloid.obj.clean.clean)
  myeloid.obj.clean.clean <- RunPCA(myeloid.obj.clean.clean, npcs = 30)
  myeloid.obj.clean.clean <- RunHarmony(myeloid.obj.clean.clean, "patient")
  myeloid.obj.clean.clean <- RunUMAP(myeloid.obj.clean.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  myeloid.obj.clean.clean <- RunUMAP(myeloid.obj.clean.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  myeloid.obj.clean.clean <- FindNeighbors(myeloid.obj.clean.clean, reduction = "harmony", dims = 1:30)
  myeloid.obj.clean.clean <- FindClusters(myeloid.obj.clean.clean, resolution = 2)

  myeloid_clean_clean_labels <- c(
    '0'  = 'Macrophage_Resident',         '1'  = 'Macrophage_Resident',
    '2'  = 'Macrophage_C1q',              '3'  = 'Macrophage_Resident',
    '4'  = 'Macrophage_Resident',         '5'  = 'Macrophage_Resident',
    '6'  = 'Macrophage_Resident',         '7'  = 'Macrophage_Resident',
    '8'  = 'Macrophage_Monocyte_Derived', '9'  = 'Macrophage_Resident',
    '10' = 'Macrophage_Resident',         '11' = 'Monocyte',
    '12' = 'Macrophage_Resident',         '13' = 'Adipocyte',
    '14' = 'Macrophage_Resident',         '15' = 'Myeloid_Proliferating',
    '16' = 'Dendritic_Cell',              '17' = 'Macrophage_Resident',
    '18' = 'Dendritic_Cell',              '19' = 'Doublet',
    '20' = 'Doublet',                     '21' = 'Macrophage_Resident',
    '22' = 'Macrophage_Resident',         '23' = 'Mast_Cell',
    '24' = 'Macrophage_Inflammatory',     '25' = 'Macrophage_TREM2',
    '26' = 'Doublet',                     '27' = 'Macrophage_Resident',
    '28' = 'Low_Quality',                 '29' = 'Doublet',
    '30' = 'Doublet',                     '31' = 'Macrophage_Resident'
  )
  myeloid.obj.clean.clean$myeloid_clean_subtype <- unname(
    myeloid_clean_clean_labels[as.character(myeloid.obj.clean.clean$seurat_clusters)])

  myeloid.obj.clean.clean <- subset(myeloid.obj.clean.clean,
    subset = myeloid_clean_subtype %in% c('Macrophage_Resident', 'Macrophage_C1q', 'Monocyte',
                                           'Dendritic_Cell', 'Myeloid_Proliferating',
                                           'Macrophage_Inflammatory', 'Macrophage_Monocyte_Derived',
                                           'Mast_Cell', 'Macrophage_TREM2'))

  myeloid.ref <- JoinLayers(myeloid.ref)
  myeloid.anchors <- FindTransferAnchors(reference = myeloid.ref, query = myeloid.obj.clean.clean,
                                         dims = 1:30, reference.reduction = 'pca')
  myeloid.predictions <- TransferData(anchorset = myeloid.anchors,
                                      refdata = myeloid.ref$Subnames, dims = 1:30)
  myeloid.obj.clean.clean <- AddMetaData(myeloid.obj.clean.clean, metadata = myeloid.predictions)

  saveRDS(myeloid.obj.clean.clean, .cache_myeloid)
  rm(myeloid.obj, myeloid.obj.clean, myeloid.ref, ref); gc()
}
idx <- match(colnames(myeloid.obj.clean.clean), colnames(xenium.obj))
.keep <- !is.na(idx)
xenium.obj$subnames[idx[.keep]] <- myeloid.obj.clean.clean$myeloid_clean_subtype[.keep]
rm(myeloid.obj.clean.clean); gc()

##############################################
#### NKT subclustering
##############################################
.cache_nkt <- file.path(.xcache, 'nkt_clean_clean.rds')
if (file.exists(.cache_nkt)) {
  message('Loading cached NKT subclustering...')
  nkt.obj.clean.clean <- readRDS(.cache_nkt)
} else {
  # NKT subclustering — full pipeline (from XeniumReanalysis.r lines 2682-2765)
  # tnk.obj.clean (T/NK cells misclassified as Myeloid by RCTD) is merged in when available.
  # If myeloid was loaded from cache it won't exist; we proceed with NKT-labeled cells only.
  nkt.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet == "NKT")
  if (exists('tnk.obj.clean')) nkt.obj <- merge(nkt.obj, tnk.obj.clean)
  nkt.obj <- JoinLayers(nkt.obj)
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(nkt.obj@meta.data), value = TRUE)
  if (length(umap_cols)) nkt.obj@meta.data[umap_cols] <- NULL
  nkt.obj[["pca"]] <- NULL; nkt.obj[["harmony"]] <- NULL
  nkt.obj[["umap"]] <- NULL; nkt.obj[["umap_orig"]] <- NULL

  nkt.obj <- FindVariableFeatures(nkt.obj)
  nkt.obj <- ScaleData(nkt.obj)
  nkt.obj <- RunPCA(nkt.obj, npcs = 30)
  nkt.obj <- RunHarmony(nkt.obj, "patient")
  nkt.obj <- RunUMAP(nkt.obj, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  nkt.obj <- RunUMAP(nkt.obj, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  nkt.obj <- FindNeighbors(nkt.obj, reduction = "harmony", dims = 1:30)
  nkt.obj <- FindClusters(nkt.obj, resolution = 0.5)

  nkt_labels <- c(
    '0'  = 'CD8_T',      '1'  = 'Fibroblast',  '2'  = 'Cardiomyocyte',
    '3'  = 'CD4_T',      '4'  = 'NK',           '5'  = 'Endothelial',
    '6'  = 'Myeloid',    '7'  = 'Neutrophil',   '8'  = 'Unknown',
    '9'  = 'Pericyte',   '10' = 'Pericyte',     '11' = 'Basophil'
  )
  nkt.obj$nkt_subtype <- unname(nkt_labels[as.character(nkt.obj$seurat_clusters)])

  nkt.obj.clean <- subset(nkt.obj, subset = nkt_subtype %in% c('CD8_T', 'CD4_T', 'NK'))
  nkt.obj.clean[["pca"]] <- NULL; nkt.obj.clean[["harmony"]] <- NULL
  nkt.obj.clean[["umap"]] <- NULL; nkt.obj.clean[["umap_orig"]] <- NULL
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(nkt.obj.clean@meta.data), value = TRUE)
  if (length(umap_cols)) nkt.obj.clean@meta.data[umap_cols] <- NULL

  nkt.obj.clean <- JoinLayers(nkt.obj.clean)
  nkt.obj.clean <- FindVariableFeatures(nkt.obj.clean)
  nkt.obj.clean <- ScaleData(nkt.obj.clean)
  nkt.obj.clean <- RunPCA(nkt.obj.clean, npcs = 30)
  nkt.obj.clean <- RunHarmony(nkt.obj.clean, "patient")
  nkt.obj.clean <- RunUMAP(nkt.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  nkt.obj.clean <- RunUMAP(nkt.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

  nkt.obj.clean.clean <- nkt.obj.clean
  saveRDS(nkt.obj.clean.clean, .cache_nkt)
  if (exists('tnk.obj.clean')) rm(tnk.obj.clean)
  rm(nkt.obj, nkt.obj.clean); gc()
}
idx <- match(colnames(nkt.obj.clean.clean), colnames(xenium.obj))
.keep <- !is.na(idx)
xenium.obj$subnames[idx[.keep]] <- nkt.obj.clean.clean$nkt_subtype[.keep]
rm(nkt.obj.clean.clean); gc()

##############################################
#### Adipo subclustering
##############################################
.cache_adipo <- file.path(.xcache, 'adipo_clean_clean.rds')
if (file.exists(.cache_adipo)) {
  message('Loading cached Adipo subclustering...')
  adipo.obj.clean.clean <- readRDS(.cache_adipo)
} else {
  adipo.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet == 'Adipo')
  adipo.obj <- JoinLayers(adipo.obj)
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(adipo.obj@meta.data), value = TRUE)
  if (length(umap_cols)) adipo.obj@meta.data[umap_cols] <- NULL
  adipo.obj[["pca"]] <- NULL; adipo.obj[["harmony"]] <- NULL
  adipo.obj[["umap"]] <- NULL; adipo.obj[["umap_orig"]] <- NULL

  adipo.obj <- FindVariableFeatures(adipo.obj)
  adipo.obj <- ScaleData(adipo.obj)
  adipo.obj <- RunPCA(adipo.obj, npcs = 30)
  adipo.obj <- RunHarmony(adipo.obj, "patient")
  adipo.obj <- RunUMAP(adipo.obj, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  adipo.obj <- RunUMAP(adipo.obj, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  adipo.obj <- FindNeighbors(adipo.obj, reduction = "harmony", dims = 1:30)
  adipo.obj <- FindClusters(adipo.obj, resolution = 0.5)

  adipo_labels <- c(
    '0'='drop_CM','1'='Adipocyte','2'='drop_FB','3'='drop_EC','4'='drop_Pericyte',
    '5'='drop_FB2','6'='drop_Myeloid','7'='Adipocyte2','8'='Adipo_Fibro',
    '9'='drop_CM2','10'='drop_Pericyte2','11'='drop_FB3','12'='drop_VSMC',
    '13'='Adipo_Vascular','14'='Adipo_Immune','15'='drop_TNK','16'='drop_Neuron',
    '17'='drop_EC2','18'='drop_CM3','19'='drop_Proliferating')
  adipo.obj$adipo_subtype <- unname(adipo_labels[as.character(adipo.obj$seurat_clusters)])

  keep_cells <- !is.na(adipo.obj$adipo_subtype) & !grepl('^drop_', adipo.obj$adipo_subtype)
  adipo.obj.clean <- subset(adipo.obj, cells = colnames(adipo.obj)[keep_cells])
  adipo.obj.clean[["pca"]] <- NULL; adipo.obj.clean[["harmony"]] <- NULL
  adipo.obj.clean[["umap"]] <- NULL; adipo.obj.clean[["umap_orig"]] <- NULL
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(adipo.obj.clean@meta.data), value = TRUE)
  if (length(umap_cols)) adipo.obj.clean@meta.data[umap_cols] <- NULL

  adipo.obj.clean <- JoinLayers(adipo.obj.clean)
  adipo.obj.clean <- FindVariableFeatures(adipo.obj.clean)
  adipo.obj.clean <- ScaleData(adipo.obj.clean)
  adipo.obj.clean <- RunPCA(adipo.obj.clean, npcs = 30)
  adipo.obj.clean <- RunHarmony(adipo.obj.clean, "patient")
  adipo.obj.clean <- RunUMAP(adipo.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  adipo.obj.clean <- RunUMAP(adipo.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  adipo.obj.clean <- FindNeighbors(adipo.obj.clean, reduction = "harmony", dims = 1:30)
  adipo.obj.clean <- FindClusters(adipo.obj.clean, resolution = 0.5)

  adipo_clean_labels <- c(
    '0'='Adipo_Mature','1'='Adipo_Quiescent','2'='Adipo_Sparse','3'='Adipo_Sparse2',
    '4'='drop_FB','5'='drop_Mixed','6'='drop_Macrophage','7'='drop_Mural',
    '8'='drop_EC','9'='Adipo_Lipogenic','10'='drop_LAM')
  adipo.obj.clean$adipo_subtype <- unname(adipo_clean_labels[as.character(adipo.obj.clean$seurat_clusters)])

  keep_clean2 <- !is.na(adipo.obj.clean$adipo_subtype) & !grepl('^drop_', adipo.obj.clean$adipo_subtype)
  adipo.obj.clean.clean <- subset(adipo.obj.clean, cells = colnames(adipo.obj.clean)[keep_clean2])
  adipo.obj.clean.clean[["pca"]] <- NULL; adipo.obj.clean.clean[["harmony"]] <- NULL
  adipo.obj.clean.clean[["umap"]] <- NULL; adipo.obj.clean.clean[["umap_orig"]] <- NULL
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(adipo.obj.clean.clean@meta.data), value = TRUE)
  if (length(umap_cols)) adipo.obj.clean.clean@meta.data[umap_cols] <- NULL

  adipo.obj.clean.clean <- JoinLayers(adipo.obj.clean.clean)
  adipo.obj.clean.clean <- FindVariableFeatures(adipo.obj.clean.clean)
  adipo.obj.clean.clean <- ScaleData(adipo.obj.clean.clean)
  adipo.obj.clean.clean <- RunPCA(adipo.obj.clean.clean, npcs = 30)
  adipo.obj.clean.clean <- RunHarmony(adipo.obj.clean.clean, "patient")
  adipo.obj.clean.clean <- RunUMAP(adipo.obj.clean.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  adipo.obj.clean.clean <- RunUMAP(adipo.obj.clean.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

  saveRDS(adipo.obj.clean.clean, .cache_adipo)
  rm(adipo.obj, adipo.obj.clean); gc()
}
idx <- match(colnames(adipo.obj.clean.clean), colnames(xenium.obj))
.keep <- !is.na(idx)
xenium.obj$subnames[idx[.keep]] <- adipo.obj.clean.clean$adipo_subtype[.keep]
rm(adipo.obj.clean.clean); gc()

##############################################
#### Mural cell subclustering (SM + PC)
##############################################
.cache_mural <- file.path(.xcache, 'mural_clean_clean.rds')
if (file.exists(.cache_mural)) {
  message('Loading cached Mural subclustering...')
  mural.obj.clean <- readRDS(.cache_mural)
} else {
  mural.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet %in% c('SM', 'PC'))
  mural.obj <- JoinLayers(mural.obj)
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(mural.obj@meta.data), value = TRUE)
  if (length(umap_cols)) mural.obj@meta.data[umap_cols] <- NULL
  mural.obj[["pca"]] <- NULL; mural.obj[["harmony"]] <- NULL
  mural.obj[["umap"]] <- NULL; mural.obj[["umap_orig"]] <- NULL

  mural.obj <- FindVariableFeatures(mural.obj)
  mural.obj <- ScaleData(mural.obj)
  mural.obj <- RunPCA(mural.obj, npcs = 30)
  mural.obj <- RunHarmony(mural.obj, "patient")
  mural.obj <- RunUMAP(mural.obj, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  mural.obj <- RunUMAP(mural.obj, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  mural.obj <- FindNeighbors(mural.obj, reduction = "harmony", dims = 1:30)
  mural.obj <- FindClusters(mural.obj, resolution = 0.5)

  mural_labels <- c(
    '0'='Pericyte_KCNJ8','1'='drop_FB','2'='VSMC_Arterial','3'='drop_CM',
    '4'='drop_EC','5'='drop_Neutrophil','6'='Pericyte','7'='VSMC',
    '8'='drop_Macrophage','9'='Pericyte_Activated','10'='Mural_Unspecified',
    '11'='Pericyte2','12'='drop_TNK','13'='drop_Unknown','14'='drop_Adipocyte',
    '15'='drop_Neuron','16'='drop_LEC','17'='drop_Proliferating')
  mural.obj$mural_subtype <- unname(mural_labels[as.character(mural.obj$seurat_clusters)])

  keep_cells <- !is.na(mural.obj$mural_subtype) & !grepl('^drop_', mural.obj$mural_subtype)
  mural.obj.clean <- subset(mural.obj, cells = colnames(mural.obj)[keep_cells])
  mural.obj.clean[["pca"]] <- NULL; mural.obj.clean[["harmony"]] <- NULL
  mural.obj.clean[["umap"]] <- NULL; mural.obj.clean[["umap_orig"]] <- NULL
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(mural.obj.clean@meta.data), value = TRUE)
  if (length(umap_cols)) mural.obj.clean@meta.data[umap_cols] <- NULL

  mural.obj.clean <- JoinLayers(mural.obj.clean)
  mural.obj.clean <- FindVariableFeatures(mural.obj.clean)
  mural.obj.clean <- ScaleData(mural.obj.clean)
  mural.obj.clean <- RunPCA(mural.obj.clean, npcs = 30)
  mural.obj.clean <- RunHarmony(mural.obj.clean, "patient")
  mural.obj.clean <- RunUMAP(mural.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  mural.obj.clean <- RunUMAP(mural.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  mural.obj.clean <- FindNeighbors(mural.obj.clean, reduction = "harmony", dims = 1:30)
  mural.obj.clean <- FindClusters(mural.obj.clean, resolution = 0.5)

  mural_clean_labels <- c(
    '0'='Pericyte_KCNJ8','1'='Pericyte','2'='VSMC_Arterial','3'='Pericyte_THBS4',
    '4'='VSMC_Synthetic','5'='VSMC','6'='Pericyte_Mixed','7'='VSMC_Inflamed',
    '8'='VSMC_Arterial2','9'='Pericyte_Quiescent','10'='Mural_Unspecified',
    '11'='VSMC_Transitional')
  mural.obj.clean$mural_subtype <- unname(mural_clean_labels[as.character(mural.obj.clean$seurat_clusters)])

  saveRDS(mural.obj.clean, .cache_mural)
  rm(mural.obj); gc()
}
idx <- match(colnames(mural.obj.clean), colnames(xenium.obj))
.keep <- !is.na(idx)
xenium.obj$subnames[idx[.keep]] <- mural.obj.clean$mural_subtype[.keep]
rm(mural.obj.clean); gc()

##############################################
#### Neuron subclustering
##############################################
.cache_neuron <- file.path(.xcache, 'neuron_clean_clean.rds')
if (file.exists(.cache_neuron)) {
  message('Loading cached Neuron subclustering...')
  neuron.obj.clean.clean <- readRDS(.cache_neuron)
} else {
  neuron.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet == 'Neuron')
  neuron.obj <- JoinLayers(neuron.obj)
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(neuron.obj@meta.data), value = TRUE)
  if (length(umap_cols)) neuron.obj@meta.data[umap_cols] <- NULL
  neuron.obj[["pca"]] <- NULL; neuron.obj[["harmony"]] <- NULL
  neuron.obj[["umap"]] <- NULL; neuron.obj[["umap_orig"]] <- NULL

  neuron.obj <- FindVariableFeatures(neuron.obj)
  neuron.obj <- ScaleData(neuron.obj)
  neuron.obj <- RunPCA(neuron.obj, npcs = 30)
  neuron.obj <- RunHarmony(neuron.obj, "patient")
  neuron.obj <- RunUMAP(neuron.obj, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  neuron.obj <- RunUMAP(neuron.obj, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  neuron.obj <- FindNeighbors(neuron.obj, reduction = "harmony", dims = 1:30)
  neuron.obj <- FindClusters(neuron.obj, resolution = 0.5)

  neuron_labels <- c(
    '0'='Schwann_Myelinating','1'='drop_CM','2'='Neuron_Autonomic',
    '3'='Schwann_Nonmyelinating','4'='drop_EC_Art','5'='drop_Neutrophil',
    '6'='drop_Macrophage','7'='drop_Pericyte','8'='drop_Pericyte2',
    '9'='drop_FB','10'='drop_TNK','11'='drop_Basophil','12'='drop_LowConf',
    '13'='drop_Unknown','14'='drop_Adipocyte','15'='drop_Proliferating',
    '16'='drop_LEC')
  neuron.obj$neuron_subtype <- unname(neuron_labels[as.character(neuron.obj$seurat_clusters)])

  neuron.obj.clean <- subset(neuron.obj,
    subset = neuron_subtype %in% c('Schwann_Myelinating','Neuron_Autonomic','Schwann_Nonmyelinating'))
  neuron.obj.clean[["pca"]] <- NULL; neuron.obj.clean[["harmony"]] <- NULL
  neuron.obj.clean[["umap"]] <- NULL; neuron.obj.clean[["umap_orig"]] <- NULL
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(neuron.obj.clean@meta.data), value = TRUE)
  if (length(umap_cols)) neuron.obj.clean@meta.data[umap_cols] <- NULL

  neuron.obj.clean <- JoinLayers(neuron.obj.clean)
  neuron.obj.clean <- FindVariableFeatures(neuron.obj.clean)
  neuron.obj.clean <- ScaleData(neuron.obj.clean)
  neuron.obj.clean <- RunPCA(neuron.obj.clean, npcs = 30)
  neuron.obj.clean <- RunHarmony(neuron.obj.clean, "patient")
  neuron.obj.clean <- RunUMAP(neuron.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  neuron.obj.clean <- RunUMAP(neuron.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  neuron.obj.clean <- FindNeighbors(neuron.obj.clean, reduction = "harmony", dims = 1:30)
  neuron.obj.clean <- FindClusters(neuron.obj.clean, resolution = 1.5)

  neuron_clean_labels <- c(
    '0'='drop_Mesenchymal','1'='drop_FB','2'='drop_LowConf',
    '3'='Schwann_Nonmyelinating','4'='drop_Ambig','5'='Neuron_Cholinergic')
  # Map known clusters; any cluster number not in the label map gets NA (will be dropped)
  neuron.obj.clean$neuron_subtype <- unname(neuron_clean_labels[as.character(neuron.obj.clean$seurat_clusters)])

  neuron.obj.clean.clean <- subset(neuron.obj.clean,
    subset = neuron_subtype %in% c('Schwann_Nonmyelinating','Neuron_Cholinergic'))

  saveRDS(neuron.obj.clean.clean, .cache_neuron)
  rm(neuron.obj, neuron.obj.clean); gc()
}
idx <- match(colnames(neuron.obj.clean.clean), colnames(xenium.obj))
.keep <- !is.na(idx)
xenium.obj$subnames[idx[.keep]] <- neuron.obj.clean.clean$neuron_subtype[.keep]
rm(neuron.obj.clean.clean); gc()

##############################################
#### Epi subclustering
##############################################
.cache_epi <- file.path(.xcache, 'epi_clean_clean.rds')
if (file.exists(.cache_epi)) {
  message('Loading cached Epi subclustering...')
  epi.obj.clean <- readRDS(.cache_epi)
} else {
  epi.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet == 'Epi')
  epi.obj <- JoinLayers(epi.obj)
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(epi.obj@meta.data), value = TRUE)
  if (length(umap_cols)) epi.obj@meta.data[umap_cols] <- NULL
  epi.obj[["pca"]] <- NULL; epi.obj[["harmony"]] <- NULL
  epi.obj[["umap"]] <- NULL; epi.obj[["umap_orig"]] <- NULL

  epi.obj <- FindVariableFeatures(epi.obj)
  epi.obj <- ScaleData(epi.obj)
  epi.obj <- RunPCA(epi.obj, npcs = 30)
  epi.obj <- RunHarmony(epi.obj, "patient")
  epi.obj <- RunUMAP(epi.obj, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  epi.obj <- RunUMAP(epi.obj, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  epi.obj <- FindNeighbors(epi.obj, reduction = "harmony", dims = 1:30)
  epi.obj <- FindClusters(epi.obj, resolution = 0.5)

  epi_labels <- c(
    '0'='Epi_Mesothelial','1'='drop_EC','2'='drop_CM','3'='drop_Myeloid',
    '4'='Epi_EPDCf','5'='drop_CM2','6'='drop_Pericyte','7'='Epi_ALDH1A2',
    '8'='drop_Pericyte2','9'='drop_VSMC','10'='drop_TNK','11'='drop_Adipocyte',
    '12'='drop_Mural','13'='Epi_Neural','14'='Epi_Mesothelial2','15'='drop_Neuron',
    '16'='drop_Proliferating','17'='drop_EC2','18'='drop_CM3')
  epi.obj$epi_subtype <- unname(epi_labels[as.character(epi.obj$seurat_clusters)])

  keep_cells <- !is.na(epi.obj$epi_subtype) & !grepl('^drop_', epi.obj$epi_subtype)
  epi.obj.clean <- subset(epi.obj, cells = colnames(epi.obj)[keep_cells])
  epi.obj.clean[["pca"]] <- NULL; epi.obj.clean[["harmony"]] <- NULL
  epi.obj.clean[["umap"]] <- NULL; epi.obj.clean[["umap_orig"]] <- NULL
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(epi.obj.clean@meta.data), value = TRUE)
  if (length(umap_cols)) epi.obj.clean@meta.data[umap_cols] <- NULL

  epi.obj.clean <- JoinLayers(epi.obj.clean)
  epi.obj.clean <- FindVariableFeatures(epi.obj.clean)
  epi.obj.clean <- ScaleData(epi.obj.clean)
  epi.obj.clean <- RunPCA(epi.obj.clean, npcs = 30)
  epi.obj.clean <- RunHarmony(epi.obj.clean, "patient")
  epi.obj.clean <- RunUMAP(epi.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  epi.obj.clean <- RunUMAP(epi.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')
  epi.obj.clean <- FindNeighbors(epi.obj.clean, reduction = "harmony", dims = 1:30)
  epi.obj.clean <- FindClusters(epi.obj.clean, resolution = 0.5)

  epi_clean_labels <- c(
    '0'='Epi_EPDCf','1'='Epi_Activated','2'='Epi_Quiescent','3'='Epi_ALDH1A2',
    '4'='Epi_Ciliated','5'='Epi_Metabolic','6'='drop_Mixed','7'='Epi_Mesothelial')
  epi.obj.clean$epi_subtype <- unname(epi_clean_labels[as.character(epi.obj.clean$seurat_clusters)])

  keep_clean <- !is.na(epi.obj.clean$epi_subtype) & !grepl('^drop_', epi.obj.clean$epi_subtype)
  epi.obj.clean <- subset(epi.obj.clean, cells = colnames(epi.obj.clean)[keep_clean])
  epi.obj.clean[["pca"]] <- NULL; epi.obj.clean[["harmony"]] <- NULL
  epi.obj.clean[["umap"]] <- NULL; epi.obj.clean[["umap_orig"]] <- NULL
  umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                    colnames(epi.obj.clean@meta.data), value = TRUE)
  if (length(umap_cols)) epi.obj.clean@meta.data[umap_cols] <- NULL

  epi.obj.clean <- JoinLayers(epi.obj.clean)
  epi.obj.clean <- FindVariableFeatures(epi.obj.clean)
  epi.obj.clean <- ScaleData(epi.obj.clean)
  epi.obj.clean <- RunPCA(epi.obj.clean, npcs = 30)
  epi.obj.clean <- RunHarmony(epi.obj.clean, "patient")
  epi.obj.clean <- RunUMAP(epi.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
  epi.obj.clean <- RunUMAP(epi.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

  saveRDS(epi.obj.clean, .cache_epi)
  rm(epi.obj); gc()
}
idx <- match(colnames(epi.obj.clean), colnames(xenium.obj))
.keep <- !is.na(idx)
xenium.obj$subnames[idx[.keep]] <- epi.obj.clean$epi_subtype[.keep]
rm(epi.obj.clean); gc()

##############################################
##############################################
#### Pseudobulk DESeq2 by subtype
##############################################
##############################################

library(DESeq2)
dir.create('./output/Figure_3/Xenium', showWarnings = FALSE, recursive = TRUE)

# Generic helper: pseudobulk DESeq2 for any subtype object
.run_pseudobulk_deseq2 <- function(obj, subtype_col, out_csv, assay_name = 'Xenium') {
  pseudo_meta <- unique(obj@meta.data[, c('patient', 'group', subtype_col)])
  pseudo_meta$patient         <- as.character(pseudo_meta$patient)
  pseudo_meta[[subtype_col]]  <- as.character(pseudo_meta[[subtype_col]])
  pseudo_meta$pseudobulk_id   <- paste0(pseudo_meta$patient, '..', pseudo_meta[[subtype_col]])
  pseudo_meta$group           <- factor(pseudo_meta$group, levels = c('NF', 'pRV', 'RVF'))
  rownames(pseudo_meta)       <- pseudo_meta$pseudobulk_id

  obj$pseudobulk_id <- paste0(
    as.character(obj$patient), '..', as.character(obj@meta.data[[subtype_col]])
  )

  pseudo_counts <- AggregateExpression(
    obj, assays = assay_name, group.by = 'pseudobulk_id',
    slot = 'counts', return.seurat = FALSE
  )[[assay_name]]
  colnames(pseudo_counts) <- gsub('-', '_', sub('^g', '', colnames(pseudo_counts)))

  cell_types    <- unique(pseudo_meta[[subtype_col]])
  deseq_results <- list()

  for (ct in cell_types) {
    ct_samples <- pseudo_meta$pseudobulk_id[pseudo_meta[[subtype_col]] == ct]
    ct_samples <- ct_samples[ct_samples %in% colnames(pseudo_counts)]
    if (length(ct_samples) < 2) next
    ct_counts  <- pseudo_counts[, ct_samples, drop = FALSE]
    ct_meta    <- pseudo_meta[ct_samples, , drop = FALSE]
    if (length(unique(ct_meta$group)) < 2) next
    keep       <- rowSums(ct_counts >= 5) >= 2
    ct_counts  <- ct_counts[keep, , drop = FALSE]
    if (nrow(ct_counts) < 10) next
    dds <- DESeqDataSetFromMatrix(countData = round(ct_counts), colData = ct_meta, design = ~ group)
    dds <- tryCatch(DESeq(dds), error = function(e) NULL)
    if (is.null(dds)) next
    for (comp in list(c('RVF', 'NF'), c('pRV', 'NF'), c('RVF', 'pRV'))) {
      comp_label <- paste0(comp[1], '_vs_', comp[2])
      res <- tryCatch(
        results(dds, contrast = c('group', comp[1], comp[2])) %>%
          as.data.frame() %>%
          tibble::rownames_to_column('gene') %>%
          dplyr::mutate(cell_type = ct, comparison = comp_label) %>%
          dplyr::arrange(padj),
        error = function(e) NULL
      )
      if (!is.null(res)) deseq_results[[paste0(ct, '_', comp_label)]] <- res
    }
  }

  all_deseq <- dplyr::bind_rows(deseq_results)
  write.csv(all_deseq, out_csv, row.names = FALSE)
  message('Saved: ', out_csv)
  invisible(all_deseq)
}

#### Pseudobulk DESeq2 by CM subtype
{
  .deseq_csv <- file.path(.x_out, 'cm_pseudobulk_deseq2.csv')
  if (!file.exists(.deseq_csv)) {
    obj <- readRDS(file.path(.xcache, 'cm_clean_clean.rds'))
    .run_pseudobulk_deseq2(obj, 'cm_subtype', .deseq_csv)
    rm(obj)
  } else { message('Using cached CM DESeq2: ', .deseq_csv) }
  gc()
}

#### Pseudobulk DESeq2 by FB subtype
{
  .deseq_csv <- file.path(.x_out, 'fb_pseudobulk_deseq2.csv')
  if (!file.exists(.deseq_csv)) {
    obj <- readRDS(file.path(.xcache, 'fb_clean_clean.rds'))
    .run_pseudobulk_deseq2(obj, 'fb_subtype', .deseq_csv)
    rm(obj)
  } else { message('Using cached FB DESeq2: ', .deseq_csv) }
  gc()
}

#### Pseudobulk DESeq2 by EC subtype
{
  .deseq_csv <- file.path(.x_out, 'ec_pseudobulk_deseq2.csv')
  if (!file.exists(.deseq_csv)) {
    obj <- readRDS(file.path(.xcache, 'ec_clean_clean.rds'))
    .run_pseudobulk_deseq2(obj, 'ec_subtype', .deseq_csv)
    rm(obj)
  } else { message('Using cached EC DESeq2: ', .deseq_csv) }
  gc()
}

#### Pseudobulk DESeq2 by Myeloid subtype
{
  .deseq_csv <- file.path(.x_out, 'myeloid_pseudobulk_deseq2.csv')
  if (!file.exists(.deseq_csv)) {
    obj <- readRDS(file.path(.xcache, 'myeloid_clean_clean.rds'))
    col <- if ('myeloid_clean_subtype' %in% colnames(obj@meta.data)) 'myeloid_clean_subtype' else 'myeloid_subtype'
    .run_pseudobulk_deseq2(obj, col, .deseq_csv)
    rm(obj)
  } else { message('Using cached Myeloid DESeq2: ', .deseq_csv) }
  gc()
}

#### Pseudobulk DESeq2 by NKT subtype
{
  .deseq_csv <- file.path(.x_out, 'nkt_pseudobulk_deseq2.csv')
  if (!file.exists(.deseq_csv)) {
    obj <- readRDS(file.path(.xcache, 'nkt_clean_clean.rds'))
    col <- grep('subtype', colnames(obj@meta.data), value = TRUE)[1]
    .run_pseudobulk_deseq2(obj, col, .deseq_csv)
    rm(obj)
  } else { message('Using cached NKT DESeq2: ', .deseq_csv) }
  gc()
}

#### Pseudobulk DESeq2 by Mural subtype
{
  .deseq_csv <- file.path(.x_out, 'mural_pseudobulk_deseq2.csv')
  if (!file.exists(.deseq_csv)) {
    obj <- readRDS(file.path(.xcache, 'mural_clean_clean.rds'))
    .run_pseudobulk_deseq2(obj, 'mural_subtype', .deseq_csv)
    rm(obj)
  } else { message('Using cached Mural DESeq2: ', .deseq_csv) }
  gc()
}

#### Pseudobulk DESeq2 by Neuron subtype
{
  .deseq_csv <- file.path(.x_out, 'neuron_pseudobulk_deseq2.csv')
  if (!file.exists(.deseq_csv)) {
    obj <- readRDS(file.path(.xcache, 'neuron_clean_clean.rds'))
    .run_pseudobulk_deseq2(obj, 'neuron_subtype', .deseq_csv)
    rm(obj)
  } else { message('Using cached Neuron DESeq2: ', .deseq_csv) }
  gc()
}

#### Pseudobulk DESeq2 by Epi subtype
{
  .deseq_csv <- file.path(.x_out, 'epi_pseudobulk_deseq2.csv')
  if (!file.exists(.deseq_csv)) {
    obj <- readRDS(file.path(.xcache, 'epi_clean_clean.rds'))
    .run_pseudobulk_deseq2(obj, 'epi_subtype', .deseq_csv)
    rm(obj)
  } else { message('Using cached Epi DESeq2: ', .deseq_csv) }
  gc()
}

#### Pseudobulk DESeq2 by Adipo subtype
{
  .deseq_csv <- file.path(.x_out, 'adipo_pseudobulk_deseq2.csv')
  if (!file.exists(.deseq_csv)) {
    obj <- readRDS(file.path(.xcache, 'adipo_clean_clean.rds'))
    .run_pseudobulk_deseq2(obj, 'adipo_subtype', .deseq_csv)
    rm(obj)
  } else { message('Using cached Adipo DESeq2: ', .deseq_csv) }
  gc()
}

##############################################
##############################################
#### Merge subclustering labels → xenium_obj_subclustered
##############################################
##############################################

.xsc_path <- file.path(.xcache, 'xenium_obj_subclustered.rds')

## Restrict xenium.obj to the real Xenium panel genes BEFORE any normalization /
## pseudobulk runs downstream. The cached subclustered object (and the upstream
## imputed object) carry imputed features from snRNA reference projection; if
## those are left in, DESeq2 size-factor estimation is computed on a hybrid
## measured+imputed matrix and the resulting log2FC / padj are unreliable.
.filter_to_panel <- function(obj) {
  .panel_csv <- './dependencies/shared/panel_genes.csv'
  if (!file.exists(.panel_csv)) {
    warning('Panel CSV missing — skipping filter; results may contain imputed features.')
    return(obj)
  }
  .pg <- read.csv(.panel_csv)
  .panel <- unique(as.character(.pg[[ if ('x' %in% names(.pg)) 'x' else ncol(.pg) ]]))
  keep <- intersect(.panel, rownames(obj))
  n_before <- nrow(obj)
  ## Seurat's subset(features=) uses NSE and (in v5) goes through method dispatch
  ## that may parse-eval feature names. Try direct row-indexing first; fall back
  ## to subset() if that fails. tryCatch surfaces the *actual* error rather than
  ## R's confusing "unexpected ')' in: ..." parse-error wrapping.
  obj <- tryCatch(
    obj[keep, ],
    error = function(e1) {
      message('[.filter_to_panel] obj[keep, ] failed: ', conditionMessage(e1))
      message('[.filter_to_panel] falling back to subset(obj, features = keep) ...')
      tryCatch(
        subset(obj, features = keep),
        error = function(e2) {
          message('[.filter_to_panel] subset() also failed: ', conditionMessage(e2))
          message('[.filter_to_panel] returning unfiltered obj — downstream pseudobulks may include imputed features.')
          obj
        }
      )
    }
  )
  message(sprintf('Filtered xenium.obj features %d -> %d (panel-only)',
                  n_before, nrow(obj)))
  obj
}

if (file.exists(.xsc_path)) {
  message('Loading cached xenium_obj_subclustered...')
  xenium.obj <- readRDS(.xsc_path)
  ## Repair the cached FOV class (missing coords_x_orientation / misc slots)
  ## up-front. This is the documented fix the validObject error itself
  ## suggests, and it keeps FOVs intact so BuildNicheAssay (Panels E/F/G)
  ## still works while feature/cell subsets (Panels B/C/D) stop crashing.
  xenium.obj <- UpdateSeuratObject(xenium.obj)
  xenium.obj <- .filter_to_panel(xenium.obj)
} else {
  xenium.obj$cell_types_subclustering <- 'Unassigned'

  subcluster_list <- list(
    list(file = file.path(.xcache, 'myeloid_clean_clean.rds'), cols = c('myeloid_clean_subtype', 'myeloid_subtype')),
    list(file = file.path(.xcache, 'nkt_clean_clean.rds'),     cols = c('nkt_subtype')),
    list(file = file.path(.xcache, 'ec_clean_clean.rds'),      cols = c('ec_subtype')),
    list(file = file.path(.xcache, 'fb_clean_clean.rds'),      cols = c('fb_subtype')),
    list(file = file.path(.xcache, 'cm_clean_clean.rds'),      cols = c('cm_subtype')),
    list(file = file.path(.xcache, 'neuron_clean_clean.rds'),  cols = c('neuron_subtype')),
    list(file = file.path(.xcache, 'mural_clean_clean.rds'),   cols = c('mural_subtype')),
    list(file = file.path(.xcache, 'epi_clean_clean.rds'),     cols = c('epi_subtype')),
    list(file = file.path(.xcache, 'adipo_clean_clean.rds'),   cols = c('adipo_subtype'))
  )

  for (entry in subcluster_list) {
    if (!file.exists(entry$file)) { cat('Skipping:', entry$file, '\n'); next }
    obj <- readRDS(entry$file)
    col <- intersect(entry$cols, colnames(obj@meta.data))[1]
    if (is.na(col)) { cat('No subtype column in:', entry$file, '\n'); next }
    sub_labels <- obj@meta.data[[col]]
    names(sub_labels) <- colnames(obj)
    shared_cells <- intersect(names(sub_labels), colnames(xenium.obj))
    xenium.obj$cell_types_subclustering[match(shared_cells, colnames(xenium.obj))] <- sub_labels[shared_cells]
    cat(col, ':', length(shared_cells), 'cells transferred\n')
    rm(obj); gc()
  }

  xenium.obj$cell_types_subclustering[is.na(xenium.obj$cell_types_subclustering)] <- 'Unassigned'
  xenium.obj$cell_types_manual <- xenium.obj$cell_types_subclustering
  cat('Unique subtypes:', length(unique(xenium.obj$cell_types_subclustering)), '\n')

  xenium.obj <- .filter_to_panel(xenium.obj)
  saveRDS(xenium.obj, file.path(.x_out, 'xenium_obj_subclustered.rds'))
}
xenium.obj$subnames <- xenium.obj$cell_types_manual

##############################################
##############################################
#### Figure 3A
##############################################
##############################################

library(Seurat)
library(hdWGCNA)
library(ggeasy)


source('./helper_scripts/spatial_functions.R')

M1 <- readRDS(.xenium_corrected_path)
M1$names <- M1$cell_type_rctd_doublet

p_3A <- PlotEmbedding(M1, group.by='cell_type_rctd_doublet', reduction='umap_orig',
                      point_size=.05, plot_under=TRUE, plot_theme=umap_theme()+NoLegend(),
                      raster_dpi=400, raster_scale=0.5)
save_figure(p_3A, 'Figure_3A.pdf', width=5, height=5)

p_3A <- PlotEmbedding(M1, group.by='patient', reduction='umap_orig',
                      point_size=.05, plot_under=TRUE, plot_theme=umap_theme()+NoLegend(),
                      raster_dpi=400, raster_scale=0.5)
save_figure(p_3A, 'Figure_3_supp_panel_A_patient.pdf', width=5, height=5)

##############################################
##############################################
#### Figure 3B-C
#### Marker genes and cross platform comparison
##############################################
##############################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(ggrastr)
  library(patchwork)
})

## Patient → group mapping (from Xenium_metadata.csv)
xen_meta_path <- './dependencies/shared/Xenium/Xenium_metadata.csv'
xen_meta      <- read.csv(xen_meta_path, row.names = 1)
patient_group <- xen_meta %>%
  dplyr::select(patient, group) %>%
  dplyr::distinct()
patient_group$patient <- as.character(patient_group$patient)

## snRNA pseudobulk DESeq2 (lineage level)
.cache_sn_deseq_c <- file.path(.xcache, 'sn_pseudobulk_lineage_deseq2.rds')
if (file.exists(.cache_sn_deseq_c)) {
  message('Loading cached snRNA lineage DESeq2 for Panel B...')
  deseq_sn_c <- readRDS(.cache_sn_deseq_c)
} else {
  sn_for_c   <- readRDS('./dependencies/shared/snRV_ref.rds')
  deseq_sn_c <- run_pseudobulk_deseq(sn_for_c, ident_col = 'Names',
                                      contrasts = list(c('RVF','NF'), c('pRV','NF'), c('RVF','pRV')),
                                      min_cells = 10)
  saveRDS(deseq_sn_c, .cache_sn_deseq_c)
  write.csv(deseq_sn_c, sub('\\.rds$', '.csv', .cache_sn_deseq_c), row.names = FALSE)
  rm(sn_for_c); gc()
}

## Xenium pseudobulk DESeq2 (lineage level, using subclustered object + cell_type_rctd_doublet)
.cache_xen_deseq_c <- file.path(.xcache, 'xenium_pseudobulk_lineage_deseq2.rds')
if (file.exists(.cache_xen_deseq_c)) {
  message('Loading cached Xenium lineage DESeq2 for Panel B...')
  deseq_xen_c <- readRDS(.cache_xen_deseq_c)
} else {
  xen_c <- if (exists('xenium.obj')) xenium.obj else readRDS(.xsc_path)
  ## Strip @images before subset — Seurat v5's subset.Seurat traverses FOVs and
  ## validObject() on the cached FOV class fails (missing coords_x_orientation/
  ## misc slots). The pseudobulk panels don't need spatial coordinates.
  xen_c@images <- list()
  xen_c <- subset(xen_c, subset = cell_types_subclustering != 'Unassigned')
  xen_c$Names <- xen_c$cell_type_rctd_doublet
  deseq_xen_c <- tryCatch(
    run_pseudobulk_deseq(xen_c, ident_col = 'Names',
                         contrasts = list(c('RVF','NF'), c('pRV','NF'), c('RVF','pRV')),
                         min_cells = 10, assay = 'Xenium'),
    error = function(e) {
      message('Xenium pseudobulk DESeq2 failed: ', conditionMessage(e)); NULL
    })
  if (!exists('xenium.obj')) rm(xen_c)
  gc()
  if (!is.null(deseq_xen_c)) {
    saveRDS(deseq_xen_c, .cache_xen_deseq_c)
    write.csv(deseq_xen_c, sub('\\.rds$', '.csv', .cache_xen_deseq_c), row.names = FALSE)
  }
}

## snRNA pseudobulk DESeq2 (sublineage level — Subnames_manual)
.cache_sn_deseq_sub <- file.path(.xcache, 'sn_pseudobulk_sublineage_deseq2.rds')
if (file.exists(.cache_sn_deseq_sub)) {
  message('Loading cached snRNA sublineage DESeq2...')
  deseq_sn_sub <- readRDS(.cache_sn_deseq_sub)
} else {
  sn_for_sub   <- readRDS('./dependencies/shared/snRV_ref.rds')
  deseq_sn_sub <- run_pseudobulk_deseq(sn_for_sub, ident_col = 'Subnames_manual',
                                        contrasts = list(c('RVF','NF'), c('pRV','NF'), c('RVF','pRV')),
                                        min_cells = 10)
  saveRDS(deseq_sn_sub, .cache_sn_deseq_sub)
  write.csv(deseq_sn_sub, sub('\\.rds$', '.csv', .cache_sn_deseq_sub), row.names = FALSE)
  rm(sn_for_sub); gc()
}

## Xenium pseudobulk DESeq2 (sublineage level — cell_types_subclustering, Unassigned removed)
.cache_xen_deseq_sub <- file.path(.xcache, 'xenium_pseudobulk_sublineage_deseq2.rds')
if (file.exists(.cache_xen_deseq_sub)) {
  message('Loading cached Xenium sublineage DESeq2...')
  deseq_xen_sub <- readRDS(.cache_xen_deseq_sub)
} else {
  xen_sub <- if (exists('xenium.obj')) xenium.obj else readRDS(.xsc_path)
  xen_sub@images <- list()  # see note above — drop FOVs for v5 subset compat
  xen_sub <- subset(xen_sub, subset = cell_types_subclustering != 'Unassigned')
  deseq_xen_sub <- tryCatch(
    run_pseudobulk_deseq(xen_sub, ident_col = 'cell_types_subclustering',
                         contrasts = list(c('RVF','NF'), c('pRV','NF'), c('RVF','pRV')),
                         min_cells = 10, assay = 'Xenium'),
    error = function(e) {
      message('Xenium sublineage pseudobulk DESeq2 failed: ', conditionMessage(e)); NULL
    })
  if (!exists('xenium.obj')) rm(xen_sub)
  gc()
  if (!is.null(deseq_xen_sub)) {
    saveRDS(deseq_xen_sub, .cache_xen_deseq_sub)
    write.csv(deseq_xen_sub, sub('\\.rds$', '.csv', .cache_xen_deseq_sub), row.names = FALSE)
  }
}

## Build plot
if (!is.null(deseq_xen_c)) {
  conc_c <- deseq_sn_c %>%
    dplyr::rename(l2fc_sn = log2FoldChange, padj_sn = padj) %>%
    dplyr::inner_join(
      deseq_xen_c %>% dplyr::rename(l2fc_xen = log2FoldChange, padj_xen = padj),
      by = c('gene', 'subtype', 'contrast')) %>%
    dplyr::filter(!is.na(l2fc_sn), !is.na(l2fc_xen)) %>%
    dplyr::mutate(
      l2fc_sn    = pmax(pmin(l2fc_sn,  10),  -10),
      l2fc_xen   = pmax(pmin(l2fc_xen, 7.5), -7.5),
      sig_both   = !is.na(padj_sn) & padj_sn < 0.05 & !is.na(padj_xen) & padj_xen < 0.05,
      sig_either = (!is.na(padj_sn) & padj_sn < 0.05) | (!is.na(padj_xen) & padj_xen < 0.05),
      lineage    = dplyr::if_else(subtype %in% c('SM', 'PC'), 'Mural', as.character(subtype)),
      concordant = sign(l2fc_sn) == sign(l2fc_xen)) %>%
    dplyr::filter(sig_either)

  r_labels <- conc_c %>%
    dplyr::group_by(contrast) %>%
    dplyr::summarise(r = round(cor(l2fc_sn, l2fc_xen, method = 'spearman', use = 'complete.obs'), 3),
                     .groups = 'drop') %>%
    dplyr::mutate(label = paste0('\u03c1 = ', r))

  p_3C <- ggplot(conc_c, aes(x = l2fc_sn, y = l2fc_xen,
                              colour = lineage, shape = concordant)) +
    ggrastr::rasterise(
      geom_point(data = subset(conc_c, !sig_both), size = PS$scatter_pt * 0.6, alpha = 0.5),
      dpi = 200) +
    geom_point(data = subset(conc_c, sig_both), size = PS$scatter_pt * 1.5, alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', linewidth = PS$geom_lw,
                colour = 'grey50') +
    geom_vline(xintercept = 0, linewidth = PS$geom_lw, colour = 'grey70') +
    geom_hline(yintercept = 0, linewidth = PS$geom_lw, colour = 'grey70') +
    scale_colour_lineage(name = 'Lineage') +
    scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 4),
                       labels = c(`TRUE` = 'Concordant', `FALSE` = 'Discordant'),
                       name = NULL) +
    geom_text(data = r_labels, aes(label = label), inherit.aes = FALSE,
              x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3,
              size = PS$text_mm, family = FONT_FAMILY) +
    facet_wrap(~ contrast, ncol = 2) +
    labs(x = 'snRNA-seq log\u2082FC', y = 'Xenium log\u2082FC',
         title = 'Cross-platform concordance') +
    theme_v52(COMP_W) +
    theme(legend.position = 'right', legend.key.size = PS$legend_key)
} else {
  p_3C <- ggplot() +
    annotate('text', x = 0.5, y = 0.5,
             label = '[Panel B: Xenium pseudobulk DESeq2 failed\n(check run_pseudobulk_deseq on xenium.obj)]',
             size = PS$text_mm, family = FONT_FAMILY, colour = 'grey50', hjust = 0.5, vjust = 0.5) +
    theme_void()
}

save_figure(p_3C, 'Figure_3C.pdf', width = 6, height = 5)


#### Figure 3B (new) — Xenium cluster marker dotplot

# FindAllMarkers on cell_type_rctd_doublet (Unassigned cells excluded)
.cache_xen_markers <- file.path(.xcache, 'xenium_corrected_findallmarkers.rds')
if (file.exists(.cache_xen_markers)) {
  message('Loading cached Xenium FindAllMarkers...')
  xen_markers <- readRDS(.cache_xen_markers)
} else {
  message('Running FindAllMarkers on Xenium corrected object...')
  ## image-stripped copy: subset.Seurat traverses FOVs and validObject() fails
  ## on the cached FOV class. FindAllMarkers needs no spatial coords.
  .xs_tmp <- xenium.obj; .xs_tmp@images <- list()
  .xen_sub <- subset(.xs_tmp, subset = cell_types_subclustering != 'Unassigned')
  rm(.xs_tmp)
  Idents(.xen_sub) <- 'cell_type_rctd_doublet'
  xen_markers <- FindAllMarkers(.xen_sub, only.pos = TRUE, min.pct = 0.25,
                                logfc.threshold = 0.25)
  rm(.xen_sub); gc()
  saveRDS(xen_markers, .cache_xen_markers)
}

# Biologically curated markers — all in xenium_genes.txt panel.
# Cross-referenced against FindAllMarkers (xenium_corrected_findallmarkers.rds):
# where FAM top hits were biologically implausible (e.g. TNNI3 in Myeloid,
# MYBPC3 in Epi), canonical markers were used instead.
# Y-axis: alphabetical, Adipo at bottom. Genes ordered left→right to match,
# so marker blocks run diagonal bottom-left (Adipo) to top-right (SM).
.ct_levels <- c('Adipo', 'CM', 'EC', 'Endo', 'Epi', 'FB', 'LEC', 'Myeloid', 'NKT', 'Neuron', 'PC', 'SM')

.dot_genes <- c(
  # Adipo
  'ADIPOQ',
  # CM 
  'TNNT2',
  # EC 
  'VWF',
  # Endo 
  'TMEM100','BMX',
  # Epi 
  'WT1','UPK3B',
  # FB 
  'DCN',
  # LEC 
  'PROX1',
  # Myeloid 
  'CD68',
  # NKT 
  'CD3E',
  # Neuron 
  'SCN7A',
  # PC 
  'PDGFRB',
  # SM 
  'MYH11'
)

## image-stripped copy (Panel 3B dotplot needs no spatial coords; avoids the
## Seurat v5 FOV-validObject crash in subset.Seurat).
.xs_tmp <- xenium.obj; .xs_tmp@images <- list()
.xen_sub <- subset(
  .xs_tmp,
  subset = cell_types_subclustering != 'Unassigned'
)
rm(.xs_tmp)
Idents(.xen_sub) <- factor(.xen_sub$cell_type_rctd_doublet, levels = .ct_levels)
p_3B_dot <- DotPlot(.xen_sub,
                    features  = .dot_genes,
                    col.min   = 0, col.max = 2,
                    dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = PS$axis_text),
        axis.text.y = element_text(size = PS$axis_text),
        legend.position  = 'bottom',
        legend.direction = 'horizontal',
        legend.box       = 'horizontal')
save_figure(p_3B_dot, 'Figure_3B.pdf', width = 6, height = 4.5)
rm(.xen_sub); gc()


##############################################
##############################################
#### Figure 3D
#### All 9 Xenium tissue slices — cell-type spatial map
#### Patients tiled within each disease-state panel using offset_coords()
##############################################
##############################################

# offset_coords is defined once in the helpers section above (row_gap_um = 10000 default).

# Use xenium.obj if still in memory (has subnames from subclustering),
# otherwise load the cached subclustered object.
.src_3b <- if (exists('xenium.obj')) xenium.obj else readRDS(.xsc_path)

meta_3b <- data.frame(
  x         = .src_3b$x_centroid,
  y         = -.src_3b$y_centroid,   # negate Y so tissue appears right-side up
  group     = .src_3b$group,
  patient   = .src_3b$patient,
  cell_type = .src_3b$cell_type_rctd_doublet,
  stringsAsFactors = FALSE
)
rm(.src_3b); gc(verbose = FALSE)

meta_3b$group <- factor(meta_3b$group, levels = c('NF', 'pRV', 'RVF'))
meta_3b <- offset_coords(meta_3b, row_gap_um = 10000)  # <-- µm offset for top tissue in 3-patient groups

# Colour palette — map $names abbreviations to canonical lineage colours.
# lineage_pal uses full names; cell types absent from lineage_pal (Endo, PC)
# borrow colours from lineage types not present in Xenium (B, Mast).
.ct_pal <- c(
  CM      = unname(lineage_pal['CM']),
  FB      = unname(lineage_pal['FB']),
  EC      = unname(lineage_pal['EC']),
  Myeloid = unname(lineage_pal['Myeloid']),
  NKT     = unname(lineage_pal['NKT']),
  SM      = unname(lineage_pal['Mural']),      # vascular smooth muscle → Mural colour
  PC      = unname(lineage_pal['Mast']),       # pericyte → borrows Mast colour
  Adipo   = unname(lineage_pal['Adipocyte']),
  Endo    = unname(lineage_pal['B']),          # endocardial → borrows B colour
  LEC     = unname(lineage_pal['Lymphatic']),
  Neuron  = unname(lineage_pal['Neuronal']),
  Epi     = unname(lineage_pal['Epicardial'])
)

# ── Figure 3D: individual tissue tile PDFs ─────────────────────────────────
# Each patient is its own rasterised PDF for independent Illustrator assembly.
# Physical scale: 3000 µm per output inch (consistent across all tiles).

suppressPackageStartupMessages(library(cowplot))

.ums_per_inch <- 3000
.tile_dir     <- file.path(V52_FIG_DIR, 'Figure_3C_tiles')
dir.create(.tile_dir, showWarnings = FALSE, recursive = TRUE)

# Legend PDF
.p_leg <- ggplot(meta_3b, aes(x = x, y = y, color = cell_type)) +
  geom_point(size = 0) +
  scale_color_manual(values = .ct_pal, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 5), ncol = 2)) +
  theme_void() +
  theme(legend.position = 'center',
        legend.text     = element_text(size = 14),
        legend.key.size = unit(0.55, 'cm'))
ggsave(file.path(.tile_dir, 'Figure_3D_legend.pdf'),
       plot   = ggdraw(get_legend(.p_leg)),
       width  = 3.5, height = 3, device = cairo_pdf)
message('Saved: Figure_3D_legend.pdf')

# Scale bar PDF — 2 mm at chosen scale
.p_sb <- ggplot() +
  annotate('segment', x = 0, xend = 1, y = 0, yend = 0,
           linewidth = 0.8, color = 'black') +
  annotate('text', x = 0.5, y = 0.2, label = '2 mm', size = 3, hjust = 0.5) +
  xlim(-0.05, 1.05) + ylim(-0.3, 0.5) +
  theme_void()
ggsave(file.path(.tile_dir, 'Figure_3D_scalebar.pdf'),
       plot   = .p_sb,
       width  = 2000 / .ums_per_inch,   # 2000 µm at scale
       height = 0.3, device = cairo_pdf)
message('Saved: Figure_3D_scalebar.pdf')

# Per-patient tile PDFs
for (.pt in unique(meta_3b$patient)) {
  .df   <- meta_3b[meta_3b$patient == .pt, ]
  .grp  <- as.character(.df$group[1])
  .xr   <- range(.df$x)
  .yr   <- range(.df$y)
  .w_in <- diff(.xr) / .ums_per_inch
  .h_in <- diff(.yr) / .ums_per_inch
  .p    <- ggplot(.df, aes(x = x, y = y, color = cell_type)) +
    ggrastr::rasterise(geom_point(size = 0.35, stroke = 0), dpi = 400) +
    scale_color_manual(values = .ct_pal, name = NULL) +
    coord_fixed(xlim = .xr, ylim = .yr) +
    theme_void() +
    theme(legend.position = 'none',
          plot.margin     = margin(0, 0, 0, 0))
  .fname <- sprintf('Figure_3D_%s_%s.pdf', .grp, .pt)
  ggsave(file.path(.tile_dir, .fname),
         plot = .p, width = .w_in, height = .h_in, device = cairo_pdf)
  message('Saved tile: ', .fname, ' (', round(.w_in, 2), ' × ', round(.h_in, 2), ' in)')
}

# ── Composite Figure_3D.pdf — each tissue placed via draw_plot (movable) ──
# Per-patient origin-normalized data (each tissue starts at (0,0))
.per_pt <- split(meta_3b, meta_3b$patient)
.per_pt <- lapply(.per_pt, function(d) {
  d$x <- d$x - min(d$x); d$y <- d$y - min(d$y); d
})

# Per-patient ggplot — each rasterised independently → movable grob
.pt_info  <- list(); .pt_plots <- list()
for (.pt in names(.per_pt)) {
  d <- .per_pt[[.pt]]
  .pt_info[[.pt]]  <- list(group = as.character(d$group[1]),
                            w_um  = diff(range(d$x)),
                            h_um  = diff(range(d$y)))
  .pt_plots[[.pt]] <- ggplot(d, aes(x = x, y = y, color = cell_type)) +
    ggrastr::rasterise(geom_point(size = 0.35, stroke = 0), dpi = 400) +
    scale_color_manual(values = .ct_pal, name = NULL) +
    coord_fixed(expand = FALSE) +
    theme_void() +
    theme(legend.position = 'none', plot.margin = margin(0, 0, 0, 0))
}

# Per-group layout (matches offset_coords): 2 largest side-by-side on bottom,
# smallest centred above.
.gap_um     <- 500       # horizontal gap within bottom row
.row_gap_um <- 1500      # vertical gap between bottom row and top tissue
.grp_gap_um <- 5000      # vertical gap between groups
.grp_levels <- c('NF', 'pRV', 'RVF')

.layout_group <- function(.g) {
  pts <- names(.pt_info)[sapply(.pt_info, function(x) x$group == .g)]
  bbs <- .pt_info[pts]
  pts <- pts[order(-sapply(bbs, function(b) b$w_um * b$h_um))]
  bbs <- bbs[pts]
  b1 <- bbs[[1]]; b2 <- bbs[[2]]; b3 <- bbs[[3]]
  bot_w <- b1$w_um + .gap_um + b2$w_um
  bot_h <- max(b1$h_um, b2$h_um)
  tiles <- list(
    list(x = 0,                     y = 0,                   w = b1$w_um, h = b1$h_um),
    list(x = b1$w_um + .gap_um,     y = 0,                   w = b2$w_um, h = b2$h_um),
    list(x = (bot_w - b3$w_um) / 2, y = bot_h + .row_gap_um, w = b3$w_um, h = b3$h_um)
  )
  names(tiles) <- pts
  list(tiles  = tiles,
       canvas = list(w = bot_w, h = bot_h + .row_gap_um + b3$h_um),
       top_pt = pts[3], bot_h = bot_h)
}
.layouts <- setNames(lapply(.grp_levels, .layout_group), .grp_levels)

.canvas_w_um <- max(sapply(.layouts, function(l) l$canvas$w))
.canvas_h_um <- sum(sapply(.layouts, function(l) l$canvas$h)) +
                (length(.grp_levels) - 1) * .grp_gap_um

# Group Y origins (bottom-up: RVF at y=0, NF on top)
.grp_y <- numeric(length(.grp_levels)); names(.grp_y) <- .grp_levels
.cur_y <- 0
for (.g in rev(.grp_levels)) {
  .grp_y[.g] <- .cur_y
  .cur_y     <- .cur_y + .layouts[[.g]]$canvas$h + .grp_gap_um
}

# Build composite: each tissue = independent draw_plot call (movable)
.canvas <- ggdraw()
for (.g in .grp_levels) {
  lay    <- .layouts[[.g]]
  grp_yo <- .grp_y[.g]
  grp_xo <- (.canvas_w_um - lay$canvas$w) / 2
  for (.pt in names(lay$tiles)) {
    t <- lay$tiles[[.pt]]
    .canvas <- .canvas + draw_plot(
      .pt_plots[[.pt]],
      x      = (grp_xo + t$x) / .canvas_w_um,
      y      = (grp_yo + t$y) / .canvas_h_um,
      width  = t$w / .canvas_w_um,
      height = t$h / .canvas_h_um
    )
  }
  .canvas <- .canvas + draw_label(
    .g,
    x = 0.01,
    y = (grp_yo + lay$canvas$h) / .canvas_h_um,
    hjust = 0, vjust = 1,
    fontfamily = FONT_FAMILY, fontface = 'bold', size = 22
  )
}

# Legend + scale bar → RVF top-right white space
.rvf_lay  <- .layouts[['RVF']]
.rvf_yo   <- .grp_y[['RVF']]
.rvf_xo   <- (.canvas_w_um - .rvf_lay$canvas$w) / 2
.top_t    <- .rvf_lay$tiles[[.rvf_lay$top_pt]]
.leg_x_um <- .rvf_xo + .top_t$x + .top_t$w + 1000
.leg_y_um <- .rvf_yo + .rvf_lay$bot_h + .row_gap_um
.leg_w_um <- (.rvf_xo + .rvf_lay$canvas$w) - .leg_x_um
.leg_h_um <- .top_t$h * 0.70
.canvas <- .canvas + draw_plot(
  get_legend(.p_leg),
  x      = .leg_x_um / .canvas_w_um,
  y      = (.leg_y_um + .top_t$h - .leg_h_um) / .canvas_h_um,
  width  = .leg_w_um / .canvas_w_um,
  height = .leg_h_um / .canvas_h_um
)

# Scale bar below the legend
.sb_y_um  <- .leg_y_um + .top_t$h - .leg_h_um - 1500
.sb_x0_um <- .leg_x_um + (.leg_w_um - 2000) / 2
.sb_x1_um <- .sb_x0_um + 2000
.canvas <- .canvas +
  draw_line(x = c(.sb_x0_um, .sb_x1_um) / .canvas_w_um,
            y = c(.sb_y_um,  .sb_y_um ) / .canvas_h_um,
            size = 0.8, color = 'black') +
  draw_label('2 mm',
             x = ((.sb_x0_um + .sb_x1_um) / 2) / .canvas_w_um,
             y = (.sb_y_um - 800) / .canvas_h_um,
             hjust = 0.5, vjust = 1, size = 14,
             fontfamily = FONT_FAMILY)

.fig_w_in <- .canvas_w_um / .ums_per_inch
.fig_h_in <- .canvas_h_um / .ums_per_inch
save_figure(.canvas, 'Figure_3D.pdf',
            width = .fig_w_in, height = .fig_h_in)
message('Assembled Figure_3D.pdf: ', round(.fig_w_in, 2),
        ' × ', round(.fig_h_in, 2), ' in (9 tissues independently movable)')

rm(meta_3b, .ct_pal, .tile_dir, .ums_per_inch,
   .p_leg, .p_sb, .pt, .df, .grp, .xr, .yr, .w_in, .h_in, .p, .fname,
   .per_pt, .pt_info, .pt_plots,
   .gap_um, .row_gap_um, .grp_gap_um, .grp_levels, .layouts,
   .canvas_w_um, .canvas_h_um, .grp_y, .cur_y, .canvas,
   .rvf_lay, .rvf_yo, .rvf_xo, .top_t,
   .leg_x_um, .leg_y_um, .leg_w_um, .leg_h_um,
   .sb_y_um, .sb_x0_um, .sb_x1_um,
   .fig_w_in, .fig_h_in, .layout_group)
gc(verbose = FALSE)




##############################################
##############################################
#### Niche computation — dependencies for Panels E, F, G
#### BuildNicheAssay + MiniBatchKmeans (k=15)
##############################################
##############################################

# Run niche computation on xenium.obj restricted to the gene panel used by M1.
.shared_genes <- intersect(rownames(M1), rownames(xenium.obj))
message('Niche computation: using xenium.obj subset to ', length(.shared_genes),
        ' shared genes (M1 panel).')
M1 <- subset(xenium.obj, features = .shared_genes)
M1 <- UpdateSeuratObject(M1)
M1$names <- if ('cell_type_rctd_doublet' %in% colnames(M1@meta.data)) {
  as.character(M1$cell_type_rctd_doublet)
} else {
  as.character(M1$cell_types_manual)
}

# Use the manual subclustered annotations carried by xenium.obj; fall back to
# $names if they are missing. Then collapse CM subtypes (CM_* → CM) to prevent
# patient-specific CM subtype variation from dominating niches.
if (!('cell_types_manual' %in% colnames(M1@meta.data))) {
  M1$cell_types_manual <- M1$names
}
M1$cell_types_manual    <- as.character(M1$cell_types_manual)
M1$cell_types_cm_collapsed <- M1$cell_types_manual
cm_mask <- grepl('^CM', M1$cell_types_cm_collapsed) & M1$cell_types_cm_collapsed != 'CM'
M1$cell_types_cm_collapsed[cm_mask] <- 'CM'
message('Collapsed ', sum(cm_mask), ' CM-subtype cells into CM; ',
        length(unique(M1$cell_types_cm_collapsed)), ' unique types after collapse.')

# Exclude Unassigned cells from niche assay (same as XeniumReanalysis.r)
# Capture clean per-patient M1 barcodes BEFORE SplitObject so we can restore
# them — Seurat 5 SplitObject on layered Xenium data renames cells via
# make.unique (e.g. "1_1" → "1_1_1_1"), which would propagate into the niche
# assay and break the M1[['niche_broad']] assignment downstream.
.m1_keep_mask <- M1$cell_types_cm_collapsed != 'Unassigned'
.m1_patient   <- as.character(M1$patient)
.m1_cells_pre <- split(colnames(M1)[.m1_keep_mask], .m1_patient[.m1_keep_mask])

M1_list <- SplitObject(subset(M1, cell_types_cm_collapsed != 'Unassigned'), split.by = 'patient')

# Restore original (clean) M1 barcodes per patient. Cell ORDER is preserved
# by subset + SplitObject, so position-wise rename is safe.
for (i in seq_along(M1_list)) {
  pid <- names(M1_list)[i]
  expected <- .m1_cells_pre[[pid]]
  if (ncol(M1_list[[i]]) != length(expected)) {
    stop(sprintf('Patient %s: SplitObject yielded %d cells but M1 (filtered) has %d.',
                 pid, ncol(M1_list[[i]]), length(expected)))
  }
  M1_list[[i]] <- RenameCells(M1_list[[i]], new.names = expected)
}

# BuildNicheAssay with collapsed CM labels (niches.k=10, neighbors.k=100)
.cache_niche_assay <- './output/Figure_3/fig3_niche_assay_cache_cm_xenobj_v4.rds'
if (file.exists(.cache_niche_assay)) {
  message('Loading cached BuildNicheAssay (collapsed CM) results...')
  M1_list <- readRDS(.cache_niche_assay)
} else {
  for ( i in seq_along(M1_list) ){ message('BuildNicheAssay patient ', i)
    object = M1_list[[i]]
    fov = Images(object)[1]
    object <- BuildNicheAssay(object = object, fov = fov,
                              group.by = 'cell_types_cm_collapsed',
                              niches.k = 10, neighbors.k = 100)
    M1_list[[i]] = object
  }
  saveRDS(M1_list, .cache_niche_assay)
}


# Concatenate neighbors across all samples
col.unique <- lapply(seq_along(M1_list), function(i) {
  colnames(t(M1_list[[i]][['niche']]@counts))
})
shared_feats <- Reduce(intersect, col.unique)

DAT.counts <- lapply(seq_along(M1_list), function(i) {
  t(M1_list[[i]][['niche']]@counts)[, shared_feats]
})
DAT.counts <- do.call('rbind', DAT.counts)

DAT.data <- lapply(seq_along(M1_list), function(i) {
  t(M1_list[[i]][['niche']]@data)[, shared_feats]
})
DAT.data <- do.call('rbind', DAT.data)

DAT.scale.data <- lapply(seq_along(M1_list), function(i) {
  t(M1_list[[i]][['niche']]@scale.data)[, shared_feats]
})
DAT.scale.data <- do.call('rbind', DAT.scale.data)

# Pad niche matrices with zeros for Unassigned cells so dimensions match M1.
all_cells     <- colnames(M1)
niche_cells   <- rownames(DAT.counts)
missing_cells <- setdiff(all_cells, niche_cells)

# Sanity check: niche barcodes must be a subset of M1.
.bad <- setdiff(niche_cells, all_cells)
if (length(.bad) > 0) {
  stop(sprintf('%d niche cells are not in M1 (e.g. %s). M1 cell sample: %s. ',
               length(.bad), paste(head(.bad, 3), collapse = ', '),
               paste(head(all_cells, 3), collapse = ', ')),
       'JoinLayers should have prevented this — check that M1 was JoinLayers-ed before SplitObject.')
}

if (length(missing_cells) > 0) {
  zero_pad <- matrix(0, nrow = length(missing_cells), ncol = ncol(DAT.counts),
                     dimnames = list(missing_cells, colnames(DAT.counts)))
  DAT.counts.full <- rbind(DAT.counts, zero_pad)[all_cells, ]
  DAT.data.full   <- rbind(DAT.data,   zero_pad)[all_cells, ]
  DAT.scale.full  <- rbind(DAT.scale.data, zero_pad)[all_cells, ]
} else {
  DAT.counts.full <- DAT.counts[all_cells, ]
  DAT.data.full   <- DAT.data[all_cells, ]
  DAT.scale.full  <- DAT.scale.data[all_cells, ]
}

niche.assay              <- CreateAssayObject(counts = t(DAT.counts.full))
niche.assay@data         <- t(DAT.data.full)
niche.assay@scale.data   <- t(DAT.scale.full)

M1[['niche_broad']]      <- niche.assay
DefaultAssay(M1)         <- 'niche_broad'

niches.k.range = 12:16



# CLR transform niche proportions to amplify differences in rare cell type fractions
# (raw proportions are dominated by the majority cell type, collapsing kmeans clusters).
# Drop the CM column before CLR + kmeans: CMs are the heart's background tissue and
# add no spatial discrimination — niches are defined by what surrounds the CMs
# (stromal/immune/vascular composition). Removing CM frees centroids that would
# otherwise collapse onto CM-density variation, giving rare niches proper resolution.
.cm_col_idx <- which(colnames(DAT.data) == 'CM')
if (length(.cm_col_idx) != 1) stop("Expected exactly one 'CM' column in DAT.data; found ", length(.cm_col_idx))
DAT.clr <- DAT.data[, -.cm_col_idx, drop = FALSE] + 1e-6  # small pseudocount to avoid log(0)
DAT.clr <- log(DAT.clr) - rowMeans(log(DAT.clr))
DAT.clr <- scale(DAT.clr)  # z-score after CLR for balanced kmeans distances

.cache_niche_kmeans <- './output/Figure_3/fig3_niche_kmeans_cache_cm_xenobj_clr_noCM.rds'
if (file.exists(.cache_niche_kmeans)) {
  message('Loading cached MiniBatchKmeans niche clustering (CLR)...')
  res.clusters <- readRDS(.cache_niche_kmeans)
} else {
  res.clusters = data.frame(rownames = rownames(DAT.clr))

  for ( k in niches.k.range ){ message("k=", k)
      newCol = paste0("kmeans_", k)
      km_mb = ClusterR::MiniBatchKmeans(
          "data" = DAT.clr
          , "clusters" = k
          , "batch_size" = 20
          , "num_init" = 20
          , "max_iters" = 100
          , "init_fraction" = 0.2
          , "initializer" = "kmeans++"
          , "early_stop_iter" = 10
          , "verbose" = F
          , "CENTROIDS" = NULL
          , "tol" = 1e-04
          , "tol_optimal_init" = 0.3
          , "seed" = 1
      )
      res.clusters[,newCol] = ClusterR::predict_MBatchKMeans(
          "data" = DAT.clr
          , "CENTROIDS" = km_mb$centroids
      )
      res.clusters[,newCol] = as.factor( res.clusters[,newCol] )
  }
  saveRDS(res.clusters, .cache_niche_kmeans)
}


# Update object — promote barcode column to row names before subsetting
rownames(res.clusters) <- res.clusters$rownames
res.clusters.ordered <- res.clusters[colnames(M1), setdiff(colnames(res.clusters), 'rownames')]


write.table(res.clusters.ordered,'./output/Figure_3/Xenium/Niche_bulk_clusters.csv',sep=',')
# Use the freshly-computed clusters (the shared/ file may have stale barcodes from a prior run)

for(i in colnames(res.clusters.ordered))
  {
    M1[[i]] = 0
    M1[[i]][rownames(res.clusters.ordered),] = res.clusters.ordered[,i]

  }

#Subset to assigned cells
M2 <- subset(M1,cell_types_cm_collapsed != 'Unassigned')

niche.patient <- table(M2$kmeans_15,M2$patient)
niche.patient <- t(t(niche.patient) / colSums(niche.patient))
niche.patient <- data.frame(niche.patient)  # Var1=cluster, Var2=patient, Freq
## Dynamic patient -> disease group (robust to any patient/cluster count;
## the old hardcoded 9-patient x replicate(15) vector broke when kmeans
## produced a different number of clusters/patients).
.pg <- unique(data.frame(patient = as.character(M2$patient),
                         group   = as.character(M2$group)))
.pg <- setNames(.pg$group, .pg$patient)
niche.patient$niche <- factor(.pg[as.character(niche.patient$Var2)],
                              levels = c('NF','pRV','RVF'))

ggplot(niche.patient,aes(Var1,Freq,color = niche))+geom_boxplot(linewidth = PS$geom_lw) + theme_classic() + scale_colour_disease() + scale_y_touch()



niche.names <- table(M2$kmeans_15, M2$cell_types_cm_collapsed)
niche.names <- niche.names[rownames(niche.names) != '0', ]   # drop Unassigned placeholder
niche.names <- niche.names / rowSums(niche.names)

niche.names_t <- table(M2$kmeans_15, M2$cell_types_cm_collapsed)
niche.names_t <- niche.names_t[rownames(niche.names_t) != '0', ]
niche.names_t <- t(t(niche.names_t) / colSums(niche.names_t))

write.csv(niche.names,   './output/Figure_3/Xenium/niche_kmeans15_composition_rowprop.csv')
write.csv(niche.names_t, './output/Figure_3/Xenium/niche_kmeans15_composition_colprop.csv')

Idents(M2) <- "kmeans_15"


library(scCustomize)

# Niche labels: k=15 kmeans on CLR-transformed collapsed-CM niche assay.
# Cluster '0' is the Unassigned placeholder for cells padded with zero counts
# (cells outside niche.assay's host set). Names follow XeniumReanalysis.r:7889.
# Assigned by auditing niche_kmeans15_composition_rowprop.csv (top fractions
# per cluster); duplicate labels collapse niches with shared dominant biology.
## DATA-DRIVEN niche assignment.
## MiniBatchKmeans cluster *numbering* is not stable across runs (upstream
## Seurat 5 / RCTD / UpdateSeuratObject / cell-order changes shift centroids),
## so a hardcoded cluster# -> name map silently scrambles every label and the
## dominant myocardium ends up mislabeled. Assign each cluster from its actual
## cell-type composition (niche.names rowprop) + a size anchor so the heart's
## background tissue is always 'Myocardium' (published Panel G: Myo ~= 100%).
.gfrac <- function(fr, pat) { ix <- grep(pat, names(fr)); if (length(ix)==0) 0 else sum(fr[ix], na.rm=TRUE) }
.get1  <- function(fr, nm)  if (nm %in% names(fr)) as.numeric(fr[nm]) else 0
## size_frac = this cluster's share of all assigned cells. Per published
## Panel G the myocardial subtypes (Stressed/Fibro/Cap/Compact/Remod-myo)
## are RARE focal niches; the bulk CM tissue is plain 'Myocardium' (~100%).
## So a CM-dominant cluster only earns a subtype label if it is SMALL
## (size_frac < .myo_subtype_max); large CM clusters are background Myo.
.myo_subtype_max <- 0.05
.assign_niche <- function(fr, size_frac = 0) {
  cm <- .get1(fr,'CM'); cap <- .get1(fr,'Capillary_EC'); art <- .get1(fr,'Arterial_EC')
  ven <- .get1(fr,'Venous_EC'); vsmc <- .gfrac(fr,'VSMC'); adipo <- .gfrac(fr,'^Adipo')
  epi <- .get1(fr,'Epi'); macq <- .get1(fr,'Mac_C1q'); fb <- .gfrac(fr,'^FB')
  fb_nox4 <- .get1(fr,'FB_NOX4'); fb_matri <- .get1(fr,'FB_Matrifibrocyte')
  fb_act <- .get1(fr,'FB_Activated'); fb_homeo <- .get1(fr,'FB_Homeostatic')
  fb_stress <- .get1(fr,'FB_Stress')
  ## genuinely non-myocardial niches (size-independent)
  if (vsmc >= 0.15)                              return('Arterial')
  if (vsmc >= 0.06 && art >= 0.02 && cap < 0.04) return('Peri-arteriole')
  if (adipo >= 0.10 && adipo >= cm) {
    infl <- (macq >= 0.04 || fb_stress >= 0.04); vasc <- (ven >= 0.08)
    if (infl && vasc) return('Inflamed adipose-venous stroma')
    if (vasc)         return('Adipose-vascular niche')
    if (infl)         return('Inflamed adipose stroma')
    return('Adipose stroma')
  }
  if (fb_homeo >= 0.12 && epi >= 0.05)           return('Sub-epicardial niche')
  if (cm < 0.30 && fb >= 0.10)                   return('Fibrotic remodeling stroma')
  ## CM-dominant: large = background Myo; small = a focal subtype niche
  if (size_frac < .myo_subtype_max) {
    if (cap      >= 0.11)                        return('Capillary-rich myocardium')
    if (fb_nox4  >= 0.04)                        return('Stressed myocardium')
    if (fb_matri >= 0.05)                        return('Fibrotic myocardium')
    if (fb_act   >= 0.04)                        return('Fibrotic remodeling myocardium')
    if (cm >= 0.55 && cap < 0.045)               return('Compact myocardium')
  }
  'Myocardium'
}
.clust_ids   <- rownames(niche.names)
.clust_sizes <- table(M2$kmeans_15)[.clust_ids]
.nz          <- setdiff(names(.clust_sizes), '0')
.tot_cells   <- sum(.clust_sizes[.nz])
niche_labels_cm15 <- setNames(
  vapply(.clust_ids, function(cl)
           .assign_niche(niche.names[cl, ],
                         size_frac = as.numeric(.clust_sizes[cl]) / .tot_cells),
         character(1)),
  .clust_ids)
niche_labels_cm15['0'] <- 'Unassigned'
## Anchor: single largest cluster is always the heart background = 'Myocardium'.
.biggest <- names(which.max(.clust_sizes[.nz]))
if (length(.biggest) == 1) niche_labels_cm15[.biggest] <- 'Myocardium'
message('Data-driven niche labels (cluster -> niche, n cells):')
for (cl in names(sort(.clust_sizes, decreasing = TRUE))) {
  if (cl == '0') next
  message(sprintf('  k%-3s n=%-8d -> %s', cl, .clust_sizes[cl], niche_labels_cm15[cl]))
}
Idents(M1) <- 'kmeans_15'


new.cluster.ids <- niche_labels_cm15[levels(M1)]
names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)
M1$niche_manual <- M1@active.ident

## Export per-cell barcode -> niche label (with kmeans cluster, patient, group)
.bc_niche <- data.frame(
  barcode      = colnames(M1),
  kmeans_15    = as.character(M1$kmeans_15),
  niche        = as.character(M1$niche_manual),
  patient      = if ('patient' %in% colnames(M1@meta.data)) as.character(M1$patient) else NA,
  group        = if ('group'   %in% colnames(M1@meta.data)) as.character(M1$group)   else NA,
  stringsAsFactors = FALSE)
write.csv(.bc_niche, './output/Figure_3/Xenium/niche_barcode_labels.csv',
          row.names = FALSE)
message('Wrote ', nrow(.bc_niche), ' barcode->niche rows to ',
        './output/Figure_3/Xenium/niche_barcode_labels.csv')

write.table(M1@meta.data,'./output/Figure_3/Xenium/metadata.csv',sep=',')


niche.celltype <- table(M1$niche_manual,M1$names)
write.table(niche.celltype,'./output/Figure_3/Xenium/niche_celltype.csv',sep=',')

##############################################
##############################################
#### Niche annotation (XeniumReanalysis lines 7490-8007)
##############################################
##############################################

library(reshape2)
library(ggpubr)

# --- 1. Heatmap of niche composition (row-normalized, collapsed to unique labels) ---
# Aggregate kmeans_15 rows that share the same niche label (e.g. clusters 3+10 → Peri-arterial)
niche_label_for_agg <- niche_labels_cm15[rownames(niche.names)]
niche.names_plot <- rowsum(niche.names, group = niche_label_for_agg)
# Drop 'Unassigned' (cluster 0 padding) row before normalizing so heatmap only shows real niches
niche.names_plot <- niche.names_plot[rownames(niche.names_plot) != 'Unassigned', , drop = FALSE]
# Hardcoded niche order (highest frequency first) — shared between Figure 3F
# heatmap y-axis and the composition bar chart x-axis.
.niche_freq_order <- c(
  'Myocardium',
  'Fibrotic remodeling stroma',
  'Capillary-rich myocardium',
  'Arterial',
  'Fibrotic remodeling myocardium',
  'Fibrotic myocardium',
  'Adipose-vascular niche',
  'Stressed myocardium',
  'Sub-epicardial niche',
  'Adipose stroma',
  'Inflamed adipose stroma',
  'Peri-arteriole',
  'Inflamed adipose-venous stroma',
  'Compact myocardium'
)
niche.names_plot <- niche.names_plot / rowSums(niche.names_plot)

niche_df <- melt(niche.names_plot)
colnames(niche_df) <- c("niche", "cell_type", "fraction")
niche_df$niche <- factor(niche_df$niche, levels = rev(.niche_freq_order))  # most frequent at top

## Supplementary view of niche composition: heatmap of fraction-of-niche
## per (cell_type × niche).  Final Panel F is the stacked bar produced below.
p_3F_heatmap <- ggplot(niche_df, aes(x = cell_type, y = niche, fill = fraction)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "darkred", midpoint = 0.15, name = "Fraction") +
  geom_text(aes(label = ifelse(fraction > 0.05, sprintf("%.0f%%", fraction * 100), "")), size = 2.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 10)) +
  labs(x = NULL, y = "Niche", title = "Cell type composition per niche (collapsed CM, 12 niches)")
save_figure(p_3F_heatmap, 'Figure_3_supp_panel_F_heatmap.pdf', width = 12, height = 6)

# --- 2. Assign niche labels (same mapping as niche_labels_cm15) ---
# Padded (cluster '0') cells get 'Unassigned' from the lookup; convert to NA
# so all downstream `!is.na()` masks/filters drop them from analysis & plots.
M1$niche_label <- unname(niche_labels_cm15[as.character(M1$kmeans_15)])
M1$niche_label[M1$niche_label == 'Unassigned'] <- NA

#######################################
#############  FIGURE 3G  #############
####### Niche frequency boxplots ######
####### NF / pRV / RVF, log-scale,    #
####### Kruskal-Wallis per niche      #
####### → Figure_3G.pdf               #
#######################################

# --- 3. Niche proportions by disease state ---
pairwise_comparisons <- list(c('NF', 'pRV'), c('NF', 'RVF'), c('pRV', 'RVF'))

niche_meta <- M1@meta.data %>%
  dplyr::filter(!is.na(niche_label)) %>%
  dplyr::select(patient, group, niche_label)

niche_totals <- niche_meta %>%
  dplyr::group_by(patient, group) %>%
  dplyr::summarise(total = n(), .groups = 'drop')

niche_prop <- niche_meta %>%
  dplyr::group_by(patient, group, niche_label) %>%
  dplyr::summarise(n = n(), .groups = 'drop') %>%
  dplyr::left_join(niche_totals, by = c('patient', 'group')) %>%
  dplyr::mutate(proportion = n / total) %>%
  tidyr::complete(tidyr::nesting(patient, group), niche_label, fill = list(n = 0, proportion = 0)) %>%
  dplyr::mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF'))) %>%
  dplyr::filter(niche_label != 'Healthy myocardium')   # drop dominant niche from boxplot

# Compute Kruskal-Wallis per niche for annotation
niche_kw_annot <- niche_prop %>%
  dplyr::group_by(niche_label) %>%
  dplyr::summarise(kw_p   = kruskal.test(proportion ~ group)$p.value, .groups = 'drop') %>%
  dplyr::mutate(kw_padj = p.adjust(kw_p, 'BH'),
                sig      = sig_stars(kw_padj))

# Order niches by median proportion (descending = most abundant at top)
niche_order <- niche_prop %>%
  dplyr::group_by(niche_label) %>%
  dplyr::summarise(med = median(proportion), .groups = 'drop') %>%
  dplyr::arrange(med) %>%
  dplyr::pull(niche_label)

niche_prop$niche_label <- factor(niche_prop$niche_label, levels = niche_order)
niche_kw_annot$niche_label <- factor(niche_kw_annot$niche_label, levels = niche_order)

# Floor zero proportions (from tidyr::complete fill) at half the smallest nonzero
# so they render on a log axis instead of being silently dropped.
.min_nonzero_prop <- min(niche_prop$proportion[niche_prop$proportion > 0], na.rm = TRUE)
.prop_pseudo      <- .min_nonzero_prop / 2
niche_prop$proportion <- pmax(niche_prop$proportion, .prop_pseudo)

x_max <- max(niche_prop$proportion, na.rm = TRUE) * 1.5  # log-scale: multiplicative offset

# Explicit decade breaks covering the data range (default log10 breaks skip decades).
.log_min_3G <- floor(log10(min(niche_prop$proportion, na.rm = TRUE)))
.log_max_3G <- ceiling(log10(max(niche_prop$proportion, na.rm = TRUE)))
.breaks_3G  <- 10^seq(.log_min_3G, .log_max_3G)
# Clean labels: drop trailing zeros so 0.1%, 1%, 10%, 100% instead of 0.10%, 1.00%, etc.
.pct_label_3G <- function(x) {
  paste0(format(x * 100, scientific = FALSE, trim = TRUE, drop0trailing = TRUE), '%')
}

p_3G <- ggplot(niche_prop,
               aes(x = proportion, y = niche_label, colour = group)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA,
               fill = 'white', linewidth = PS$geom_lw, width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.75),
             size = PS$scatter_pt * 0.8, alpha = 0.8) +
  geom_text(data = niche_kw_annot,
            aes(x = x_max, y = niche_label, label = sig),
            inherit.aes = FALSE,
            hjust = 0, size = PS$text_mm * 1.1, family = FONT_FAMILY, colour = 'black') +
  scale_colour_disease(name = 'Disease') +
  scale_x_log10(breaks = .breaks_3G, labels = .pct_label_3G,
                expand = expansion(mult = c(0.02, 0.15))) +
  theme_v52(COMP_W) +
  theme(legend.position = c(0.85, 0.15)) +
  labs(x = 'Proportion of cells (log scale)', y = NULL,
       title = 'Spatial niche distribution by disease state')
save_figure(p_3G, 'Figure_3G.pdf', width = 5.25, height = 8)

# KW stats across all niches
niche_stats <- niche_prop %>%
  dplyr::group_by(niche_label) %>%
  dplyr::summarise(kw_p = kruskal.test(proportion ~ group)$p.value, .groups = 'drop') %>%
  dplyr::mutate(kw_padj = p.adjust(kw_p, method = 'BH'),
                sig = sprintf('p=%.3f', kw_padj))
cat('\nNiche Kruskal-Wallis results (BH-adjusted):\n')
print(as.data.frame(niche_stats %>% dplyr::arrange(kw_padj)), row.names = FALSE)

#######################################
#############  FIGURE 3E  #############
####### Spatial niche maps  ###########
####### (BuildNicheAssay, 14 niches)  #
####### → per-patient tiles +         #
####### composite Figure_3E.pdf       #
#######################################

# --- 4. Spatial niche map (same assembly pattern as Figure 3C) ---
# Each patient tissue rasterised independently and placed via cowplot::draw_plot
# so every tile is a movable element in the assembled Figure_3E.pdf.

suppressPackageStartupMessages(library(cowplot))

.ums_per_inch_D <- 3000
.pt_size_D      <- 0.35
.tile_dir_D     <- file.path(V52_FIG_DIR, 'Figure_3E_tiles')
dir.create(.tile_dir_D, showWarnings = FALSE, recursive = TRUE)

# niche_cells_mask (on M1) is reused by the downstream niche-composition section.
niche_cells_mask <- !is.na(M1$niche_label)

# Transfer niche_label from M1 → xenium.obj via shared cell barcodes for spatial plot.
.shared_cells <- intersect(colnames(xenium.obj), colnames(M1))
.nl_vec <- setNames(rep(NA_character_, ncol(xenium.obj)), colnames(xenium.obj))
.nl_vec[.shared_cells] <- as.character(M1@meta.data[.shared_cells, 'niche_label'])
xenium.obj$niche_label <- .nl_vec

.xn_mask <- !is.na(xenium.obj$niche_label)
meta_niche <- data.frame(
  x = xenium.obj$x_centroid[.xn_mask],
  y = -xenium.obj$y_centroid[.xn_mask],
  group   = xenium.obj$group[.xn_mask],
  patient = xenium.obj$patient[.xn_mask],
  niche   = xenium.obj$niche_label[.xn_mask]
)
meta_niche$group <- factor(meta_niche$group, levels = c('NF', 'pRV', 'RVF'))

# Niche colour palette (consistent across tiles + composite + legend)
meta_niche$niche <- factor(as.character(meta_niche$niche))
.niche_levels <- levels(meta_niche$niche)
.niche_pal    <- setNames(scales::hue_pal()(length(.niche_levels)), .niche_levels)

# Legend grob — build from a synthetic 1-row-per-niche data frame so every
# niche renders a visible key (sampling from meta_niche missed 3 niches).
.leg_src_df <- data.frame(
  x     = seq_along(.niche_levels),
  y     = 0,
  niche = factor(.niche_levels, levels = .niche_levels)
)
.p_leg_D <- ggplot(.leg_src_df, aes(x = x, y = y, color = niche)) +
  geom_point(size = 3) +
  scale_color_manual(values = .niche_pal, name = 'Niche', drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 5), ncol = 1)) +
  theme(legend.title    = element_text(size = 14, face = 'bold'),
        legend.text     = element_text(size = 12),
        legend.key.size = unit(0.55, 'cm'),
        legend.key      = element_rect(fill = 'white', colour = NA))
.leg_grob_D <- get_legend(.p_leg_D)
ggsave(file.path(.tile_dir_D, 'Figure_3E_legend.pdf'),
       plot = ggdraw(.leg_grob_D),
       width = 4.5, height = 4.5, device = cairo_pdf)
message('Saved: Figure_3E_legend.pdf (', length(.niche_levels), ' niches)')

# Scale bar standalone PDF
.p_sb_D_solo <- ggplot() +
  annotate('segment', x = 0, xend = 1, y = 0, yend = 0,
           linewidth = 0.8, color = 'black') +
  annotate('text', x = 0.5, y = 0.2, label = '2 mm', size = 3, hjust = 0.5) +
  xlim(-0.05, 1.05) + ylim(-0.3, 0.5) + theme_void()
ggsave(file.path(.tile_dir_D, 'Figure_3E_scalebar.pdf'),
       plot = .p_sb_D_solo,
       width = 2000 / .ums_per_inch_D, height = 0.3, device = cairo_pdf)

# Per-patient origin-normalized data
.per_pt_D <- split(meta_niche, meta_niche$patient)
.per_pt_D <- lapply(.per_pt_D, function(d) {
  d$x <- d$x - min(d$x); d$y <- d$y - min(d$y); d
})

.pt_info_D <- list(); .pt_plots_D <- list()
for (.pt_D in names(.per_pt_D)) {
  d <- .per_pt_D[[.pt_D]]
  .pt_info_D[[.pt_D]]  <- list(group = as.character(d$group[1]),
                                w_um  = diff(range(d$x)),
                                h_um  = diff(range(d$y)))
  .pt_plots_D[[.pt_D]] <- ggplot(d, aes(x = x, y = y, color = niche)) +
    ggrastr::rasterise(geom_point(size = .pt_size_D, stroke = 0), dpi = 400) +
    scale_color_manual(values = .niche_pal, guide = 'none') +
    coord_fixed(expand = FALSE) +
    theme_void() +
    theme(legend.position = 'none', plot.margin = margin(0, 0, 0, 0))
}

# Individual tile PDFs
for (.pt_D in names(.pt_plots_D)) {
  .w_in_D  <- .pt_info_D[[.pt_D]]$w_um / .ums_per_inch_D
  .h_in_D  <- .pt_info_D[[.pt_D]]$h_um / .ums_per_inch_D
  .grp_D   <- .pt_info_D[[.pt_D]]$group
  .fname_D <- sprintf('Figure_3E_%s_%s.pdf', .grp_D, .pt_D)
  ggsave(file.path(.tile_dir_D, .fname_D),
         plot = .pt_plots_D[[.pt_D]],
         width = .w_in_D, height = .h_in_D, device = cairo_pdf)
  message('Saved tile: ', .fname_D, ' (', round(.w_in_D, 2), ' × ', round(.h_in_D, 2), ' in)')
}

# Composite assembly (same pattern as 3C)
.gap_um_D     <- 500
.row_gap_um_D <- 1500
.grp_gap_um_D <- 5000
.grp_levels_D <- c('NF', 'pRV', 'RVF')

.layout_group_D <- function(.g) {
  pts <- names(.pt_info_D)[sapply(.pt_info_D, function(x) x$group == .g)]
  bbs <- .pt_info_D[pts]
  pts <- pts[order(-sapply(bbs, function(b) b$w_um * b$h_um))]
  bbs <- bbs[pts]
  b1 <- bbs[[1]]; b2 <- bbs[[2]]; b3 <- bbs[[3]]
  bot_w <- b1$w_um + .gap_um_D + b2$w_um
  bot_h <- max(b1$h_um, b2$h_um)
  tiles <- list(
    list(x = 0,                     y = 0,                     w = b1$w_um, h = b1$h_um),
    list(x = b1$w_um + .gap_um_D,   y = 0,                     w = b2$w_um, h = b2$h_um),
    list(x = (bot_w - b3$w_um) / 2, y = bot_h + .row_gap_um_D, w = b3$w_um, h = b3$h_um)
  )
  names(tiles) <- pts
  list(tiles  = tiles,
       canvas = list(w = bot_w, h = bot_h + .row_gap_um_D + b3$h_um),
       top_pt = pts[3], bot_h = bot_h)
}
.layouts_D <- setNames(lapply(.grp_levels_D, .layout_group_D), .grp_levels_D)

.canvas_w_um_D <- max(sapply(.layouts_D, function(l) l$canvas$w))
.canvas_h_um_D <- sum(sapply(.layouts_D, function(l) l$canvas$h)) +
                  (length(.grp_levels_D) - 1) * .grp_gap_um_D

.grp_y_D <- numeric(length(.grp_levels_D)); names(.grp_y_D) <- .grp_levels_D
.cur_y_D <- 0
for (.g_D in rev(.grp_levels_D)) {
  .grp_y_D[.g_D] <- .cur_y_D
  .cur_y_D       <- .cur_y_D + .layouts_D[[.g_D]]$canvas$h + .grp_gap_um_D
}

.canvas_D <- ggdraw()
for (.g_D in .grp_levels_D) {
  lay    <- .layouts_D[[.g_D]]
  grp_yo <- .grp_y_D[.g_D]
  grp_xo <- (.canvas_w_um_D - lay$canvas$w) / 2
  for (.pt_D in names(lay$tiles)) {
    t <- lay$tiles[[.pt_D]]
    .canvas_D <- .canvas_D + draw_plot(
      .pt_plots_D[[.pt_D]],
      x      = (grp_xo + t$x) / .canvas_w_um_D,
      y      = (grp_yo + t$y) / .canvas_h_um_D,
      width  = t$w / .canvas_w_um_D,
      height = t$h / .canvas_h_um_D
    )
  }
  .canvas_D <- .canvas_D + draw_label(
    .g_D,
    x = 0.01,
    y = (grp_yo + lay$canvas$h) / .canvas_h_um_D,
    hjust = 0, vjust = 1,
    fontfamily = FONT_FAMILY, fontface = 'bold', size = 22
  )
}

# Legend + scale bar → RVF top-right white space
.rvf_lay_D  <- .layouts_D[['RVF']]
.rvf_yo_D   <- .grp_y_D[['RVF']]
.rvf_xo_D   <- (.canvas_w_um_D - .rvf_lay_D$canvas$w) / 2
.top_t_D    <- .rvf_lay_D$tiles[[.rvf_lay_D$top_pt]]
.leg_x_um_D <- .rvf_xo_D + .top_t_D$x + .top_t_D$w + 1000
.leg_y_um_D <- .rvf_yo_D + .rvf_lay_D$bot_h + .row_gap_um_D
.leg_w_um_D <- (.rvf_xo_D + .rvf_lay_D$canvas$w) - .leg_x_um_D
.leg_h_um_D <- .top_t_D$h * 0.95
.canvas_D <- .canvas_D + draw_plot(
  .leg_grob_D,
  x      = .leg_x_um_D / .canvas_w_um_D,
  y      = (.leg_y_um_D + .top_t_D$h - .leg_h_um_D) / .canvas_h_um_D,
  width  = .leg_w_um_D / .canvas_w_um_D,
  height = .leg_h_um_D / .canvas_h_um_D
)

.sb_y_um_D  <- .leg_y_um_D - 800
.sb_x0_um_D <- .leg_x_um_D + (.leg_w_um_D - 2000) / 2
.sb_x1_um_D <- .sb_x0_um_D + 2000
.canvas_D <- .canvas_D +
  draw_line(x = c(.sb_x0_um_D, .sb_x1_um_D) / .canvas_w_um_D,
            y = c(.sb_y_um_D,  .sb_y_um_D ) / .canvas_h_um_D,
            size = 0.8, color = 'black') +
  draw_label('2 mm',
             x = ((.sb_x0_um_D + .sb_x1_um_D) / 2) / .canvas_w_um_D,
             y = (.sb_y_um_D - 800) / .canvas_h_um_D,
             hjust = 0.5, vjust = 1, size = 14,
             fontfamily = FONT_FAMILY)

.fig_w_in_D <- .canvas_w_um_D / .ums_per_inch_D
.fig_h_in_D <- .canvas_h_um_D / .ums_per_inch_D
save_figure(.canvas_D, 'Figure_3E.pdf',
            width = .fig_w_in_D, height = .fig_h_in_D)
message('Assembled Figure_3E.pdf: ', round(.fig_w_in_D, 2),
        ' × ', round(.fig_h_in_D, 2), ' in (9 niche tissues independently movable)')

# --- 5. Per-niche highlight maps (one PDF per niche: target red, others gray) ---
# Reuses .per_pt_D, .pt_info_D, .layouts_D, .grp_y_D, and canvas dims from
# section 4. Target-niche points are plotted on top of gray background points.
.highlight_dir_D <- file.path(V52_FIG_DIR, 'Figure_3E_highlights')
dir.create(.highlight_dir_D, showWarnings = FALSE, recursive = TRUE)

for (.target_niche in .niche_levels) {
  .pt_plots_H <- list()
  for (.pt_D in names(.per_pt_D)) {
    d <- .per_pt_D[[.pt_D]]
    d$highlight <- ifelse(as.character(d$niche) == .target_niche, 'target', 'other')
    # Draw 'other' first so 'target' points render on top
    d <- d[order(d$highlight == 'target'), ]
    d$highlight <- factor(d$highlight, levels = c('other', 'target'))
    .pt_plots_H[[.pt_D]] <- ggplot(d, aes(x = x, y = y, color = highlight)) +
      ggrastr::rasterise(geom_point(size = .pt_size_D, stroke = 0), dpi = 400) +
      scale_color_manual(values = c('other' = 'gray85', 'target' = 'red'),
                         guide = 'none', drop = FALSE) +
      coord_fixed(expand = FALSE) +
      theme_void() +
      theme(legend.position = 'none', plot.margin = margin(0, 0, 0, 0))
  }

  .canvas_H <- ggdraw()
  for (.g_D in .grp_levels_D) {
    lay    <- .layouts_D[[.g_D]]
    grp_yo <- .grp_y_D[.g_D]
    grp_xo <- (.canvas_w_um_D - lay$canvas$w) / 2
    for (.pt_D in names(lay$tiles)) {
      t <- lay$tiles[[.pt_D]]
      .canvas_H <- .canvas_H + draw_plot(
        .pt_plots_H[[.pt_D]],
        x      = (grp_xo + t$x) / .canvas_w_um_D,
        y      = (grp_yo + t$y) / .canvas_h_um_D,
        width  = t$w / .canvas_w_um_D,
        height = t$h / .canvas_h_um_D
      )
    }
    .canvas_H <- .canvas_H + draw_label(
      .g_D,
      x = 0.01,
      y = (grp_yo + lay$canvas$h) / .canvas_h_um_D,
      hjust = 0, vjust = 1,
      fontfamily = FONT_FAMILY, fontface = 'bold', size = 22
    )
  }

  # Niche name as red title at top of canvas
  .canvas_H <- .canvas_H + draw_label(
    .target_niche,
    x = 0.5, y = 0.995,
    hjust = 0.5, vjust = 1,
    fontfamily = FONT_FAMILY, fontface = 'bold', size = 24, colour = 'red'
  )

  .fname_H <- sprintf('Figure_3E_highlight_%s.pdf',
                      gsub('[^A-Za-z0-9]+', '_', .target_niche))
  ggsave(file.path(.highlight_dir_D, .fname_H),
         plot = .canvas_H,
         width = .fig_w_in_D, height = .fig_h_in_D, device = cairo_pdf)
  message('Saved highlight: ', .fname_H)
}
message('Assembled ', length(.niche_levels),
        ' per-niche highlight maps in Figure_3E_highlights/')

rm(.ums_per_inch_D, .pt_size_D, .tile_dir_D, .p_leg_D, .p_sb_D_solo,
   .niche_levels, .niche_pal,
   .per_pt_D, .pt_info_D, .pt_plots_D,
   .pt_D, .grp_D, .w_in_D, .h_in_D, .fname_D,
   .gap_um_D, .row_gap_um_D, .grp_gap_um_D, .grp_levels_D, .layouts_D,
   .canvas_w_um_D, .canvas_h_um_D, .grp_y_D, .cur_y_D, .canvas_D,
   .rvf_lay_D, .rvf_yo_D, .rvf_xo_D, .top_t_D,
   .leg_x_um_D, .leg_y_um_D, .leg_w_um_D, .leg_h_um_D,
   .sb_y_um_D, .sb_x0_um_D, .sb_x1_um_D,
   .fig_w_in_D, .fig_h_in_D, .layout_group_D, .g_D)
gc(verbose = FALSE)

#######################################
#############  FIGURE 3F  #############
####### Stacked bar of cell-type ######
####### composition per niche  ########
####### → Figure_3F.pdf               #
#######################################

# --- 5. Panel F: stacked bar of cell-type composition per niche -------------
.niche_ct_csv <- './output/Figure_3/Xenium/niche_celltype.csv'
.nc <- read.csv(.niche_ct_csv, check.names = FALSE, row.names = 1)
# Drop 'Unassigned' (cluster 0 padding) row before plotting so composition bar shows only real niches
.nc <- .nc[rownames(.nc) != 'Unassigned', , drop = FALSE]
.nc_long <- reshape2::melt(as.matrix(.nc), varnames = c('niche', 'cell_type'),
                           value.name = 'n')
.nc_long <- .nc_long %>%
  dplyr::group_by(niche) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

# Order niches by frequency — same order as Figure 3F (highest count leftmost)
.nc_long$niche <- factor(.nc_long$niche, levels = .niche_freq_order)

# Order cell types by overall abundance (largest segment at bottom)
.ct_order <- names(sort(colSums(.nc), decreasing = TRUE))
.nc_long$cell_type <- factor(.nc_long$cell_type, levels = rev(.ct_order))

# 12-color palette for coarse cell types
.ct_pal_E <- setNames(
  c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFD92F',
    '#A65628','#F781BF','#66C2A5','#8DA0CB','#17BECF','#B3B3B3'),
  c('CM','EC','FB','SM','Myeloid','Adipo',
    'NKT','Epi','PC','LEC','Endo','Neuron')
)

p_niche_comp_bar <- ggplot(.nc_long, aes(x = niche, y = prop, fill = cell_type)) +
  geom_bar(stat = 'identity', width = 0.8) +
  scale_fill_manual(values = .ct_pal_E, breaks = .ct_order, name = 'Cell type') +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(accuracy = 1)) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(face = 'bold'),
        plot.margin = margin(8, 8, 8, 8)) +
  labs(x = NULL, y = 'Fraction of niche')

save_figure(p_niche_comp_bar, 'Figure_3F.pdf',
            width = 10, height = 5)

##############################################
##############################################
#### Cell type frequency by disease state (XeniumReanalysis lines 8009-8110)
##############################################
##############################################

if ('cell_types_manual' %in% colnames(xenium.obj@meta.data)) {
  ct_freq <- data.frame(
    patient = xenium.obj$patient,
    group = xenium.obj$group,
    cell_type = xenium.obj$cell_types_manual
  )
  ct_freq <- ct_freq[ct_freq$cell_type != 'Unassigned', ]
  ct_freq$group <- factor(ct_freq$group, levels = c('NF', 'pRV', 'RVF'))

  ct_prop <- ct_freq %>%
    dplyr::group_by(patient, group, cell_type) %>%
    dplyr::summarise(n = n(), .groups = 'drop') %>%
    dplyr::group_by(patient) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()

  ct_order <- ct_prop %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(mean_prop = mean(prop), .groups = 'drop') %>%
    dplyr::arrange(dplyr::desc(mean_prop)) %>%
    dplyr::pull(cell_type)
  ct_prop$cell_type <- factor(ct_prop$cell_type, levels = ct_order)

  p_ct_freq <- ggplot(ct_prop, aes(x = cell_type, y = prop * 100, color = group)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.4) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
                size = 1.2, alpha = 0.7) +
    scale_color_manual(values = c(NF = 'steelblue', pRV = 'darkorange', RVF = 'firebrick')) +
    coord_flip() + theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(size = 7)) +
    labs(x = NULL, y = 'Proportion of cells (%)', color = 'Group',
         title = 'Cell type frequency by disease state')

  save_figure(p_ct_freq, 'Figure_3_celltype_frequency_by_group.pdf', width = 10, height = 14)

  p_ct_stacked <- ggplot(ct_prop %>% dplyr::mutate(cell_type = factor(cell_type, levels = rev(ct_order))),
      aes(x = patient, y = prop * 100, fill = cell_type)) +
    geom_bar(stat = 'identity', width = 0.8) +
    facet_wrap(~ group, scales = 'free_x', nrow = 1) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.text = element_text(size = 6), legend.key.size = unit(0.3, 'cm'),
          strip.text = element_text(face = 'bold')) +
    labs(x = NULL, y = 'Proportion of cells (%)', fill = 'Cell type',
         title = 'Cell type composition per patient')

  save_figure(p_ct_stacked, 'Figure_3_celltype_composition_per_patient.pdf', width = 12, height = 7)

  # Top disease-variable cell types
  pairwise_comparisons <- list(c('NF','pRV'), c('pRV','RVF'), c('NF','RVF'))
  ct_signal <- ct_prop %>%
    dplyr::group_by(cell_type, group) %>%
    dplyr::summarise(group_mean = mean(prop), group_sd = sd(prop), .groups = 'drop') %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(between_cv = sd(group_mean) / mean(group_mean),
                     within_cv_mean = mean(group_sd / (group_mean + 1e-6)), .groups = 'drop') %>%
    dplyr::mutate(signal_ratio = between_cv / (within_cv_mean + 0.01)) %>%
    dplyr::arrange(dplyr::desc(signal_ratio))

  top_signal_cts <- head(ct_signal$cell_type, 15)
  ct_prop_top <- ct_prop %>% dplyr::filter(cell_type %in% top_signal_cts)
  ct_prop_top$cell_type <- factor(ct_prop_top$cell_type, levels = top_signal_cts)

  p_ct_top <- ggplot(ct_prop_top, aes(x = group, y = prop * 100, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, size = 1.5) +
    stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test',
                       label = 'p.signif', size = 3) +
    stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 2.5) +
    scale_fill_manual(values = c(NF = 'steelblue', pRV = 'darkorange', RVF = 'firebrick')) +
    facet_wrap(~ cell_type, scales = 'free_y', ncol = 5) +
    theme_bw(base_size = 10) +
    theme(legend.position = 'none', strip.text = element_text(face = 'bold', size = 8),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = 'Proportion of cells (%)',
         title = 'Top disease-variable cell types (ranked by signal-to-noise ratio)')

  save_figure(p_ct_top, 'Figure_3_celltype_top_disease_variable.pdf', width = 14, height = 10)
}

##############################################
##############################################
#### Supplementary: Niche cluster frequency statistics
##############################################
##############################################

niche.counts <- table(M1$niche_manual)
niche.counts <- niche.counts[names(niche.counts) != 'Unassigned']   # drop padded cells

niche.patient <- table(M1$niche_manual,M1$patient)
niche.patient <- niche.patient[rownames(niche.patient) != 'Unassigned', ]
niche.patient <- t(t(niche.patient) / colSums(niche.patient))
niche.patient <- data.frame(niche.patient)  # Var1=niche, Var2=patient, Freq
## Dynamic patient -> disease group (replaces fragile hardcoded vector).
.pg2 <- unique(data.frame(patient = as.character(M1$patient),
                          group   = as.character(M1$group)))
.pg2 <- setNames(.pg2$group, .pg2$patient)
niche.patient$niche <- factor(.pg2[as.character(niche.patient$Var2)],
                               levels = c('NF','pRV','RVF'))
niche.patient$Var1 <- factor(niche.patient$Var1, levels=rev(names(sort(niche.counts))))

library(ggpubr)
p_niche_freq_stats <- ggboxplot(niche.patient,x="niche",y="Freq",fill="niche",group="niche", size = PS$geom_lw)+
  theme_v52(COMP_W) +
  scale_fill_disease() +
  facet_wrap(~Var1,ncol=7) +
  scale_y_touch() +
  stat_compare_means(aes(group=niche),method="anova", size = PS$text_mm, symnum.args = pub_signif_args)
save_figure(p_niche_freq_stats, 'Figure_3_niche_clust_freq_stats.pdf', width=12.5, height=15)




message('Figure 3 (v53) per-panel PDFs (A–G + supplementary) written to ', V52_FIG_DIR)

