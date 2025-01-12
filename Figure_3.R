
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
Crop_custom <- 
  function(x = NULL, 
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
      if (crop_molecules && 
          grep("molecules", names(fov)) != 0) {
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



### SLIDE 1

# Load the Xenium data
data <- ReadXenium('~/Downloads/hdWGCNA_TOM/Xenium/proseg-output-XETG00217__0038213__Region_1__20241206__182124',
	outs = c("matrix", "microns"), type = c("centroids", "segmentations"))


segmentations.data <- list(
  centroids = CreateCentroids(data$centroids),
  segmentation = CreateSegmentation(data$segmentations))
coords <- CreateFOV(
  coords = segmentations.data, 
  type = c("segmentation", "centroids"), 
  molecules = data$microns, 
  assay = "Xenium")
xenium.obj <- CreateSeuratObject(
  counts = data$matrix[["Gene Expression"]], 
  assay = "Xenium")
xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
xenium.obj[["fov"]] <- coords

# remove cells with 0 counts
xenium.obj1 <- subset(xenium.obj, subset = nCount_Xenium > 0)

region1 <- read.table('~/Downloads/hdWGCNA_TOM/Xenium/output-XETG00217__0038213__Region_1__20241206__182124/Region1_coordinates.csv',sep=',')
region2 <- read.table('~/Downloads/hdWGCNA_TOM/Xenium/output-XETG00217__0038213__Region_1__20241206__182124/Region2_coordinates.csv',sep=',')
region3 <- read.table('~/Downloads/hdWGCNA_TOM/Xenium/output-XETG00217__0038213__Region_1__20241206__182124/Region3_coordinates.csv',sep=',')

names1 <- rep("seg1", length(region1$V1))
names2 <- rep("seg2", length(region2$V1))
names3 <- rep("seg3", length(region3$V1))

reg1 <- data.frame(x = region1$V1, y = region1$V2, cell = names1)
seg1 <- CreateSegmentation(reg1)

reg2 <- data.frame(x = region2$V1, y = region2$V2, cell = names2)
seg2 <- CreateSegmentation(reg2)

reg3 <- data.frame(x = region3$V1, y = region3$V2, cell = names3)
seg3 <- CreateSegmentation(reg3)

xenium.obj1[["seg1"]] <- Overlay(xenium.obj1[["fov"]], seg1)
xenium.obj1[["seg2"]] <- Overlay(xenium.obj1[["fov"]], seg2)
xenium.obj1[["seg3"]] <- Overlay(xenium.obj1[["fov"]], seg3)

cells_1697 <- Cells(xenium.obj1[['seg1']])
cells_1691 <- Cells(xenium.obj1[['seg2']])
cells_1618 <- Cells(xenium.obj1[['seg3']])

xenium.1697 <- subset(xenium.obj1,cells = cells_1697)
xenium.1691 <- subset(xenium.obj1,cells = cells_1691)
xenium.1618 <- subset(xenium.obj1,cells = cells_1618)

xenium.1697$patient <- '1697'
xenium.1691$patient <- '1691'
xenium.1618$patient <- '1618'

xenium.1697$group <- 'NF'
xenium.1691$group <- 'NF'
xenium.1618$group <- 'pRV'

saveRDS(xenium.1697,'~/Downloads/hdWGCNA_TOM/Xenium/xenium_1697.rds')
saveRDS(xenium.1691,'~/Downloads/hdWGCNA_TOM/Xenium/xenium_1691.rds')
saveRDS(xenium.1618,'~/Downloads/hdWGCNA_TOM/Xenium/xenium_1618.rds')



### Sample 2

# Load the Xenium data
data <- ReadXenium('~/Downloads/hdWGCNA_TOM/Xenium/proseg-output-XETG00217__0038216__Region_1__20241206__182124',
	outs = c("matrix", "microns"), type = c("centroids", "segmentations"))


segmentations.data <- list(
  centroids = CreateCentroids(data$centroids),
  segmentation = CreateSegmentation(data$segmentations))
coords <- CreateFOV(
  coords = segmentations.data, 
  type = c("segmentation", "centroids"), 
  molecules = data$microns, 
  assay = "Xenium")
xenium.obj <- CreateSeuratObject(
  counts = data$matrix[["Gene Expression"]], 
  assay = "Xenium")
xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
xenium.obj[["fov"]] <- coords

# remove cells with 0 counts
xenium.obj2 <- subset(xenium.obj, subset = nCount_Xenium > 0)

region1 <- read.table('~/Downloads/hdWGCNA_TOM/Xenium/output-XETG00217__0038216__Region_1__20241206__182124/Region1_coordinates.csv',sep=',')
region2 <- read.table('~/Downloads/hdWGCNA_TOM/Xenium/output-XETG00217__0038216__Region_1__20241206__182124/Region2_coordinates.csv',sep=',')

names1 <- rep("seg1", length(region1$V1))
names2 <- rep("seg2", length(region2$V1))

reg1 <- data.frame(x = region1$V1, y = region1$V2, cell = names1)
seg1 <- CreateSegmentation(reg1)

reg2 <- data.frame(x = region2$V1, y = region2$V2, cell = names2)
seg2 <- CreateSegmentation(reg2)

xenium.obj2[["seg1"]] <- Overlay(xenium.obj2[["fov"]], seg1)
xenium.obj2[["seg2"]] <- Overlay(xenium.obj2[["fov"]], seg2)

cells_1567 <- Cells(xenium.obj2[['seg1']])
cells_1692 <- Cells(xenium.obj2[['seg2']])

xenium.1567 <- subset(xenium.obj2,cells = cells_1567)
xenium.1692 <- subset(xenium.obj2,cells = cells_1692)

xenium.1567$patient <- '1567'
xenium.1692$patient <- '1692'

xenium.1567$group <- 'pRV'
xenium.1692$group <- 'pRV'

saveRDS(xenium.1567,'~/Downloads/hdWGCNA_TOM/Xenium/xenium_1567.rds')
saveRDS(xenium.1692,'~/Downloads/hdWGCNA_TOM/Xenium/xenium_1692.rds')


### Sample 3

# Load the Xenium data
data <- ReadXenium('~/Downloads/hdWGCNA_TOM/Xenium/proseg-output-XETG00217__0038290__Region_1__20241212__142808',
	outs = c("matrix", "microns"), type = c("centroids", "segmentations"))


segmentations.data <- list(
  centroids = CreateCentroids(data$centroids),
  segmentation = CreateSegmentation(data$segmentations))
coords <- CreateFOV(
  coords = segmentations.data, 
  type = c("segmentation", "centroids"), 
  molecules = data$microns, 
  assay = "Xenium")
xenium.obj <- CreateSeuratObject(
  counts = data$matrix[["Gene Expression"]], 
  assay = "Xenium")
xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
xenium.obj[["fov"]] <- coords

# remove cells with 0 counts
xenium.obj3 <- subset(xenium.obj, subset = nCount_Xenium > 0)

cells_1561 <- Cells(xenium.obj3[['fov']])

xenium.1561 <- subset(xenium.obj3,cells = cells_1561)

xenium.1561$patient <- '1561'

xenium.1561$group <- 'NF'

saveRDS(xenium.1561,'~/Downloads/hdWGCNA_TOM/Xenium/xenium_1561.rds')



### Sample 4

# Load the Xenium data
data <- ReadXenium('~/Downloads/hdWGCNA_TOM/Xenium/proseg-output-XETG00217__0038290__Region_2__20241212__142808',
	outs = c("matrix", "microns"), type = c("centroids", "segmentations"))


segmentations.data <- list(
  centroids = CreateCentroids(data$centroids),
  segmentation = CreateSegmentation(data$segmentations))
coords <- CreateFOV(
  coords = segmentations.data, 
  type = c("segmentation", "centroids"), 
  molecules = data$microns, 
  assay = "Xenium")
xenium.obj <- CreateSeuratObject(
  counts = data$matrix[["Gene Expression"]], 
  assay = "Xenium")
xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
xenium.obj[["fov"]] <- coords

# remove cells with 0 counts
xenium.obj4 <- subset(xenium.obj, subset = nCount_Xenium > 0)

cells_1343 <- Cells(xenium.obj4[['fov']])

xenium.1343 <- subset(xenium.obj4,cells = cells_1343)

xenium.1343$patient <- '1343'

xenium.1343$group <- 'RVF'

saveRDS(xenium.1343,'~/Downloads/hdWGCNA_TOM/Xenium/xenium_1343.rds')


### Sample 5

# Load the Xenium data
data <- ReadXenium('~/Downloads/hdWGCNA_TOM/Xenium/proseg-output-XETG00217__0038291__Region_1__20241212__142808',
	outs = c("matrix", "microns"), type = c("centroids", "segmentations"))


segmentations.data <- list(
  centroids = CreateCentroids(data$centroids),
  segmentation = CreateSegmentation(data$segmentations))
coords <- CreateFOV(
  coords = segmentations.data, 
  type = c("segmentation", "centroids"), 
  molecules = data$microns, 
  assay = "Xenium")
xenium.obj <- CreateSeuratObject(
  counts = data$matrix[["Gene Expression"]], 
  assay = "Xenium")
xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
xenium.obj[["fov"]] <- coords

# remove cells with 0 counts
xenium.obj5 <- subset(xenium.obj, subset = nCount_Xenium > 0)

cells_1467 <- Cells(xenium.obj5[['fov']])

xenium.1467 <- subset(xenium.obj5,cells = cells_1467)

xenium.1467$patient <- '1467'

xenium.1467$group <- 'RVF'

saveRDS(xenium.1467,'~/Downloads/hdWGCNA_TOM/Xenium/xenium_1467.rds')



### Sample 6

# Load the Xenium data
data <- ReadXenium('~/Downloads/hdWGCNA_TOM/Xenium/proseg-output-XETG00217__0038291__Region_2__20241212__142808',
	outs = c("matrix", "microns"), type = c("centroids", "segmentations"))


segmentations.data <- list(
  centroids = CreateCentroids(data$centroids),
  segmentation = CreateSegmentation(data$segmentations))
coords <- CreateFOV(
  coords = segmentations.data, 
  type = c("segmentation", "centroids"), 
  molecules = data$microns, 
  assay = "Xenium")
xenium.obj <- CreateSeuratObject(
  counts = data$matrix[["Gene Expression"]], 
  assay = "Xenium")
xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
xenium.obj[["fov"]] <- coords

# remove cells with 0 counts
xenium.obj6 <- subset(xenium.obj, subset = nCount_Xenium > 0)

cells_1632 <- Cells(xenium.obj6[['fov']])


xenium.1632 <- subset(xenium.obj6,cells = cells_1632)

xenium.1632$patient <- '1632'

xenium.1632$group <- 'RVF'

saveRDS(xenium.1632,'~/Downloads/hdWGCNA_TOM/Xenium/xenium_1632.rds')


##############################################
##############################################
#### MERGING AND PREPROCESSING
##############################################
##############################################

xenium.1343 <- readRDS('~/Downloads/hdWGCNA_TOM/Xenium/xenium_1343.rds')
xenium.1467 <- readRDS('~/Downloads/hdWGCNA_TOM/Xenium/xenium_1467.rds')
xenium.1561 <- readRDS('~/Downloads/hdWGCNA_TOM/Xenium/xenium_1561.rds')
xenium.1567 <- readRDS('~/Downloads/hdWGCNA_TOM/Xenium/xenium_1567.rds')
xenium.1618 <- readRDS('~/Downloads/hdWGCNA_TOM/Xenium/xenium_1618.rds')
xenium.1632 <- readRDS('~/Downloads/hdWGCNA_TOM/Xenium/xenium_1632.rds')
xenium.1691 <- readRDS('~/Downloads/hdWGCNA_TOM/Xenium/xenium_1691.rds')
xenium.1692 <- readRDS('~/Downloads/hdWGCNA_TOM/Xenium/xenium_1692.rds')
xenium.1697 <- readRDS('~/Downloads/hdWGCNA_TOM/Xenium/xenium_1697.rds')




xenium.obj <- merge(xenium.1343,c(xenium.1467,xenium.1561,xenium.1567,xenium.1618,xenium.1632,xenium.1691,xenium.1692,xenium.1697))
rm(xenium.1343)
rm(xenium.1467)
rm(xenium.1561)
rm(xenium.1567)
rm(xenium.1618)
rm(xenium.1632)
rm(xenium.1691)
rm(xenium.1692)
rm(xenium.1697)
gc()


xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunHarmony(xenium.obj,'patient')
xenium.obj <- RunUMAP(xenium.obj, reduction = "harmony", dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "harmony", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.3)

#saveRDS(xenium.obj,'~/Downloads/hdWGCNA_TOM/Xenium/collated_xenium_data.rds')

DimPlot(xenium.obj,label=T)

###### A lot of doublets which need to be removed
###### Cluster at high resolution to identify doublet clusters
xenium.obj <- FindClusters(xenium.obj, resolution = 1)

# Markers
# Adipo - ADIPOQ
# CM - TNNT2
# EC - VWF/PECAM1
# Endo - LEPR/BMP6
# Epi - WT1
# FB - DCN
# LEC - PROX1
# Myeloid - LYVE1/F13A1
# Neuron - SCN7A
# NKT - CD2
# PC - PDGFRB
# SM - MYH11

DotPlot(xenium.obj,c('ADIPOQ','TNNT2','VWF','LEPR','WT1','DCN','PROX1','F13A1','SCN7A','CD2','PDGFRB','MYH11'))

# 0 EC
# 1 FB
# 2 DOUBLET - EC/CM
# 3 CM
# 4 EC
# 5 FB
# 6 CM
# 7 CM
# 8 CM
# 9 DOUBLET - FB/CM
# 10 DOUBLET - FB/CM
# 11 Myelod
# 12 CM
# 13 SM/PC
# 14 CM
# 15 CM
# 16 DOUBLET - FB/Myeloid/CM
# 17 Myeloid
# 18 Adipo
# 19 CM
# 20 CM
# 21 FB
# 22 FB
# 23 DOUBLET - CM/FB
# 24 FB
# 25 T cell
# 26 FB
# 27 FB
# 28 LEC
# 29 FB
# 30 DOUBLET - Myeloid/CM
# 31 CM
# 32 DOUBLET - CM/FB
# 33 Myeloid
# 34 CM
# 35 CM
# 36 CM
# 37 CM
# 38 CM
# 39 DOUBLET - CM/SM/PC
# 40 CM
# 41 CM
# 42 CM
# 43 CM

# DOUBLET clusters = 2, 9, 10, 16, 23, 30, 32, 33, 39
# Note we keep cluster 2 for now to clarify if some of it may be biological
xenium.obj <- subset(xenium.obj,idents = c(9,10,16,23,30,32,33,39), invert = TRUE)

# Clust 2 seems to be composed off some true ECs expressing CM transcrits and then missegmented CMs
clust2 <- subset(xenium.obj,idents=2)
FeaturePlot(clust2,'VWF',label=T)

# Identify true ECs
clust2 <- FindNeighbors(clust2, reduction = "harmony", dims = 1:30)
clust2 <- FindClusters(clust2, resolution = 1)
FeaturePlot(clust2,'VWF',label=T)
# We want to keep 0, 1, 2, 4, 12, and 13
# Get rid of 3, 5, 6, 7, 8, 9, 10, 11
clust2_del<- subset(clust2,idents=c(3,5,6,7,8,9,10,11))
to_remove <- Cells(clust2_del)

xenium.obj <- subset(xenium.obj,cells = to_remove, invert = TRUE)

# Re-transform and re-cluster
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunHarmony(xenium.obj,'patient')
xenium.obj <- RunUMAP(xenium.obj, reduction = "harmony", dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "harmony", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 1)

#saveRDS(xenium.obj,'~/Downloads/hdWGCNA_TOM/Xenium/collated_xenium_data_clean.rds')



#### Final ROUND OF DOUBLET REMOVAL
#### Include some manual removal here
DimPlot(xenium.obj,label=T)


#36 is a doublet
xenium.obj <- subset(xenium.obj,idents = c(36), invert = TRUE)

cells.located <- CellSelector(plot = DimPlot(xenium.obj,label=T))
xenium.obj <- subset(xenium.obj,cells = cells.located, invert = TRUE)

xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunHarmony(xenium.obj,'patient')
xenium.obj <- RunUMAP(xenium.obj, reduction = "harmony", dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "harmony", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.5)


##############################################
##############################################
#### COARSE CELL TYPE ASSIGNMENT
##############################################
##############################################


new.cluster.ids <- c("FB","EC","CM","CM","EC",
	"CM","Myeloid","PC","CM","SM",
	"Adipo","CM","FB","CM","FB",
	"NKT","FB","LEC","CM","CM",
	"CM","CM","CM","CM","CM")
names(new.cluster.ids) <- levels(xenium.obj)
xenium.obj <- RenameIdents(xenium.obj, new.cluster.ids)
xenium.obj$names <- xenium.obj@active.ident

saveRDS(xenium.obj,'~/Downloads/hdWGCNA_TOM/Xenium/collated_xenium_data_FINAL.rds')

VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

xenium.obj$subnames <- as.character(xenium.obj$names)
xenium.obj$markernames <- as.character(xenium.obj$names)


ImageDimPlot(xenium.obj, fov = "fov.2",split.by='ident',cols = "polychrome", axes = TRUE)

##############################################
##############################################
#### FINE-GRAINED CELL TYPE ASSIGNMENT
##############################################
##############################################


### Cardiomyocytes

cm <- subset(xenium.obj,idents=c('CM'))

cm <- SCTransform(cm, assay = "Xenium")
cm <- RunPCA(cm, npcs = 30, features = rownames(cm))
cm <- RunHarmony(cm,'patient')
cm <- RunUMAP(cm, reduction = "harmony", dims = 1:30)
cm <- FindNeighbors(cm, reduction = "harmony", dims = 1:30)
cm <- FindClusters(cm, resolution = 0.5)

#Rename slots to recorrect UMI
slot(object = cm@assays$SCT@SCTModel.list[[1]], name="umi.assay") <- "Xenium"
slot(object = cm@assays$SCT@SCTModel.list[[2]], name="umi.assay") <- "Xenium"
slot(object = cm@assays$SCT@SCTModel.list[[3]], name="umi.assay") <- "Xenium"
slot(object = cm@assays$SCT@SCTModel.list[[4]], name="umi.assay") <- "Xenium"
slot(object = cm@assays$SCT@SCTModel.list[[5]], name="umi.assay") <- "Xenium"
slot(object = cm@assays$SCT@SCTModel.list[[6]], name="umi.assay") <- "Xenium"
slot(object = cm@assays$SCT@SCTModel.list[[7]], name="umi.assay") <- "Xenium"
slot(object = cm@assays$SCT@SCTModel.list[[8]], name="umi.assay") <- "Xenium"
slot(object = cm@assays$SCT@SCTModel.list[[9]], name="umi.assay") <- "Xenium"
SCTResults(object=cm, slot="umi.assay")
cm<-PrepSCTFindMarkers(cm)

# Find cluster markers
m0 <- FindAllMarkers(cm,recorrect_umi=TRUE) 
m1 <- subset(m0,p_val_adj < 0.05 & avg_log2FC > 0)
m2 <- split(m1, m1$cluster)
m3 <- lapply(m2,function(x){paste(unlist(na.omit(x[1:5,]$gene)),collapse='_')})


cm$subnames <- paste0('Cm_',cm@active.ident)
cm$markernames <- paste0('Cm_',unlist(m3[cm@active.ident]))

#saveRDS(cm,'~/Downloads/hdWGCNA_TOM/Xenium/cm.rds')

#Assign back to parent object
idx <- match(colnames(cm),colnames(xenium.obj))
xenium.obj$subnames[idx] <- cm$subnames
xenium.obj$markernames[idx] <- cm$markernames


### Fibroblasts


fb <- subset(xenium.obj,idents=c('FB'))

fb <- SCTransform(fb, assay = "Xenium")
fb <- RunPCA(fb, npcs = 30, features = rownames(fb))
fb <- RunHarmony(fb,'patient')
fb <- RunUMAP(fb, reduction = "harmony", dims = 1:30)
fb <- FindNeighbors(fb, reduction = "harmony", dims = 1:30)
fb <- FindClusters(fb, resolution = 0.5)

#Rename slots to recorrect UMI
slot(object = fb@assays$SCT@SCTModel.list[[1]], name="umi.assay") <- "Xenium"
slot(object = fb@assays$SCT@SCTModel.list[[2]], name="umi.assay") <- "Xenium"
slot(object = fb@assays$SCT@SCTModel.list[[3]], name="umi.assay") <- "Xenium"
slot(object = fb@assays$SCT@SCTModel.list[[4]], name="umi.assay") <- "Xenium"
slot(object = fb@assays$SCT@SCTModel.list[[5]], name="umi.assay") <- "Xenium"
slot(object = fb@assays$SCT@SCTModel.list[[6]], name="umi.assay") <- "Xenium"
slot(object = fb@assays$SCT@SCTModel.list[[7]], name="umi.assay") <- "Xenium"
slot(object = fb@assays$SCT@SCTModel.list[[8]], name="umi.assay") <- "Xenium"
slot(object = fb@assays$SCT@SCTModel.list[[9]], name="umi.assay") <- "Xenium"
SCTResults(object=fb, slot="umi.assay")
fb<-PrepSCTFindMarkers(fb)

# Find cluster markers
m0 <- FindAllMarkers(fb,recorrect_umi=TRUE) 
m1 <- subset(m0,p_val_adj < 0.05 & avg_log2FC > 0)
m2 <- split(m1, m1$cluster)
m3 <- lapply(m2,function(x){paste(unlist(na.omit(x[1:5,]$gene)),collapse='_')})


fb$subnames <- paste0('Fb_',fb@active.ident)
fb$markernames <- paste0('Fb_',unlist(m3[fb@active.ident]))

#saveRDS(fb,'~/Downloads/hdWGCNA_TOM/Xenium/fb.rds')

#Assign back to parent object
idx <- match(colnames(fb),colnames(xenium.obj))
xenium.obj$subnames[idx] <- fb$subnames
xenium.obj$markernames[idx] <- fb$markernames

### EC

ec <- subset(xenium.obj,idents=c('EC'))

ec <- SCTransform(ec, assay = "Xenium")
ec <- RunPCA(ec, npcs = 30, features = rownames(ec))
ec <- RunHarmony(ec,'patient')
ec <- RunUMAP(ec, reduction = "harmony", dims = 1:30)
ec <- FindNeighbors(ec, reduction = "harmony", dims = 1:30)
ec <- FindClusters(ec, resolution = 0.5)

#Rename slots to recorrect UMI
slot(object = ec@assays$SCT@SCTModel.list[[1]], name="umi.assay") <- "Xenium"
slot(object = ec@assays$SCT@SCTModel.list[[2]], name="umi.assay") <- "Xenium"
slot(object = ec@assays$SCT@SCTModel.list[[3]], name="umi.assay") <- "Xenium"
slot(object = ec@assays$SCT@SCTModel.list[[4]], name="umi.assay") <- "Xenium"
slot(object = ec@assays$SCT@SCTModel.list[[5]], name="umi.assay") <- "Xenium"
slot(object = ec@assays$SCT@SCTModel.list[[6]], name="umi.assay") <- "Xenium"
slot(object = ec@assays$SCT@SCTModel.list[[7]], name="umi.assay") <- "Xenium"
slot(object = ec@assays$SCT@SCTModel.list[[8]], name="umi.assay") <- "Xenium"
slot(object = ec@assays$SCT@SCTModel.list[[9]], name="umi.assay") <- "Xenium"
SCTResults(object=ec, slot="umi.assay")
ec<-PrepSCTFindMarkers(ec)

# Find cluster markers
m0 <- FindAllMarkers(ec,recorrect_umi=TRUE) 
m1 <- subset(m0,p_val_adj < 0.05 & avg_log2FC > 0)
m2 <- split(m1, m1$cluster)
m3 <- lapply(m2,function(x){paste(unlist(na.omit(x[1:5,]$gene)),collapse='_')})


ec$subnames <- paste0('Ec_',ec@active.ident)
ec$markernames <- paste0('Ec_',unlist(m3[ec@active.ident]))

#saveRDS(ec,'~/Downloads/hdWGCNA_TOM/Xenium/ec.rds')

#Assign back to parent object
idx <- match(colnames(ec),colnames(xenium.obj))
xenium.obj$subnames[idx] <- ec$subnames
xenium.obj$markernames[idx] <- ec$markernames

### Myeloid

myeloid <- subset(xenium.obj,idents=c('Myeloid'))

myeloid <- SCTransform(myeloid, assay = "Xenium")
myeloid <- RunPCA(myeloid, npcs = 30, features = rownames(myeloid))
myeloid <- RunHarmony(myeloid,'patient')
myeloid <- RunUMAP(myeloid, reduction = "harmony", dims = 1:30)
myeloid <- FindNeighbors(myeloid, reduction = "harmony", dims = 1:30)
myeloid <- FindClusters(myeloid, resolution = 0.5)

#Rename slots to recorrect UMI
slot(object = myeloid@assays$SCT@SCTModel.list[[1]], name="umi.assay") <- "Xenium"
slot(object = myeloid@assays$SCT@SCTModel.list[[2]], name="umi.assay") <- "Xenium"
slot(object = myeloid@assays$SCT@SCTModel.list[[3]], name="umi.assay") <- "Xenium"
slot(object = myeloid@assays$SCT@SCTModel.list[[4]], name="umi.assay") <- "Xenium"
slot(object = myeloid@assays$SCT@SCTModel.list[[5]], name="umi.assay") <- "Xenium"
slot(object = myeloid@assays$SCT@SCTModel.list[[6]], name="umi.assay") <- "Xenium"
slot(object = myeloid@assays$SCT@SCTModel.list[[7]], name="umi.assay") <- "Xenium"
slot(object = myeloid@assays$SCT@SCTModel.list[[8]], name="umi.assay") <- "Xenium"
slot(object = myeloid@assays$SCT@SCTModel.list[[9]], name="umi.assay") <- "Xenium"
SCTResults(object=myeloid, slot="umi.assay")
myeloid<-PrepSCTFindMarkers(myeloid)

# Find cluster markers
m0 <- FindAllMarkers(myeloid,recorrect_umi=TRUE) 
m1 <- subset(m0,p_val_adj < 0.05 & avg_log2FC > 0)
m2 <- split(m1, m1$cluster)
m3 <- lapply(m2,function(x){paste(unlist(na.omit(x[1:5,]$gene)),collapse='_')})


myeloid$subnames <- paste0('Myeloid_',myeloid@active.ident)
myeloid$markernames <- paste0('Myeloid_',unlist(m3[myeloid@active.ident]))

#saveRDS(myeloid,'~/Downloads/hdWGCNA_TOM/Xenium/myeloid.rds')

#Assign back to parent object
idx <- match(colnames(myeloid),colnames(xenium.obj))
xenium.obj$subnames[idx] <- myeloid$subnames
xenium.obj$markernames[idx] <- myeloid$markernames

### NKT

nkt <- subset(xenium.obj,idents=c('NKT'))

nkt <- SCTransform(nkt, assay = "Xenium")
nkt <- RunPCA(nkt, npcs = 30, features = rownames(nkt))
nkt <- RunHarmony(nkt,'patient')
nkt <- RunUMAP(nkt, reduction = "harmony", dims = 1:30)
nkt <- FindNeighbors(nkt, reduction = "harmony", dims = 1:30)
nkt <- FindClusters(nkt, resolution = 0.5)

#Rename slots to recorrect UMI
slot(object = nkt@assays$SCT@SCTModel.list[[1]], name="umi.assay") <- "Xenium"
slot(object = nkt@assays$SCT@SCTModel.list[[2]], name="umi.assay") <- "Xenium"
slot(object = nkt@assays$SCT@SCTModel.list[[3]], name="umi.assay") <- "Xenium"
slot(object = nkt@assays$SCT@SCTModel.list[[4]], name="umi.assay") <- "Xenium"
slot(object = nkt@assays$SCT@SCTModel.list[[5]], name="umi.assay") <- "Xenium"
slot(object = nkt@assays$SCT@SCTModel.list[[6]], name="umi.assay") <- "Xenium"
slot(object = nkt@assays$SCT@SCTModel.list[[7]], name="umi.assay") <- "Xenium"
slot(object = nkt@assays$SCT@SCTModel.list[[8]], name="umi.assay") <- "Xenium"
slot(object = nkt@assays$SCT@SCTModel.list[[9]], name="umi.assay") <- "Xenium"
SCTResults(object=nkt, slot="umi.assay")
nkt<-PrepSCTFindMarkers(nkt)

# Find cluster markers
m0 <- FindAllMarkers(nkt,recorrect_umi=TRUE) 
m1 <- subset(m0,p_val_adj < 0.05 & avg_log2FC > 0)
m2 <- split(m1, m1$cluster)
m3 <- lapply(m2,function(x){paste(unlist(na.omit(x[1:5,]$gene)),collapse='_')})


nkt$subnames <- paste0('NKT_',nkt@active.ident)
nkt$markernames <- paste0('NKT_',unlist(m3[nkt@active.ident]))

#saveRDS(nkt,'~/Downloads/hdWGCNA_TOM/Xenium/nkt.rds')

#Assign back to parent object
idx <- match(colnames(nkt),colnames(xenium.obj))
xenium.obj$subnames[idx] <- nkt$subnames
xenium.obj$markernames[idx] <- nkt$markernames

### Adipo

adipo <- subset(xenium.obj,idents=c('Adipo'))

adipo <- SCTransform(adipo, assay = "Xenium")
adipo <- RunPCA(adipo, npcs = 30, features = rownames(adipo))
adipo <- RunHarmony(adipo,'patient')
adipo <- RunUMAP(adipo, reduction = "harmony", dims = 1:30)
adipo <- FindNeighbors(adipo, reduction = "harmony", dims = 1:30)
adipo <- FindClusters(adipo, resolution = 0.5)

#Rename slots to recorrect UMI
slot(object = adipo@assays$SCT@SCTModel.list[[1]], name="umi.assay") <- "Xenium"
slot(object = adipo@assays$SCT@SCTModel.list[[2]], name="umi.assay") <- "Xenium"
slot(object = adipo@assays$SCT@SCTModel.list[[3]], name="umi.assay") <- "Xenium"
slot(object = adipo@assays$SCT@SCTModel.list[[4]], name="umi.assay") <- "Xenium"
slot(object = adipo@assays$SCT@SCTModel.list[[5]], name="umi.assay") <- "Xenium"
slot(object = adipo@assays$SCT@SCTModel.list[[6]], name="umi.assay") <- "Xenium"
slot(object = adipo@assays$SCT@SCTModel.list[[7]], name="umi.assay") <- "Xenium"
slot(object = adipo@assays$SCT@SCTModel.list[[8]], name="umi.assay") <- "Xenium"
slot(object = adipo@assays$SCT@SCTModel.list[[9]], name="umi.assay") <- "Xenium"
SCTResults(object=adipo, slot="umi.assay")
adipo<-PrepSCTFindMarkers(adipo)

# Find cluster markers
m0 <- FindAllMarkers(adipo,recorrect_umi=TRUE) 
m1 <- subset(m0,p_val_adj < 0.05 & avg_log2FC > 0)
m2 <- split(m1, m1$cluster)
m3 <- lapply(m2,function(x){paste(unlist(na.omit(x[1:5,]$gene)),collapse='_')})


adipo$subnames <- paste0('Adipo_',adipo@active.ident)
adipo$markernames <- paste0('Adipo_',unlist(m3[adipo@active.ident]))

#saveRDS(adipo,'~/Downloads/hdWGCNA_TOM/Xenium/Adipo.rds')

#Assign back to parent object
idx <- match(colnames(adipo),colnames(xenium.obj))
xenium.obj$subnames[idx] <- adipo$subnames
xenium.obj$markernames[idx] <- adipo$markernames


### PC

pc <- subset(xenium.obj,idents=c('PC'))

pc <- SCTransform(pc, assay = "Xenium")
pc <- RunPCA(pc, npcs = 30, features = rownames(pc))
pc <- RunHarmony(pc,'patient')
pc <- RunUMAP(pc, reduction = "harmony", dims = 1:30)
pc <- FindNeighbors(pc, reduction = "harmony", dims = 1:30)
pc <- FindClusters(pc, resolution = 0.5)

#Rename slots to recorrect UMI
slot(object = pc@assays$SCT@SCTModel.list[[1]], name="umi.assay") <- "Xenium"
slot(object = pc@assays$SCT@SCTModel.list[[2]], name="umi.assay") <- "Xenium"
slot(object = pc@assays$SCT@SCTModel.list[[3]], name="umi.assay") <- "Xenium"
slot(object = pc@assays$SCT@SCTModel.list[[4]], name="umi.assay") <- "Xenium"
slot(object = pc@assays$SCT@SCTModel.list[[5]], name="umi.assay") <- "Xenium"
slot(object = pc@assays$SCT@SCTModel.list[[6]], name="umi.assay") <- "Xenium"
slot(object = pc@assays$SCT@SCTModel.list[[7]], name="umi.assay") <- "Xenium"
slot(object = pc@assays$SCT@SCTModel.list[[8]], name="umi.assay") <- "Xenium"
slot(object = pc@assays$SCT@SCTModel.list[[9]], name="umi.assay") <- "Xenium"
SCTResults(object=pc, slot="umi.assay")
pc<-PrepSCTFindMarkers(pc)

# Find cluster markers
m0 <- FindAllMarkers(pc,recorrect_umi=TRUE) 
m1 <- subset(m0,p_val_adj < 0.05 & avg_log2FC > 0)
m2 <- split(m1, m1$cluster)
m3 <- lapply(m2,function(x){paste(unlist(na.omit(x[1:5,]$gene)),collapse='_')})


pc$subnames <- paste0('Pc_',pc@active.ident)
pc$markernames <- paste0('Pc_',unlist(m3[pc@active.ident]))

#saveRDS(pc,'~/Downloads/hdWGCNA_TOM/Xenium/Pc.rds')

#Assign back to parent object
idx <- match(colnames(pc),colnames(xenium.obj))
xenium.obj$subnames[idx] <- pc$subnames
xenium.obj$markernames[idx] <- pc$markernames


### SM

sm <- subset(xenium.obj,idents=c('SM'))

sm <- SCTransform(sm, assay = "Xenium")
sm <- RunPCA(sm, npcs = 30, features = rownames(sm))
sm <- RunHarmony(sm,'patient')
sm <- RunUMAP(sm, reduction = "harmony", dims = 1:30)
sm <- FindNeighbors(sm, reduction = "harmony", dims = 1:30)
sm <- FindClusters(sm, resolution = 0.5)

#Rename slots to recorrect UMI
slot(object = sm@assays$SCT@SCTModel.list[[1]], name="umi.assay") <- "Xenium"
slot(object = sm@assays$SCT@SCTModel.list[[2]], name="umi.assay") <- "Xenium"
slot(object = sm@assays$SCT@SCTModel.list[[3]], name="umi.assay") <- "Xenium"
slot(object = sm@assays$SCT@SCTModel.list[[4]], name="umi.assay") <- "Xenium"
slot(object = sm@assays$SCT@SCTModel.list[[5]], name="umi.assay") <- "Xenium"
slot(object = sm@assays$SCT@SCTModel.list[[6]], name="umi.assay") <- "Xenium"
slot(object = sm@assays$SCT@SCTModel.list[[7]], name="umi.assay") <- "Xenium"
slot(object = sm@assays$SCT@SCTModel.list[[8]], name="umi.assay") <- "Xenium"
slot(object = sm@assays$SCT@SCTModel.list[[9]], name="umi.assay") <- "Xenium"
SCTResults(object=sm, slot="umi.assay")
sm<-PrepSCTFindMarkers(sm)

# Find cluster markers
m0 <- FindAllMarkers(sm,recorrect_umi=TRUE) 
m1 <- subset(m0,p_val_adj < 0.05 & avg_log2FC > 0)
m2 <- split(m1, m1$cluster)
m3 <- lapply(m2,function(x){paste(unlist(na.omit(x[1:5,]$gene)),collapse='_')})


sm$subnames <- paste0('Pc_',sm@active.ident)
sm$markernames <- paste0('Pc_',unlist(m3[sm@active.ident]))

#saveRDS(sm,'~/Downloads/hdWGCNA_TOM/Xenium/SM.rds')

#Assign back to parent object
idx <- match(colnames(sm),colnames(xenium.obj))
xenium.obj$subnames[idx] <- sm$subnames
xenium.obj$markernames[idx] <- sm$markernames




### LEC

lec <- subset(xenium.obj,idents=c('LEC'))

lec <- SCTransform(lec, assay = "Xenium")
lec <- RunPCA(lec, npcs = 30, features = rownames(lec))
lec <- RunHarmony(lec,'patient')
lec <- RunUMAP(lec, reduction = "harmony", dims = 1:30)
lec <- FindNeighbors(lec, reduction = "harmony", dims = 1:30)
lec <- FindClusters(lec, resolution = 0.5)

#Rename slots to recorrect UMI
slot(object = lec@assays$SCT@SCTModel.list[[1]], name="umi.assay") <- "Xenium"
slot(object = lec@assays$SCT@SCTModel.list[[2]], name="umi.assay") <- "Xenium"
slot(object = lec@assays$SCT@SCTModel.list[[3]], name="umi.assay") <- "Xenium"
slot(object = lec@assays$SCT@SCTModel.list[[4]], name="umi.assay") <- "Xenium"
slot(object = lec@assays$SCT@SCTModel.list[[5]], name="umi.assay") <- "Xenium"
slot(object = lec@assays$SCT@SCTModel.list[[6]], name="umi.assay") <- "Xenium"
slot(object = lec@assays$SCT@SCTModel.list[[7]], name="umi.assay") <- "Xenium"
slot(object = lec@assays$SCT@SCTModel.list[[8]], name="umi.assay") <- "Xenium"
slot(object = lec@assays$SCT@SCTModel.list[[9]], name="umi.assay") <- "Xenium"
SCTResults(object=lec, slot="umi.assay")
lec<-PrepSCTFindMarkers(lec)

# Find cluster markers
m0 <- FindAllMarkers(lec,recorrect_umi=TRUE) 
m1 <- subset(m0,p_val_adj < 0.05 & avg_log2FC > 0)
m2 <- split(m1, m1$cluster)
m3 <- lapply(m2,function(x){paste(unlist(na.omit(x[1:5,]$gene)),collapse='_')})


lec$subnames <- paste0('Pc_',lec@active.ident)
lec$markernames <- paste0('Pc_',unlist(m3[lec@active.ident]))

#saveRDS(lec,'~/Downloads/hdWGCNA_TOM/Xenium/LEC.rds')

#Assign back to parent object
idx <- match(colnames(lec),colnames(xenium.obj))
xenium.obj$subnames[idx] <- lec$subnames
xenium.obj$markernames[idx] <- lec$markernames


#saveRDS(xenium.obj,'~/Downloads/hdWGCNA_TOM/Xenium/collated_xenium_data_FINAL.rds')



##############################################
##############################################
#### Figure XA
##############################################
##############################################

library(Seurat)
library(hdWGCNA)
library(ggeasy)


source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')

M1<-readRDS('~/Downloads/hdWGCNA_TOM/Xenium/collated_xenium_data_FINAL.rds')

pdf(paste0('~/Downloads/hdWGCNA_TOM/Xenium/', 'Xenium_snUMAP.pdf'), width=5, height=5)
PlotEmbedding(M1,group.by='names',point_size=.05,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()

##############################################
##############################################
#### Figure XB
##############################################
##############################################

M1<-readRDS('~/Downloads/hdWGCNA_TOM/Xenium/collated_xenium_data_FINAL.rds')

pdf(paste0('~/Downloads/hdWGCNA_TOM/Xenium/', 'Xenium_1697_celltype.pdf'), width=4, height=5)
ImageDimPlot(M1, fov = "fov.2",split.by='names',cols = "polychrome", size =0.35,
	axes = FALSE, dark.background = FALSE) + NoLegend()
dev.off()



##############################################
##############################################
#### Figure XC
##############################################
##############################################


# Recalculate just for broad celltypes

# Calculate neighbors

for ( i in seq_along(M1_list) ){ message(i)
	object = M1_list[[i]]
	fov = fov_list[[i]]
	object <- BuildNicheAssay(object = object, fov = fov, group.by = "names",niches.k = 10, neighbors.k = 30)
	M1_list[[i]] = object
}


# Concatenate neighbors across all samples

col.unique = lapply( seq_along(M1_list), function(i){
    colnames(t( M1_list[[i]][['niche']]@counts ))
}  )
shared_feats <- Reduce(intersect,col.unique)


DAT.counts = lapply( seq_along(M1_list), function(i){
    t( M1_list[[i]][['niche']]@counts )[,shared_feats]
}  )
DAT.counts <- do.call("rbind", DAT.counts)

DAT.data = lapply( seq_along(M1_list), function(i){
    t( M1_list[[i]][['niche']]@data )[,shared_feats]
}  )
DAT.data <- do.call("rbind", DAT.data)

DAT.scale.data = lapply( seq_along(M1_list), function(i){
    t( M1_list[[i]][['niche']]@scale.data )[,shared_feats]
}  )
DAT.scale.data <- do.call("rbind", DAT.scale.data)

niche.assay <- CreateAssayObject(counts = t(DAT.counts))
niche.assay@data <- t(DAT.data)
niche.assay@scale.data <- t(DAT.scale.data)

M1[['niche_broad']] <- niche.assay
DefaultAssay(M1) <- 'niche_broad'

niches.k.range = 2:30

res.clusters = data.frame(row.names = rownames(DAT.scale.data))

for ( k in niches.k.range ){ message("k=", k)
    # new column name
    newCol = paste0("kmeans_", k)
    # get centroids
    km_mb = ClusterR::MiniBatchKmeans(
        "data" = DAT.scale.data
        , "clusters" = k # the number of clusters
        , "batch_size" = 20 # the size of the mini batches
        , "num_init" = 20 # number of times the algorithm will be run with different centroid seeds
        , "max_iters" = 100 # the maximum number of clustering iterations. 
        , "init_fraction" = 0.2 # percentage of data to use for the initialization centroids (applies if initializer is kmeans++ or optimal_init). Should be a float number between 0.0 and 1.0.
        , "initializer" = "kmeans++" # the method of initialization. One of, optimal_init, quantile_init, kmeans++ and random. See details for more information
        , "early_stop_iter" = 10 # continue that many iterations after calculation of the best within-cluster-sum-of-squared-error
        , "verbose" = F
        , "CENTROIDS" = NULL
        , "tol" = 1e-04
        , "tol_optimal_init" = 0.3
        , "seed" = 1
    )
    
    # use centroids to get clusters

    res.clusters[,newCol] = ClusterR::predict_MBatchKMeans( # This function takes the data and the output centroids and returns the clusters.
        "data" = DAT.scale.data
        , "CENTROIDS" = km_mb$centroids
    )
    res.clusters[,newCol] = as.factor( res.clusters[,newCol] ) # change clusters to factors
    
}


# Update object
res.clusters.ordered <- res.clusters[colnames(M1),]
for(i in colnames(res.clusters.ordered))
	{
		M1[[i]] = 0
		M1[[i]][rownames(res.clusters.ordered),] = res.clusters.ordered[,i]

	}

write.table(res.clusters.ordered,'~/Downloads/hdWGCNA_TOM/Xenium/Niche_bulk_clusters.csv',sep=',')

niche.patient <- table(M1$kmeans_15,M1$patient)
niche.patient <- t(t(niche.patient) / colSums(niche.patient))

disease <- c('RVF','RVF','NF','pRV','pRV','RVF','NF','pRV','NF')
disease <- c(t(replicate(15,disease)))

niche.patient <- data.frame(niche = disease,niche.patient)

niche.patient$niche <- factor(niche.patient$niche, levels=c('NF','pRV','RVF'))

ggplot(niche.patient,aes(Var1,Freq,color = niche))+geom_boxplot() + theme_classic()



niche.names <- table(M1$kmeans_15,M1$names)
niche.names <- niche.names / rowSums(niche.names)

niche.names_t <- table(M1$kmeans_15,M1$names)
niche.names_t <- t(t(niche.names_t) / colSums(niche.names_t))

#1 - SM, FB, EC
#2 - FB, LEC, EC
#3 - FB, Myeloid with T cells and absent CMs
#4 - CM, FB
#5 - EC, CM, FB with T cells
#6 - CM, EC, PC
#7 - EC, FB, CM with T cells
#8 - FB, CM, EC
#9 - CM
#10 - Adipo, FB
#11 - CM, EC
#12 - EC, CM, FB
#13 - SM, EC
#14 - CM, EC, and PC
#15 - CM, EC with T cells

#1 and 13 are SM/EC
#2 are LEC
#3 are Myeloid
#4 and 8 are CM/FB
#5 and 6 and 7 and 8 and 12 are CM/FB/EC
#9 are CM
#10 are adipo
#11 and 14 and 15 are CM/EC

library(scCustomize)

# CMs are in 7, 9, 14, 15, 4, 5, 6, 8, 11, and 12
ImageDimPlot(M1, fov = "fov.2", group.by = "kmeans_15", 
	size = .5, dark.background = F, cells = WhichCells(M1,expression = kmeans_15 %in%  c('7','9','14','15')))

ImageDimPlot(M1, fov = "fov.2", group.by = "names", 
	size = .5, dark.background = F, cells = WhichCells(M1,expression = kmeans_15 %in%  c('15')))


niche_names <- c()
#15 has less adipose and FBs and more CMs and ECs nearby, these vessels are RGCC rich
FindMarkers(M1,ident.1= '15',ident.2=c('7'),only.pos=T)
FindMarkers(M2,ident.1 = '15',ident.2=c('7'),recorrect_umi=F,only.pos=T)
ImageDimPlot(M1, fov = "fov.2", group.by = "kmeans_15", 
	size = .5, dark.background = F, cells = WhichCells(M1,expression = kmeans_15 %in%  c('7','15','1')))
DefaultAssay(M1) <- 'SCT'
M2 <- subset(M1,names=='EC')
M2 <- SetIdent(M2,value="kmeans_15")
niche_names['15'] <- 'Perivascular CMs'


#4
ImageDimPlot(M1, fov = "fov.2", group.by = "kmeans_15", 
	size = .5, dark.background = F, cells = WhichCells(M1,expression = kmeans_15 %in%  c('5','1','7','15')))
FindMarkers(M1,ident.1 = '4',ident.2=c('9','14','15','5'),recorrect_umi=F,only.pos=T)
FindMarkers(M1,ident.1 = '4',ident.2=c('5'),recorrect_umi=F,only.pos=T)
FindMarkers(M1,ident.1 = '5',ident.2=c('4'),recorrect_umi=F,only.pos=T)
M3 <- subset(M1,names=='NKT')
M4 <- subset(M1,names=='Myeloid')
FindMarkers(M3,ident.1 = '5',ident.2=c('4'),recorrect_umi=F,only.pos=T)
FindMarkers(M3,ident.1 = '4',ident.2=c('5'),recorrect_umi=F,only.pos=T)
FindMarkers(M4,ident.1 = '5',ident.2=c('4'),recorrect_umi=F,only.pos=T)
FindMarkers(M4,ident.1 = '4',ident.2=c('5'),recorrect_umi=F,only.pos=T)


niche_names['4'] <- 'Isolated CMs'

#5
niche_names['5'] <- 'Peri-immune CMs'


#6
ImageDimPlot(M1, fov = "fov.2", group.by = "kmeans_15", 
	size = .5, dark.background = F, 
	cells = WhichCells(M1,expression = kmeans_15 %in%  c('4','6')),molecules='THBS4')
M5 <- subset(M1,names == 'CM')
FindMarkers(M1,ident.1 = '6',ident.2='8',recorrect_umi=F,only.pos=T)

niche_names['6'] <- 'Stromal CMs/THBS4 FBs'


#8 NPPA and NPPB
FindMarkers(M1,ident.1 = '8',recorrect_umi=F,only.pos=T)
niche_names['8'] <- 'Stroma'
ImageDimPlot(M1, fov = "fov.2", group.by = "kmeans_15", 
	size = .5, dark.background = F, 
	cells = WhichCells(M1,expression = kmeans_15 %in%  c('8')))


#11
FindMarkers(M1,ident.1 = '11',ident.2='4',recorrect_umi=F,only.pos=T)
ImageDimPlot(M1, fov = "fov.2", group.by = "kmeans_15", 
	size = .5, dark.background = F, 
	cells = WhichCells(M1,expression = kmeans_15 %in%  c('11','4')))
niche_names['11'] <- 'Isolated CMs'


#12
FindMarkers(M1,ident.1 = '12',ident.2='11',recorrect_umi=F,only.pos=T)

ImageDimPlot(M1, fov = "fov.2", group.by = "kmeans_15", 
	size = .5, dark.background = F, 
	cells = WhichCells(M1,expression = kmeans_15 %in%  c('12','2')))
niche_names['12'] <- 'Peri-lymphatic/venous CMs'


#9 is endocardial CMs
marks.9 <- FindMarkers(M1,ident.1= '9',ident.2=c('7','14','15'))
ImageDimPlot(M1, fov = "fov.8", group.by = "names", 
	size = .5, dark.background = F, molecules = 'NPPB',
	cells = WhichCells(M1,expression = kmeans_15 %in%  c('9')))
ImageDimPlot(M1, fov = "fov.7", group.by = "names", 
	size = .5, dark.background = F,molecules = 'NPPB')
M1 <- SetIdent(M1,value="kmeans_15")

FindMarkers(M1,ident.1 = '9',ident.2=c('7','15','14'),recorrect_umi=T,only.pos=T)
niche_names['9'] <- 'Endocardial CMs'


#14 is surrounded by FBs, adipocytes, ECs, and myeloid
marks.14 <- FindMarkers(M1,ident.1= '14',ident.2=c('7','9','15'))
FindMarkers(M1,ident.1 = '14',ident.2=c('7'),recorrect_umi=T,only.pos=T)
ImageDimPlot(M1, fov = "fov.7", group.by = "names", 
	size = .5, dark.background = F,
	cells = WhichCells(M1,expression = kmeans_15 %in%  c('14')))
niche_names['14'] <- 'Stromal CMs'



#1 Vasculature (smooth muscle poor)
FindMarkers(M2,ident.1 = '1',ident.2=c('13'),recorrect_umi=F,only.pos=T)
ImageDimPlot(M1, fov = "fov.7", group.by = "names", 
	size = .5, dark.background = F,
	cells = WhichCells(M1,expression = kmeans_15 %in%  c('1')))
niche_names['1'] <- 'Arterioles'

#13 Macrovasculature  (smooth muscle rich)
FindMarkers(M1,ident.1 = '13',ident.2=c('1'),recorrect_umi=F,only.pos=T)
ImageDimPlot(M1, fov = "fov.4", group.by = "names", 
	size = .5, dark.background = F,
	cells = WhichCells(M1,expression = kmeans_15 %in%  c('13')))
niche_names['13'] <- 'Macrovasculature'

#2 LEC and peri-lymphatic stroma
ImageDimPlot(M1, fov = "fov.6", group.by = "names", 
	size = .5, dark.background = F,
	cells = WhichCells(M1,expression = kmeans_15 %in%  c('2')))
ImageDimPlot(M1, fov = "fov.6", group.by = "kmeans_15", 
	size = .5, dark.background = F,
	cells = WhichCells(M1,expression = kmeans_15 %in%  c('2','3','7','10','13')))
FindMarkers(M1,ident.1 = '2',ident.2=c('3'),recorrect_umi=F,only.pos=T)
niche_names['2'] <- 'Lymphatics and peri-lymphatic stroma'


#3 Adipose-rich Stroma
FindMarkers(M1,ident.1 = '3',ident.2=c('2'),recorrect_umi=F,only.pos=T)
niche_names['3'] <- 'Adipose-rich stroma'
ImageDimPlot(M1, fov = "fov.6", group.by = "names", 
	size = .5, dark.background = F,
	cells = WhichCells(M1,expression = kmeans_15 %in%  c('3')))

#7 Capillary-rich
# More RGCC and CD36 than 1 or 13
FindMarkers(M2,ident.1 = '7',ident.2=c('1'),recorrect_umi=F,only.pos=T)
FindMarkers(M2,ident.1 = '7',ident.2=c('13'),recorrect_umi=F,only.pos=T)
niche_names['7'] <- 'Capillaries'


#1 is more arterial
FindMarkers(M2,ident.1 = '1',ident.2=c('7'),recorrect_umi=F,only.pos=T)

ImageDimPlot(M1, fov = "fov.7", group.by = "kmeans_15", 
	size = .5, dark.background = F,
	cells = WhichCells(M1,expression = kmeans_15 %in%  c('7','15')))

#10 Adipose
niche_names['10'] <- 'Adipose'
FindMarkers(M1,ident.1 = '10',ident.2=c('3'),recorrect_umi=F,only.pos=T)
FindMarkers(M1,ident.1 = '3',ident.2=c('10'),recorrect_umi=F,only.pos=T)

ImageDimPlot(M1, fov = "fov.5", group.by = "kmeans_15", 
	size = .5, dark.background = F,
	cells = WhichCells(M1,expression = kmeans_15 %in%  c('15','10','3')))


new.cluster.ids<-niche_names[levels(M1)]
names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)
M1$niche_manual <- M1@active.ident

pdf('~/Downloads/hdWGCNA_TOM/Xenium/niche_fov2.pdf',width=5,height=5)
ImageDimPlot(M1, fov = "fov.2", group.by = "niche_manual", 
	size = .25, dark.background = F)
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/Xenium/niche_fov9.pdf',width=5,height=5)
ImageDimPlot(M1, fov = "fov.9", group.by = "niche_manual", 
	size = .25, dark.background = F)
dev.off()


write.table(M1@meta.data,'~/Downloads/hdWGCNA_TOM/Xenium/metadata.csv',sep=',')

##############################################
##############################################
#### Figure XD
##############################################
##############################################


niche.patient <- table(M1$niche_manual,M1$patient)
niche.patient <- t(t(niche.patient) / colSums(niche.patient))

disease <- c('RVF','RVF','NF','pRV','pRV','RVF','NF','pRV','NF')
disease <- c(t(replicate(14,disease)))

niche.patient <- data.frame(niche = disease,niche.patient)

niche.patient$niche <- factor(niche.patient$niche, levels=c('NF','pRV','RVF'))
niche.patient$Var1 <- factor(niche.patient$Var1, levels=rev(names(sort(niche.counts))))

niche.counts <- table(M1$niche_manual)


pdf('~/Downloads/hdWGCNA_TOM/Xenium/niche_counts.pdf',width=12,height=5)
ggplot(niche.patient,aes(Var1,Freq,color = niche))+geom_boxplot() + theme_classic()
dev.off()

library(ggpubr)
pdf('~/Downloads/hdWGCNA_TOM/Xenium/niche_clust_freq_stats.pdf',width=12.5,height=15)
p <- ggboxplot(niche.patient,x="niche",y="Freq",fill="niche",group="niche")+
	theme_classic() + 
	theme(axis.text.x=element_text(size=16),
	axis.text.y=element_text(size=16),
	axis.title.x=element_text(size=16),
	axis.title.y=element_text(size=16),
	legend.title=element_text(size=16),
	legend.text=element_text(size=16),
	text=element_text(color='black'),
	axis.text=element_text(color='black')) + 
	facet_wrap(~Var1,ncol=7) + 
	#stat_compare_means(aes(group=group),comparisons=my_comparisons,method="t.test",ref.group="NF")+
	stat_compare_means(aes(group=niche),method="anova")
p
dev.off()



##############################################
##############################################
#### Figure XE
##############################################
##############################################


M2 <- subset(M1,niche_manual=='Arterioles')
M2 <- SetIdent(M2,value="group")
FindMarkers(M2,ident.1='RVF',ident.2='pRV',recorrect_umi=F)


M2 <- subset(M1,niche_manual=='Capillaries')
M2 <- SetIdent(M2,value="group")
FindMarkers(M2,ident.1='RVF',ident.2='pRV',recorrect_umi=F)


M2 <- subset(M1,niche_manual=='Peri-immune CMs')
M2 <- SetIdent(M2,value="group")
FindMarkers(M2,ident.1='RVF',ident.2='pRV',recorrect_umi=F)

M1 <- SetIdent(M1,value="group")
markers.RVF <- FindMarkers(M1,ident.1='RVF',ident.2='pRV',recorrect_umi=T,assay='Xenium')
RVF_vs_pRV_25_up <- rownames(subset(markers.RVF,avg_log2FC > 0 & p_val_adj < 0.05))[1:25]
RVF_vs_pRV_25_down <- rownames(subset(markers.RVF,avg_log2FC < 0 & p_val_adj < 0.05))[1:25]
markers.pRV <- FindMarkers(M1,ident.1='pRV',ident.2='NF',recorrect_umi=T)
pRV_vs_NF_25_up <- rownames(subset(markers.pRV,avg_log2FC > 0 & p_val_adj < 0.05))[1:200]
pRV_vs_NF_25_down <- rownames(subset(markers.pRV,avg_log2FC < 0 & p_val_adj < 0.05))[1:200]

up_int <- intersect(RVF_vs_pRV_25_up,pRV_vs_NF_25_up)

down_int <- intersect(RVF_vs_pRV_25_down,pRV_vs_NF_25_down)


to_plot <- c(RVF_vs_pRV_25_up,RVF_vs_pRV_25_down)

library(EnhancedVolcano)
EnhancedVolcano(markers.RVF,lab = rownames(markers.RVF), x = 'avg_log2FC',y = 'p_val_adj',pCutoff = 0.05,FCcutoff=0.25)

M1.aggreg <- AggregateExpression(object = M1, group.by = c('patient'),assays='Xenium')


#Pseudobulk
#M1.aggreg <- AggregateExpression(object = M1, group.by = c('patient'), return.seurat = T, assays='Xenium')
mycol <- colorpanel(1000,"blue","white","red")
#disease <- c('RVF','RVF','NF','pRV','pRV','RVF','NF','pRV','NF')
#M1.aggreg$group <- disease
#M1.aggreg <- SetIdent(M1.aggreg,value="group")
#markers.RVF <- FindMarkers(M1.aggreg,ident.1='RVF',ident.2=c('pRV'),test.use = "DESeq2")
#RVF_vs_pRVF_25_up <- rownames(subset(markers.RVF,avg_log2FC > 0))[1:25]
#RVF_vs_pRVF_25_down <- rownames(subset(markers.RVF,avg_log2FC < 0))[1:25]
#to_plot <- c(RVF_vs_pRVF_25_up,RVF_vs_pRVF_25_down)





M1.aggreg <- AggregateExpression(object = M1, group.by = c('niche_manual'),return.seurat=T,assay='Xenium')
M1.aggreg <- FindVariableFeatures(M1.aggreg,selection.method = "vst")
to_plot <- VariableFeatures(M1.aggreg)[1:50]
M1.aggreg <- ScaleData(M1.aggreg)
data <- GetAssayData(M1.aggreg,layer='scale.data')

pdf('~/Downloads/hdWGCNA_TOM/Xenium/xenium_heatmap_niche.pdf',width=4,height=10)

heatmap.2(as.matrix(data[to_plot,rev(names(sort(niche.counts)))]), scale="row",
   labRow=to_plot, 
   col=mycol, margin=c(6,6),trace="none", density.info="none", lhei=c(1,10,3),lwid=c(1,10), dendrogram='none',breaks = seq(-4, 4, length.out = 1001),
   Rowv = FALSE, Colv=FALSE,srtCol=90,lmat = rbind(c(0,3),c(2,1),c(0,4)))

dev.off()










# #' Construct an assay for spatial niche analysis
# #'
# #' This function will construct a new assay where each feature is a
# #' cell label The values represents the sum of a particular cell label
# #' neighboring a given cell.
# #'
# #' @param list.object list of Seurat objects to do clustering on 
# #' @param list.fov list of fov names to use for grabbing cells to cluster from list.object.  Should be the same length as list.object
# #' @param group.by Cell classifications to count in spatial neighborhood
# #' @param assay Name for spatial neighborhoods assay
# #' @param cluster.name Name of output clusters
# #' @param neighbors.k Number of neighbors to consider for each cell
# #' @param niches.k.range Number of clusters to return based on the niche assay.  provide a range
# #' @param batch_size Number of mini-batches for ClusterR::MiniBatchKmeans
# #' @param num_init  # number of times the algorithm will be run with different centroid seeds for ClusterR::MiniBatchKmeans
# #'
# #' @importFrom stats kmeans
# #' @return Seurat object containing a new assay
# #' @concept clustering
# #' @export
# #'
# # BuildNicheAssay.multiple_FOVs.MiniBatchKmeans <- function(
# #         list.object,
# #         list.fov,
# #         group.by,
# #         assay = "niche",
# #         cluster.name = "niches",
# #         neighbors.k = 20,
# #         niches.k.range = 2:30 , 
# #         batch_size = 20, 
# #         num_init = 20
# # ) {
# #     # check for fov in sample set
# #     # remove if not found in object
# #     remove = NULL # init list of indices to remove
# #     for ( i in seq_along(list.object) ){ # message(i)
# #         # get object and fov for each object
# #         object = list.object[[i]]
# #         fov = list.fov[[i]]
        
# #         if( !(fov %in%  names(object@images)) ){
# #             warning( "fov is not found in the i-th object.  Removing the object from the list.object and list.fov.  i =", i)
# #            remove = c(remove, i)
# #         }
# #     }
# #     for (i in rev(remove) ){
# #         list.object[[i]] = NULL
# #         list.fov[[i]] = NULL
# #     }
    
    
# #     for ( i in seq_along(list.object) ){ message(i)
# #         # get object and fov for each object
# #         object = list.object[[i]]
# #         fov = list.fov[[i]]
        
# #         # initialize an empty cells x groups binary matrix
# #         cells <- Cells( object[[fov]] )
# #         group.labels <- unlist(object[[group.by]][cells, ] )
# #         groups <- sort( unique(group.labels) )
# #         cell.type.mtx <- matrix(
# #             "data" = 0
# #             , "nrow" = length(cells)
# #             , "ncol" = length(groups)
# #         )
# #         rownames(cell.type.mtx) <- cells
# #         colnames(cell.type.mtx) <- groups

# #         # populate the binary matrix 
# #         cells.idx <- seq_along(cells)
# #         group.idx <- match(group.labels, groups)
# #         cell.type.mtx[cbind(cells.idx, group.idx)] <- 1
        
# #         # find neighbors based on tissue position
# #         coords <- Seurat::GetTissueCoordinates( object[[fov]], "which" = "centroids" )
# #         rownames(coords) <- coords[["cell"]]
# #         coords <- as.matrix(coords[ , c("x", "y")])
# #         neighbors <- Seurat::FindNeighbors(
# #             "object" = coords
# #             , "k.param" = neighbors.k # Defines k for the k-nearest neighbor algorithm
# #             , "compute.SNN" = F
# #         )
        
# #         # create niche assay
# #         sum.mtx <- as.matrix( neighbors[["nn"]] %*% cell.type.mtx )
# #         niche.assay <- CreateAssayObject( "counts" = t(sum.mtx) )
# #         object[[assay]] <- niche.assay
# #         DefaultAssay(object) <- assay
        
# #         # scale data 
# #         object <- ScaleData(object)
        
# #         # return edited object to list
# #         list.object[[i]] = object
        
        
# #     }
    
# #     # get aggregate data for ClusterR::MiniBatchKmeans
# #     # columns = features
# #     # rows = cells
# #     # cells = values
# #     DAT = lapply( seq_along(list.object), function(i){
# #         t( list.object[[i]][[assay]]@scale.data )
# #     }  )
# #     DAT <- do.call("rbind", DAT)
    
    
    
# #     res.clusters = data.frame(row.names = rownames(DAT))
    
# #     for ( k in niches.k.range ){ message("k=", k)
# #         # new column name
# #         newCol = paste0("kmeans_", k)
# #         # get centroids
# #         km_mb = ClusterR::MiniBatchKmeans(
# #             "data" = DAT
# #             , "clusters" = k # the number of clusters
# #             , "batch_size" = batch_size # the size of the mini batches
# #             , "num_init" = num_init # number of times the algorithm will be run with different centroid seeds
# #             , "max_iters" = 100 # the maximum number of clustering iterations. 
# #             , "init_fraction" = 0.2 # percentage of data to use for the initialization centroids (applies if initializer is kmeans++ or optimal_init). Should be a float number between 0.0 and 1.0.
# #             , "initializer" = "kmeans++" # the method of initialization. One of, optimal_init, quantile_init, kmeans++ and random. See details for more information
# #             , "early_stop_iter" = 10 # continue that many iterations after calculation of the best within-cluster-sum-of-squared-error
# #             , "verbose" = F
# #             , "CENTROIDS" = NULL
# #             , "tol" = 1e-04
# #             , "tol_optimal_init" = 0.3
# #             , "seed" = 1
# #         )
        
# #         # use centroids to get clusters
  
# #         res.clusters[,newCol] = ClusterR::predict_MBatchKMeans( # This function takes the data and the output centroids and returns the clusters.
# #             "data" = DAT
# #             , "CENTROIDS" = km_mb$centroids
# #         )
# #         res.clusters[,newCol] = as.factor( res.clusters[,newCol] ) # change clusters to factors
        
# #     }
    
# #     # get clusters back onto the objects
# #     colnames(res.clusters) = paste0(cluster.name,".", colnames(res.clusters))
# #     for ( i in seq_along(list.object) ){ message(i)
# #         # get object and fov for each object
# #         object = list.object[[i]]
        
# #         # get clusters in correct cell row order into metadata of object
# #         object[[]] = res.clusters[rownames(object[[]]),]
        
# #         # return edited object to list
# #         list.object[[i]] = object
# #     }

    
# #     return(list.object)
# # }

# # library(magrittr)
# # library(dplyr)

# # M1_list <- SplitObject(M1, split.by = "patient")
# # fov_list <- c('fov','fov.2','fov.3','fov.4','fov.5','fov.6','fov.7','fov.8','fov.9')
# # NicheAssay <- BuildNicheAssay.multiple_FOVs.MiniBatchKmeans(list.object = M1_list,
# # 	list.fov = fov_list, group.by = 'markernames', assay = 'niche',
# # 	cluster.name = 'niches', neighbors.k = 50, niches.k.range = 2:30,
# # 	batch_size = 20, num_init = 20)


# # for ( i in seq_along(M1_list) ){ message(i)
# #     # get object and fov for each object
# #     object = M1_list[[i]]
# #     fov = fov_list[[i]]
    
# #     # initialize an empty cells x groups binary matrix
# #     cells <- Cells( object[[fov]] )
# #     group.labels <- unlist(object[['markernames']][cells, ] )
# #     groups <- sort( unique(group.labels) )
# #     cell.type.mtx <- matrix(
# #         "data" = 0
# #         , "nrow" = length(cells)
# #         , "ncol" = length(groups)
# #     )
# #     rownames(cell.type.mtx) <- cells
# #     colnames(cell.type.mtx) <- groups

# #     # populate the binary matrix 
# #     cells.idx <- seq_along(cells)
# #     group.idx <- match(group.labels, groups)
# #     cell.type.mtx[cbind(cells.idx, group.idx)] <- 1
    
# #     # find neighbors based on tissue position
# #     coords <- Seurat::GetTissueCoordinates( object[[fov]], "which" = "centroids" )
# #     rownames(coords) <- coords[["cell"]]
# #     coords <- as.matrix(coords[ , c("x", "y")])
# #     neighbors <- Seurat::FindNeighbors(
# #         "object" = coords
# #         , "k.param" = 50 # Defines k for the k-nearest neighbor algorithm
# #         , "compute.SNN" = F
# #     )
    
# #     # create niche assay
# #     sum.mtx <- as.matrix( neighbors[["nn"]] %*% cell.type.mtx )
# #     niche.assay <- CreateAssayObject( "counts" = t(sum.mtx) )
# #     object[['niche']] <- niche.assay
# #     DefaultAssay(object) <- 'niche'
    
# #     # scale data 
# #     object <- ScaleData(object)
    
# #     # return edited object to list
# #     M1_list[[i]] = object
        
# # }










# # get aggregate data for ClusterR::MiniBatchKmeans
# # columns = features
# # rows = cells
# # cells = values
# niches.k.range = 10:12

# res.clusters = data.frame(row.names = rownames(DAT.scale.data))

# for ( k in niches.k.range ){ message("k=", k)
#     # new column name
#     newCol = paste0("kmeans_", k)
#     # get centroids
#     km_mb = ClusterR::MiniBatchKmeans(
#         "data" = DAT.scale.data
#         , "clusters" = k # the number of clusters
#         , "batch_size" = 20 # the size of the mini batches
#         , "num_init" = 20 # number of times the algorithm will be run with different centroid seeds
#         , "max_iters" = 100 # the maximum number of clustering iterations. 
#         , "init_fraction" = 0.2 # percentage of data to use for the initialization centroids (applies if initializer is kmeans++ or optimal_init). Should be a float number between 0.0 and 1.0.
#         , "initializer" = "kmeans++" # the method of initialization. One of, optimal_init, quantile_init, kmeans++ and random. See details for more information
#         , "early_stop_iter" = 10 # continue that many iterations after calculation of the best within-cluster-sum-of-squared-error
#         , "verbose" = F
#         , "CENTROIDS" = NULL
#         , "tol" = 1e-04
#         , "tol_optimal_init" = 0.3
#         , "seed" = 1
#     )
    
#     # use centroids to get clusters

#     res.clusters[,newCol] = ClusterR::predict_MBatchKMeans( # This function takes the data and the output centroids and returns the clusters.
#         "data" = DAT.scale.data
#         , "CENTROIDS" = km_mb$centroids
#     )
#     res.clusters[,newCol] = as.factor( res.clusters[,newCol] ) # change clusters to factors
    
# }


# # Update object
# res.clusters.ordered <- res.clusters[colnames(M1),]
# for(i in colnames(res.clusters.ordered))
# 	{
# 		M1[[i]] = 0
# 		M1[[i]][rownames(res.clusters.ordered),] = res.clusters.ordered[,i]

# 	}

# write.table(res.clusters.ordered,'~/Downloads/hdWGCNA_TOM/Xenium/Niche_clusters.csv',sep=',')

# DimPlot(M1,group.by='kmeans_2',label=T)


# ImageDimPlot(M1, fov = "fov.3", group.by = "kmeans_2", size = .5, dark.background = F)

# ImageDimPlot(M1, fov = "fov.2", group.by = "kmeans_30", size = .5, dark.background = F)






# M1 <- BuildNicheAssay(object = M1, fov = "fov", group.by = "markernames",niches.k = 10, neighbors.k = 200)

# ImageDimPlot(M1, fov = "fov.2", group.by = "kmeans_30", size = .5, dark.background = F)










# ImageDimPlot(M1, fov = "fov.2", group.by = "niches", size = .5, dark.background = F)

# DimPlot(subset(M1,patient='1697'),group.by='niches',label=T)

# table(M1$names, M1$niches)




# DimPlot(xenium.obj,group.by='markernames',label=T) + NoLegend()

# xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov.2", group.by = "markernames",niches.k = 10, neighbors.k = 30)


# celltype.plot <- ImageDimPlot(xenium.obj, fov = "fov.2", group.by = "markernames", size = 1.5, cols = "polychrome",
#     dark.background = F) + ggtitle("Cell type") + NoLegend()
# niche.plot <- ImageDimPlot(xenium.obj, fov = "fov.2", group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches") +
#     scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))
# celltype.plot | niche.plot

# niche_comp <- table(xenium.obj$markernames, xenium.obj$niches)

# t(t(table(xenium.obj$names, xenium.obj$niches)) / colSums(table(xenium.obj$names, xenium.obj$niches)))


# ImageDimPlot(xenium.obj, fov = "fov.2",split.by='subnames',cols = "polychrome", axes = TRUE)


# head(subset(m0,p_val_adj < 0.05 & cluster == 28))



# xenium.1697 <- FindSpatiallyVariableFeatures(xenium.1697, fov = 'fov',assay = "SCT", features = rownames(xenium.1697),
#     selection.method = "moransi")


# top.features = rownames(
#   dplyr::slice_min(
#     xenium.1697[["SCT"]]@meta.features,
#     moransi.spatially.variable.rank,
#     n = 20
#   )
# )


# ImageFeaturePlot(xenium.1697, features = top.features)


# cm <- SetIdent(cm,value='group')
# a <- FindMarkers(cm,ident.1='RVF',ident.2='pRV')
# library(EnhancedVolcano)

# EnhancedVolcano(a,lab = rownames(a), x = 'avg_log2FC',y = 'p_val_adj',pCutoff = 0.05,FCcutoff=0.25)

# #ImageDimPlot(xenium.obj, fov = "fov", molecules = c("TREM2", "LYVE1", "CD4"), nmols = 20000)


# ImageFeaturePlot(xenium.obj, fov = "fov.2", features = c("NPPA"), max.cutoff = c(3), size = 0.75, cols = c("white", "red"))
# ImageFeaturePlot(xenium.obj, fov = "fov.2", features = c("SMYD2"), max.cutoff = c(3), size = 0.75, cols = c("white", "red"))
# ImageFeaturePlot(xenium.obj, fov = "fov.2", features = c("IDH3A"), max.cutoff = c(3), size = 0.75, cols = c("white", "red"))





# FeaturePlot(xenium.obj, features = c("PECAM1","LYVE1","DCN","TNNT2"))

# m0 <- FindAllMarkers(xenium.obj,recorrect_umi=FALSE) #Ideally would recorrect UMI here


# #0, 3, 7, 9, 10, 13, 19, 20,21, 22, 24, 25
# #20, 21, 22, 24, 25 are outliers
# cm <- subset(xenium.obj,idents=c('0','3','7','9','10','13','19'))

# cm <- RunPCA(cm, npcs = 30, features = rownames(cm))
# cm <- RunHarmony(cm,'dataset')
# cm <- RunUMAP(cm, dims = 1:30,reduction = "harmony")
# cm <- FindNeighbors(cm, reduction = "harmony", dims = 1:30)
# cm <- FindClusters(cm, resolution = 0.3)

# DimPlot(cm,label=T)

# m1 <- FindAllMarkers(cm,recorrect_umi=FALSE)

# head(subset(m1,avg_log2FC > 0 & p_val_adj < 0.05 & cluster == 1)$gene,n=5)
# head(subset(m1,avg_log2FC < 0 & p_val_adj < 0.05 & cluster == 1)$gene,n=5)
# head(subset(m0,p_val_adj < 0.05 & cluster == 0))

# #0 Low - DST MYO18B TNNT1 VEGFA SORBS2
# #1 DST MYO18B VEGFA TNNT1 SORBS2 FHL2
# #2 MYH6 GPX3 PALLD; Low - TNNT1, CYP2J2 SORBS2
# #3 Low - FHL2 SORBS2 MYBPC3 CDH2 VEGFA
# #4 NPPA NPPB
# #5 SORBS2 CYP2J2 TNNT1
# #6 SMYD2
# #7 IDH3A
# #8 MYPN
# #9 TMEM65
# #10 Low MYH2, ACTC1, MYH7, CDH2, S100A1, ANKRD1
# #11 MYO18B
# #12 PDE3A
# #13 CYP2J2
# #14 SORBS2
# #15 FB?
# #16 TNNI3
# #17 TIMP3

# #NPPA/B ACTA1 TMEM65 ACTC1, SORBS2, SMYD2, IDH3A, VEGFB, RGCC, TIMP3, MYPN,CNBP, DST





# #6 is EC


# clust_freq <- 100 * t(table(xenium.obj@active.ident,xenium.obj$group)) / colSums(table(xenium.obj@active.ident,xenium.obj$group))

# clust_freq <- 100 * t(table(xenium.obj@active.ident,xenium.obj$patient)) / colSums(table(xenium.obj@active.ident,xenium.obj$patient))

# df <- clust_freq[,3]
# df <- as.data.frame(df)
# df$disease <- c('RVF','RVF','NF','pRV','pRV','RVF','NF','pRV','NF')


# ggplot(df,aes(disease,df))+geom_boxplot()

# xenium.1697 <- FindSpatiallyVariableFeatures(xenium.1697, assay = "SCT", features = rownames(xenium.1697),
#     selection.method = "moransi")
# top.features <- head(SpatiallyVariableFeatures(xenium.1697, selection.method = "moransi"), 6)
# SpatialFeaturePlot(xenium.1697, features = top.features, ncol = 3, alpha = c(0.1, 1))



# #Unique are 7,8,0,2,12

# FeaturePlot(cm, features = c("NPPB"))


# #TMEM65 #ATC1 #MYH7 cluster

# cm <- subset(xenium.obj,idents=c('0','3','6','7','8','9','11','12','14','15','16'))


# myeloid <- subset(xenium.obj,idents=c('5'))

# myeloid <- RunPCA(myeloid, npcs = 30, features = rownames(myeloid))
# myeloid <- RunUMAP(myeloid, dims = 1:30)
# myeloid <- FindNeighbors(myeloid, reduction = "pca", dims = 1:30)
# myeloid <- FindClusters(myeloid, resolution = 0.3)

# DimPlot(myeloid,label=T)

# m1 <- FindAllMarkers(myeloid)
# head(subset(m1,avg_log2FC > 0 & p_val_adj < 0.05 & cluster == 3)$gene,n=5)

# FeaturePlot(myeloid, features = c("LYVE1"))


# m3 <- FindAllMarkers(xenium.obj)
# head(subset(m3,avg_log2FC > 0 & p_val_adj < 0.05 & cluster == 3)$gene,n=5) 

# #10 MYH11 - SM
# #1 - EC CD34, RGCC, EFNB2 (cap, art)
# # ART markers - EFNB2, HEY1, FBLN5, maybe MECOM FeaturePlot(xenium.obj, features = c("MECOM","HEY1","FBLN5","EFNB2"), label = T)
# # Venous - NR2F2, EPHB4
# # Cap - RGCC, MFSD2A
# #4 - venous EC NR2F2
# #17 - lymphatic EC PROX1
# #19 Adipo
# #2,13 FB DCN
# #18 - NKT

# ec <- subset(xenium.obj,idents=c('1','4'))

# ec <- RunPCA(ec, npcs = 30, features = rownames(ec))
# ec <- RunUMAP(ec, dims = 1:30)
# ec <- FindNeighbors(ec, reduction = "pca", dims = 1:30)
# ec <- FindClusters(ec, resolution = 0.3)

# DimPlot(ec,label=T)

# m2 <- FindAllMarkers(ec)
# head(subset(m2,avg_log2FC > 0 & p_val_adj < 0.05 & cluster == 3)$gene,n=5) 

# FeaturePlot(ec, features = c("NR2F2"))



