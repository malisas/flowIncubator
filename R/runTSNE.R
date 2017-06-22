#' run tSNE from (R pkg 'Rtsne') on a gatingSet
#' Will sample the minimal number of cells available in all samples to generate balanced cell counts
#' 
#' IMPORTANT: Requires a valid GatingSet with cytokine gates downstream of a parent gate
#' Also expects that pData(gs) contains at least columns: 'name', 'ptid' so we can identify cells later
#' 
#' @param gs a GatingSet object, properly gated data with annotation in its pData
#' @param parentGate a \code{string} describing the gate upstream of the cytokine gates (eg. "CD4", "cd8+", etc...)
#' @param degreeFilterGates a \code{vector} of \code{strings} describing the marker gates immediately downstream of parentGate, eg: "CD4/IL2", "CD4/IFNg". Each entry is assumed to corresspond to one marker, which is used to degree filter the data and input into tsne.
#' @param otherMarkers a \code{vector} of \code{strings} describing other markers (aside from those corresponding to degreeFilterNodes) which are input into tsne but not used to degree filter
#' @param notRunMarkers a \code{vector} of \code{strings} describing markers which are NOT used for the run itself but are imported as columns in the final tsne output matrix
#' @param gateMarkerMap named list of gate names to marker names, eg. list("CD4/IL2" = "IL2","CD4/IFNg" = "IFNg"). Must contain entries for parentGate and degreeFilterGates. If entries are provided for otherMarkers and notRunMarkers, their corresponding xx_Bool columns will appear in the final tsne output matrix.
#' @param swap boolean for whether marker and gate names (from markerMap above) should be swapped. Passed onto getSingleCellExpression()
#' @param groupBy a \code{vector} of \code{strings} describing columns of the \code{gatingSet}'s phenoData, same number of cells will be sampled from each group
#' @param degreeFilter keep cells of this degree and higher, useful when tSNE takes too long to run
#' @param seed since tSNE is random, need a random seed so we can reproduce results
#' @param theta parameter to be passed to the \code{Rtsne} function
#' @param numThreads if > 1, uses multicore t-SNE instead of regular t-SNE
#' @param ... other parameters to be passed to the \code{Rtsne} function
#' @return a \code{matrix} of X and Y coordinates
#' @import flowWorkspace
#' @import data.table
#' @import plyr
#' @import magrittr
#' @import Rtsne
#' @import Rtsne.multicore
#' TODO: maybe allow clone = FALSE option to save memory...?
runTSNE <- function (gs=NULL, parentGate = NULL, degreeFilterGates = c(), otherMarkers = c(), notRunMarkers = c(), gateMarkerMap = NULL, swap = FALSE,
                     groupBy = c(), degreeFilter = 0, seed = 999, theta = 0.9, numThreads = 1, ...) {
  
  if (is.null(gs)) stop ("required gs is missing ! STOPPING....")
  if (is.null(parentGate)) stop ("required parentGate is missing ! STOPPING....")
  if (is.null(gateMarkerMap)) stop ("required gateMarkerMap is missing ! STOPPING....")
  
  # The GatingSet gets modified below, so we clone it in order to avoid modifying the original object
  gsClone <- clone(gs)
  set.seed(seed)
  meta_cols <- colnames(pd)
  meta_cols <- c(meta_cols, "degree", notRunMarkers,
                 unlist(lapply(unlist(gateMarkerMap, use.names = FALSE), function(x) { paste0(x, "_Bool") })))

  for(group in groupBy) {
    cat("getting total cell counts from parent gate", parentGate, 
        "\n")
    parent_count <- unlist(lapply(gsClone, function(gh) getTotal(gh, 
                                                            parentGate)))
    parent_count = ldply(parent_count)
    setnames(parent_count, c("name", parentGate))
    pd <- as.data.table(pData(gsClone))
    pd <- merge(pd, parent_count, by = "name")
    nTcells <- min(pd[, sum(get(parentGate)), by = group][, 
                                                            V1])
    cat("after grouping by '", group, "', all groups will have at least", 
        nTcells, "cells.\n")
    pd[, {
      # Sample from the total number of events
      # in the parent population. 
      totalEvents <- sum(get(parentGate))
      gInd <- 1:totalEvents
      # sample, without replacement, a vector of cell indices for the entire group, length totalEvents
      gInd <- sample.int(totalEvents, size = nTcells)
      gInd.logical <- rep(F, totalEvents)
      gInd.logical[gInd] <- T
      sn.factor <- unlist(sapply(name, function(sn) rep(sn, 
                                                        .SD[name == sn, get(parentGate)])))
      ind.vec <- split(gInd.logical, sn.factor)
      for (sn in name) {
        thisInd <- ind.vec[[sn]]
        gh <- gsClone[[sn]]
        updateIndices(gh, parentGate, thisInd)
      }
    }, by = group] 
    cat("subsampling complete ! recomputing... \n")
    nodes <- getChildren(gsClone[[1]], parentGate, path = 2)
    for (node in nodes) recompute(gsClone, node)
  }

  cat("generating event masks \n")

  # Make sure degreeFilterGates exist as children of the parent node. This can probably be checked better.
  childNodeOptions <- c(getChildren(gsClone[[1]], parentGate, path = 2), getChildren(gsClone[[1]], parentGate, path = 1),
                        getChildren(gsClone[[1]], parentGate, path = "full"), getChildren(gsClone[[1]], parentGate, path = "auto"))
  for(node in c(degreeFilterGates)) { stopifnot(node %in% childNodeOptions & node %in% names(gateMarkerMap)) }
  # Also make sure parentGate is in gateMarkerMap
  stopifnot(parentGate %in% names(gateMarkerMap))
  
  # First we want to obtain the single cell fluorescence data for ALL cells that are members of the parentGate.
  # Pass the parentGate argument to the nodes argument to ensure this happens (all members of parentGate are retrieved)
  res <- getSingleCellExpression(gsClone, nodes = c(parentGate, degreeFilterGates), other.markers = c(otherMarkers, notRunMarkers), 
                                 map = gateMarkerMap, threshold = FALSE, swap=swap)
  # Then remove the parentGate column from res
  for(i in seq_along(res)) { res[[i]] <- res[[i]][,!colnames(res[[i]]) %in% gateMarkerMap[[parentGate]]] }
  
  cat("\n There are", sum(unlist(lapply(res, function(x) { nrow(x) } ))), "rows after subsampling cells")

  # Same as res, but all data below gating thresholds are set to 0. Used to create "degree" and "xx_Bool" columns.
  # Therefore "degree" and "xx_Bool" columns are based solely on the degreeFilterGates and associated markers provided in gateMarkerMap.
  # First we obtain the "boolean" information for all the nodes/gates listed in gateMarkerMap
  res_mask <- getSingleCellExpression(gsClone, nodes = names(gateMarkerMap),
                                      map = gateMarkerMap, 
                                      threshold = T, swap=swap)
  # Remove the parentGate column from res_mask
  for(i in seq_along(res_mask)) { res_mask[[i]] <- res_mask[[i]][,!colnames(res_mask[[i]]) %in% gateMarkerMap[[parentGate]]] }
  
  # Find the column indices of res_mask which corresspond to degreeFilterGates markers:
  degreeFilterCols <- which(colnames(res_mask[[1]]) %in% gateMarkerMap[degreeFilterGates])
  # And their corresponding names
  degreeFilterColNamesOrdered <- colnames(res_mask[[1]])[degreeFilterCols]
  
  # Concatenate all the res matrices together into one big tsne-friendly matrix, after
  # using the res_mask list of matrices to create additional "degree" and "xx_Bool" columns
  res_collapse <- ldply(names(res), function(sn) { # each element of res corresponds to a sample
    message(".", appendLF = F)
    # message(sn)
    mat <- as.data.frame(res[[sn]])
    if (nrow(mat) > 0) {
      mat_mask <- res_mask[[sn]]
      mat_mask[mat_mask > 0] <- 1
      
      # Create the degree column for the run markers:
      if(length(degreeFilterCols)) {
        mat_mask_degreeFilterCols <- mat_mask[,degreeFilterCols]
        mat[, "degree"] <- if(class(mat_mask_degreeFilterCols) == "numeric") {
          mat_mask_degreeFilterCols
        } else {
          rowSums(mat_mask_degreeFilterCols)
        }
      } else {
        mat[, "degree"] <- rep(0, nrow(mat))
      }
      
      # Create the "xx_Bool" columns for all markers in gateMarkerMap
      colnames(mat_mask) <- unlist(lapply(colnames(res_mask[[1]]), function(x) { paste0(x, "_Bool") }))
      
      pd <- pData(gsClone[[sn]])
      rownames(pd) <- NULL
      # Return the combined matrices
      cbind(mat, mat_mask, pd)
    } 
    else NULL
  })

  cat("\n degreeFilter Markers:", unlist(gateMarkerMap[degreeFilterGates], use.names = FALSE))
  cat("\n other markers:", as.character(otherMarkers))
  cat("\n input matrix has", nrow(res_collapse), "rows...")
  res_collapse <- subset(res_collapse, degree >= degreeFilter)
  cat("\n input matrix has", nrow(res_collapse), "rows after filtering for cells of degree >=", 
      degreeFilter)
  input_mat <- as.matrix(res_collapse[, !names(res_collapse) %in% 
                                        meta_cols])
  cat("\n starting tSNE run at ", date(), " with ", numThreads, " threads\n")
  system.time(tsne_out <- if(numThreads > 1) {
    Rtsne.multicore::Rtsne.multicore(X = input_mat, check_duplicates = FALSE, num_threads = numThreads, 
                                  ...)
  } else {
    Rtsne(input_mat, check_duplicates = FALSE, 
                    ...)
  })
  
  dat <- tsne_out$Y
  colnames(dat) <- c("x", "y")
  dat <- cbind(dat, res_collapse)
  dat <- data.table(dat)
  cat("completed tSNE run at", date(), "!\n")
  return(dat)
}