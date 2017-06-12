#' run tSNE from (R pkg 'Rtsne') on a gatingSet
#' Will sample the minimal number of cells available in all samples to generate balanced cell counts
#' 
#' IMPORTANT: Requires a valid gatingSet with cytokine gates downstream of a parent gate
#' Also expects that pData(gs) contains at least columns: 'name', 'ptid' so we can identify cells later
#' 
#' @param gs a GatingSet object, properly gated data with annotation in its pData
#' @param parentGate a \code{string} describing the gate upstream of the cytokine gates (eg. "CD4", "cd8+", etc...)
#' @param parentGateMarker a \code{string} describing the marker which corresponds to parentGate
#' @param degreeFilterMarkers a \code{vector} of \code{strings} describing the marker gates immediately downstream of parentGate, eg: "IL2", "IFNg"
#' @param otherMarkers the remaining markers of the data
#' @param markerMap named list of degreeFilterMarkers marker names to gate names, eg. list("CD4/IL2" = "IL2","CD4/IFNg" = "IFNg")
#' @param notRunMarkers a \code{vector} of \code{strings} describing markers which are NOT used for the run itself but are imported as columns in the final tsne output matrix
#' @param notRunMarkersMap named list of notRunMarkers marker names to gate names, eg. list("CD3/CD45RA" = "CD45RA","CD3/CD4" = "CD4"). If provided, a "poly_NotRunMarkers" column will be created.
#' @param swap boolean for whether marker and gate names (from markerMap above) should be swapped. Passed onto getSingleCellExpression()
#' @param groupBy columns of the \code{gatingSet}'s phenoData, same number of cells will be sampled from each group
#' @param degreeFilter keep cells of this degree and higher, useful when tSNE takes too long to run
#' @param seed since tSNE is random, need a random seed so we can reproduce results
#' @param theta parameter to be passed to the \code{Rtsne} function
#' @param numThreads if > 1, uses the multicore t-SNE instead of regular t-SNE
#' @param ... other parameters to be passed to the \code{Rtsne} function
#' @return a \code{matrix} of X and Y coordinates
#' @import flowWorkspace
#' @import data.table
#' @import plyr
#' @import Rtsne
#' @import Rtsne.multicore
#' TODO: the `degreeFilterMarkers` parameter is not actually required, just the markerMap.
#' Also, this version of the function allows degreeFilter == 0, but requires the parentGate to be a marginal 1D gate which maps to a specific marker
#' (due to the way getSingleCellExpression works)
runTSNE <- function (gs, parentGate, parentGateMarker, degreeFilterMarkers, otherMarkers, markerMap, notRunMarkers = c(), notRunMarkersMap = list(), swap = FALSE,
                     groupBy, degreeFilter = 0, seed = 999, theta = 0.9, numThreads = 1, ...) {
  
  if (is.null(markerMap)) stop ("required markerMap is missing ! STOPPING....")
  
  set.seed(seed)
  pd <- as.data.table(pData(gs))
  meta_cols <- colnames(pd)
  meta_cols <- c(meta_cols, "degree", "poly", "poly_notRunMarkers", notRunMarkers)
  cat("getting total cell counts from parent gate", parentGate, 
      "\n")
  parent_count <- unlist(lapply(gs, function(gh) getTotal(gh, 
                                                          parentGate)))
  parent_count = ldply(parent_count)
  setnames(parent_count, c("name", parentGate))
  pd <- merge(pd, parent_count, by = "name")
  nTcells <- min(pd[, sum(get(parentGate)), by = groupBy][, 
                                                          V1])
  cat("after grouping by '", groupBy, "', all groups will have at least", 
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
      gh <- gs[[sn]]
      updateIndices(gh, parentGate, thisInd)
    }
  }, by = groupBy] 
  cat("subsampling complete ! recomputing... \n")
  nodes <- getChildren(gs[[1]], parentGate, path = 2)
  for (node in nodes) recompute(gs, node)
  
  cat("generating event masks \n")

  # Obtain the list of nodes of interest from the markerMap.
  userNodes <- names(markerMap)
  # Make sure they all exist as children of the parent node.
  childNodeOptions <- c(nodes, getChildren(gs[[1]], parentGate, path = 1), getChildren(gs[[1]], parentGate, path = "full"), getChildren(gs[[1]], parentGate, path = "auto"))
  for(userNode in userNodes) { stopifnot(userNode %in%  childNodeOptions) }
  # For now we want to obtain the single cell fluorescence data for ALL cells that are members of the parentGate.
  # We pass the parentGate argument to the nodes argument to ensure this happens.
  # TODO: Figure out a method to get all cells w/o requiring parentGate to be a marginal 1D node with an associated marker
  markerMap[[parentGate]] <- parentGateMarker
  res <- getSingleCellExpression(gs, nodes = c(userNodes, parentGate), other.markers = c(otherMarkers, notRunMarkers), 
                                 map = markerMap, threshold = FALSE, swap=swap)
  # Then remove the parentGate column from res
  for(i in seq_along(res)) { res[[i]] <- res[[i]][,!colnames(res[[i]]) %in% parentGate] }
  
  cat("\n There are", sum(unlist(lapply(res, function(x) { nrow(x) } ))), "rows after subsampling cells")

  # Same as res, but all data below gating thresholds are set to 0. Used to create "degree", "poly", and "poly_notRunMarkers" columns.
  # Therefore "degree" and "poly" columns are based solely on the degreeFilterMarkers listed under nodes and degreeFilterMarkers.
  res_mask <- getSingleCellExpression(gs, nodes = unique(if(length(notRunMarkersMap)) { c(userNodes, parentGate, notRunMarkers[notRunMarkers %in% notRunMarkersMap]) } else { c(userNodes, parentGate) }),
                                      map = c(markerMap, notRunMarkersMap), 
                                      threshold = T, swap=swap)
  # Remove the parentGate column from res_mask
  for(i in seq_along(res_mask)) { res_mask[[i]] <- res_mask[[i]][,!colnames(res_mask[[i]]) %in% parentGate] }
  
  # Concatenate all the res matrices together into one big tsne-friendly matrix, after
  # using the res_mask list of matrices to create additional columns "poly" and "degree"
  res_collapse <- ldply(names(res), function(sn) { # each element of res corresponds to a sample
    message(".", appendLF = F)
    # message(sn)
    mat <- as.data.frame(res[[sn]])
    if (nrow(mat) > 0) {
      nColsUserNodes <- length(userNodes[!userNodes %in% c(parentGate)])
      # Create the degree and poly columns for the run markers:
      mat_mask <- res_mask[[sn]][,1:nColsUserNodes]
      mat_mask[mat_mask > 0] <- 1
      mat[, "degree"] <- rowSums(mat_mask)
      mat[, "poly"] <- Reduce(paste0, as.list(as.data.frame(mat_mask)))
      pd <- pData(gs[[sn]])
      rownames(pd) <- NULL
      if(length(notRunMarkersMap)) {
        # Create the poly_notRunMarkers column for the notRunMarkers:
        mat_mask_notRun <- res_mask[[sn]][,(nColsUserNodes+1):ncol(res_mask[[sn]])]
        mat_mask_notRun[mat_mask_notRun > 0] <- 1
        mat[, "poly_notRunMarkers"] <- Reduce(paste0, as.list(as.data.frame(mat_mask_notRun))) 
      }
      # Return the combined matrices
      cbind(mat, pd)
    } 
    else NULL
  })
  
  cat("\n degreeFilterMarkers:", degreeFilterMarkers)
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