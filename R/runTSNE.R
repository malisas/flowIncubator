createTsneInputMatrix <- function(gs=NULL, parentGate = NULL, degreeFilterGates = c(), otherMarkers = c(), notRunMarkers = c(), gateMarkerMap = NULL, swap = FALSE,
                                  groupBy = c(), degreeFilter = 0, seed = 999, theta = 0.9, cloneGs = TRUE) {
  if (is.null(gs)) stop ("required gs is missing ! STOPPING....")
  if (is.null(parentGate)) stop ("required parentGate is missing ! STOPPING....")
  if (is.null(gateMarkerMap)) stop ("required gateMarkerMap is missing ! STOPPING....")
  if (length(groupBy) > 2) stop ("groupBy length can be at most 2")
  
  # The GatingSet gets modified below, so we clone it in order to avoid modifying the original object
  gsClone <- if(class(gs) == "character") {
    load_gs(gs)
  } else {
    if (cloneGs) {
      clone(gs)
    } else {
      gs
      }
  } 
  if (!all(groupBy %in% colnames(pData(gsClone)))) stop("all groupBy values must be columns of gs metadata")
  
  set.seed(seed)
  pd <- as.data.table(pData(gsClone))
  meta_cols <- colnames(pd)
  meta_cols <- c(meta_cols, "degree", notRunMarkers,
                 unlist(lapply(unlist(gateMarkerMap, use.names = FALSE), function(x) { paste0(x, "_Bool") })))
  
  # If desired, sample the cells
  cat("getting total cell counts from parent gate", parentGate, 
      "\n")
  parent_count <- unlist(lapply(gsClone, function(gh) getTotal(gh, 
                                                               parentGate)))
  parent_count = ldply(parent_count)
  setnames(parent_count, c("name", parentGate))
  pd <- merge(pd, parent_count, by = "name")
  
  if (length(groupBy) == 1) {
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
  } else if (length(groupBy) == 2) {
    # Calculate the "minimum maximum group size" which lets the size of all groupBy[[1]] categories be equal
    # and all groupBy[[2]] category sizes equal *within* a groupBy[[1]] category
    # Usually, this means: Find maximum group size with equal # cells per sample within a group
    pdAgg <- aggregate(pd[,get(parentGate)], by=list(pd[,get(groupBy[[1]])], pd[,get(groupBy[[2]])]), sum)
    # pdAgg now has 3 columns, Group.1, Group.2, and x (the sum of all sample cell counts for that condition)
    finalGroup1Size <- min((pdAgg %>%
                              group_by(Group.1) %>%
                              summarize(maxGroupSize = min(x) * length(x)))
                           $maxGroupSize)
    
    cat("after grouping by '", groupBy[[1]], "' and '", groupBy[[2]], "', all '", groupBy[[1]], "', groups will have ~", 
        finalGroup1Size, "cells.\n")
    
    # Define the number of cells that should be sampled per groupBy[[2]] category, dependent on groupBy[[1]] category
    finalGroup2Sizes <- pdAgg %>%
      group_by(Group.1) %>%
      summarize(Group.2.Size = finalGroup1Size %/% length(Group.2), Group.2.Size.Remainder = finalGroup1Size %% length(Group.2)) # dividing by a non-divisor could result in unequal group sizes
    
    pd2 <- merge(pd, finalGroup2Sizes[,c("Group.1", "Group.2.Size")], by.x = groupBy[[1]], by.y = "Group.1")
    pd2$nameTmp <- pd2$name # in the likely case that groupBy[[2]] is "name", store it in an extra column so it doesn't get erased in the next step(?)
    
    # When calculating Group.2.Size, the divisor might not have divided evenly into the dividend, resulting in a remainder and slightly unequal final Group.1.Sizes.
    # If that's the case, this is where we distribute the remaining cells across samples within the Group.1 category/ies in which the remainder occured.
    # For each Group 1 value, distribute the remainder if applicable
    pd2 <- pd2[, {
      mySD <- copy(.SD)
      remainder <- finalGroup2Sizes[which(finalGroup2Sizes$Group.1 == .BY[[1]]), c("Group.2.Size.Remainder")][[1]]
      remainderOriginal <- remainder
      # Obtain rows where the current number of cells to be sampled is smaller than the pool of parentGate cells available for that row.
      availableRows <- which(mySD[,get(parentGate)] > mySD[,c("Group.2.Size")])
      while(remainder > 0 && length(availableRows) > 0) {
        numRowsToSample <- min(remainder, length(availableRows))
        # Sample without replacement "numRowsToSample" rows from availableRows
        rowsToIncrement <- sample(availableRows, numRowsToSample, replace=F)
        mySD[rowsToIncrement, Group.2.Size := mySD[rowsToIncrement, c("Group.2.Size")] + 1 ]
        # .SD[rowsToIncrement, c("Group.2.Size")] <- .SD[rowsToIncrement, c("Group.2.Size")] + 1

        # Prepare variables for the next loop
        remainder <- remainder - numRowsToSample
        availableRows <- which(mySD[,get(parentGate)] > mySD[,c("Group.2.Size")])
      }
      if(remainder > 0 && length(availableRows) == 0) {
        print(paste0("Alert: ", remainder, " cells unable to be sampled from Group ", groupBy[[1]], " = ", .BY[[1]], ". "))
      } else if(remainder == 0 && remainderOriginal > 0) {
        print(paste0(remainderOriginal, " cells successfully distributed across ", groupBy[[1]], " = ", .BY[[1]], " samples. "))
      }
      mySD
    }, by=c(groupBy[1])]
    
    # After the previous step, there should be equal number of cells for each Group.1 value:
    cellsPerGroup1 <- pd2[,{sum(.SD[,c("Group.2.Size")])}, by=c(groupBy[1])]
    if(!(length(unique(cellsPerGroup1$V1)) == 1)) {
      print("Unequal number of cells to be samples for each Group.1 value:")
      print(cellsPerGroup1)
    }
    
    pd2[, {
      nTcells <- .SD$Group.2.Size[[1]]
      
      # Sample from the total number of events
      # in the parent population. 
      totalEvents <- sum(get(parentGate))
      gInd <- 1:totalEvents
      # sample, without replacement, a vector of cell indices for the entire group, length totalEvents
      gInd <- sample.int(totalEvents, size = nTcells)
      gInd.logical <- rep(F, totalEvents)
      gInd.logical[gInd] <- T
      sn.factor <- unlist(sapply(nameTmp, function(sn) rep(sn, 
                                                           .SD[nameTmp == sn, get(parentGate)])))
      ind.vec <- split(gInd.logical, sn.factor)
      for (sn in nameTmp) {
        thisInd <- ind.vec[[sn]]
        gh <- gsClone[[sn]]
        updateIndices(gh, parentGate, thisInd)
      }
    }, by=groupBy]
    print(pd2)
  }
  
  cat("subsampling complete ! recomputing... \n")
  nodes <- getChildren(gsClone[[1]], parentGate, path = 2)
  for (node in nodes) recompute(gsClone, node)
  # Question: Why doesn't this include the parentGate? parentGate counts seem correct after, but I'm unsure why.
  
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
  
  return(list(input_mat = input_mat, res_collapse = res_collapse))
}

#' Run tSNE from (R pkg 'Rtsne' or 'Rtsne.multicore') on a GatingSet
#' 
#' If groupBy is empty, no sampling is done. If length is 1, equal number of cells are sampled from each category in the column.
#' If length is 2, equal number of cells are sampled from each category in the first column. Additionally, within each category,
#' equal numbers of cells are sampled from the second column (usually the Sample names column).
#' 
#' 
#' IMPORTANT: Requires a valid GatingSet with cytokine gates downstream of a parent gate
#' Also expects that pData(gs) contains at least columns: 'name', 'ptid' so we can identify cells later
#' 
#' @param gs a GatingSet object or path to a GatingSet directory, properly gated data with annotation in its pData
#' @param parentGate a \code{string} describing the gate upstream of the cytokine gates (eg. "CD4", "cd8+", etc...)
#' @param degreeFilterGates a \code{vector} of \code{strings} describing the marker gates immediately downstream of parentGate, eg: "CD4/IL2", "CD4/IFNg". Each entry is assumed to corresspond to one marker, which is used to degree filter the data and input into tsne.
#' @param otherMarkers a \code{vector} of \code{strings} describing other markers (aside from those corresponding to degreeFilterNodes) which are input into tsne but not used to degree filter
#' @param notRunMarkers a \code{vector} of \code{strings} describing markers which are NOT used for the run itself but are imported as columns in the final tsne output matrix
#' @param gateMarkerMap named list of gate names to marker names, eg. list("CD4/IL2" = "IL2","CD4/IFNg" = "IFNg"). Must contain entries for parentGate and degreeFilterGates. If entries are provided for otherMarkers and notRunMarkers, their corresponding xx_Bool columns will appear in the final tsne output matrix.
#' @param swap boolean for whether marker and gate names (from markerMap above) should be swapped. Passed onto getSingleCellExpression()
#' @param groupBy a \code{vector} of \code{strings} describing columns of the \code{gatingSet}'s phenoData. Affects sampling (see function description).
#' @param degreeFilter keep cells of this degree and higher, useful when tSNE takes too long to run
#' @param seed since tSNE is random, need a random seed so we can reproduce results
#' @param theta parameter to be passed to the \code{Rtsne} function
#' @param numThreads if > 1, uses multicore t-SNE instead of regular t-SNE
#' @param cloneGs boolean for whether to to clone the GatingSet. Cloning is safer when gs is a GatingSet object because the GatingSet may get modified in this function. Cloning takes up time and memory, however, so if you are not using the in-memory-GatingSet again you can choose FALSE.
#' @param ... other parameters to be passed to the \code{Rtsne} function
#' @return a \code{matrix} of X and Y coordinates
#' @export
#' @import flowWorkspace
#' @import data.table
#' @import plyr
#' @import magrittr
#' @import Rtsne
#' @import Rtsne.multicore
runTSNE <- function (gs=NULL, parentGate = NULL, degreeFilterGates = c(), otherMarkers = c(), notRunMarkers = c(), gateMarkerMap = NULL, swap = FALSE,
                     groupBy = c(), degreeFilter = 0, seed = 999, theta = 0.9, numThreads = 1, cloneGs = TRUE, ...) {
  data4tsne <- createTsneInputMatrix(gs = gs,
                                     parentGate = parentGate,
                                     degreeFilterGates = degreeFilterGates,
                                     otherMarkers = otherMarkers,
                                     notRunMarkers = notRunMarkers,
                                     gateMarkerMap = gateMarkerMap,
                                     swap = swap,
                                     groupBy = groupBy,
                                     degreeFilter = degreeFilter,
                                     seed = seed,
                                     theta = theta,
                                     cloneGs = cloneGs)

  cat("\n starting tSNE run at ", date(), " with ", numThreads, " threads\n")
  system.time(tsne_out <- if(numThreads > 1) {
    Rtsne.multicore::Rtsne.multicore(X = data4tsne$input_mat, check_duplicates = FALSE, num_threads = numThreads, 
                                  ...)
  } else {
    Rtsne(data4tsne$input_mat, check_duplicates = FALSE, 
                    ...)
  })
  
  dat <- tsne_out$Y
  colnames(dat) <- c("x", "y")
  dat <- cbind(dat, data4tsne$res_collapse)
  dat <- data.table(dat)
  cat("completed tSNE run at", date(), "!\n")
  return(dat)
}

#' A function to run One-SENSE. Has an additional dimensionMarkers parameter. For now, assuming people won't run more than 3 dimensions.
#' 
#' Note lack of numThreads argument. Multicore tsne doesn't work for One-SENSE (i.e. dims = 1). Since the input matrix dimensions will be smaller, it seems
#' that multicore tsne might not provide much of a speed-up anyway. See: https://github.com/DmitryUlyanov/Multicore-TSNE#what-to-expect
#' @param dimensionMarkers a list of character vectors, containing the markers which are to be run for each One-SENSE dimension
#' @export
#' @import flowWorkspace
#' @import data.table
#' @import plyr
#' @import magrittr
#' @import Rtsne
runOneSense <- function (gs=NULL, parentGate = NULL, degreeFilterGates = c(), otherMarkers = c(), notRunMarkers = c(), gateMarkerMap = NULL, swap = FALSE,
                     groupBy = c(), degreeFilter = 0, seed = 999, theta = 0.9, cloneGs = TRUE, dimensionMarkers = list(), ...) {
  if(length(dimensionMarkers) > 3) stop("There is a maximum of 3 output dimensions (feel free to modify code if you want more)")
  data4tsne <- createTsneInputMatrix(gs = gs,
                                     parentGate = parentGate,
                                     degreeFilterGates = degreeFilterGates,
                                     otherMarkers = otherMarkers,
                                     notRunMarkers = notRunMarkers,
                                     gateMarkerMap = gateMarkerMap,
                                     swap = swap,
                                     groupBy = groupBy,
                                     degreeFilter = degreeFilter,
                                     seed = seed,
                                     theta = theta,
                                     cloneGs = cloneGs)
  
  # Make sure that each marker in dimensionMarkers exists in colnames(data4tsne$input_mat). This corressponds to degreeFilterGates markers and otherMarkers
  lapply(dimensionMarkers, function(v) { if(!all(v %in% colnames(data4tsne$input_mat))) stop("all dimensionMarkers markers must exist in input_mat columns")})
  
  numThreads <- 1 # For some reason, Multicore tsne crashes when dims = 1
  cat("\n starting One-SENSE run at ", date(), " with ", numThreads, " threads\n")
  dat <- lapply(dimensionMarkers, function(d) {
    cat("\nStarting run for the dimension with the following markers: ", d, "\n")
    input_mat_sub <- data4tsne$input_mat[,d]
    system.time(tsne_out <- Rtsne(input_mat_sub, check_duplicates = FALSE, dims = 1, ...))
    tsne_out$Y
  })
  dat <- do.call(cbind, dat)
  
  colNames <- c("x", "y", "z") # Assuming people won't request more than 3 dimension categories
  colnames(dat) <- colNames[1:ncol(dat)]
  dat <- cbind(dat, data4tsne$res_collapse)
  dat <- data.table(dat)
  cat("completed One-SENSE run at", date(), "!\n")
  return(dat)
}
