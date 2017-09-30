#' Stratified Boxplots of the Proportion of the Boolean/Existing Subset seen
#' 
#' Creates a stratified boxplot of the per-patient proportions of the given subset.
#' Iff booleanSubset is provided, adds the given booleanSubset to the gating tree first.
#'
#' @param gsOrGsListOrPath Either a GatingSet, a GatingSetList, or path to one of these objects on disk
#' @param existingSubset (existingSubset or booleanSubset must be provided) The name of the existing subset of interest, e.g. "CD3/GMM"
#' @param booleanSubset (existingSubset or booleanSubset must be provided) The booleanSubset (a combination of existing gates) in string format, e.g. "8+/GMM+&!8+/GAMMADELTA"
#' @param parentGate The gate under which the booleansubset is added
#' @param parentGateForProportionCalc The parent gate which is to be used in proportion calculations (e.g. "8+" or "3+"). Defaults to parentGate
#' @param groupBy optional string. What variable to stratify the boxplot by
#' @param outdir Where to save the boxplot and stats, if desired.
#' @param overrideGate If this is set to TRUE, booleanSubset will override any previous node of the same name under the same parentGate.
#' @param ylimits optional numeric vector
#' @param sampleIDCol (optional) Name of the GatingSet metadata column which contains individual sample identifiers. If sampleIDCol is provided, outliers will be labeled with this value.
#' @return boxplot and stats, and optionally saves the plot to outdir
#' @import coin
#' @import data.table
#' @import flowWorkspace
#' @import grDevices
#' @import svglite
#' @export boxplot.subset.proportion.stats
#' @keywords Box Plot BooleanSubset Proportion
#' @examples
#' \dontrun{
#' boxplot.subset.proportion.stats(gsOrGsListOrPath="/Path/To/GatingSet",
#'                                 booleanSubset="GMM&GAMMADELTA",
#'                                 parentGate="CD3",
#'                                 groupBy="HIVStatus",
#'                                 outdir="Path/To/Plot/Directory",
#'                                 overrideGate=FALSE,
#'                                 ylimits=NULL)
#' }
boxplot.subset.proportion.stats <- function(gsOrGsListOrPath,
                                            existingSubset=NULL,
                                            booleanSubset=NULL,
                                            parentGate,
                                            parentGateForProportionCalc=parentGate,
                                            groupBy=NULL,
                                            outdir=NULL,
                                            overrideGate=FALSE,
                                            ylimits=NULL,
                                            baseSize=19,
                                            sampleIDCol=NULL
) {
  try(if(missing(gsOrGsListOrPath) || (is.null(booleanSubset) && is.null(existingSubset)) || missing(parentGate)) stop("Required arguments missing.") )
  
  gs <- if(class(gsOrGsListOrPath) == "GatingSet" || class(gsOrGsListOrPath) == "GatingSetList") {
    gsOrGsListOrPath
  } else {
    try(if(!(class(gsOrGsListOrPath) == "character")) stop("gsOrGsListOrPath must be either a GatingSet, a GatingSetList, or the path to the folder containing one of these objects on disk"))
    # Load the saved GatingSetList or GatingSet:
    loadGSListOrGS <- function (gsOrGsListOrPath) {
      out <- try(flowWorkspace::load_gslist(gsOrGsListOrPath))
      if (class(out) == "try-error") {
        cat("Caught an error during flowWorkspace::load_gslist, trying flowWorkspace::load_gs.\n")
        out <- flowWorkspace::load_gs(gsOrGsListOrPath)
      }
      out
    }
    loadGSListOrGS(gsOrGsListOrPath)
  }
  
  subsetInfo <- if(!(is.null(booleanSubset))) {
    message(paste0("booleanSubset ", booleanSubset, " given"))
    
    # Check if booleanSubset already exists under parentGate
    booleanSubsetName <- paste0(parentGate, ":", booleanSubset)
    booleanSubsetName <- gsub("/", ":", booleanSubsetName)
    if(booleanSubsetName %in% getNodes(gs, path="auto")) {
      # If the gate already exists, decide what to do with it.
      if(overrideGate) {
        message(paste0("Gate ", booleanSubsetName, " already exists. Deleting old gate and adding new gate..."))
        call <- substitute(flowWorkspace::booleanFilter(v), list(v = as.symbol(booleanSubset)))
        g <- eval(call)
        flowWorkspace::Rm(booleanSubsetName, gs)
        flowWorkspace::add(gs, g, parent = parentGate, name=booleanSubsetName)
        flowWorkspace::recompute(gs, booleanSubsetName)
      } else {
        message(paste0("Gate ", booleanSubsetName, " already exists. Keeping old gate..."))
      }
    } else {
      # If the gate doesn't exist, add it.
      message(paste0("Adding new gate ", booleanSubsetName, " to GatingSet/GatingSetList..."))
      call <- substitute(flowWorkspace::booleanFilter(v), list(v = as.symbol(booleanSubset)))
      g <- eval(call)
      flowWorkspace::add(gs, g, parent = parentGate, name=booleanSubsetName)
      flowWorkspace::recompute(gs, booleanSubsetName)
    }
    
    list("subsetName" = booleanSubsetName, "popStats" = flowWorkspace::getPopStats(gs, flowJo=FALSE, subpopulations=c(booleanSubsetName)))
  } else if(!is.null(existingSubset)) {
    message(paste0("existingSubset ", existingSubset, " given"))
    if(! (existingSubset %in% getChildren(gs[[1]], parentGateForProportionCalc, path="auto"))) {
      stop(paste0("existingSubset ", existingSubset, " not a child of ", parentGateForProportionCalc))
      }

    list("subsetName" = existingSubset, "popStats" = flowWorkspace::getPopStats(gs, flowJo=FALSE, subpopulations=c(existingSubset)))
  }

  booleanSubsetPopStats <- subsetInfo$popStats
  parentGateForProportionCalcPopStats <- flowWorkspace::getPopStats(gs, flowJo=FALSE, subpopulations=c(parentGateForProportionCalc))
  mergedPopStats <- merge(booleanSubsetPopStats[,c("name", "Count")], parentGateForProportionCalcPopStats[,c("name", "Count")], by="name", suffixes=c(".boolSubset", ".parentForProportion"))
  mergedPopStats$Proportion <- mergedPopStats$Count.boolSubset / mergedPopStats$Count.parentForProportion
  
  # Function that takes in vector of data and a coefficient,
  # returns boolean vector if a certain point is an outlier or not
  # From: https://stackoverflow.com/a/33525306/6282213
  check_outlier <- function(v, coef=1.5){
    quantiles <- quantile(v,probs=c(0.25,0.75))
    IQR <- quantiles[2]-quantiles[1]
    res <- v < (quantiles[1]-coef*IQR)|v > (quantiles[2]+coef*IQR)
    return(res)
  }
  pData4Plot <- pData(gs)
  pData4Plot$name <- rownames(pData4Plot)
  mergedPopStats <- merge(mergedPopStats, pData4Plot, by="name")
    
  mergedPopStats$GroupTmp <- as.factor(mergedPopStats[, get(groupBy)])
  subsetNameForPlot <- subsetInfo$subsetName
  myPlot <- if(is.null(groupBy)) {
    # Make a simple un-stratified boxplot of the proportion of subsetNameForPlot cells

    plottitle <- paste0("Boxplot of Proportion of\n", subsetNameForPlot, " Cells\nof total ", parentGateForProportionCalc, " Cells")
    p <- ggplot2::ggplot(mergedPopStats, ggplot2::aes(x=factor(0), y = Proportion)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_jitter() +
      ggplot2::theme_set(ggplot2::theme_gray(base_size = baseSize)) +
      ggplot2::theme(plot.title=ggplot2::element_text(vjust=-0.8, hjust=0.5)) +
      ggplot2::labs(x="All Samples", y=paste0("Proportion of ", parentGateForProportionCalc, " Cells"),
                    title=plottitle)
    if(!is.null(sampleIDCol)) {
      mergedPopStats[,outlier:=check_outlier(Proportion)]
      mergedPopStats[,label:=ifelse(outlier, get(sampleIDCol),"")]
      p <- p + ggplot2::geom_text(ggplot2::aes(label=label),vjust=-0.2)
    }

    if(!is.null(ylimits)) {
      p <- p + ggplot2::coord_cartesian(ylim=ylimits) # adjust visible data for y axis but keep points
    }
    p
  } else {
    # Make a stratified boxplot of the proportion of subsetNameForPlot cells
    test <- coin::wilcox_test(Proportion ~ GroupTmp, data=mergedPopStats)

    plottitle <- paste0("Boxplot of Proportion of\n", subsetNameForPlot, " Cells\nof total ", parentGateForProportionCalc, " Cells")
    subtitle <- paste0("Stratified by ", groupBy, "\np = ", signif(coin::pvalue(test), 4), ",  Z = ", signif(coin::statistic(test), 4))
    p <- ggplot2::ggplot(mergedPopStats, ggplot2::aes_string(x = "GroupTmp", y = "Proportion")) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_jitter() +
      ggplot2::theme_set(ggplot2::theme_gray(base_size = baseSize)) +
      ggplot2::theme(plot.title=ggplot2::element_text(vjust=-0.8, hjust=0.5)) +
      ggplot2::labs(x=groupBy, y=paste0("Proportion of ", parentGateForProportionCalc, " Cells"),
                    title=plottitle, subtitle=subtitle)
    if(!is.null(sampleIDCol)) {
      mergedPopStats[,outlier:=check_outlier(Proportion), by=GroupTmp]
      mergedPopStats[,label:=ifelse(outlier, get(sampleIDCol),"")]
      p <- p + ggplot2::geom_text(ggplot2::aes(label=label),vjust=-0.2)
    }
    if(!is.null(ylimits)) {
      p <- p + ggplot2::coord_cartesian(ylim=ylimits) # adjust visible data for y axis but keep points
    }
    p
  }
  # if(is.null(outdir)) {
  #   list(plot = myPlot, data = mergedPopStats)
  # } else {
  #   filename <- if(is.null(groupBy)) {
  #     paste0("BoxplotProportion_", subsetNameForPlot, "_of_", parentGateForProportionCalc, ".png")
  #   } else {
  #     paste0("BoxplotProportion_", subsetNameForPlot, "_of_", parentGateForProportionCalc ,"_by_", groupBy, ".png")
  #   }
  #   message(paste0("Saving plot to ", file.path(outdir, filename)))
  #   ggplot2::ggsave(filename=file.path(outdir, filename),
  #                   plot=myPlot,
  #                   width=6.66, height=7.85)
  # }
  if(!(is.null(outdir))) {
    parentGateForProportionCalc4File <- gsub("/", ":", parentGateForProportionCalc)
    filename <- if(is.null(groupBy)) {
      paste0("BoxplotProportion_", subsetNameForPlot, "_of_", parentGateForProportionCalc4File, ".png")
    } else {
      paste0("BoxplotProportion_", subsetNameForPlot, "_of_", parentGateForProportionCalc4File ,"_by_", groupBy, ".png")
    }
    message(paste0("Saving plot to ", file.path(outdir, filename)))
    ggplot2::ggsave(filename=file.path(outdir, filename),
                    plot=myPlot,
                    width=8, height=7.85)
  }
  list(plot = myPlot, data = mergedPopStats)
}