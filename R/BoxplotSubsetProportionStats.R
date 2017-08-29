#' Stratified Boxplots of the Proportion of the Boolean Subset seen
#' 
#' Adds the given booleanSubset to the gating tree and creates a stratified boxplot of the per-patient proportions of this subset.
#'
#' @param gsOrGsListOrPath Either a GatingSet, a GatingSetList, or path to one of these objects on disk
#' @param booleanSubset The booleanSubset (a combination of existing gates) in string format, e.g. "8+/GMM+&!8+/GAMMADELTA"
#' @param parentGate The gate under which the booleansubset is added
#' @param parentGateForProportionCalc The parent gate which is to be used in proportion calculations (e.g. "8+" or "3+"). Defaults to parentGate
#' @param groupBy optional string. What variable to stratify the boxplot by
#' @param outdir Where to save the boxplot and stats, if desired
#' @param overrideGate If this is set to TRUE, booleanSubset will override any previous node of the same name under the same parentGate.
#' @param ylimits optional numeric vector
#' @return boxplot and stats, unless outdir is specified
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
                                            booleanSubset,
                                            parentGate,
                                            parentGateForProportionCalc=parentGate,
                                            groupBy=NULL,
                                            outdir=NULL,
                                            overrideGate=FALSE,
                                            ylimits=NULL
) {
  try(if(missing(gsOrGsListOrPath) || missing(booleanSubset) || missing(parentGate)) stop("Required arguments missing.") )
  
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
  
  booleanSubsetPopStats <- flowWorkspace::getPopStats(gs, flowJo=FALSE, subpopulations=c(booleanSubsetName))
  parentGateForProportionCalcPopStats <- flowWorkspace::getPopStats(gs, flowJo=FALSE, subpopulations=c(parentGateForProportionCalc))
  mergedPopStats <- merge(booleanSubsetPopStats[,c("name", "Count")], parentGateForProportionCalcPopStats[,c("name", "Count")], by="name", suffixes=c(".boolSubset", ".parentForProportion"))
  mergedPopStats$Proportion <- mergedPopStats$Count.boolSubset / mergedPopStats$Count.parentForProportion
  
  myPlot <- if(is.null(groupBy)) {
    # Make a simple un-stratified boxplot of the proportion of booleanSubsetName cells
    plottitle <- paste0("Boxplot of Proportion of\n", booleanSubsetName, " Cells\nof total ", parentGateForProportionCalc, " Cells")
    p <- ggplot2::ggplot(mergedPopStats, ggplot2::aes(x=factor(0), y = Proportion)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_jitter() +
      ggplot2::theme_set(ggplot2::theme_gray(base_size = 19)) +
      ggplot2::theme(plot.title=ggplot2::element_text(vjust=-0.8, hjust=0.5)) +
      ggplot2::labs(x="All Samples", y=paste0("Proportion of ", parentGateForProportionCalc, " Cells"),
                    title=plottitle)
    if(!is.null(ylimits)) {
      p <- p + ggplot2::coord_cartesian(ylim=ylimits) # adjust visible data for y axis but keep points
    }
    p
  } else {
    pData4Plot <- pData(gs)
    pData4Plot$name <- rownames(pData4Plot)
    mergedPopStatsMeta <- merge(mergedPopStats, pData4Plot, by="name")
    mergedPopStatsMeta[, groupBy] <- as.factor(mergedPopStatsMeta[, get(groupBy)])
    test <- coin::wilcox_test(Proportion ~ get(groupBy), data=mergedPopStatsMeta)
    
    plottitle <- paste0("Boxplot of Proportion of\n", booleanSubsetName, " Cells\nof total ", parentGateForProportionCalc, " Cells")
    subtitle <- paste0("Stratified by ", groupBy, "\np = ", signif(coin::pvalue(test), 4), ",  Z = ", signif(coin::statistic(test), 4))
    p <- ggplot2::ggplot(mergedPopStatsMeta, ggplot2::aes_string(x = get("groupBy"), y = "Proportion")) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_jitter() +
      ggplot2::theme_set(ggplot2::theme_gray(base_size = 19)) +
      ggplot2::theme(plot.title=ggplot2::element_text(vjust=-0.8, hjust=0.5)) +
      ggplot2::labs(x=groupBy, y=paste0("Proportion of ", parentGateForProportionCalc, " Cells"),
                    title=plottitle, subtitle=subtitle)
    if(!is.null(ylimits)) {
      p <- p + ggplot2::coord_cartesian(ylim=ylimits) # adjust visible data for y axis but keep points
    }
    p
  }
  if(is.null(outdir)) {
    myPlot
  } else {
    filename <- paste0("BoxplotProportion_", booleanSubsetName, "_of_", parentGateForProportionCalc ,"_by_", groupBy, ".png")
    message(paste0("Saving plot to ", file.path(outdir, filename)))
    ggplot2::ggsave(filename=file.path(outdir, filename),
                    plot=myPlot,
                    width=6.66, height=7.85)
  }
}