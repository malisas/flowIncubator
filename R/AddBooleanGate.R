#' Adds boolean gate to the GatingSet or GatingSetList
#' 
#' @param gs Either a GatingSet, a GatingSetList, or path to one of these objects on disk
#' @param booleanSubset The booleanSubset (a combination of existing gates) in string format, e.g. "8+/GMM+&!8+/GAMMADELTA"
#' @param parentGate The gate under which the booleansubset is added
#' @param booleanGateName optional. What to call the new gate
#' @import flowWorkspace
#' @export addBooleanGate
addBooleanGate <- function(gs,
               booleanSubset,
               parentGate,
               overrideGate=FALSE,
               booleanGateName=NULL) {
  try(if(missing(gs) || missing(booleanSubset) || missing(parentGate)) stop("Required arguments missing.") )
  # Check if booleanSubset already exists under parentGate
  booleanSubsetName <- if(is.null(booleanGateName)) {
    gsub("/", ":", paste0(parentGate, ":", booleanSubset))
  } else {
    booleanGateName
  }
  message("Checking for gate ", booleanSubsetName)
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
}