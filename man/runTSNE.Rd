% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runTSNE.R
\name{runTSNE}
\alias{runTSNE}
\title{run tSNE from (R pkg 'Rtsne') on a gatingSet
Will sample the minimal number of cells available in all samples to generate balanced cell counts}
\usage{
runTSNE(gs, parentGate, cytokines, otherMarkers, markerMap, swap = FALSE,
  groupBy, degreeFilter = 0, seed = 999, theta = 0.9, ...)
}
\arguments{
\item{gs}{a GatingSet object, properly gated data with annotation in its pData}

\item{parentGate}{a \code{string} describing the gate upstream of the cytokine gates (eg. "CD4", "cd8+", etc...)}

\item{otherMarkers}{the remaining markers of the data}

\item{markerMap}{named list of marker names to gate names, eg. list("CD4/IL2" = "IL2","CD4/IFNg" = "IFNg")}

\item{swap}{boolean for whether marker and gate names (from markerMap above) should be swapped. Passed onto getSingleCellExpression()}

\item{groupBy}{columns of the \code{gatingSet}'s phenoData, same number of cells will be sampled from each group}

\item{degreeFilter}{keep cells of this degree and higher, useful when tSNE takes too long to run}

\item{seed}{since tSNE is random, need a random seed so we can reproduce results}

\item{theta}{parameter to be passed to the \code{Rtsne} function}

\item{...}{other parameters to be passed to the \code{Rtsne} function}

\item{cytokine}{a \code{vector} of \code{strings} describing the cytokine gates immediately downstream of parentGate, eg: "IL2", "IFNg"}
}
\value{
a \code{matrix} of X and Y coordinates
}
\description{
IMPORTANT: Requires a valid gatingSet with cytokine gates downstream of a parent gate
Also expects that pData(gs) contains at least columns: 'name', 'ptid' so we can identify cells later
}

