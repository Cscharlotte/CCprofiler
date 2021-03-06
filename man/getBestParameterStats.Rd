% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getBestFeatureParameters.R
\name{getBestParameterStats}
\alias{getBestParameterStats}
\title{Get the best performing ParameterStats}
\usage{
getBestParameterStats(grid_search_stats, FDR = 0.1)
}
\arguments{
\item{grid_search_stats}{Table of grid search statistics.
(obtained from \code{\link{estimateGridSearchDecoyFDR}}).}

\item{FDR}{Numeric, maximum FDR that should be considered. Default = \code{0.1}.}
}
\value{
data.table one row with best stats across grid search.
}
\description{
Pick the parameter set in a grid search stats output that performs
best in terms of identifications while staying within a specified FDR cutoff.
}
\examples{

## NOT RUN
gridStats # see function \\code{\\link{estimateGridSearchDecoyFDR}} to see how to generate this object.

## Pick best parameter set
getBestParameterStats(gridStats)

}

