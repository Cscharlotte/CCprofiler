% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getBestFeatureParameters.R
\name{getBestParameterStats_constraintFDR}
\alias{getBestParameterStats_constraintFDR}
\title{getBestParameterStats_constraintFDR}
\usage{
getBestParameterStats_constraintFDR(complex_feature_grid_filtered,
  complex_hypotheses, n_subsets = 3, FDR_cutoff = 0.1,
  grid_search_params = c("corr", "window", "rt_height", "smoothing_length",
  "min_feature_completeness", "min_hypothesis_completeness", "min_subunits",
  "min_peak_corr", "min_monomer_distance_factor"))
}
\arguments{
\item{complex_feature_grid_filtered}{list with all features from grid search.}

\item{complex_hypotheses}{data.table with complex hypotheses including decoys.}

\item{n_subsets}{Numeric, number of subsets used for cross-calculation of FDR, default = 3.}

\item{FDR_cutoff}{Numeric, maximum FDR that should be considered, default = 0.1.}
}
\value{
data.table one row with best stats across grid search.
}
\description{
getBestParameterStats_constraintFDR. Do x-fold cross-FDR estimation.
Works only for complexFeatures at the moment.
}
\examples{
## Generate example data
complexFeaturesGrid <- list(exampleComplexFeatures)
complexHypotheses <- exampleComplexHypotheses
complexFeaturesGridFiltered <- filterGridSearchResults(complexFeaturesGrid,
                                                       peak_corr_cutoffs = c(0.5,0.75,0.9),
                                                       feature_completeness_cutoffs = c(0,0.5,1),
                                                       hypothesis_completeness_cutoffs = c(0.5,1),
                                                       n_subunits_cutoffs =c(2,3,4),
                                                       monomer_distance_cutoffs = c(1,2),
                                                       remove_decoys=FALSE)

## Pick the best parameter set with a cross checked FDR constraint of 0.1
getBestParameterStats_constraintFDR(complexFeaturesGridFiltered,
                                    complexHypotheses,
                                    FDR = 0.1,
                                    grid_search_params =c("min_feature_completeness",
                                                          "min_hypothesis_completeness",
                                                          "min_subunits",
                                                          "min_peak_corr",
                                                          "min_monomer_distance_factor"))

}

