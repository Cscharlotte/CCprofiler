% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findProteinFeatures.R
\name{findProteinFeatures}
\alias{findProteinFeatures}
\title{Protein feature detection}
\usage{
findProteinFeatures(traces, corr_cutoff = 0.95, window_size = 12,
  parallelized = FALSE, n_cores = 1, collapse_method = "apex_only",
  perturb_cutoff = "5\%", rt_height = 5, smoothing_length = 9,
  useRandomDecoyModel = TRUE)
}
\arguments{
\item{traces}{An object of class traces (type "peptide").}

\item{corr_cutoff}{Numeric, the correlation value for chromatograms above which
peptides are considered to be coeluting, default=0.95.}

\item{window_size}{Numeric, size of the window (in fractions), default=12}

\item{parallelized}{Logical, wether the computation should be done in parallel, default=FALSE}

\item{n_cores}{Integer, the number of cores to use for parallel processing
(only applies if parallelized is TRUE), default=1}

\item{collapse_method}{Method for collapsing multiple features into one feature: 
\itemize{
\item "apex_only": collapses by apex
\item "apex_network": collapses by apex and connected network cluster
}
Default="apex_only"}

\item{perturb_cutoff}{Numeric, the quantile to use in estimating the perturbation level, default="5%".
Intensity values that are zero are replaced with random values that are
below the specified quantile of the input values. Alternatively a
cutoff value can be specified as an upper limit for perturbation values.
This is nescessary for correlation calculation.}

\item{rt_height}{Numeric, RT cutoff for collapsing features, default is 5}

\item{smoothing_length}{Numeric, smoothing length of Savitzky-Golay filter, default is 7}

\item{useRandomDecoyModel}{Logical, wether random peptide protein associations should be used as decoy model, default = TRUE}
}
\value{
A data.table containing protein features.
}
\description{
Run the sliding window algorithm to find protein features.
}
\examples{
## Load example data
peptideTraces <- examplePeptideTracesFiltered
## Subset traces for shorter processing time
testProteins = unique(peptideTraces$trace_annotation$protein_id)[1:5]
peptideTracesSubset = subset(peptideTraces,trace_subset_ids = testProteins, trace_subset_type = "protein_id")

## Perform co-elution signal detection
proteinFeatures <- findProteinFeatures(traces=peptideTracesSubset)

## Inspect complex features
head(proteinFeatures,n=3)
}

