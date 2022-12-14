% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R4scSHARP.R
\name{run_scsorter}
\alias{run_scsorter}
\title{Runs scSorter}
\usage{
run_scsorter(data, markers = NULL, ref = NULL)
}
\arguments{
\item{data}{Gene expression matrix}

\item{markers}{Marker genes for cell classification}

\item{ref}{Reference dataset for predictions.}
}
\value{
Returns scSorter predictions as an R data frame
}
\description{
Runs scSorter
}
