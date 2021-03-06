% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlations.R
\name{pairwise_correlation_counts}
\alias{pairwise_correlation_counts}
\title{pairwise correlation counts}
\usage{
pairwise_correlation_counts(keep_matrix)
}
\arguments{
\item{keep_matrix}{\code{logical} matrix}
}
\value{
matrix
}
\description{
given the \code{keep} entry from \code{pairwise_correlation}, get the number
of things that would have been used for each of the pairwise comparisons. Note
that pairwise comparisons of the rows are used.
}

