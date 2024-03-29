% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cell-type-specific-methylation.r
\name{meffil.cell.type.specific.methylation}
\alias{meffil.cell.type.specific.methylation}
\title{Reduce methylation profiles to most cell-type specific sites}
\usage{
meffil.cell.type.specific.methylation(
  beta,
  cell.types,
  number.sites = 50,
  verbose = F
)
}
\arguments{
\item{beta}{Numeric matrix (values = 0..1; rows = CpG sites; columns = samples).}

\item{cell.types}{Name of cell type for each column of beta.}

\item{number.sites}{For each cell type, the number of sites less methylated and the number
more methylated than other cell types to include in the reduced methylation profiles.}
}
\value{
Numeric matrix (values = 0..1; rows = CpG sites; columns = cell types)
with \code{number.sites} CpG sites per cell type more methylated than other cell types
and the same number less methylated.  Values are the mean CpG site methylation levels
of all original samples of the same cell type.#
}
\description{
Reduce methylation profiles to most cell-type specific sites
}
