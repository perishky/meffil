% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/handle-outliers.r
\name{meffil.handle.outliers}
\alias{meffil.handle.outliers}
\title{Handle outliers in a methylation matrix}
\usage{
meffil.handle.outliers(beta, winsorize.pct = 0.05, outlier.iqr.factor = NA)
}
\arguments{
\item{beta}{Methylation matrix (rows=CpG sites, columns=samples, values=methylation levels).}

\item{winsorize.pct}{Apply all regression models to methylation levels
winsorized to the given level. Set to NA to avoid winsorizing (Default: 0.05).}

\item{outlier.iqr.factor}{For each CpG site, prior to fitting regression models,
set methylation levels less than
\code{Q1 - outlier.iqr.factor * IQR} or more than
\code{Q3 + outlier.iqr.factor * IQR} to NA. Here IQR is the inter-quartile
range of the methylation levels at the CpG site, i.e. Q3-Q1.
Set to NA to skip this step (Default: NA).}
}
\value{
\code{beta} after winsorizing and outliers set to NA.
}
\description{
Handle outliers in a methylation matrix
}
