% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-pc-fit.r
\name{meffil.plot.pc.fit}
\alias{meffil.plot.pc.fit}
\title{Number of control matrix principal components}
\usage{
meffil.plot.pc.fit(
  qc.objects,
  fixed.effects = NULL,
  random.effects = NULL,
  n.cross = 10,
  name = "autosomal.ii"
)
}
\arguments{
\item{qc.objects}{A list of outputs from \code{\link{meffil.create.qc.object}()}.}

\item{fixed.effects}{Names of columns in samplesheet that should be included as fixed effects
along with control matrix principal components (Default: NULL).}

\item{random.effects}{Names of columns in samplesheet that should be included as random effects
(Default: NULL).}

\item{number.pcs}{Number of principal components to include in the design matrix (Default: all).}
}
\value{
A list containing a data frame with the mean squared residuals for different numbers of principal components
and a plot of these residuals.
}
\description{
Fits probe intensities to principal components of the microarray control matrix
and calculates the resulting mean squared residuals for different
numbers of principal components.
}
