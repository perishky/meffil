% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/plot-ewas.r
\name{meffil.ewas.cpg.plot}
\alias{meffil.ewas.cpg.plot}
\title{Scatter plots for a CpG site in an EWAS}
\usage{
meffil.ewas.cpg.plot(ewas.object, cpg, title = cpg, beta = NULL)
}
\arguments{
\item{ewas.object}{Return object from \code{\link{meffil.ewas()}}.}

\item{cpg}{CpG site to plot.}

\item{title}{Title of the plot (Default: \code{cpg}).}

\item{beta}{Optional matrix of methylation levels used to create the \code{ewas.object} (Default: NULL).}

\item{\code{\link{ggplot}}}{object showing the scatterplots of DNA methylation vs the variable of interest
in the EWAS.  Each plot corresponds to a covariate set.
Methylation levels are in fact residuals from fitting a model with DNA methylation and the covariates.}
}
\description{
Scatter plots for a CpG site in an EWAS
}
