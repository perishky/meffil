% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc-report.r
\name{meffil.plot.cell.counts}
\alias{meffil.plot.cell.counts}
\title{Cell count estimate quality plot}
\usage{
meffil.plot.cell.counts(qc.objects)
}
\arguments{
\item{qc.objects}{Output from \code{\link{meffil.qc}()}.}

\item{reference}{Object describing methylation profiles of purified cell populations
obtained from \code{\link{meffil.add.cell.type.reference}()}.}
}
\value{
Two \link{ggplot2} boxplot objects:
\itemize{
\item \code{betas} Contains one box per sample or reference cell type
representing the distribution of methylation levels for the CpG sites used to
estimate cell counts.
\item \code{counts} Contains one box per reference cell type
representing the distribution of cell count estimates across the samples.
}
}
\description{
Cell count estimate quality plot
}
