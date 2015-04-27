#' Quantile normalize sample
#'
#' Normalize quantiles of a sample to supplied quantiles (Infinium HumanMethylation450 BeadChip).
#'
#' @param object An object created by \code{\link{meffil.compute.normalization.object}()}.
#' @param mu (Optional) Methylated and unmethylated intensities for a sample
#' created by \code{\link{meffil.rg.to.mu}}.  If this is not supplied,
#' then \code{object} is used to create such an object from the corresponding IDAT files.
#' The intensities will have been background corrected and dye bias corrected.
#' @param subsets Named list of probe sets (see \code{quantiles}).
#' @param quantiles Named list of numeric sequences defining the target distributions
#' of probes in corresponding sets specified by \code{subsets}.
#' For a set of probes \code{subset[[name]]}, \code{quantiles[[name]]$M}
#' is an ordered sequence
#' of numbers defining the quantiles of the target distribution
#' the methylated signal of the probes.  Similarly
#' \code{quantiles[[name]]$U} defines the target distribution of their
#' unmethylated signal distribution.
#' @param probes Probe annotation used to construct the control matrix
#' (Default: \code{\link{meffil.probe.info}()}).
#' @param verbose If \code{TRUE}, then status messages are printed during execution (Default: \code{FALSE}).
#' @return \code{mu} modified so that for probe set \code{subsets[[name]]},
#' \code{mu$M[subsets[[name]]]} has the distribution defined by \code{quantiles[[name]]$M}
#' and \code{mu$U[subsets[[name]]]} has the distribution defined by \code{quantiles[[name]]$U}.
#'
#' @export
meffil.quantile.normalize.sample <- function(object, mu=NULL, subsets, quantiles,
                                             dye.intensity=5000,
                                             probes=meffil.probe.info(), verbose=F) {
    stopifnot(is.normalization.object(object))
    stopifnot(length(subsets) == length(quantiles))
    stopifnot(all(names(subsets) %in% names(quantiles)))

    if (is.null(mu)) {
        rg <- meffil.read.rg(object$basename, probes, verbose=verbose)
        rg <- meffil.background.correct(rg, probes, verbose=verbose)
        rg <- meffil.dye.bias.correct(rg, dye.intensity, probes, verbose=verbose)
        mu <- meffil.rg.to.mu(rg, probes, verbose=verbose)
    }

    for (subset.name in names(subsets)) {
        subset <- subsets[[subset.name]]
        for (target in names(mu)) {
            data <- matrix(mu[[target]][subset])
            full.quantiles <- quantiles[[subset.name]][[target]]
            full.quantiles <- approx(1:length(full.quantiles), full.quantiles, 1:length(data))$y
            mu[[target]][subset] <- preprocessCore::normalize.quantiles.use.target(data,full.quantiles)
        }
    }
    mu
}
