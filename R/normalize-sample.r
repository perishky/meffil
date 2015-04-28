#' Normalize Infinium HumanMethylation450 BeadChip
#'
#' Normalize sample methylation data using normalized quantiles.
#'
#' @param object A list element from output of \code{\link{meffil.normalize.objects}()}.
#' @param probes Probe annotation used to construct the control matrix
#' (Default: \code{\link{meffil.probe.info}()}).
#' @param verbose If \code{TRUE}, then status messages are printed during execution (Default: \code{FALSE}).
#' @return List containing normalized methylated and unmethylated signals.
#'
#' @examples
#'
#' path <- ...
#' basenames <- meffil.basenames(path)
#' norm.objects <- lapply(basenames, function(basename) {
#'   meffil.compute.normalization.object(basename)
#' })
#' norm.objects <- meffil.normalize.objects(norm.objects, number.pcs=2)
#' mu1 <- mefill.normalize.sample(norm.objects[[1]])
#' beta1 <- mefill.get.beta(mu1)
#'
#' @export
meffil.normalize.sample <- function(object, mu=NULL, probes=meffil.probe.info(), verbose=F) {
    stopifnot(is.normalization.object(object))
    stopifnot("norm" %in% names(object))
    
    probe.names <- unique(na.omit(probes$name))
    probe.names <- probe.names[which(substring(probe.names,1,2) %in% c("cg","ch"))]

    if (is.null(mu)) {
        rg <- meffil.read.rg(object$basename, probes, verbose=verbose)
        rg <- meffil.background.correct(rg, probes, verbose=verbose)
        rg <- meffil.dye.bias.correct(rg, object$reference.intensity, probes, verbose=verbose)
        mu <- meffil.rg.to.mu(rg, probes, verbose=verbose)
        
        mu$M <- mu$M[probe.names]
        mu$U <- mu$U[probe.names]
    }

    msg("Normalizing methylated and unmethylated signals.", verbose=verbose)
    probe.subsets <- get.quantile.probe.subsets(probes)
    for (name in names(object$norm)) {
        for (target in names(object$norm[[name]])) {
            probe.idx <- which(names(mu[[target]]) %in% probe.subsets[[name]])
            if (length(probe.idx) > 0) {
                orig.signal <- mu[[target]][probe.idx]
                norm.target <- compute.quantiles.target(object$norm[[name]][[target]])
                norm.signal <- preprocessCore::normalize.quantiles.use.target(matrix(orig.signal),
                                                                              norm.target)
                mu[[target]][probe.idx] <- norm.signal
            }
        }
    }
    mu
}

compute.quantiles.target <- function(quantiles) {
    n <- length(quantiles)
    unlist(lapply(1:(n-1), function(j) {
        start <- quantiles[j]
        end <- quantiles[j+1]
        seq(start,end,(end-start)/n)[-n]
    }))
}
