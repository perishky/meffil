#' Normalize Infinium HumanMethylation450 BeadChip
#'
#' Normalize sample methylation data using normalized quantiles.
#'
#' @param norm.object An element of \code{\link{meffil.normalize.quantiles}()}.
#' @param verbose If \code{TRUE}, then status messages are printed during execution (Default: \code{FALSE}).
#' @return List containing normalized methylated and unmethylated signals.
#'
#' @export
meffil.normalize.sample <- function(norm.object, verbose=F) {
    stopifnot(is.normalized.object(norm.object))

    probe.names <- meffil.get.sites()

    rg <- read.rg(norm.object$basename, verbose=verbose)
    rg <- background.correct(rg, verbose=verbose)
    rg <- dye.bias.correct(rg, norm.object$reference.intensity, verbose=verbose)
    mu <- rg.to.mu(rg)
        
    mu$M <- mu$M[probe.names]
    mu$U <- mu$U[probe.names]

    msg("Normalizing methylated and unmethylated signals.", verbose=verbose)
    probe.subsets <- get.quantile.probe.subsets()
    for (name in names(norm.object$norm)) {
        for (target in names(norm.object$norm[[name]])) {
            probe.idx <- which(names(mu[[target]]) %in% probe.subsets[[name]])
            if (length(probe.idx) > 0) {
                orig.signal <- mu[[target]][probe.idx]
                norm.target <- compute.quantiles.target(norm.object$norm[[name]][[target]])
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
