#' Normalize Infinium HumanMethylation450 BeadChip
#'
#' Normalize sample methylation data using normalized quantiles.
#'
#' @param norm.object An element of \code{\link{meffil.normalize.quantiles}()}.
#' @param remove.poor.signal Set methylation values for poorly detected probes
#' to missing  (Default: \code{FALSE}). Poor signal was 
#' identified during QC by \code{\link{meffil.qc}()}
#' as signal that failed to pass the 
#' detection p-value threshold (\code{detection.threshold})
#' or bead threshold (\code{bead.threshold}).
#' @param verbose If \code{TRUE}, then status messages are printed during execution (Default: \code{FALSE}).
#' @return List containing normalized methylated and unmethylated signals.
#'
#' @export
meffil.normalize.sample <- function(norm.object, remove.poor.signal=F, verbose=F) {
    stopifnot(is.normalized.object(norm.object))

    ## begin backwards compatibility
    if (any(c("chrX","chrY") %in% names(norm.object$quantiles))) {
        idx <- which(names(norm.object$quantiles) %in% c("chrX","chrY"))
        names(norm.object$quantiles)[idx] <- tolower(names(norm.object$quantiles)[idx])
    }
    if (is.null(norm.object$featureset)) {
        norm.object$chip <- norm.object$featureset <- "450k"
    }
    ## end backwards compatibility
    
    probes <- meffil.probe.info(norm.object$chip, norm.object$featureset)
   
    rg <- read.rg(norm.object$basename, verbose=verbose)
    rg <- background.correct(rg, probes, verbose=verbose)
    rg <- dye.bias.correct(rg, probes, norm.object$reference.intensity, verbose=verbose)
    mu <- rg.to.mu(rg, probes)

    sites <- meffil.get.sites(norm.object$featureset)
    mu$M <- mu$M[sites]
    mu$U <- mu$U[sites]

    msg("Normalizing methylated and unmethylated signals.", verbose=verbose)
    probe.subsets <- get.quantile.site.subsets(norm.object$featureset)
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

    if (remove.poor.signal) {
        bad.probes <- union(
            names(norm.object$bad.probes.beadnum),
            names(norm.object$bad.probes.detectionp))
        for (target in setdiff(names(mu), c("class","version"))) {
            probe.idx <- which(names(mu[[target]]) %in% bad.probes)
            mu[[target]][probe.idx] <- NA
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
