
#' Normalize Infinium HumanMethylation450 BeadChip
#'
#' Normalize sample methylation data using normalized quantiles.
#'
#' @param object A list element from output of \code{\link{meffil.normalize.objects}()}.
#' @param probes Probe annotation used to construct the control matrix
#' (Default: \code{\link{meffil.probe.info}()}).
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
meffil.normalize.sample <- function(object, probes=meffil.probe.info()) {
    stopifnot(is.normalization.object(object))

    probe.names <- unique(na.omit(probes$name))

    rg <- meffil.read.rg(object$basename, probes)
    rg.correct <- meffil.background.correct(rg, probes)
    rg.correct <- meffil.dye.bias.correct(rg.correct, object$reference.intensity, probes)
    mu <- meffil.rg.to.mu(rg.correct, probes)

    mu$M <- mu$M[probe.names]
    mu$U <- mu$U[probe.names]

    msg("Normalizing methylated and unmethylated signals.")
    probe.subsets <- get.quantile.probe.subsets(probes)
    for (name in names(object$norm)) {
        for (target in names(object$norm[[name]])) {
            probe.idx <- which(names(mu[[target]]) %in% probe.subsets[[name]][[target]])
            orig.signal <- mu[[target]][probe.idx]
            norm.target <- compute.quantiles.target(object$norm[[name]][[target]])
            norm.signal <- preprocessCore::normalize.quantiles.use.target(matrix(orig.signal),
                                                                          norm.target)
            mu[[target]][probe.idx] <- norm.signal
        }
    }
    mu
}
