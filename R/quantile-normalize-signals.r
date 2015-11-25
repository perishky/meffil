quantile.normalize.signals <- function(mu, subsets, quantiles, verbose=F) {
    stopifnot(is.mu(mu))
    stopifnot(length(subsets) == length(quantiles))
    stopifnot(all(names(subsets) %in% names(quantiles)))

    for (subset.name in names(subsets)) {
        subset <- subsets[[subset.name]]
        for (target in c("M","U")) {
            data <- matrix(mu[[target]][subset])
            full.quantiles <- quantiles[[subset.name]][[target]]
            full.quantiles <- approx(1:length(full.quantiles), full.quantiles, 1:length(data))$y
            mu[[target]][subset] <- preprocessCore::normalize.quantiles.use.target(data,full.quantiles)
        }
    }
    mu
}
