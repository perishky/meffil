quantile.normalize.betas <- function(beta, subsets, quantiles, verbose=F) {
    stopifnot(is.matrix(beta))
    stopifnot(length(subsets) == length(quantiles))
    stopifnot(all(names(subsets) %in% names(quantiles)))
    
    for (subset.name in names(subsets)) {
        subset <- intersect(subsets[[subset.name]], rownames(beta))
        if (length(subset) == 0)
            stop(paste("subset", subset.name, "and beta matrix have no features in common"))
        full.quantiles <- quantiles[[subset.name]]$beta
        full.quantiles <- approx(1:length(full.quantiles), full.quantiles, 1:length(subset))$y
        beta[subset,] <- preprocessCore::normalize.quantiles.use.target(beta[subset,,drop=F],
                                                                        full.quantiles)
    }
    beta
}
