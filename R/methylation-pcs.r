#' Compute principal components of a methylation matrix.
#'
#' @param  normalized.beta Output from \code{meffil.normalize.samples}
#' @param  probe.range Default = 5000. How many probes to be used in calculating PCs
#' @param  autosomal Default = TRUE. If true, remove probes on sex chromosomes.  
#' @param  verbose=T Print progress messages?
#' @return the principal components of \code{normalized.beta}.
#'
#' @export
meffil.methylation.pcs <- function(normalized.beta, 
                                   probe.range=5000, autosomal=T, verbose=F) {
    if (autosomal) {
        featureset <- guess.featureset(rownames(normalized.beta))
        autosomal.sites <- meffil.get.autosomal.sites(featureset)
        autosomal.sites <- intersect(autosomal.sites, rownames(normalized.beta))
        normalized.beta <- normalized.beta[autosomal.sites,]
    }
        
    msg("Calculating variances", verbose=verbose)
    
    var.sites <- meffil.most.variable.cpgs(normalized.beta, n=probe.range)
    var.idx <- match(var.sites, rownames(normalized.beta))
    
    msg("Calculating beta PCs", verbose=verbose)
    prcomp(t(meffil:::impute.matrix(normalized.beta[var.idx,], margin=1)))$x
}
