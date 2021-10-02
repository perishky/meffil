#' Most variable CpG sites
#'
#' Returns the most variable CpG sites (rows) in the methylation matrix.
#'
#' @param  beta Output from \code{\link{meffil.normalize.samples}()},
#' either a matrix or a GDS filename.
#' @param n Number of CpG sites to return.
#' @param sites Subset of CpG sites to consider (row names of beta) (Default: NULL).
#' @param samples Subset of samples to consider (column names of beta) (Default: NULL).
#' @param  autosomal If true, remove probes on sex chromosomes (Default: TRUE).
#' @param winsorize.pct Apply to methylation levels
#' winsorized to the given level. Set to NA to avoid winsorizing (Default: NA).
#' @param outlier.iqr.factor Apply to methylation after setting,
#' for each CpG site, values less than
#' \code{Q1 - outlier.iqr.factor * IQR} or more than 
#' \code{Q3 + outlier.iqr.factor * IQR} to NA. Here IQR is the inter-quartile
#' range of the methylation levels at the CpG site, i.e. Q3-Q1.
#' Set to NA to skip this step (Default: NA).
#' @return The \code{n} CpG site identifiers (rownames of \code{x}) with the greatest variance in \code{x}.
#' 
#' @export
meffil.most.variable.cpgs <- function(beta, n=1000, sites=NULL, samples=NULL, autosomal=T, winsorize.pct=NA, outlier.iqr.factor=NA) {
    stopifnot(n > 0)
    if (is.matrix(beta)) {    
        stopifnot(!is.null(rownames(beta)))
        stopifnot(!is.null(colnames(beta)))
        all.sites <- rownames(beta)
        all.samples <- colnames(beta)
    } else {
        beta.dims <- meffil:::retrieve.gds.dims(beta)
        all.sites <- beta.dims[[1]]
        all.samples <- beta.dims[[2]]
    }
    
    if (is.null(sites)) sites <- all.sites
    else sites <- intersect(sites, all.sites)
    stopifnot(length(sites)>0)

    if (autosomal) {
        featureset <- meffil:::guess.featureset(all.sites)
        autosomal.sites <- meffil.get.autosomal.sites(featureset)
        sites <- intersect(autosomal.sites, sites)
    }
    
    if (is.null(samples)) samples <- all.samples
    else samples <- intersect(samples, all.samples)
    stopifnot(length(samples)>0)
    
    if (is.matrix(beta)) {
        beta <- beta[sites,samples,drop=F]
        beta <- meffil.handle.outliers(beta, winsorize.pct, outlier.iqr.factor)
        vars <- matrixStats::rowVars(beta, na.rm=T)
    } else {               
        var.fun <- function(x,winsorize.pct, outlier.iqr.factor) var(x, na.rm=T)
        if (is.numeric(winsorize.pct) || is.numeric(outlier.iqr.factor))
            var.fun <- function(x, winsorize.pct, outlier.iqr.factor) {
                x <- meffil::meffil.handle.outliers(x, winsorize.pct, outlier.iqr.factor)
                var(x, na.rm=T)
            }

        vars <- meffil:::lapply.gds(
            beta,
            margin=1,
            sites=sites,
            samples=samples,
            type="double",
            FUN=var.fun, 
            winsorize.pct=winsorize.pct,
            outlier.iqr.factor=outlier.iqr.factor)
    }
    n <- min(length(sites), n)
    sites[order(vars, decreasing=TRUE)[1:n]]
}


