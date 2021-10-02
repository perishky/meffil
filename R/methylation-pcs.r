#' Compute principal components of a methylation matrix.
#'
#' @param beta Output from \code{\link{meffil.normalize.samples}()},
#' either a matrix or a GDS filename.
#' @param probe.range Default = 50000. How many probes to be used in calculating PCs.
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
#' @param  full.obj Default = FALSE.  If true, then return the full `prcomp` object rather than just the PCs.
#' @param  verbose=T Print progress messages?
#' @return the principal components of \code{normalized.beta}.
#'
#' @export
meffil.methylation.pcs <- function (beta, probe.range = 50000, sites=NULL, samples=NULL, autosomal=T, winsorize.pct=NA, outlier.iqr.factor=NA, full.obj=F, verbose = F) {
    meffil:::msg("Calculating CpG variance", verbose = verbose)
    
    var.sites <- meffil.most.variable.cpgs(beta, probe.range, sites, samples, autosomal, winsorize.pct, outlier.iqr.factor)

    if (is.matrix(beta)) {
        beta <- beta[var.sites,,drop=F]
        if (!is.null(samples))
            beta <- beta[,samples,drop=F]
    } else {
        beta <- meffil:::retrieve.gds.methylation(beta, var.sites, samples) 
    }
    beta <- meffil.handle.outliers(beta, winsorize.pct, outlier.iqr.factor)
    beta <- meffil:::impute.matrix(beta, margin=1)
    
    meffil:::msg("Calculating beta PCs", verbose = verbose)
    ret <- prcomp(t(beta))$x
    if (!full.obj)
        ret <- ret$x
    ret
}

