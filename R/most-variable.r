#' Most variable CpG sites
#'
#' Returns the most variable CpG sites (rows) in the methylation matrix.
#'
#' @param beta Methylation matrix (rows = CpG sites, columns = samples).
#' @param n Number of CpG sites to return.
#' @return The \code{n} CpG site identifiers (rownames of \code{x}) with the greatest variance in \code{x}.
#' 
#' @export
meffil.most.variable.cpgs <- function(beta, n=1000) {
    stopifnot(!is.null(rownames(beta)))
    var <- matrixStats::rowVars(beta, na.rm=T)
    rownames(beta)[order(var, decreasing=TRUE)[1:n]]
}
