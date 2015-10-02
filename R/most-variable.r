#' Most variable CpG sites
#'
#' Returns the most variable CpG sites (rows) in the methylation matrix.
#'
#' @param x Methylation matrix (rows = CpG sites, columns = samples).
#' @param n Number of CpG sites to return.
#' @param autosomal If \code{TRUE}, then consider only autosomal CpG sites (Default: TRUE).
#' @return The \code{n} CpG site identifiers (rownames of \code{x}) with the greatest variance in \code{x}.
#' 
#' @export
meffil.most.variable.cpgs <- function(x, n=1000, autosomal=T) {
    stopifnot(!is.null(rownames(x)))
    if (autosomal)
        x <- x[which(rownames(x) %in% meffil.get.autosomal.sites()),,drop=F]
    var <- matrixStats::rowVars(x, na.rm=T)
    rownames(x)[order(var, decreasing=TRUE)[1:n]]
}
