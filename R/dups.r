#' Collapse duplicate probes
#'
#' Collapse duplicated probes by replacing them with a summary.
#'
#' @param beta Methylation matrix returned by
#' \code{\link{meffil.normalize.samples}()}.
#' @param dup.fun Function to collapse duplicate probes
#' (Default: median).
#' @return The input matrix with duplicated probes
#' (i.e. row names identical after stripping everything
#' after the "_" character) replaced by summaries defined by
#' \code{dup.fun}.
#'
#' @export
meffil.collapse.dups <- function(beta, dup.fun=function(x) median(x,na.rm=T)) {
    stopifnot(is.matrix(beta))
    dups <- identify.dups(rownames(beta))
    if (length(dups) > 0)
        collapse.dups(beta, dups, dup.fun)
    else
        beta
}

identify.dups <- function(sites) {
    is.dup <- grepl("_", sites)
    sites[is.dup] <- sub("_.*", "", sites[is.dup])
    is.dup <- is.dup | sites %in% sites[is.dup]
    if (sum(is.dup) > 0)
        split(which(is.dup), sites[is.dup])
    else
        NULL
}

collapse.dups <- function(beta, dups, dup.fun=function(x) median(x,na.rm=T)) {
    if (is.vector(beta)) {
        beta.nodup <- beta[-unlist(dups)]
        beta.undup <- sapply(dups, function(idx) dup.fun(beta[idx]))
        c(beta.nodup, beta.undup)
    }
    else {
        stopifnot(is.matrix(beta))
        colFUN <- function(beta) apply(beta,2,dup.fun)
        beta.nodup <- beta[-unlist(dups),,drop=F]
        beta.undup <- sapply(dups, function(idx) colFUN(beta[idx,,drop=F]))
        rbind(beta.nodup, t(beta.undup))
    }
}
