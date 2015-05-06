#' Normalize Infinium HumanMethylation450 BeadChips
#'
#' Normalize a set of samples using their normalized quality control objects.
#'
#' @param norm.objects The list or sublist of \code{\link{meffil.normalize.quantiles}()}.
#' @param beta If \code{TRUE} (default), then the function returns
#' the normalized matrix of methylation levels; otherwise, it returns
#' the normalized matrices of methylated and unmethylated signals.
#' @param pseudo Value to add to the denominator to make the methylation
#' estimate more stable when calculating methylation levels (Default: 100).
#' @param verbose If \code{TRUE}, then detailed status messages are printed during execution (Default: \code{FALSE}).
#' @param temp.dir Option to choose alternative temporary directory, otherwise uses result from \code{tempdir()}
#' @param cpglist.remove Optional list of CpGs to exclude from final output
#' @param ... Arguments passed to \code{\link[parallel]{mclapply}()}
#' except for \code{ret.bytes}.
#' @return Matrix of normalized methylation levels if \code{beta} is \code{TRUE};
#' otherwise matrices of normalized methylated and unmethylated signals.
#' Matrices returned have one column per sample and one row per CpG site.
#' Methylation levels are values between 0 and 1
#' equal to methylated signal/(methylated + unmethylated signal + pseudo).
#'
#' @export
meffil.normalize.samples <- function(norm.objects, beta=T, pseudo=100, verbose=F,
                                     temp.dir=tempdir(), cpglist.remove=NULL, ...) {
    if (beta) 
        normalize.betas(norm.objects, pseudo, verbose, temp.dir, cpglist.remove, ...)
    else 
        normalize.signals(norm.objects, pseudo, verbose, temp.dir, cpglist.remove, ...)
}


normalize.betas <- function(norm.objects, pseudo=100, verbose=F,
                            temp.dir=tempdir(), cpglist.remove=NULL, ...) {    
    stopifnot(length(norm.objects) >= 2)
    stopifnot(all(sapply(norm.objects, is.normalized.object)))
    
    basenames <- sapply(norm.objects, function(object) object$basename)
    probe.names <- get.genomic.probes()
    n.sites <- length(probe.names)
    
    example <- rep(NA_real_, n.sites)
    ret.bytes <- as.integer(object.size(example))
    
    ret <- meffil.mclapply(norm.objects, function(object) {
        msg("Normalizing", object$basename, which(basenames == object$basename),
            verbose=verbose)
        ret <- meffil.normalize.sample(object, verbose=verbose)
        ret <- get.beta(ret$M, ret$U, pseudo)
        unname(ret[probe.names])
    }, ret.bytes=ret.bytes, temp.dir=temp.dir, ...)
    ret <- do.call(cbind, ret)
    colnames(ret) <- names(norm.objects)
    rownames(ret) <- probe.names
    if(!is.null(cpglist.remove)) {
        misscpgs <- cpglist.remove[! cpglist.remove %in% rownames(ret)]
        if(length(misscpgs) > 0)
            warning("The following CPGs were not found: ", paste(misscpgs, collapse=", "))
        ret <- ret[! rownames(ret) %in% cpglist.remove, ]
    }
    ret
}

normalize.signals <- function(norm.objects, pseudo=100, verbose=T,
                              temp.dir=tempdir(), cpglist.remove=NULL, ...) {    
    stopifnot(length(norm.objects) >= 2)
    stopifnot(all(sapply(norm.objects, is.normalized.object)))

    basenames <- sapply(norm.objects, function(object) object$basename)
    probe.names <- get.genomic.probes()
    n.sites <- length(probe.names)

    example <- list(rep(NA_real_, n.sites), rep(NA_real_, n.sites))
    ret.bytes <- as.integer(object.size(example))
    
    ret <- meffil.mclapply(norm.objects, function(object) {
        msg("Normalizing", object$basename, which(basenames == object$basename),
            verbose=verbose)
        ret <- meffil.normalize.sample(object, verbose=verbose)
        ret$M <- unname(ret$M[probe.names])
        ret$U <- unname(ret$U[probe.names])
        ret
    }, ret.bytes=ret.bytes, temp.dir=temp.dir, ...)
    names(ret) <- sapply(norm.objects, function(object) object$basename)
    ret <- list(M=sapply(ret, function(x) x$M),
                U=sapply(ret, function(x) x$U))
    rownames(ret$M) <- rownames(ret$U) <- probe.names
    colnames(ret$M) <- colnames(ret$U) <- names(norm.objects)
    if(!is.null(cpglist.remove)) {
        misscpgs <- cpglist.remove[! cpglist.remove %in% rownames(ret$M)]
        if(length(misscpgs) > 0)
            warning("The following CPGs were not found: ", paste(misscpgs, collapse=", "))
        ret$M <- ret$M[! rownames(ret$M) %in% cpglist.remove, ]
        ret$U <- ret$U[! rownames(ret$U) %in% cpglist.remove, ]
    }
    ret
}
