#' Normalize Infinium HumanMethylation450 BeadChips
#'
#' Normalize a set of samples using their normalized quality control objects.
#'
#' @param norm.objects The list or sublist of \code{\link{meffil.normalize.quantiles}()}.
#' @param pseudo Value to add to the denominator to make the methylation
#' estimate more stable when calculating methylation levels (Default: 100).
#' @param just.beta If \code{TRUE}, then return just the normalized methylation levels; otherwise,
#' return the normalized methylated and unmethylated matrices (Default: TRUE).
#' @param cpglist.remove Optional list of CpGs to exclude from final output
#' @param verbose If \code{TRUE}, then detailed status messages are printed during execution (Default: \code{FALSE}).
#' @param gds.filename If not \code{NULL} (default), then saves the output to a GDS (Genomic Data Structure).
#' This is for cases where the output is too large to fit into main memory.
#' The GDS option assumes that argument \code{just.beta == TRUE}.
#' @param ... Arguments passed to \code{\link[parallel]{mclapply}()}.
#' @return If \code{just.beta == TRUE}, the normalized matrix of 
#' methylation levels between between 0 and 1
#' equal to methylated signal/(methylated + unmethylated signal + pseudo).
#' Otherwise, a list containing two matrices, the normalized methylated and unmethylated signals.
#' If \code{gds.filename} is not \code{NULL}, then the output is saved to the GDS file
#' rather than retained in memory and returned to the caller.
#' The library 'gdsfmt' must be installed.
#' 
#' @export
meffil.normalize.samples <- function(norm.objects, 
                                     pseudo=100,
                                     just.beta=T,
                                     cpglist.remove=NULL,
                                     max.bytes=2^30-1, ## maximum number of bytes for mclapply
                                     gds.filename=NULL,
                                     verbose=F,
                                     ...) {
    stopifnot(length(norm.objects) >= 2)
    stopifnot(all(sapply(norm.objects, is.normalized.object)))
    stopifnot(is.null(gds.filename) || just.beta)

    if (length(unique(sapply(norm.objects, function(object) object$featureset))) > 1)
        stop(paste("Heterogeneous microarray formats included without a common featureset.",
                   "Need to set the 'featureset' argument when creating QC objects."))

    featureset <- norm.objects[[1]]$featureset
    
    if (is.null(featureset)) { ## backwards compatibility
        featureset <- "450k"
    }   

    sites <- meffil.get.sites(featureset)
    if(!is.null(cpglist.remove))
        sites <- setdiff(sites, cpglist.remove)

    if (is.null(gds.filename)) {
        ret <- mcsapply.safe(
            norm.objects,
            FUN=function(object) {
                mu <- meffil.normalize.sample(object, verbose=verbose)
                m.idx <- match(sites, names(mu$M))
                u.idx <- match(sites, names(mu$U))
                if (just.beta)
                    get.beta(unname(mu$M[m.idx]), unname(mu$U[u.idx]), pseudo)
                else
                    c(unname(mu$M[u.idx]), unname(mu$U[u.idx]))
            },
            ...,
            max.bytes=max.bytes)

        if (!just.beta) {
            ret <- list(M=ret[1:length(sites),],
                        U=ret[(length(sites)+1):nrow(ret),])
            dimnames(ret$M) <- dimnames(ret$U) <- list(sites, names(norm.objects))
        }
        else
            dimnames(ret) <- list(sites, names(norm.objects))
        ret
    } else {
        require(gdsfmt)
        mcsapply.gds(
            norm.objects,
            FUN=function(object) {
                mu <- meffil.normalize.sample(object, verbose=verbose)
                m.idx <- match(sites, names(mu$M))
                u.idx <- match(sites, names(mu$U))
                ret <- get.beta(unname(mu$M[m.idx]), unname(mu$U[u.idx]), pseudo)
                names(ret) <- sites
                ret
            },
            ...,
            gds.filename=gds.filename,
            storage="float64",
            max.bytes=max.bytes)
        gds.filename
    }
}

#' We use \code{\link[parallel]{mclapply}()} to reduce running time by taking advantage of the fact
#' that each sample can be normalized independently of the others.
#' Unfortuantely \code{\link[parallel]{mclapply}()} has two important limitations.
#' The size in memory of the list returned may be at most 2Gb otherwise
#' \code{\link[parallel]{mclapply}()} fails with the following error message:
#'    Error in sendMaster(try(lapply(X = S, FUN = FUN, ...), silent = TRUE)) :
#'    long vectors not supported ...
#' 
#' A non-elegant solution to this problem is to guess the size of each element
#' in the returned list and then apply \code{\link[parallel]{mclapply}()} sequentially to a sequence
#' appropriately sized input subsets.
#' A solution for \code{lapply} is to allocate the final object (e.g. a matrix)
#' prior to calling \code{lapply} and then populate the object during
#' the call to \code{lapply} using the global assignment operator '<<-'.
#' Unfortunately this is not a solution for \code{\link[parallel]{mclapply}()}
#' because \code{\link[parallel]{mclapply}()} immediately
#' duplicates the object, applies any modifications to the duplicate
#' and then deletes it prior to completion losing all modifications.
#' I'm not sure why the duplicate is not copied onto the original.
#' This is a solution if the object is a \code{\link{bigmemory::big.matrix}}.
#' However, we tried this but encountered random errors.
#' Sometimes the function completed without incident but other times,
#' with the same data, entire columns of the output matrix would be NA,
#' implying that the meffil.normalize.sample() function failed.
#' However, no errors were generated (tested with a tryCatch).
#' It seems that mclapply and big.matrix do not play well together all the time.
#' We have replaced this with the less elegant approach implemented in mcsapply().

