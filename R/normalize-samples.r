#' Normalize Infinium HumanMethylation450 BeadChips
#'
#' Normalize a set of samples using their normalized quality control objects.
#'
#' @param norm.objects The list or sublist of \code{\link{meffil.normalize.quantiles}()}.
#' @param pseudo Value to add to the denominator to make the methylation
#' estimate more stable when calculating methylation levels (Default: 100).
#' @param cpglist.remove Optional list of CpGs to exclude from final output
#' @param just.beta Return a normalized methylation levels if TRUE (Default); otherwise return
#' normalized methylated and unmethylated signals.
#' @param filename Optional filename if the returned \code{\link[bigmemory]{big.matrix}} object(s) should
#' be file-backed rather than memory-backed (Default: NULL).  
#' @param verbose If \code{TRUE}, then detailed status messages are printed during execution (Default: \code{FALSE}).
#' @param ... Arguments passed to \code{\link[parallel]{mclapply}()}.
#' @return Either a single (\code{\link[bigmemory]{big.matrix}}) object of normalized methylation levels
#' between 0 and 1 if \code{just.beta == TRUE}; otherwise, a list of two (\code{\link[bigmemory]{big.matrix}})
#' objects "M" and "U" containing normalized methylated and unmethylated signals.
#'
#' Any matrices returned have one column per sample and one row per CpG site.
#'
#' Methylation levels are values between 0 and 1
#' equal to methylated signal/(methylated + unmethylated signal + pseudo).
#'
#' The returned \code{\link[bigmemory]{big.matrix}} object can be converted to an R matrix
#' using the sub-matrix operator, e.g. \code{beta[]}.
#' 
#' If \code{filename} is not \code{NULL}, then two files are created for each returned matrix
#' with names appended to \code{filename}: 'beta.bin' and 'beta.desc'
#' if \code{just.beta == T}, and 'M.bin', 'M.desc', 'U.bin' and 'U.desc' otherwise.
#' These matrices can be reloaded into other R sessions
#' using \code{\link[bigmemory]{attach.big.matrix}()} with the '.desc' filename as an argument
#' (e.g. \code{attach.big.matrix("file.desc", backingpath="path-to-files")}.
#'
#' @export
meffil.normalize.samples <- function(norm.objects, 
                                     pseudo=100, 
                                     cpglist.remove=NULL,
                                     just.beta=T,
                                     filename=NULL,
                                     verbose=F,
                                     ...) {
    stopifnot(length(norm.objects) >= 2)
    stopifnot(all(sapply(norm.objects, is.normalized.object)))
    
    sites <- meffil.get.sites()
    if(!is.null(cpglist.remove)) 
        sites <- setdiff(sites, cpglist.remove)

    if (just.beta) {
        if (!is.null(filename)) filename <- paste0(filename, "beta")
        beta <- create.big.matrix(norm.objects, sites, "double", filename=filename)
    }
    else {
        m.filename <- u.filename <- NULL
        if (!is.null(filename)) m.filename <- paste0(filename, "M")
        if (!is.null(filename)) u.filename <- paste0(filename, "U")
        M <- create.big.matrix(norm.objects, sites, "integer", filename=m.filename)
        U <- create.big.matrix(norm.objects, sites, "integer", filename=u.filename)
    }
        
    status <- mclapply(1:length(norm.objects), function(i) {
        msg("Normalizing", names(norm.objects)[i],
            i, "of", length(norm.objects), verbose=verbose)
        ret <- meffil.normalize.sample(norm.objects[[i]], verbose=verbose)

        if (just.beta)
            beta[,i] <- get.beta(ret$M[rownames(beta)], ret$U[rownames(beta)], pseudo)
        else {
            M[,i] <- ret$M[rownames(M)]
            U[,i] <- ret$U[rownames(U)]
        }
        TRUE
    }, ...)
    
    errors <- check.for.normalization.errors(status)
    if (!is.null(errors)) return(errors)
    
    gc()
    if (just.beta)
        beta
    else
        list(M=M,U=U)
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
#' Fortunately, this is a solution if the object is a \code{\link{bigmemory::big.matrix}},
#' and we use it here.


create.big.matrix <- function(norm.objects, sites, type, filename=NULL) {
    if (is.null(filename)) 
        bigmemory::big.matrix(init=NA,
                              ncol=length(norm.objects),
                              nrow=length(sites),
                              dimnames=list(sites,names(norm.objects)),
                              type=type)
    else
        filebacked.big.matrix(init=NA, type=type,
                              nrow=length(sites), ncol=length(norm.objects),
                              dimnames=list(sites,names(norm.objects)),
                              backingfile=paste(basename(filename), "bin", sep="."),
                              descriptorfile=paste(basename(filename), "desc", sep="."),
                              backingpath=dirname(filename),
                              binarydescriptor=T)

}

check.for.normalization.errors <- function(ret) {
    is.error <- sapply(ret, class) == "try-error"
    if (any(is.error)) {
        warning("Errors were generated by meffil.normalize.sample(), normalization failed.")
        ret[which(is.error)]
    }
    else
        NULL
}

