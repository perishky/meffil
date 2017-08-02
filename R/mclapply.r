#' mclapply without the 2Gb output memory bound.
#'
#' @param X Same as \code{\link[parallel]{mclapply}()}:
#' A vector (atomic or list) or an expressions vector.
#' Other objects (including classed objects) will be coerced by
#' \code{\link{as.list}()}.
#' @param FUN Same as \code{\link[parallel]{mclapply}()}:
#' The function to be applied to each element of 'X'.
#' @param ... Optional arguments to \code{FUN} and \code{\link[parallel]{mclapply}()}.
#' @param max.bytes The size in memory of the largest object that can
#' be returned by \code{\link[parallel]{mclapply}} (Default: 2Gb-1).
#' @return Same as \code{\link[parallel]{mclapply}()},
#' a list of the same length as 'X' and named by 'X'.
#' Element i is equal to \code{FUN(X[[i]])}.
#'
#' mclapply() has another problem besides an output memory limit.
#' If some magical process in linux called OOM (out of memory) decides
#' that a fork is using too much memory, then it simply kills it without
#' any warning or message.  In such cases, mclapply() will simply return NULL.
#' 
#' http://stackoverflow.com/questions/20674538/mclapply-returns-null-randomly
#'
#' In these functions, the output is tested for NULL return values.
#' They are assumed to be due to memory errors so the FUN argument
#' should not return NULL.
#' 
mclapply.safe <- function (X, FUN, ..., max.bytes=2^30-1) {
    stopifnot(length(X) > 0)

    cores <- options()$mc.cores
    if (is.null(cores) || cores == 1) return(mclapply(X, FUN, ...))
    
    first <- mclapply(X[1], FUN, ...)
    ret.bytes <- as.numeric(object.size(first))
    
    if (length(X) == 1) return(first)

    X <- X[-1]
    
    max.bytes <- min(max.bytes, ret.bytes * length(X))
    n.fun <- floor(max.bytes/ret.bytes)
    if (n.fun < 1)
        stop(paste("The max.bytes parameter is too small.  Try setting it > ", ret.bytes, ".", sep=""))
    
    n.mclapply <- ceiling(length(X)/n.fun)

    partitions <- partition.integer.subsequence(1,length(X),n.mclapply)
    c(first, do.call(c, lapply(1:nrow(partitions), function(i) {
        idx <- partitions[i,"start"]:partitions[i,"end"]
        ret <- mclapply(X[idx], FUN, ...)
        if (length(idx) != length(ret) || any(sapply(ret, is.null)))
            stop(paste("The operating system has decided that some forks of mclapply are using too much memory.\n",
                       "Try reducing the max.bytes parameter or the R option 'mc.cores'."))
        ret
    })))
}

#' yes, could be done with mclapply.safe but then
#' when the list of vectors is merged into a matrix
#' the same data would be using twice as much memory.
mcsapply.safe <- function (X, FUN, ..., max.bytes=2^30-1) {
    stopifnot(length(X) > 0)

    cores <- options()$mc.cores
    if (is.null(cores) || cores == 1) {
        ret <- mclapply(X, FUN, ...)
        ret <- do.call(cbind, ret)
        colnames(ret) <- names(X)
        return(ret)
    }
    
    first <- mclapply(X[1], FUN, ...)
    ret.bytes <- as.numeric(object.size(first))

    if (length(X) == 1) return(first[[1]])

    X <- X[-1]
    
    max.bytes <- min(max.bytes, ret.bytes * length(X))
    n.fun <- floor(max.bytes/ret.bytes)
    if (n.fun < 1)
        stop(paste("The max.bytes parameter is too small.  Try setting it > ", ret.bytes, ".", sep=""))

    n.mclapply <- ceiling(length(X)/n.fun)

    ret <- matrix(NA, ncol=length(X)+1, nrow=length(first[[1]]))
    rownames(ret) <- names(first[[1]])
    colnames(ret) <- c(names(first), names(X))
    ret[,1] <- first[[1]]
   
    partitions <- partition.integer.subsequence(1,length(X),n.mclapply)
    
    lapply(1:nrow(partitions), function(i) {
        idx <- partitions[i,"start"]:partitions[i,"end"]
        mc.ret <- mclapply(X[idx], FUN, ...)
        if (length(idx) != length(mc.ret) || any(sapply(mc.ret, is.null)))
            stop(paste("The operating system has decided that some forks of mclapply are using too much memory.\n",
                       "Try reducing the max.bytes parameter or the R option 'mc.cores'."))
        is.error <- sapply(mc.ret, class) == "try-error"
        idx <- idx[which(!is.error)]
        if (length(idx) > 0)
            ret[,idx+1] <<- do.call(cbind, mc.ret)
        TRUE
    })
    ret
}


partition.integer.subsequence <- function(start, end, n) {
    stopifnot(start <= end)
    stopifnot(n <= end-start+1)
    partitions <- floor(seq(start,end+1,length.out=n+1))
    cbind(start=head(partitions, n=-1),
          end=tail(partitions, n=-1) - 1)
}
