
## just like mcsapply but saves the matrix output to a GDS file
mcsapply.to.gds <- function (X, FUN, ..., gds.filename, storage="float64",  max.bytes=2^30-1) {
    stopifnot(length(X) > 0)

    gds.file <- createfn.gds(gds.filename)
    on.exit(closefn.gds(gds.file))
    
    cores <- options()$mc.cores
    if (is.null(cores) || cores <= 1 || length(X) == 1) {
        ret <- mclapply(X, FUN, ...)
        if (length(X) > 1)
            ret <- do.call(cbind, ret)
        else
            ret <- matrix(ret[[1]], ncol=1, nrow=length(ret[[1]]))
        colnames(ret) <- names(X)

        create.gds.matrix(gds.file, storage=storage,
                          n.row=nrow(ret), n.col=ncol(ret),
                          row.names=rownames(ret), col.names=colnames(ret), mat=ret)
    } else {    
        first <- mclapply(X[1], FUN, ...)[[1]]
        
        row.names <- names(first)
        col.names <- names(X)
        n.row <- length(first)
        n.col <- length(X)

        matrix.node <- create.gds.matrix(gds.file, storage=storage,
                                         n.row=n.row, n.col=n.col,
                                         row.names=row.names, col.names=col.names)
        append.gds.columns(matrix.node, mat=first, start=1)
 
        X <- X[-1]
        
        ret.bytes <- as.numeric(object.size(first))
                
        max.bytes <- min(max.bytes, ret.bytes * length(X))
        n.fun <- floor(max.bytes/ret.bytes)
        if (n.fun < 1)
            stop(paste("The max.bytes parameter is too small.  Try setting it > ", ret.bytes, ".", sep=""))
        
        n.mclapply <- ceiling(length(X)/n.fun)
               
        partitions <- partition.integer.subsequence(1,length(X),n.mclapply)
        
        for (i in 1:nrow(partitions)) {
            idx <- partitions[i,"start"]:partitions[i,"end"]
            mc.ret <- mclapply(X[idx], FUN, ...)
            if (length(idx) != length(mc.ret) || any(sapply(mc.ret, is.null)))
                stop(paste("The operating system has decided that some forks of mclapply are using too much memory.\n",
                           "Try reducing the max.bytes parameter or the R option 'mc.cores'."))
            is.error <- sapply(mc.ret, class) == "try-error"
            mc.ret <- sapply(1:length(mc.ret), function(i) {
                if (is.error[i])
                    rep(NA, n.row)
                else if (length(mc.ret[[i]]) == n.row) {
                    mc.ret[[i]]
                }
                else {
                    warning("Item ", i, " has the wrong length.")
                    rep(NA, n.row)
                }
            })
            append.gds.columns(matrix.node, mat=mc.ret, start=idx[1]+1)
        }
    }
    return(gds.filename)
}

## create GDS file storing a matrix with row names and column names
create.gds.matrix <- function(gds.file, storage, n.row, n.col, row.names=NULL, col.names=NULL, mat=NULL) {
    if (!is.null(mat)) {
        stopifnot(is.matrix(mat))
        stopifnot(nrow(mat) == n.row)
        stopifnot(ncol(mat) == n.col)
    }

    if (!is.null(row.names)) {
        stopifnot(n.row == length(row.names))
        add.gdsn(gds.file, "row.names", row.names)
    }
    if (!is.null(col.names)) {
        stopifnot(n.col == length(col.names))
        add.gdsn(gds.file, "col.names", col.names)
    }
    add.gdsn(gds.file,
             name="matrix",
             storage=storage,
             valdim=c(n.row, n.col),
             val=mat,
             replace=TRUE)
}

## add the columns in 'mat' to the GDS matrix
## starting at the given column
append.gds.columns <- function(node, mat, start) {
    if (!is.matrix(mat))
        mat <- matrix(mat, ncol=1)
    write.gdsn(node, mat, start=c(1,start), count=c(nrow(mat), ncol(mat)))
}

## load a methylation matrix from a GDS file
retrieve.gds.methylation <- function(gds.filename, sites, samples) {
    retrieve.gds.matrix(gds.filename, sites, samples)
}

## load a matrix from a GDS file
retrieve.gds.matrix <- function(gds.filename, sites, samples) {
    stopifnot(file.exists(gds.filename))
    gds.file <- openfn.gds(gds.filename)
    on.exit(closefn.gds(gds.file))        

    all.sites <- read.gdsn(index.gdsn(gds.file, "row.names"))
    if (is.null(sites)) sites <- all.sites
    else stopifnot(all(sites %in% all.sites))

    all.samples <- read.gdsn(index.gdsn(gds.file, "col.names"))
    if (is.null(samples)) samples <- all.samples
    else stopifnot(all(samples %in% all.samples))

    matrix.node <- index.gdsn(gds.file, "matrix")
    mat <- readex.gdsn(
        matrix.node,
        sel=list(all.sites %in% sites, all.samples %in% samples), 
        simplify="none")
    rownames(mat) <- all.sites[which(all.sites %in% sites)]
    colnames(mat) <- all.samples[which(all.samples %in% samples)]
    mat[sites,samples,drop=F]
}

## load the dimension names from the GDS matrix
retrieve.gds.dims <- function(gds.filename) {
    stopifnot(file.exists(gds.filename))
    gds.file <- openfn.gds(gds.filename)
    on.exit(closefn.gds(gds.file))
    list(read.gdsn(index.gdsn(gds.file, "row.names")),
         read.gdsn(index.gdsn(gds.file, "col.names")))                  
}

retrieve.gds.cpg.sites <- function(gds.filename) {
    retrieve.gds.dims(gds.filename)[[1]]
}

retrieve.gds.samples <- function(gds.filename) {
    retrieve.gds.dims(gds.filename)[[2]]
}

## apply a function to the rows or columns of the GDS matrix
lapply.gds <- function(gds.filename, margin, sites=NULL, samples=NULL, type, FUN, ...) {
    beta.dims <- retrieve.gds.dims(gds.filename)
    all.sites <- beta.dims[[1]]
    all.samples <- beta.dims[[2]]

    if (is.null(sites)) sites <- all.sites
    else sites <- intersect(sites, all.sites)
    stopifnot(length(sites)>0)
    
    if (is.null(samples)) samples <- all.samples
    else samples <- intersect(samples, all.samples)
    stopifnot(length(samples)>0)

    cores <- getOption("mc.cores", 1)
    cl <- parallel::makeCluster(cores)
    
    ret <- clusterApply.gdsn(
        cl=cl,
        gds.fn=gds.filename,
        node.name="matrix",
        margin=margin,
        selection=list(all.sites %in% sites, all.samples %in% samples),
        as.is=type,
        FUN=FUN,
        ...)
    if (margin == 1) {
        names(ret) <- all.sites[which(all.sites %in% sites)]
        ret <- ret[sites]
    }
    else {
        names(ret) <- all.samples[which(all.samples %in% samples)]
        ret <- ret[samples]
    }
    ret
}


#' Retrieve methylation or detection p-value matrix row and column names
#'
#' @param gds.filename Name of GDS file generated by \code{\link{meffil.normalize.samples}()} or \code{\link{meffil.save.detection.pvalues}()}.
#' @return A list of two vectors,
#' the first providing the row names (CpG sites)
#' and the second providing the column names (sample identifiers).
#'
#' @export
meffil.gds.dims <- function(gds.filename) {
    retrieve.gds.dims(gds.filename)
}



#' Return a vector or list of values obtained by applying a function
#' to the margins of a methylation or detection p-value matrix
#' stored in a GDS file.
#'
#' @param gds.filename Name of GDS file generated by \code{\link{meffil.normalize.samples}()}
#' @param bysite If `TRUE`, then apply function to each CpG site (row),
#' otherwise to each sample (column) (Default: TRUE).
#' @param FUN the function to be applied.
#' @param type returned value.
#' @param sites Names of CpG sites to apply to, `NULL` means all sites (Default: NULL).
#' @param samples Names of samples to apply to, `NULL` means all samples (Default: NULL).
#' @param ... 
#' @return 
#'
#' @export
meffil.gds.apply <- function(
    gds.filename,
    bysite=T,
    type=c("list", "none", "integer", "double", "character", "logical", "raw"),
    FUN,
    sites=NULL,
    samples=NULL,
    ...) {
    
    lapply.gds(gds.filename, margin=sign(!bysite)+1, sites=sites, samples=samples,
               type=type, FUN=FUN, ...)
}


#' Retrieve methylation levels from GDS file
#'
#' @param gds.filename Name of GDS file generated by \code{\link{meffil.normalize.samples}()}.
#' @param sites Names of CpG sites to load, if `NULL` then load all (Default: NULL).
#' @param samples Names of samples to load, if `NULL` then load all (Default: NULL).
#' @return Matrix of methylation levels with rows corresponding to CpG sites
#' and columns to samples.  Rows restricted \code{sites} if not \code{NULL},
#' and columns restricted to \code{samples} if not \code{NULL}.
#' 
#'
#' @export
meffil.gds.methylation <- function(gds.filename, sites=NULL, samples=NULL) {
    retrieve.gds.methylation(gds.filename, sites, samples)
}
