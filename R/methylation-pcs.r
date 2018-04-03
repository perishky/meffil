#' Compute principal components of a methylation matrix.
#'
#' @param  normalized.beta Output from \code{\link{meffil.normalize.samples}()},
#' either a matrix or a GDS filename.
#' @param  probe.range Default = 5000. How many probes to be used in calculating PCs
#' @param  autosomal Default = TRUE. If true, remove probes on sex chromosomes.  
#' @param  verbose=T Print progress messages?
#' @return the principal components of \code{normalized.beta}.
#'
#' @export
meffil.methylation.pcs <- function (normalized.beta, probe.range = 5000, autosomal = T, verbose = F) {
    if (is.matrix(normalized.beta))
        sites <- rownames(normalized.beta)
    else {
        if (is.character(normalized.beta)) {
            stopifnot(file.exists(normalized.beta))
            gds.file <- openfn.gds(gds.filename)
            on.exit(closefn.gds(gds.file))
        }
        else
            gds.file <- normalized.beta
        sites <- read.gdsn(index.gdsn(gds.file, "row.names"))
    }
    
    subset <- sites
      
    if (autosomal) {
        featureset <- meffil:::guess.featureset(sites)
        autosomal.sites <- meffil.get.autosomal.sites(featureset)
        subset <- intersect(autosomal.sites, subset)
    }
    
    meffil:::msg("Calculating variances", verbose = verbose)
    if (is.matrix(normalized.beta)) {
        var.sites <- meffil.most.variable.cpgs(normalized.beta[subset,], n = probe.range)
        
    } else {
        cores <- getOption("mc.cores", 1)
        cl <- parallel::makeCluster(cores)
        vars <- clusterApply.gdsn(cl=cl,
                                  gds.fn=gds.filename,
                                  node.name="matrix",
                                  margin=1,
                                  as.is="double",
                                  FUN=function(x) var(x, na.rm=T))
        vars <- vars[match(subset, sites)]
        var.sites <- subset[order(vars, decreasing=T)[1:probe.range]]
    }
    var.idx <- match(var.sites, sites)
    
    meffil:::msg("Calculating beta PCs", verbose = verbose)
    if (is.matrix(normalized.beta)) 
        mat <- normalized.beta[var.idx,]
    else {
        matrix.node <- index.gdsn(gds.file, "matrix")
        mat <- t(sapply(var.idx, function(idx) 
                        read.gdsn(matrix.node, start=c(idx,1), count=c(1,-1))))
    }
        
    prcomp(t(meffil:::impute.matrix(mat, margin = 1)))$x
}

