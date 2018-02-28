#' Load detection p-value matrix
#'
#' @param qc.objects A list of outputs from \code{\link{meffil.create.qc.object}()}.
#' @param verbose If \code{TRUE}, then detailed status messages are printed during execution (Default: \code{FALSE}).
#' @param ... Arguments passed to \code{\link[parallel]{mclapply}()}.
#' @return Matrix of probe detection p-values.
#' 
#' @export
meffil.load.detection.pvalues <- function(qc.objects,
                                        max.bytes=2^30-1,
                                        verbose=F,
                                        ...) {
    stopifnot(all(sapply(qc.objects, is.qc.object)))

    featuresets <- sapply(qc.objects, function(qc.object) qc.object$featureset)
    featureset <- featuresets[1]

    if (is.list(featuresets)) ## backwards compatibility
        featureset <- featuresets <- "450k"

    if (any(featuresets != featureset)) 
        stop("Multiple feature sets were used to create these QC objects:",
             paste(unique(featuresets), collapse=", "))

    feature.names <- meffil.get.features(featureset)$name
    
    if (!all(sapply(qc.objects, function(qc.object) exists.rg(qc.object$basename))))
         stop("IDAT files are not accessible for all QC objects")
                          
    ret <- mcsapply.safe(qc.objects, function(qc.object) {
        if (is.null(qc.object$featureset)) ## backwards compatibility
            qc.object$chip <- "450k" 

        rg <- read.rg(qc.object$basename, verbose=verbose)        
        probes <- meffil.probe.info(qc.object$chip)        
        pvalues <- extract.detection.pvalues(rg, probes, verbose=verbose)
        unname(pvalues[feature.names])
    }, ..., max.bytes=max.bytes)

    dimnames(ret) <- list(feature.names, names(qc.objects))
    ret
}
