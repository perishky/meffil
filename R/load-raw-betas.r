#' Load raw beta matrix
#'
#' @param qc.objects A list of outputs from \code{\link{meffil.create.qc.object}()}.
#' @param pseudo Value to add to the denominator to make the methylation
#' estimate more stable when calculating methylation levels (Default: 100).
#' @param just.beta If \code{TRUE}, then return just the methylation levels; otherwise,
#' return the methylated and unmethylated matrices (Default: TRUE).
#' @param background.correct Function for performing background correction.
#' Arguments are the same as \code{\link{meffil.background.correct}()} which is the default.
#' @param dye.bias.correct  Function for performing dye bias correction.
#' Arguments are the same as \code{\link{meffil.dye.bias.correct}()} which is the default.
#' @param verbose If \code{TRUE}, then detailed status messages are printed during execution (Default: \code{FALSE}).
#' @param ... Arguments passed to \code{\link[parallel]{mclapply}()}.
#' @return If \code{just.beta == TRUE}, the matrix of 
#' methylation levels between between 0 and 1
#' equal to methylated signal/(methylated + unmethylated signal + pseudo).
#' Otherwise, a list containing two matrices, the methylated and unmethylated signals.
#' 
#' @export
meffil.load.raw.data <- function(qc.objects,
                                 pseudo=100,
                                 just.beta=T,
                                 max.bytes=2^30-1,
                                 background.correct=meffil.background.correct,
                                 dye.bias.correct=meffil.dye.bias.correct,
                                 verbose=F,
                                 ...) {
    stopifnot(all(sapply(qc.objects, is.qc.object)))

    featuresets <- sapply(qc.objects, function(qc.object) qc.object$featureset)
    featureset <- featuresets[1]

    if (is.list(featuresets)) { ## backwards compatibility
        featureset <- featuresets <- "450k"
    }   

    if (any(featuresets != featureset)) 
        stop("Multiple feature sets were used to create these QC objects:",
             paste(unique(featuresets), collapse=", "))

    sites <- meffil.get.sites(featureset)
    
    intensity.R <- sapply(qc.objects, function(object) object$intensity.R)
    intensity.G <- sapply(qc.objects, function(object) object$intensity.G)
    valid.idx <- which(intensity.R + intensity.G > 200)
    if (length(valid.idx) == 0) {
        valid.idx <- 1:length(intensity.R)
        warning("All of the microarrays have very low intensity")
    }
    reference.idx <- valid.idx[which.min(abs(intensity.R/intensity.G-1)[valid.idx])]
    dye.intensity <- (intensity.R + intensity.G)[reference.idx]/2

    ret <- mcsapply.safe(qc.objects, function(qc.object) {
        
        if (is.null(qc.object$featureset)) { ## backwards compatibility
            qc.object$chip <- qc.object$featureset <- "450k"
        }   
        
        probes <- meffil.probe.info(qc.object$chip)
        rg <- read.rg(qc.object$basename, verbose=verbose)
        rg <- background.correct(rg, probes, verbose=verbose)
        rg <- dye.bias.correct(rg, probes, dye.intensity, verbose=verbose)
        mu <- rg.to.mu(rg, probes)
        
        m.idx <- match(sites, names(mu$M))
        u.idx <- match(sites, names(mu$U))
        if (just.beta)
            get.beta(unname(mu$M[m.idx]), unname(mu$U[u.idx]), pseudo)
        else
            c(unname(mu$M[m.idx]), unname(mu$U[u.idx]))
    }, ..., max.bytes=max.bytes)

    if (!just.beta) {
        ret <- list(M=ret[1:length(sites),],
                    U=ret[(length(sites)+1):nrow(ret),])
        dimnames(ret$M) <- dimnames(ret$U) <- list(sites, names(qc.objects))
    }
    else
        dimnames(ret) <- list(sites, names(qc.objects))
    ret                       
}
