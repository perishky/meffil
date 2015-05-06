##library(MASS) ## for huber()
##library(limma) ## normexp.signal

#' Normalize microarray quantiles
#'
#' Normalize microarray quantiles using controls extracted (Infinium HumanMethylation450 BeadChip).
#'
#' @param qc.objects A list of outputs from \code{\link{meffil.create.qc.object}()}.
#' @param number.pcs Number of control matrix principal components to adjust for (Default: 2).
#' @param verbose If \code{TRUE}, then status messages are printed during execution (Default: \code{FALSE}).
#' @return Same list as input with additional elements added for each sample
#' including normalized quantiles needed for normalizing each sample.
#'
#' @export
meffil.normalize.quantiles <- function(qc.objects,
                                    number.pcs=2, 
                                    verbose=F) {

    stopifnot(length(qc.objects) >= 2)
    stopifnot(all(sapply(qc.objects, is.qc.object)))
    stopifnot(number.pcs >= 1)

    msg("selecting dye correction reference", verbose=verbose)
    intensity.R <- sapply(qc.objects, function(object) object$intensity.R)
    intensity.G <- sapply(qc.objects, function(object) object$intensity.G)
    reference.idx <- which.min(abs(intensity.R/intensity.G-1))
    dye.intensity <- (intensity.R + intensity.G)[reference.idx]/2

    sex <- sapply(qc.objects, function(obj) obj$predicted.sex)
    sex.summary <- table(sex)
    has.both.sexes <- length(sex.summary) >= 2 & min(sex.summary) > 1

    male.idx <- which(sex == "M")
    female.idx <- which(sex == "F")

    msg("creating control matrix", verbose=verbose)
    design.matrix <- meffil.design.matrix(qc.objects, number.pcs)
    if (has.both.sexes) {
        design.male <- meffil.design.matrix(qc.objects[male.idx])
        design.female <- meffil.design.matrix(qc.objects[female.idx])
    }

    msg("normalizing quantiles", verbose=verbose)
    subset.names <- names(qc.objects[[1]]$quantiles)
    normalized.quantiles <- sapply(subset.names, function(name) {
        sapply(c("M","U"), function(target) {
            msg(name, target, verbose=verbose)
            original <- sapply(qc.objects, function(object) {
                object$quantiles[[name]][[target]] * dye.intensity/object$dye.intensity
            })

            if (is.sex.specific.subset(name) && has.both.sexes) {
                norm.male <- normalize.quantiles(original[,male.idx,drop=F], design.male)
                norm.female <- normalize.quantiles(original[,female.idx,drop=F], design.female)
                norm <- original
                norm[,male.idx] <- norm.male
                norm[,female.idx] <- norm.female
            }
            else
                norm <- normalize.quantiles(original, design.matrix)
            norm
        }, simplify=F)
    }, simplify=F)

    norm.objects <- lapply(1:length(qc.objects), function(i) {
        object <- qc.objects[[i]]
        object$reference.intensity <- dye.intensity
        subset.names <- applicable.quantile.probe.subsets(object$predicted.sex, has.both.sexes)
        object$norm <- sapply(subset.names, function(subset.name) {
            list(M=normalized.quantiles[[subset.name]]$M[,i],
                 U=normalized.quantiles[[subset.name]]$U[,i])
        },simplify=F)

        object
    })
    names(norm.objects) <- sapply(norm.objects, function(object) object$sample.name)
    norm.objects
}

is.normalized.object <- function(object) {
    is.qc.object(object) && "norm" %in% names(object)
}

normalize.quantiles <- function(quantiles, design.matrix) {
    stopifnot(is.matrix(quantiles))
    stopifnot(is.matrix(design.matrix))
    stopifnot(ncol(quantiles) == nrow(design.matrix))

    quantiles[1,] <- 0
    safe.increment <- 1000
    quantiles[nrow(quantiles),] <- quantiles[nrow(quantiles)-1,] + safe.increment
    mean.quantiles <- rowMeans(quantiles)
    fit <- lm.fit(x=design.matrix, y=t(quantiles - mean.quantiles))
    norm.quantiles <- mean.quantiles + t(residuals(fit))
    
    ## make sure that the quantiles in monotonically increasing, not usually a problem but ...
    for (i in 2:nrow(norm.quantiles))
        norm.quantiles[i,] <- apply(norm.quantiles[(i-1):i,], 2, max, na.rm=T)
    norm.quantiles
}
