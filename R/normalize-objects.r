##library(MASS) ## for huber()
##library(limma) ## normexp.signal

#' Normalize objects
#'
#' Normalize microarray quantiles using controls extracted (Infinium HumanMethylation450 BeadChip).
#'
#' @param objects A list of outputs from \code{\link{meffil.compute.normalization.object}()}.
#' @param number.pcs Number of control matrix principal components to adjust for (Default: 2).
#' @param verbose If \code{TRUE}, then status messages are printed during execution (Default: \code{FALSE}).
#' @return Same list as input with additional elements added for each sample
#' including normalized quantiles needed for normalizing each sample.
#'
#' @export
meffil.normalize.objects <- function(objects,
                                    number.pcs=2, 
                                    verbose=F) {

    stopifnot(length(objects) >= 2)
    stopifnot(all(sapply(objects, is.normalization.object)))
    stopifnot(number.pcs >= 1)

    msg("selecting dye correction reference", verbose=verbose)
    intensity.R <- sapply(objects, function(object) object$intensity.R)
    intensity.G <- sapply(objects, function(object) object$intensity.G)
    reference.idx <- which.min(abs(intensity.R/intensity.G-1))
    dye.intensity <- (intensity.R + intensity.G)[reference.idx]/2

    sex <- sapply(objects, function(obj) obj$predicted.sex)
    sex.summary <- table(sex)
    has.both.sexes <- length(sex.summary) >= 2 & min(sex.summary) > 1

    male.idx <- which(sex == "M")
    female.idx <- which(sex == "F")

    msg("creating control matrix", verbose=verbose)
    design.matrix <- meffil.design.matrix(objects, number.pcs)
    if (has.both.sexes) {
        design.male <- meffil.design.matrix(objects[male.idx],
                                            min(length(male.idx), number.pcs))
        design.female <- meffil.design.matrix(objects[female.idx],
                                              min(length(female.idx), number.pcs))
    }

    msg("normalizing quantiles", verbose=verbose)
    subset.names <- names(objects[[1]]$quantiles)
    normalized.quantiles <- sapply(subset.names, function(name) {
        sapply(c("M","U"), function(target) {
            msg(name, target, verbose=verbose)
            original <- sapply(objects, function(object) {
                object$quantiles[[name]][[target]] * dye.intensity/object$dye.intensity
            })

            if (name %in% sex.specific.quantile.probe.subsets() && has.both.sexes) {
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

    out <- lapply(1:length(objects), function(i) {
        object <- objects[[i]]
        object$reference.intensity <- dye.intensity

        subset.names <- applicable.quantile.probe.subsets(object$predicted.sex, has.both.sexes)
        object$norm <- sapply(subset.names, function(subset.name) {
            list(M=normalized.quantiles[[subset.name]]$M[,i],
                 U=normalized.quantiles[[subset.name]]$U[,i])
        },simplify=F)

        object
    })
    names(out) <- sapply(out, function(x) x$Sample_Name)
    return(out)
}

#' Infinium HumanMethylation450 BeadChip normalization design matrix
#'
#' Design matrix derived by applying principal components analysis to control probes.
#'
#' @param objects A list of outputs from \code{\link{meffil.compute.normalization.object}()}.
#' @param number.pcs Number of principal components to include in the design matrix (Default: \code{length(objects)}).
#' @return Design matrix with one column for each of the first \code{number.pcs} prinicipal
#' components.
#'
#' @export
meffil.design.matrix <- function(objects, number.pcs) {
    stopifnot(length(objects) >= 2)

    if (missing(number.pcs))
        number.pcs <- length(objects)

    stopifnot(number.pcs >= 1 && number.pcs <= length(objects))

    control.matrix <- meffil.control.matrix(objects)

    control.matrix <- impute.matrix(control.matrix)
    control.matrix <- scale(t(control.matrix))
    control.matrix[control.matrix > 3] <- 3
    control.matrix[control.matrix < -3] <- -3
    control.matrix <- t(scale(control.matrix))
    
    control.components <- prcomp(t(control.matrix))$x[,1:number.pcs,drop=F]
    design.matrix <- model.matrix(~control.components-1)
    rownames(design.matrix) <- colnames(control.matrix)
    design.matrix
}

#' Infinium HumanMethylation450 BeadChip control matrix
#'
#' Matrix containing control probe intensities from the Infinium HumanMethylation450 BeadChip.
#'
#' @param objects A list of outputs from \code{\link{meffil.compute.normalization.object}()}.
#' @return Matrix with one column per object consisting of control probe intensities and summaries.
#'
#' @export
meffil.control.matrix <- function(objects) {
    stopifnot(length(objects) >= 2)
    stopifnot(all(sapply(objects, is.normalization.object)))

    names(objects) <- sapply(objects, function(object) object$basename)
    sapply(objects, function(object) object$controls)
}


impute.matrix <- function(x, FUN=function(x) mean(x, na.rm=T)) {
    idx <- which(is.na(x), arr.ind=T)
    if (length(idx) > 0) {
        na.rows <- unique(idx[,"row"])
        v <- apply(x[na.rows,],1,FUN)
        v[which(is.na(v))] <- FUN(v) ## if any row imputation is NA ...
        x[idx] <- v[match(idx[,"row"],na.rows)]
    }
    x
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

