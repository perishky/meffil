#' Infinium HumanMethylation450 BeadChip normalization design matrix
#'
#' Design matrix derived by applying principal components analysis to control probes.
#'
#' @param qc.objects A list of outputs from \code{\link{meffil.create.qc.object}()}.
#' @param number.pcs Number of principal components to include in the design matrix (Default: all).
#' @return Design matrix with one column for each of the first \code{number.pcs} prinicipal
#' components.
#'
#' @export
meffil.design.matrix <- function(qc.objects, number.pcs) {
    control.pca <- meffil.pcs(qc.objects)
    
    if (missing(number.pcs))
        number.pcs <- ncol(control.pca$x)

    stopifnot(number.pcs >= 1 && number.pcs <= ncol(control.pca$x))

    control.components <- control.pca$x[,1:number.pcs,drop=F]
    design.matrix <- model.matrix(~control.components-1)
    rownames(design.matrix) <- rownames(control.pca$x)
    design.matrix
}


#' Calculate control probe PCs
#'
#' @param qc.objects A list of outputs from \code{\link{meffil.create.qc.object}()}.
#' @export
#' @return PCA of control probes
meffil.pcs <- function(qc.objects)
{
    stopifnot(length(qc.objects) >= 2)
    control.matrix <- meffil.control.matrix(qc.objects)
    control.matrix <- impute.matrix(control.matrix)
    control.matrix <- scale(t(control.matrix))
    control.matrix[control.matrix > 3] <- 3
    control.matrix[control.matrix < -3] <- -3
    control.matrix <- t(scale(control.matrix))    
    control.pca <- prcomp(t(control.matrix))
    return(control.pca)
}



#' Infinium HumanMethylation450 BeadChip control matrix
#'
#' Matrix containing control probe intensities from the Infinium HumanMethylation450 BeadChip.
#'
#' @param qc.objects A list of outputs from \code{\link{meffil.create.qc.object}()}.
#' @return Matrix with one column per object consisting of control probe intensities and summaries.
#'
#' @export
meffil.control.matrix <- function(qc.objects) {
    stopifnot(length(qc.objects) >= 2)
    stopifnot(all(sapply(qc.objects, is.qc.object)))

    sapply(qc.objects, function(object) object$controls)
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
