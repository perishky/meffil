#' Infinium HumanMethylation450 BeadChip normalization design matrix
#'
#' Design matrix derived by applying principal components analysis to control probes.
#'
#' @param qc.objects A list of outputs from \code{\link{meffil.create.qc.object}()}.
#' @param number.pcs Number of principal components to include in the design matrix (Default: all).
#' @param fixed.effects Names of columns in samplesheet that should be included as fixed effects
#' along with control matrix principal components (Default: NULL).
#' @param random.effects Names of columns in samplesheet that should be included as random effects
#' (Default: NULL).
#' @return Design matrix with one column for each of the first \code{number.pcs} prinicipal
#' components.
#'
#' @export
meffil.design.matrix <- function(qc.objects, number.pcs, fixed.effects=NULL, random.effects=NULL) {
    if (!is.null(fixed.effects)) {
        fixed.effects <- do.call(rbind, lapply(qc.objects, function(object) {
            object$samplesheet[,unique(fixed.effects),drop=F]
        }))
        for (name in colnames(fixed.effects))
            stopifnot(length(unique(na.omit(fixed.effects[,name]))) > 1)
    }
    if (!is.null(random.effects)) {
        random.effects <- do.call(rbind, lapply(qc.objects, function(object) {
            object$samplesheet[,unique(random.effects),drop=F]
        }))
        for (name in colnames(random.effects))
            stopifnot(length(unique(na.omit(random.effects[,name]))) > 1)
    }
    
    control.pca <- meffil.pcs(qc.objects,fixed.effects,random.effects)
    
    if (missing(number.pcs))
        number.pcs <- ncol(control.pca$x)

    stopifnot(number.pcs >= 1 && number.pcs <= ncol(control.pca$x))

    control.components <- control.pca$x[,1:number.pcs,drop=F]
    design.matrix <- model.matrix(~control.components-1)
    rownames(design.matrix) <- rownames(control.pca$x)
    list(fixed=cbind(design.matrix, fixed.effects),
         random=random.effects)
}


#' Calculate control probe PCs
#'
#' @param qc.objects A list of outputs from \code{\link{meffil.create.qc.object}()}.
#' @export
#' @return PCA of control probes
meffil.pcs <- function(qc.objects, fixed.effects=NULL, random.effects=NULL) {
    stopifnot(length(qc.objects) >= 2)
    stopifnot(is.null(fixed.effects) || nrow(fixed.effects) == length(qc.objects))
    stopifnot(is.null(random.effects) || nrow(random.effects) == length(qc.objects))
    
    control.matrix <- meffil.control.matrix(qc.objects)
    control.matrix <- impute.matrix(control.matrix, margin=2)
    control.matrix <- scale(control.matrix)
    control.matrix[control.matrix > 3] <- 3
    control.matrix[control.matrix < -3] <- -3
    
    if (!is.null(fixed.effects) || !is.null(random.effects)) {
        control.matrix <- adjust.columns(control.matrix, fixed.effects, random.effects)
        control.matrix <- scale(control.matrix)
    }
    control.pca <- prcomp(control.matrix)
    return(control.pca)
}



#' Infinium HumanMethylation450 BeadChip control matrix
#'
#' Matrix containing control probe intensities from the Infinium HumanMethylation450 BeadChip.
#'
#' @param qc.objects A list of outputs from \code{\link{meffil.create.qc.object}()}.
#' @return Matrix with one row per object consisting of control probe intensities and summaries.
#'
#' @export
meffil.control.matrix <- function(qc.objects) {
    stopifnot(length(qc.objects) >= 2)
    stopifnot(all(sapply(qc.objects, is.qc.object)))

    t(sapply(qc.objects, function(object) object$controls))
}

impute.matrix <- function(x, margin=1, fun=function(x) mean(x, na.rm=T)) {
    if (margin == 2) x <- t(x)
    
    idx <- which(is.na(x) | !is.finite(x), arr.ind=T)
    if (length(idx) > 0) {
        na.idx <- unique(idx[,"row"])
        v <- apply(x[na.idx,],margin,fun) ## v = summary for each row
        v[which(is.na(v))] <- fun(v)      ## if v[i] is NA, v[i] = fun(v)
        x[idx] <- v[match(idx[,"row"],na.idx)] ##
    }

    if (margin == 2) x <- t(x)
    x
}
