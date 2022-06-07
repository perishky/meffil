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
meffil.design.matrix <- function(qc.objects, number.pcs, fixed.effects = NULL, random.effects = NULL){
    if (!is.null(fixed.effects) && is.character(fixed.effects))
        fixed.effects <- extract.from.samplesheet(qc.objects, fixed.effects)
    if (!is.null(random.effects) && is.character(random.effects))
        random.effects <- extract.from.samplesheet(qc.objects, random.effects)
    
    pca.ret <- meffil.pcs(qc.objects, fixed.effects, random.effects)
    pca.to.design.matrix(pca.ret, number.pcs, fixed.effects, random.effects)
}

#' Calculate control probe PCs
#'
#' @param qc.objects A list of outputs from \code{\link{meffil.create.qc.object}()}.
#' @param fixed.effects Names of columns in samplesheet that should be included as fixed effects
#' along with control matrix principal components (Default: NULL).
#' @param random.effects Names of columns in samplesheet that should be included as random effects
#' (Default: NULL).
#' @return PCA of control probes
#' 
#' @export
meffil.pcs <- function(qc.objects, fixed.effects=NULL, random.effects=NULL) {
    control.matrix <- meffil.control.matrix(qc.objects, normalize=T,
                                            fixed.effects=fixed.effects,
                                            random.effects=random.effects)
    control.pca <- prcomp(control.matrix)
    return(control.pca)
}



#' Infinium HumanMethylation450 BeadChip control matrix
#'
#' Matrix containing control probe intensities from the Infinium HumanMethylation450 BeadChip.
#'
#' @param qc.objects A list of outputs from \code{\link{meffil.create.qc.object}()}.
#' @param normalize If \code{TRUE}, then control matrix is scaled and specified
#' fixed and random effects removed from the matrix.  Otherwise, the raw control matrix is returned.
#' (Default: \code{FALSE}).
#' @param fixed.effects Names of columns in samplesheet that should be included as fixed effects
#' along with control matrix principal components (Default: NULL).
#' @param random.effects Names of columns in samplesheet that should be included as random effects
#' (Default: NULL).
#' @return Matrix with one row per object consisting of control probe intensities and summaries.
#'
#' @export
meffil.control.matrix <- function(qc.objects, normalize=F,
                                  fixed.effects = NULL, random.effects = NULL) {
    if (!is.null(fixed.effects) && is.character(fixed.effects))
        fixed.effects <- extract.from.samplesheet(qc.objects, fixed.effects)
    if (!is.null(random.effects) && is.character(random.effects))
        random.effects <- extract.from.samplesheet(qc.objects, random.effects)
    
    stopifnot(all(sapply(qc.objects, meffil:::is.qc.object)))
    stopifnot(length(qc.objects) >= 2)
    stopifnot(is.null(fixed.effects) || nrow(fixed.effects) == length(qc.objects))
    stopifnot(is.null(random.effects) || nrow(random.effects) == length(qc.objects))

    control.matrix <- t(sapply(qc.objects, function(object) object$controls))

    if (normalize) {
        control.matrix <- meffil:::impute.matrix(control.matrix, margin = 2)

        ## check for zero-variance control matrix variables
        cv <- colVars(control.matrix,na.rm=T)
        zero.idx <- which(cv < 2e-16)
        
        ## normalize the control matrix
        control.matrix <- scale(control.matrix)

        ## add noise to zero variance columns
        if (length(zero.idx) > 0) {
            warning(paste("Some control matrix variables have zero variance:",
                          paste(colnames(control.matrix)[zero.idx], collapse=", ")))
            if (length(zero.idx) == ncol(control.matrix))
                stop("All control matrix variables have zero variance.")
            for (idx in zero.idx) 
                control.matrix[,idx] <- rnorm(nrow(control.matrix), mean=0, sd=4e-16)

            ## remove 'scale' attributes from control.matrix now that we have
            ## changed the sd of some columns so that prcomp which uses these attributes
            ## to identify zero-variance columns doesn't get confused!
            control.matrix <- control.matrix[,] 
        }

        control.matrix[control.matrix > 3] <- 3
        control.matrix[control.matrix < -3] <- -3
        if (!is.null(fixed.effects) || !is.null(random.effects)) {
            control.matrix <- meffil:::adjust.columns(control.matrix, fixed.effects, random.effects)
            control.matrix <- scale(control.matrix)
        }
    }
    control.matrix
}

# Given a partition of the QC objects into two sets,
# use the design matrix from one partition to predict the design
# matrix for the other partition.
# @param qc.objects A list of outputs from \code{\link{meffil.create.qc.object}()}.
# @param number.pcs Number of principal components to include in the design matrix (Default: all).
# @param new.idx Indices into \code{qc.objects} that identify the objects for which
# the design matrix will be predicted.  The initial design matrix will be derived from
# the other objects.
# @param fixed.effects Names of columns in samplesheet that should be included as fixed effects
# along with control matrix principal components (Default: NULL).
# @param random.effects Names of columns in samplesheet that should be included as random effects
# (Default: NULL).
# @return Design matrix with one column for the \code{new.idx} partition predicted from the
# first \code{number.pcs} prinicipal components of the objects. 
predict.design.matrix <- function (qc.objects, number.pcs, new.idx,
                                   fixed.effects = NULL,
                                   random.effects = NULL) {
    stopifnot(all(sapply(qc.objects, meffil:::is.qc.object)))
    stopifnot(all(new.idx %in% 1:length(qc.objects)))
    
    if (!is.null(fixed.effects) && is.character(fixed.effects))
        fixed.effects <- meffil:::extract.from.samplesheet(qc.objects, fixed.effects)
    if (!is.null(random.effects) && is.character(random.effects))
        random.effects <- meffil:::extract.from.samplesheet(qc.objects, random.effects)
    
    old.pca <- meffil.pcs(qc.objects[-new.idx],
                          fixed.effects[-new.idx,,drop = F],
                          random.effects[-new.idx,, drop = F])
    new.controls <- meffil.control.matrix(qc.objects, normalize = T,
                                          fixed.effects = fixed.effects,
                                          random.effects = random.effects)[new.idx,,drop=F]
    new.pca <- predict(old.pca, newdata = data.frame(new.controls, check.names = F))
    meffil:::pca.to.design.matrix(new.pca, number.pcs,
                                  fixed.effects[new.idx,, drop = F],
                                  random.effects[new.idx,, drop = F])
}




# From PCA applied to the raw control matrix,
# derive the corresponding design matrix using the first \code{number.pcs} principal components.
pca.to.design.matrix <- function(pca.ret, number.pcs, fixed.effects=NULL, random.effects=NULL) {
    if (! "matrix" %in% class(pca.ret))
        pca.ret <- pca.ret$x
    if (missing(number.pcs)) 
        number.pcs <- ncol(pca.ret)
    stopifnot(number.pcs >= 1 && number.pcs <= ncol(pca.ret))
    control.components <- pca.ret[, 1:number.pcs, drop = F]
    design.matrix <- model.matrix(~control.components - 1)
    rownames(design.matrix) <- rownames(pca.ret)
    list(fixed = cbind(design.matrix, fixed.effects), random = random.effects)   
}

