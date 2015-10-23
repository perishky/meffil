#' Estimate cell counts from a reference
#'
#' Estimate cell type ratios from methylation profiles of purified cell populations
#' (Infinium HumanMethylation450 BeadChip).
#'
#' @param object An object created by \code{\link{meffil.create.qc.object}()}.
#' @param verbose If \code{TRUE}, then status messages are printed during execution
#' (Default: \code{FALSE}).
#' @param cell.type.reference Character string name of the cell type reference
#' to use for estimating cell counts. Estimates are not generated if set to NULL (default).
#' See \code{\link{meffil.get.cell.type.references}()} for a list of available
#' references.  New references can be created using
#' \code{\link{meffil.create.cell.type.reference}()}. 
#' @return A list:
#' - \code{counts} Cell count estimates.
#' - \code{beta} Normalized methylation levels of sites used to differentiate
#' - \code{reference} Name of the cell type reference used.
#' between reference cell types.
#'
#' Results should be nearly identical to \code{\link[minfi]{estimateCellCounts}()}.
#' 
#' @export
meffil.estimate.cell.counts <- function(object, cell.type.reference, verbose=T) {
    stopifnot(is.qc.object(object))
    check.cell.type.reference.exists(cell.type.reference)
    
    rg <- read.rg(object$basename, verbose=verbose)
    rg <- background.correct(rg, verbose=verbose)
    rg <- dye.bias.correct(rg, object$dye.intensity, verbose=verbose)
    mu <- rg.to.mu(rg)

    estimate.cell.counts.from.mu(mu, cell.type.reference, verbose)
}

estimate.cell.counts.from.mu <- function(mu, cell.type.reference, verbose=F) {
    check.cell.type.reference.exists(cell.type.reference)

    reference.object <- get.cell.type.reference.object(cell.type.reference)

    mu <- quantile.normalize.signals(mu, reference.object$subsets, reference.object$quantiles, verbose=F)
    beta <- get.beta(mu$M, mu$U)
    beta <- beta[rownames(reference.object$beta)]
    counts <- estimate.cell.counts.from.beta(beta, reference.object$beta)
    
    list(origin="meffil.estimate.cell.counts", counts=counts, beta=beta, reference=cell.type.reference)
}    

is.cell.count.object <- function(object) {
    object$origin == "meffil.estimate.cell.counts"
}

## Input:
## beta[CpG] = methylation level of the CpG in the sample
## beta.cell.types[CpG,cell type] = methylation level of CpG in the cell type
##
## Output:
## counts[cell type] = proportion of cell type in the sample
## that minimizes ( beta - beta.cell.types * counts )^2
## subject to counts >= 0.
##
## Based on code from PMID 22568884.
estimate.cell.counts.from.beta <- function(beta, beta.cell.types) {
    stopifnot(length(beta) == nrow(beta.cell.types))
    
    ## Let r = counts[cell type],
    ##     b = beta[cpg]
    ##     bc = beta.cell.types[cpg,cell type]
    ## for a given sample.
    ## We want to find r >= 0 that minimizes (b - bc r)^2.
    ## The unknown r is quadratic in the expression so
    ## we need a quadratic program solver.
    ## 
    ## We use the function solve.QP(D,d,A,b0).
    ## It finds a vector b1 that minimizes
    ##   b1T D b1 - 2 dT b1 subject to AT b1 >= b0.
    ## We need to put (b - bc r)^2 into this form.
    ##
    ## (b - bc r)^2
    ## = (b - bc r)T (b - bc r)
    ## = bT b - 2 bT bc r - rT bcT bc r
    ## <=> finding r to minimize -2 bT bc r + rT (bcT bc) r
    ## = rT (bcT bc) r - 2 (bcT b)T r
    ## <=> solve.QP(bcT bc, bcT b, I, z) where I = identity matrix, z = 0
    require(quadprog)
    I <- diag(ncol(beta.cell.types))
    zero <- rep(0, ncol(beta.cell.types))
    cpg.idx <- which(!is.na(beta))
    bcT.x.bc <- t(beta.cell.types[cpg.idx,]) %*% beta.cell.types[cpg.idx,]
    bcT.x.b <- t(beta.cell.types[cpg.idx,]) %*% matrix(beta[cpg.idx], nrow=length(cpg.idx))
    counts <- solve.QP(bcT.x.bc, bcT.x.b, I, zero)$sol
    names(counts) <- colnames(beta.cell.types)
    counts
}

#' Estimate cell counts for a beta matrix from a reference
#'
#' Estimate cell type ratios from methylation profiles of purified cell populations
#' (Infinium HumanMethylation450 BeadChip).
#'
#' @param beta Matrix of Illumina 450K methylation levels (rows = CpG sites, columns = subjects).
#' @param verbose If \code{TRUE}, then status messages are printed during execution
#' (Default: \code{FALSE}).
#' @param cell.type.reference Character string name of the cell type reference
#' to use for estimating cell counts. Estimates are not generated if set to NULL (default).
#' See \code{\link{meffil.get.cell.type.references}()} for a list of available
#' references.  New references can be created using
#' \code{\link{meffil.create.cell.type.reference}()}. 
#' @return A matrix of cell count estimates.
#'
#' Results should be nearly identical to \code{\link[minfi]{estimateCellCounts}()}.
#' 
#' @export
meffil.estimate.cell.counts.from.betas <- function(beta, cell.type.reference, verbose=F) {
    stopifnot(is.matrix(beta))
    meffil:::check.cell.type.reference.exists(cell.type.reference)
    
    reference.object <- meffil:::get.cell.type.reference.object(cell.type.reference)
    
    beta <- quantile.normalize.betas(beta, reference.object$subsets, reference.object$quantiles, verbose=verbose)
    
    beta <- beta[intersect(rownames(beta), rownames(reference.object$beta)),]
    reference.beta <- reference.object$beta[rownames(beta),]

    t(apply(beta, 2, meffil:::estimate.cell.counts.from.beta, reference.beta))
}    

