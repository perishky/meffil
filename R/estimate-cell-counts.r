#' Estimate cell counts from a reference
#'
#' Estimate cell type ratios from methylation profiles of purified cell populations
#' (Infinium HumanMethylation450 BeadChip).
#'
#' @param object An object created by \code{\link{meffil.compute.normalization.object}()}.
#' @param mu (Optional) Methylated and unmethylated intensities for a sample
#' created by \code{\link{meffil.rg.to.mu}}.  If this is not supplied,
#' then \code{object} is used to create such an object from the corresponding IDAT files.
#' The intensities will have been background corrected and dye bias corrected.
#' @param reference Object describing methylation profiles of purified cell populations.
#' @param probes Probe annotation used to construct the control matrix
#' (Default: \code{\link{meffil.probe.info}()}).
#' @param verbose If \code{TRUE}, then status messages are printed during execution (Default: \code{FALSE}).
#' @return A list:
#' - \code{counts} Cell count estimates.
#' - \code{beta} Normalized methylation levels of sites used to differentiate
#' between reference cell types.
#'
#' Results should be nearly identical to \code{\link[minfi]{estimateCellCounts}()}.
#' 
#' @export
meffil.estimate.cell.counts <- function(object, mu=NULL, reference,  
                                        probes=meffil.probe.info(), verbose=F) {

    stopifnot(is.cell.type.reference(reference))
    
    mu <- meffil.quantile.normalize.sample(object, mu=mu,
                                           subsets=reference$subsets,
                                           quantiles=reference$quantiles,
                                           probes=probes, verbose=verbose)
    
    beta <- meffil.get.beta(mu)
    beta <- beta[rownames(reference$beta)]
    counts <- estimate.cell.counts(beta, reference$beta)

    list(counts=counts, beta=beta)
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
estimate.cell.counts <- function(beta, beta.cell.types) {
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





