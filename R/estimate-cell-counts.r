#' Estimate cell counts from a reference
#'
#' Estimate cell type ratios from methylation profiles of purified cell populations
#' (Infinium HumanMethylation450 BeadChip) using the the Houseman algorithm
#' (PMID 22568884).
#'
#' @param object An object created by \code{\link{meffil.create.qc.object}()}.
#' @param verbose If \code{TRUE}, then status messages are printed during execution
#' (Default: \code{FALSE}).
#' @param cell.type.reference Character string name of the cell type reference
#' to use for estimating cell counts. 
#' See \code{\link{meffil.list.cell.type.references}()} for a list of available
#' references.  New references can be created using
#' \code{\link{meffil.add.cell.type.reference}()}. 
#' @return A list:
#' - \code{counts} Cell count estimates.
#' - \code{beta} Normalized methylation levels of sites used to differentiate
#' - \code{reference} Name of the cell type reference used.
#' between reference cell types.
#'
#' Results should be nearly identical to \code{\link[minfi]{estimateCellCounts}()}.
#' 
#' @export
meffil.estimate.cell.counts <- function(qc.object, cell.type.reference, verbose=T) {
    stopifnot(is.qc.object(qc.object))
    stopifnot(is.character(cell.type.reference) && cell.type.reference %in% meffil.list.cell.type.references())

    reference.object <- get.cell.type.reference(cell.type.reference)

    if (is.null(qc.object$featureset)) { ## backwards compatibility
        qc.object$chip <- qc.object$featureset <- "450k"
    }   
    
    rg <- read.rg(qc.object$basename, verbose=verbose)
    probes <- meffil.probe.info(qc.object$chip, reference.object$featureset)
    rg <- background.correct(rg, probes, verbose=verbose)
    rg <- dye.bias.correct(rg, probes, qc.object$dye.intensity, verbose=verbose)
    mu <- rg.to.mu(rg, probes)

    estimate.cell.counts.from.mu(mu, cell.type.reference, verbose)
}

estimate.cell.counts.from.mu <- function(mu, cell.type.reference, verbose=F) {
    stopifnot(is.character(cell.type.reference) && cell.type.reference %in% meffil.list.cell.type.references())

    reference.object <- get.cell.type.reference(cell.type.reference)

    mu <- quantile.normalize.signals(mu, reference.object$subsets, reference.object$quantiles, verbose=F)
    beta <- get.beta(mu$M, mu$U)
    beta <- beta[rownames(reference.object$beta)]
    counts <- estimate.cell.counts.from.beta(beta, reference.object$beta)
    
    list(class="cell.counts",
         version=packageVersion("meffil"),
         counts=counts,
         beta=beta,
         reference=cell.type.reference)
}    

is.cell.count.object <- function(object)
    is.list(object) &&
    ("class" %in% names(object) && object$class == "cell.counts"
     || "origin" %in% names(object) && object$origin == "meffil.estimate.cell.counts") ## backwards compatibility

quantile.normalize.signals <- function(mu, subsets, quantiles, verbose=F) {
    stopifnot(is.mu(mu))
    stopifnot(length(subsets) == length(quantiles))
    stopifnot(all(names(subsets) %in% names(quantiles)))

    for (subset.name in names(subsets)) {
        subset <- subsets[[subset.name]]
        for (target in c("M","U")) {
            data <- matrix(mu[[target]][subset])
            full.quantiles <- quantiles[[subset.name]][[target]]
            full.quantiles <- approx(1:length(full.quantiles), full.quantiles, 1:length(data))$y
            mu[[target]][subset] <- preprocessCore::normalize.quantiles.use.target(data,full.quantiles)
        }
    }
    mu
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

#' Estimate cell counts for a methylation matrix from a reference
#'
#' @param beta Matrix of methylation levels (rows = CpG sites, columns = subjects).
#' @param verbose If \code{TRUE}, then status messages are printed during execution
#' (Default: \code{FALSE}).
#' @param cell.type.reference Character string name of the cell type reference
#' to use for estimating cell counts. 
#' See \code{\link{meffil.list.cell.type.references}()} for a list of available
#' references.  New references can be created using
#' \code{\link{meffil.add.cell.type.reference}()}. 
#' @return A matrix of cell count estimates.
#'
#' Results should be nearly identical to \code{\link[minfi]{estimateCellCounts}()}.
#' 
#' @export
meffil.estimate.cell.counts.from.betas <- function(beta, cell.type.reference, verbose=F) {
    reference.object <- get.cell.type.reference(cell.type.reference)

    stopifnot(is.matrix(beta))

    beta <- quantile.normalize.betas(beta, reference.object$subsets, reference.object$quantiles, verbose=verbose)

    cpg.sites <- intersect(rownames(beta), rownames(reference.object$beta))
    stopifnot(length(cpg.sites) > 0)
    
    beta <- beta[cpg.sites,]
    reference.beta <- reference.object$beta[cpg.sites,]

    t(apply(beta, 2, estimate.cell.counts.from.beta, reference.beta))
}    

quantile.normalize.betas <- function(beta, subsets, quantiles, verbose=F) {
    stopifnot(is.matrix(beta))
    stopifnot(length(subsets) == length(quantiles))
    stopifnot(all(names(subsets) %in% names(quantiles)))
    
    for (subset.name in names(subsets)) {
        subset <- intersect(subsets[[subset.name]], rownames(beta))
        if (length(subset) == 0)
            stop(paste("subset", subset.name, "and beta matrix have no features in common"))
        full.quantiles <- quantiles[[subset.name]]$beta
        full.quantiles <- approx(1:length(full.quantiles), full.quantiles, 1:length(subset))$y
        beta[subset,] <- preprocessCore::normalize.quantiles.use.target(beta[subset,,drop=F],
                                                                        full.quantiles)
    }
    beta
}
