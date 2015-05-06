#' Create a cell type reference object
#'
#' Create a cell type reference object for estimating cell counts
#' with the Infinium HumanMethylation450 BeadChip.
#'
#' @param name Character string providing the name used to reference to the reference.
#' @param M Matrix of methylated probe intensities (rows=CpG sites, columns=samples).
#' @param U Matrix of unmethylatd probe intensities (rows=CpG sites, columns=samples).
#' @param cell.types Vector of cell type names corresponding to sample \code{basename}s.
#' @param number.quantiles Length of numeric sequence to specify probe intensity distributions
#' (Default: 500).
#' @param number.sites Number of probes to characterise cell type methylation (Default: 50).
#' For each cell type, this number of probes with greater methylation than other cell types
#' and the same number with lesser methylation than the other cell types will be included.
#' @param verbose If \code{TRUE}, then status messages are printed during execution
#' (Default: \code{FALSE}).
#' @return A list specifying a cell type reference object that can be used by
#' \code{\link{meffil.estimate.cell.counts}()} can use to estimate cell counts.
#' The object is a list containing:
#' - \code{beta} The normalized methylation values of sites
#' differentially methylated between cell types.
#' - \code{quantiles} The average quantiles of methylated and unmethylated signals
#' of probe sets defined by \code{subsets} (see below).
#' e.g. \code{quantiles[[name]]$M} provides the quantiles (\code{number.quantiles} quantiles)
#' of the probes specified by \code{subsets[[name]]}.
#' - \code{subsets} Probes on the microarray partitioned by relationship to CpG islands,
#' either in an island, in a shore or far from an island.
#'
#' @export
meffil.create.cell.type.reference <- function(name, M, U, cell.types,
                                              number.sites=50,
                                              number.quantiles=500,
                                              subsets=NULL,
                                              verbose=F) {
    stopifnot(!cell.type.reference.exists(name))

    object <- create.cell.type.reference(M,U,cell.types,number.sites,number.quantiles,
                                         subsets,verbose)
    object$name <- name
    assign(name, object, reference.globals)
    object
}

create.cell.type.reference <- function(M, U, cell.types,
                                       number.sites=50,
                                       number.quantiles=500,
                                       subsets=NULL,
                                       verbose=F) {
    probes <- meffil.probe.info()
    
    stopifnot(nrow(M) == nrow(U))
    stopifnot(ncol(M) == ncol(M))
    stopifnot(all(rownames(M) %in% probes$name))
    stopifnot(all(rownames(M) == rownames(U)))
    stopifnot(length(cell.types) == ncol(M))
    stopifnot(number.sites > 0)
    
    autosomal.sites <- get.autosomal.probes()
    beta <- meffil.get.beta(M=M,U=U)
    beta <- beta[which(rownames(beta) %in% autosomal.sites),]
    specific.beta <- meffil.cell.type.specific.methylation(beta, cell.types, number.sites, verbose)

    if (is.null(subsets))
        subsets <- get.island.probe.subsets()
    
    ## create target quantiles only from the Type ii probes
    typeii <- get.typeii.probes()
    subsets.ii <- lapply(subsets, function(set) intersect(set, typeii))

    ## probs argument for the quantile function
    probs <- seq(0,1,length.out=number.quantiles)
    quantiles <- lapply(subsets.ii, function(set) {
        set <- intersect(set, rownames(M))
        list(M=apply(apply(M[set,], 2, quantile, probs=probs, na.rm=T), 1, mean, na.rm=T),
             U=apply(apply(U[set,], 2, quantile, probs=probs, na.rm=T), 1, mean, na.rm=T))
    })

    list(beta=specific.beta, quantiles=quantiles, subsets=subsets)
}


#' List of available cell type references
#'
#' @export
meffil.get.cell.type.references <- function() {
    ls(reference.globals)
}

#' Select cell type reference for cell count estimates
#'
#' @param name Name of cell type reference previously created by
#' \code{\link{meffil.create.cell.type.reference}()}.
#'
#' Predefined references include "blood gse35069" and "blood gse35069 complete".
#' Both use methylation profiles from Reinius et al. 2012 (PMID 25424692)
#' for purified blood cell types.
#' The first is based on
#' six cell types: CD4T, CD8T, Mono, Bcell, NK, Gran.
#' The second is based on 
#' the same cell types but with Gran replaced by Neu and Eos.
#' @export
meffil.set.current.cell.type.reference <- function(name) {
    if (!is.null(name))
        check.cell.type.reference.exists(name)
    assign("current.cell.type.reference", name, pkg.globals)
}

#' Currently selected cell type reference
#' 
#' @export
meffil.get.current.cell.type.reference <- function() {
    if (!exists("current.cell.type.reference", pkg.globals))
        NULL
    else
        get("current.cell.type.reference", pkg.globals)
}

check.cell.type.reference.exists <- function(name) {
    if (!cell.type.reference.exists(name))
        stop(paste("Cell type reference", name,
                   "has not been created with meffil.create.cell.type.reference()"))
}    

cell.type.reference.exists <- function(name) {
    exists(name, reference.globals)
}

get.cell.type.reference.object <- function(name) {
    check.cell.type.reference.exists(name)
    get(name, reference.globals)
}
