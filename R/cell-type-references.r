#' List of available cell type references
#'
#' @details Names and description of available references:
#' ```{r results="asis",echo=F}
#' references <- comment(meffil.list.cell.type.references())
#' references <- references[order(names(references))]
#' for (i in 1:length(references)) 
#'   cat("\n*", paste0('"',names(references)[i],'"'), references[[i]])
#' cat("\n")
#' #' @md
#' ```
#' 
#' @examples
#' ## obtain a list of references
#' references <- meffil.list.cell.type.references()
#' ## show descriptions for each
#' comment(references)
#' 
#' @export
meffil.list.cell.type.references <- function() {
    references <- ls(reference.globals)
    structure(
        references,
        comment=sapply(reference.globals, function(ref) ref$description))
}
 
get.cell.type.reference <- function(name) {
    stopifnot(is.character(name) && name %in% meffil.list.cell.type.references())
    get(name, reference.globals)
}

#' Create a cell type reference object
#'
#' Create a cell type reference object for estimating cell counts
#' with the Infinium HumanMethylation450 BeadChip.
#'
#' @param name Character string providing the name of the reference.
#' @param M Matrix of methylated probe intensities (rows=CpG sites, columns=samples).
#' @param U Matrix of unmethylatd probe intensities (rows=CpG sites, columns=samples).
#' @param cell.types Vector of cell type names corresponding to sample \code{basename}s.
#' @param chip Name returned by \code{\link{meffil.list.chips()}} (Default: NA).
#' @param featureset Name returned by \code{\link{meffil.list.featuresets()}} (Default: chip).
#' @param number.quantiles Length of numeric sequence to specify probe intensity distributions
#' (Default: 500).
#' @param number.sites Number of probes to characterise cell type methylation (Default: 50).
#' For each cell type, this number of probes with greater methylation than other cell types
#' and the same number with lesser methylation than the other cell types will be included.
#' @param specific.sites If not null (default), then \code{number.sites} is ignored and
#' the supplied site identifiers are used to differentiate between cell types instead
#' of those maximally different between the cell types within the reference.
#' @param object Cell type reference previously created by this function.
#' If not \code{NULL}, then this reference is added
#' and all other function arguments are ignored (Default: NULL).
#' @param description Text description of the reference (Default: NULL).
#' @param verbose If \code{TRUE}, then status messages are printed during execution
#' (Default: \code{FALSE}).
#' @return A list specifying a cell type reference object that can be used by
#' \code{\link{meffil.estimate.cell.counts}()} to estimate cell counts
#' in another dataset.
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
meffil.add.cell.type.reference <- function(name, M, U, cell.types,
                                           chip=NA,
                                           featureset=chip,
                                           number.sites=50,
                                           specific.sites=NULL,
                                           number.quantiles=500,
                                           subsets=NULL,
                                           object=NULL,
                                           description=NULL,
                                           verbose=F) {
    if (is.null(object)) {
        object <- create.cell.type.reference(M,U,cell.types,chip,featureset,
                                             number.sites,specific.sites,number.quantiles,
                                             subsets,verbose)
        object$name <- name
        object$description <- description
    }
    else
        stopifnot(is.cell.type.reference(object))
    assign(object$name, object, reference.globals)
    invisible(object)
}

create.cell.type.reference <- function(M, U, cell.types,
                                       chip=NA,
                                       featureset=chip,
                                       number.sites=50,
                                       specific.sites=NULL,
                                       number.quantiles=500,
                                       subsets=NULL,
                                       verbose=F) {
    chip <- guess.chip(M, chip)
    
    if (is.na(featureset))
        featureset <- chip

    stopifnot(is.compatible.chip(featureset, chip))
    
    stopifnot(nrow(M) == nrow(U))
    stopifnot(ncol(M) == ncol(M))
    stopifnot(all(rownames(M) == rownames(U)))
    stopifnot(length(cell.types) == ncol(M))
    stopifnot(number.sites > 0)
    
    autosomal.sites <- meffil.get.autosomal.sites(featureset)
    beta <- meffil.get.beta(M=M,U=U)
    beta <- beta[which(rownames(beta) %in% autosomal.sites),]
    if (is.null(specific.sites))
        specific.beta <- meffil.cell.type.specific.methylation(beta, cell.types,
                                                               number.sites, verbose)
    else {
        stopifnot(is.character(specific.sites))
        specific.sites <- intersect(specific.sites, rownames(beta))
        specific.beta <- meffil.cell.type.specific.methylation(beta[specific.sites,],
                                                               cell.types,
                                                               length(specific.sites),
                                                               verbose)
    }

    if (is.null(subsets))
        subsets <- get.island.site.subsets(featureset)
    
    ## create target quantiles only from the Type ii probes
    typeii <- meffil.get.typeii.sites(featureset)
    subsets.ii <- lapply(subsets, function(set) intersect(set, typeii))

    ## probs argument for the quantile function
    probs <- seq(0,1,length.out=number.quantiles)
    quantiles <- lapply(subsets.ii, function(set) {
        set <- intersect(set, rownames(M))
        list(M=apply(apply(M[set,], 2, quantile, probs=probs, na.rm=T), 1, mean, na.rm=T),
             U=apply(apply(U[set,], 2, quantile, probs=probs, na.rm=T), 1, mean, na.rm=T),
             beta=apply(apply(get.beta(M[set,],U[set,]), 2, quantile, probs=probs, na.rm=T), 1, mean, na.rm=T))
    })

    list(class="cell.type.reference",
         version=packageVersion("meffil"),
         featureset=featureset,
         beta=specific.beta,
         quantiles=quantiles,
         subsets=subsets)
}

is.cell.type.reference <- function(object)
    "class" %in% names(object) && object[["class"]] == "cell.type.reference"
