#' Create a cell type reference object
#'
#' Create a cell type reference object for estimating cell counts
#' with the Infinium HumanMethylation450 BeadChip.
#'
#' @param basename IDAT file basename (see \code{\link{meffil.basenames}()}).
#' @param cell.types Vector of cell type names corresponding to sample \code{basename}s.
#' @param probes Probe annotation used to construct the control matrix
#' (Default: \code{\link{meffil.probe.info}()}).
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
#' - \code{norm.objects} The normalization objects used to compute \code{beta}
#' using \code{\link{meffil.normalize.objects}()} corresponding to the IDAT files
#' identified by \code{basenames}.
#' - \code{quantiles} The average quantiles of methylated and unmethylated signals
#' of probe sets defined by \code{subsets} (see below).
#' e.g. \code{quantiles[[name]]$M} provides the quantiles (\code{number.quantiles} quantiles)
#' of the probes specified by \code{subsets[[name]]}.
#' - \code{subsets} Probes on the microarray partitioned by relationship to CpG islands,
#' either in an island, in a shore or far from an island.
#'
#' @export
meffil.create.cell.type.reference <- function(basenames, cell.types, 
                                              number.quantiles=500,
                                              number.sites=50,
                                              probes=meffil.probe.info(),
                                              temp.dir=NULL, verbose=F) {
    msg("cell types", paste(unique(cell.types), collapse=", "), verbose=verbose)
    
    stopifnot(length(basenames) == length(cell.types))

    norm.objects <- mclapply(basenames, meffil.compute.normalization.object, verbose=verbose)

    design.matrix <- meffil.design.matrix(norm.objects)
    variance <- apply(design.matrix, 2, var)
    number.pcs <- max(2, length(which(variance/sum(variance) > 0.01)))
    norm.objects <- meffil.normalize.objects(norm.objects, number.pcs=number.pcs)

    mu <- meffil.normalize.samples(norm.objects, beta=F, temp.dir=temp.dir, verbose=verbose)
        
    beta <- meffil.get.beta(mu)
    autosomal.sites <- probes$name[which(probes$chr %in% paste0("chr",1:22))]
    beta <- beta[which(rownames(beta) %in% autosomal.sites),]
    specific.beta <- meffil.cell.type.specific.methylation(beta, cell.types, number.sites, verbose)

    ## target probes
    subsets <- get.island.probe.subsets(probes)

    ## create target quantiles only from the Type ii probes
    typeii <- unique(probes$name[which(probes$type == "ii")])
    subsets.ii <- lapply(subsets, function(set) intersect(set, typeii))

    ## probs argument for the quantile function
    probs <- seq(0,1,length.out=number.quantiles)

    quantiles <- lapply(subsets.ii, function(set) {
        list(M=apply(apply(mu$M[set,], 2, quantile, probs=probs, na.rm=T), 1, mean, na.rm=T),
             U=apply(apply(mu$U[set,], 2, quantile, probs=probs, na.rm=T), 1, mean, na.rm=T))
    })

    list(beta=specific.beta, norm.objects=norm.objects, quantiles=quantiles, subsets=subsets)
}

is.cell.type.reference <- function(object) {
    all(c("beta","norm.objects","quantiles","subsets") %in% names(object))
}
