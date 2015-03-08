#library(illuminaio) ## for readIDAT()


#' Read Infinium HumanMethylation450 BeadChip.
#'
#' Reads Cy5 and Cy3 files for a given Infinium HumanMethylation450 BeadChip.
#'
#' @param basename IDAT file basename (see \code{\link{meffil.basenames}}).
#' @param probes Probe annotation used to construct the control matrix
#' (Default: \code{\link{meffil.probe.info}()}).
#' @param verbose If \code{TRUE}, then status messages are printed during execution (Default: \code{FALSE}).
#' @return List containing raw Cy5 ('R') and Cy3 ('G') intensities for the sample.
#'
#' @export
meffil.read.rg <- function(basename, probes=meffil.probe.info(), verbose=F) {
    rg <- list(G=read.idat(paste(basename, "_Grn.idat", sep = ""), verbose=verbose),
               R=read.idat(paste(basename, "_Red.idat", sep=""), verbose=verbose))
    rg$R <- rg$R[which(names(rg$R) %in% probes$address[which(probes$dye == "R")])]
    rg$G <- rg$G[which(names(rg$G) %in% probes$address[which(probes$dye == "G")])]

    stopifnot(length(rg$R) > 100000)
    stopifnot(length(rg$R) > 100000)

    rg
}

read.idat <- function(filename, verbose=F) {
    msg("Reading", filename, verbose=verbose)

    if (!file.exists(filename))
        stop("Filename does not exist:", filename)
    illuminaio::readIDAT(filename)$Quants[,"Mean"]
}

is.rg <- function(rg) {
    (all(c("R","G") %in% names(rg))
     && length(rg$R) >= 100000 & length(rg$G) >= 100000
     && is.vector(rg$R) && is.vector(rg$G)
     && length(names(rg$G)) == length(rg$G)
     && length(names(rg$R)) == length(rg$R))
}
