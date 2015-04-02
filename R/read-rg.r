#library(illuminaio) ## for readIDAT()


#' Read Infinium HumanMethylation450 BeadChip.
#'
#' Reads Cy5 and Cy3 files for a given Infinium HumanMethylation450 BeadChip.
#'
#' @param basename IDAT file basename (see \code{\link{meffil.basenames}}).
#' @param probes Probe annotation used to construct the control matrix
#' (Default: \code{\link{meffil.probe.info}()}).
#' @param verbose If \code{TRUE}, then status messages are printed during execution (Default: \code{FALSE}).
#' @return List containing raw Cy5 ('R') and Cy3 ('G') data including
#' the intensity mean, intensity standard deviation and number of contributing beads.
#'
#' @export
meffil.read.rg <- function(basename, probes=meffil.probe.info(), verbose=F) {
    rg <- list(G=read.idat(paste(basename, "_Grn.idat", sep = ""), verbose=verbose),
               R=read.idat(paste(basename, "_Red.idat", sep=""), verbose=verbose))

    stopifnot(is.rg(rg))
    rg
}

read.idat <- function(filename, verbose=F) {
    msg("Reading", filename, verbose=verbose)

    if (!file.exists(filename))
        stop("Filename does not exist:", filename)
    illuminaio::readIDAT(filename)$Quants
}

is.rg <- function(rg) {
    (all(c("R","G") %in% names(rg))
     && is.matrix(rg$R) && is.matrix(rg$G)
     && nrow(rg$R) >= 100000 & nrow(rg$G) >= 100000
     && ncol(rg$R) == 3 && ncol(rg$G) == 3
     && length(rownames(rg$G)) == nrow(rg$G)
     && length(rownames(rg$R)) == nrow(rg$R))
}
