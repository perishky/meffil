#library(illuminaio) ## for readIDAT()


#' Read Infinium HumanMethylation450 BeadChip.
#'
#' Reads Cy5 and Cy3 files for a given Infinium HumanMethylation450 BeadChip.
#'
#' @param basename IDAT file basename (see \code{\link{meffil.basenames}}).
#' @param verbose If \code{TRUE}, then status messages are printed during execution (Default: \code{FALSE}).
#' @return List containing raw Cy5 ('R') and Cy3 ('G') data including
#' the intensity mean, intensity standard deviation and number of contributing beads.
read.rg <- function(basename, verbose=F) {
    list(G=read.idat(paste(basename, "_Grn.idat", sep = ""), verbose=verbose),
         R=read.idat(paste(basename, "_Red.idat", sep=""), verbose=verbose),
         class="rg",
         version=packageVersion("meffil"),
         basename=basename)
}

read.idat <- function(filename, verbose=F) {
    msg("Reading", filename, verbose=verbose)

    if (!file.exists(filename)) {
        filename <- paste(filename, "gz", sep=".")
        if (!file.exists(filename))
            stop("Filename does not exist:", filename, sub(".gz", "", filename))
    }
    illuminaio::readIDAT(filename)$Quants
}

exists.rg <- function(basename, verbose=F) {
    file.exists(paste(basename, "_Grn.idat", sep = "")) && file.exists(paste(basename, "_Red.idat", sep=""))
}

is.rg <- function(rg) {
    is.list(rg) && "class" %in% names(rg) && rg$class == "rg"
}
