#' IDAT file basenames
#'
#' List IDAT file basenames in a given directory.
#'
#' @param path Directory containing the IDAT files.
#' @param recursive If \code{TRUE}, search for IDAT files in subdirectories as well
#' (Default: \code{FALSE}).
#' @return Character vector of IDAT file basenames
#' (i.e. filenames with "_Grn.idat" and "_Red.idat" removed).
#' In other words, each identifies the Cy5 and Cy3 output files corresponding to a single microarray.
#'
#' @export
meffil.basenames <- function(path,recursive=FALSE) {
    grn.files <- list.files(path, pattern = "_Grn.idat(|\\.gz)$", recursive = recursive,
                            ignore.case = TRUE, full.names = TRUE)
    red.files <- list.files(path, pattern = "_Red.idat(|\\.gz)$", recursive = recursive,
                            ignore.case = TRUE, full.names = TRUE)
    get.basenames(c(grn.files, red.files))
}

get.basenames <- function(filenames)
    unique(gsub("_(Red|Grn).idat(|\\.gz)$", "", filenames))
