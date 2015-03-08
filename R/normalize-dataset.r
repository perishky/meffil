#' Functional normalization
#'
#' Apply functional normalization to a set of Infinium HumanMethylation450 BeadChip IDAT files.
#'
#' Fortin JP, Labbe A, Lemire M, Zanke BW, Hudson TJ, Fertig EJ, Greenwood CM, Hansen KD.
#' Functional normalization of 450k methylation array data
#' improves replication in large cancer studies.
#' Genome Biol. 2014 Dec 3;15(12):503. doi: 10.1186/s13059-014-0503-2.
#' PMID: 25599564
#'
#' @param path Directory containing the idat files.  Ignored if \code{filenames} is defined.
#' @param recursive If \code{TRUE}, any idat file in \code{path} or a subdirectory is included
#' in the normalization; otherwise, only those in the immediate directory are included
#' (Default: \code{FALSE}).
#' @param filenames Optional character vector
#' listing the idat files to include in the normalization. Filenames may omit
#' the "_Grn.idat"/"_Red.idat" suffix.
#' @param number.pcs Number of control matrix principal components to adjust for (Default: 2).
#' @param sex Optional character vector assigning a sex label ("M" or "F") to each sample.
#' @param beta If \code{TRUE} (default), then the function returns
#' the normalized matrix of methylation levels; otherwise, it returns
#' the normalized matrices of methylated and unmethylated signals.
#' @param pseudo Value to add to the denominator to make the methylation
#' estimate more stable when calculating methylation levels (Default: 100).
#' @param probes Probe annotation used to construct the control matrix
#' (Default: \code{\link{meffil.probe.info}()}).
#' @param verbose If \code{TRUE}, then status messages are printed during execution (Default: \code{FALSE}).
#' @param ... Arguments passed to \code{\link[parallel]{mclapply}()}.
#' @return Matrix of normalized methylation levels if \code{beta} is \code{TRUE};
#' otherwise matrices of normalized methylated and unmethylated signals.
#' Matrices returned have one column per sample and one row per CpG site.
#' Methylation levels are values between 0 and 1
#' equal to methylated signal/(methylated + unmethylated signal + pseudo).
#'
#' @return Normalized beta matrix (rows = CG sites, columns = samples, values = 0..1).
#'
#' @export
meffil.normalize.dataset <- function(path, recursive=F, filenames, number.pcs=2,
                                     sex=NULL,
                                     beta=T, pseudo=100,
                                     probes=meffil.probe.info(),
                                     verbose=F,
                                     ...) {

    stopifnot(!missing(filenames) || !missing(path))

    msg("Collecting idat files ...")
    if (missing(filenames))
        basenames <- meffil.basenames(path, recursive)
    else
        basenames <- get.basenames(filenames)

    stopifnot(length(basenames) > 1)

    msg("Computing normalization objects ...")
    norm.objects <- mclapply(basenames,
                             meffil.compute.normalization.object,
                             probes=probes, verbose=verbose, ...)

    msg("Normalizing objects ...")
    norm.objects <- meffil.normalize.objects(norm.objects,
                                             number.pcs=number.pcs,
                                             sex=sex,
                                             verbose=verbose)

    msg("Normalizing methylation data from normalized objects ...")
    norm.objects <- meffil.normalize.samples(norm.objects,
                                             beta=beta, pseudo=pseudo,
                                             probes=probes,
                                             verbose=verbose,
                                             ...)
    msg("Finished.")
    norm.objects
}
