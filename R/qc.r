#' Perform QC on HumanMethylation450 idat files
#'
#' Read in control matrices for each sample.
#' Perform background correction and R/G dye bias correction. 
#' Predict sex 
#' @param samplesheet Data frame containing IDAT file and sample info (see \code{\link{meffil.read.samplesheet}} pr \code{\link{meffil.create.samplesheet}}).
#' @param number.quantiles Number of quantiles to compute for probe subset (Default: 500).
#' @param dye.intensity Reference intensity for scaling each color channel (Default: 5000).
#' @param detection.threshold Default value = 0.01.
#' All probes above this detection threshold detected.
#' @param bead.threshold Default value = 3.
#' All probes with less than this number of beads detected.
#' @param sex.cutoff Sex prediction cutoff. Default value = -2.
#' @param chip Name returned by \code{\link{meffil.list.chips()}} (Default: \code{NA}).
#' @param featureset Name returned by \code{\link{meffil.list.featuresets()}} (Default: \code{chip}).
#' @param cell.type.reference Character string name of the cell type reference
#' to use for estimating cell counts. Estimates are not generated if set to NA (default).
#' See \code{\link{meffil.list.cell.type.references}()} for a list of available
#' references.  New references can be created using
#' \code{\link{meffil.create.cell.type.reference}()}. 
#' @param verbose If \code{TRUE}, then status messages are printed during execution (Default: \code{FALSE}).
#' @return List containing control probe information, probe summaries
#' and quantiles.
#'
#' @export
meffil.qc <- function(samplesheet, number.quantiles=500, dye.intensity=5000,
                      detection.threshold=0.01, bead.threshold=3, sex.cutoff=-2,
                      chip=NA,
                      featureset=chip,
                      cell.type.reference=NA,
                      max.bytes=2^30-1, ## maximum number of bytes that can be returned by mclapply
                      verbose=F, ...) {
    check.samplesheet(samplesheet)

    stopifnot(is.na(featureset) || is.character(featureset) && featureset %in% meffil.list.featuresets())
    stopifnot(is.na(chip) || is.character(chip) && chip %in% meffil.list.chips())
    if (!is.na(featureset) && !is.na(chip))
        stopifnot(is.compatible.chip(featureset, chip))

    samplesheet.row <- lapply(1:nrow(samplesheet), function(i) samplesheet[i,,drop=F])
    names(samplesheet.row) <- samplesheet$Sample_Name

    qc.objects <- mclapply.safe(
        samplesheet.row, 
        meffil.create.qc.object,
        number.quantiles=number.quantiles,
        dye.intensity=dye.intensity,
        verbose=verbose,
        detection.threshold=detection.threshold,
        bead.threshold=bead.threshold,
        sex.cutoff=sex.cutoff,
        featureset=featureset,
        chip=chip,
        cell.type.reference=cell.type.reference,
        ...,
        max.bytes=max.bytes)

    names(qc.objects) <- sapply(qc.objects, function(x) {
        if ("samplesheet" %in% names(x))
            if ("Sample_Name" %in% names(x$samplesheet))
                return(x$samplesheet$Sample_Name)
        "Error"
    })
    
    is.error <- sapply(qc.objects, class) == "try-error"
    if (any(is.error))
        warning("Errors were encountered, do this: which(sapply(qc.objects, class) == 'try-error')")
    
    if (length(unique(sapply(qc.objects[which(!is.error)], function(object) object$featureset))) > 1)
        warning("Heterogeneous microarray formats included without setting 'featureset'.")
    
    return(qc.objects)
}

