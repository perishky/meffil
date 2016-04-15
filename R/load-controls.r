#' Load control probes
#'
#' @param samplesheet Sample info (see \code{\link{meffil.read.samplesheet}} or \code{\link{meffil.create.samplesheet}}).
#' @param verbose (Default: FALSE).
#' @param chip Name returned by \code{\link{meffil.list.chips()}} (Default: \code{NA}).
#' @param featureset Name returned by \code{\link{meffil.list.featuresets()}} (Default: \code{chip}).
#' @param ... Arguments to mclapply.
#' @return List containing two elements: probes and values.
#' The probes item is a data frame describing the control probes.
#' The values item is a matrix providing the intensities of the
#' control probes for each samples (rows=probes, columns=samples).
#' 
#' @export
meffil.load.controls <- function(samplesheet,
                                 chip=NA,
                                 featureset=chip,                                 
                                 verbose=F, ...) {
    meffil:::check.samplesheet(samplesheet)
    
    stopifnot(is.na(featureset) || is.character(featureset) && featureset %in% meffil.list.featuresets())
    stopifnot(is.na(chip) || is.character(chip) && chip %in% meffil.list.chips())
    if (!is.na(featureset) && !is.na(chip))
        stopifnot(is.compatible.chip(featureset, chip))

    probes.R <- probes.G <- NULL
    values <- meffil:::mclapply.safe(1:nrow(samplesheet), function(i) {
        msg(samplesheet$Sample_Name[i], samplesheet$Basename[i])
        rg <- meffil:::read.rg(samplesheet$Basename[i], verbose=verbose)
        this.chip <- meffil:::guess.chip(rg, chip)
        if (is.na(chip))
            chip <<- this.chip
        if (is.na(featureset))
            featureset <- chip
        if (this.chip != chip || is.null(probes.R)) {
            probes <- meffil.probe.info(this.chip, featureset)
            probes <- probes[which(probes$type == "control"),]
            rownames(probes) <- with(probes, paste(target, dye, name, sep="."))
            probes.R <<- probes[which(probes$dye == "R"),]
            probes.G <<- probes[which(probes$dye == "G"),]
            probes.R <<- probes.R[order(rownames(probes.R)),]
            probes.G <<- probes.G[order(rownames(probes.G)),]
        }
        c(rg$R[match(probes.R$address, rownames(rg$R)), "Mean"],
          rg$G[match(probes.G$address, rownames(rg$G)), "Mean"])
    }, ...)
    is.error <- sapply(values, class) == "try-error"
    if (any(is.error))
        stop(values[which(is.error)[1]])
    names(values) <- samplesheet$Sample_Name
    values <- do.call(cbind, values)
    probes <- rbind(probes.R, probes.G)
    probes$address <- NULL
    rownames(values) <- rownames(probes)
    list(probes=probes, values=values)
}
