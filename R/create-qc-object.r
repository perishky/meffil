#' Quality control object
#'
#' Create a quality control object for a given Infinium HumanMethylation450 BeadChip.
#'
#' @param samplesheet.row Row from the data frame containing IDAT file and sample info (see \code{\link{meffil.read.samplesheet}} or \code{\link{meffil.create.samplesheet}}).
#' @param number.quantiles Number of quantiles to compute for probe subset (Default: 500).
#' @param dye.intensity Reference intensity for scaling each color channel (Default: 5000).
#' @param verbose If \code{TRUE}, then status messages are printed during execution (Default: \code{FALSE}).
#' @param detection.threshold Default value = 0.01.
#' All probes above this detection threshold detected.
#' @param bead.threshold Default value = 3.
#' All probes with less than this number of beads detected.
#' @param sex.cutoff Sex prediction cutoff. Default value = -2.
#' @param chip Name returned by \code{\link{meffil.list.chips()}} (Default: NA).
#' @param featureset Name returned by \code{\link{meffil.list.featuresets()}} (Default: \code{chip}).
#' @param cell.type.reference Character string name of the cell type reference
#' to use for estimating cell counts. Estimates are not generated if set to NA (default).
#' See \code{\link{meffil.list.cell.type.references}()} for a list of available
#' references.  New references can be created using
#' \code{\link{meffil.add.cell.type.reference}()}. 
#' @return List containing control probe information, probe summaries
#' and quantiles.  We call this a "QC object".
#'
#' @export
meffil.create.qc.object <- function(samplesheet.row,
                                    number.quantiles=500,
                                    dye.intensity=5000,
                                    verbose=F,
                                    detection.threshold=0.01,
                                    bead.threshold=3,
                                    sex.cutoff=-2,
                                    chip=NA,
                                    featureset=chip,
                                    cell.type.reference=NA) {
    stopifnot(number.quantiles >= 100)
    stopifnot(dye.intensity >= 100)
    stopifnot(samplesheet.row$Sex %in% c(NA, "F", "M"))

    rg <- read.rg(samplesheet.row$Basename, verbose=verbose)

    chip <- guess.chip(rg, chip)
    
    if (is.na(featureset))
        featureset <- chip
    
    probes <- meffil.probe.info(chip, featureset)

    detectionp <- extract.detection.pvalues(rg, probes, verbose=verbose)
    bad.probes.detectionp <- detectionp[which(detectionp > detection.threshold)]
    
    beadnum <- extract.beadnum(rg, probes, verbose=verbose)
    bad.probes.beadnum <- beadnum[which(beadnum < bead.threshold)]
    
    snp.betas <- extract.snp.betas(rg, probes, verbose=verbose)

    controls <- extract.controls(rg, probes, verbose=verbose)

    tryCatch({
        rg <- background.correct(rg, probes, verbose=verbose)
    }, error = function(e) {
        save.image(paste(samplesheet.row$Basename, "rda", sep="."))
        stop(e)
    })

    intensity.R <- calculate.intensity.R(rg, probes)
    intensity.G <- calculate.intensity.G(rg, probes)

    rg <- dye.bias.correct(rg, probes, dye.intensity, verbose=verbose)

    mu <- rg.to.mu(rg, probes)
    
    sites.x <- meffil.get.x.sites(featureset)
    x.signal <- median(log(mu$M[sites.x] + mu$U[sites.x], 2), na.rm=T)

    sites.y <- meffil.get.y.sites(featureset)
    y.signal <- median(log(mu$M[sites.y] + mu$U[sites.y], 2), na.rm=T)

    probs <- seq(0,1,length.out=number.quantiles)

    quantiles <- lapply(get.quantile.site.subsets(featureset), function(sets) {
        list(M=unname(quantile(mu$M[sets], probs=probs,na.rm=T)),
             U=unname(quantile(mu$U[sets], probs=probs,na.rm=T)))
    })

    msg("predicting sex", verbose=verbose)
    xy.diff <- y.signal-x.signal
    predicted.sex <- ifelse(xy.diff < sex.cutoff, "F","M")

    cell.counts <- NULL
    if (!is.na(cell.type.reference))
        cell.counts <- estimate.cell.counts.from.mu(mu, cell.type.reference, verbose)
    
    list(class="qc.object",
         version=packageVersion("meffil"),
         sample.name=samplesheet.row$Sample_Name,
         basename=samplesheet.row$Basename,
         controls=controls,
         featureset=featureset,
         chip=chip,
         quantiles=quantiles,
         dye.intensity=dye.intensity,
         intensity.R=intensity.R,
         intensity.G=intensity.G,
         x.signal=x.signal,
         y.signal=y.signal,
         xy.diff=xy.diff,
         sex=samplesheet.row$Sex,
         sex.cutoff=sex.cutoff,
         predicted.sex=predicted.sex,
         median.m.signal=median(mu$M,na.rm=T),
         median.u.signal=median(mu$U,na.rm=T),
         bad.probes.detectionp=bad.probes.detectionp,
         bad.probes.beadnum=bad.probes.beadnum,
         bad.probes.detectionp.threshold=detection.threshold,
         bad.probes.beadnum.threshold=bead.threshold,
         snp.betas=snp.betas,
         cell.counts=cell.counts,
         samplesheet=samplesheet.row
         )
}

is.qc.object <- function(object) {
    is.list(object) &&
        ("class" %in% names(object) && object$class %in% "qc.object"
         || "origin" %in% names(object) && object$origin == "meffil.create.qc.object") ## backward compatibility
}




