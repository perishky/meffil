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
                                    sex.cutoff=-2) {
    stopifnot(number.quantiles >= 100)
    stopifnot(dye.intensity >= 100)
    stopifnot(samplesheet.row$Sex %in% c(NA, "F", "M"))

    rg <- read.rg(samplesheet.row$Basename, verbose=verbose)

    bad.probes.detectionp <- identify.bad.probes.detectionp(rg, detection.threshold, verbose=verbose)

    bad.probes.beadnum <- identify.bad.probes.beadnum(rg, bead.threshold, verbose=verbose)

    snp.betas <- extract.snp.probe.betas(rg, verbose=verbose)

    controls <- extract.controls(rg, verbose=verbose)

    rg <- background.correct(rg, verbose=verbose)

    intensity.R <- calculate.intensity.R(rg)
    intensity.G <- calculate.intensity.G(rg)

    rg <- dye.bias.correct(rg, dye.intensity, verbose=verbose)

    mu <- rg.to.mu(rg)
    
    probes.x <- get.x.probes()
    x.signal <- median(log(mu$M[probes.x] + mu$U[probes.x], 2), na.rm=T)

    probes.y <- get.y.probes()
    y.signal <- median(log(mu$M[probes.y] + mu$U[probes.y], 2), na.rm=T)

    probs <- seq(0,1,length.out=number.quantiles)

    quantiles <- lapply(get.quantile.probe.subsets(), function(sets) {
        list(M=unname(quantile(mu$M[sets], probs=probs,na.rm=T)),
             U=unname(quantile(mu$U[sets], probs=probs,na.rm=T)))
    })

    msg("predicting sex", verbose=verbose)
    xy.diff <- y.signal-x.signal
    predicted.sex <- ifelse(xy.diff < sex.cutoff, "F","M")

    if (!is.null(meffil.get.current.cell.type.reference()))
        cell.counts <- estimate.cell.counts.from.mu(mu, verbose)
    else
        cell.counts <- NULL
    
    list(origin="meffil.create.qc.object",
         sample.name=samplesheet.row$Sample_Name,
         basename=samplesheet.row$Basename,
         controls=controls,
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
    "origin" %in% names(object) && object$origin == "meffil.create.qc.object"
}


