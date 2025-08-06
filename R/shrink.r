#' Shrink quality control object
#' 
#' Remove any information not absolutely needed for normalizing quantiles
#'
#' @param qc.object An output from \code{\link{meffil.create.qc.object}()}.
#' @param keep.vars Samplesheet variables to retain if the object will be
#' used to normalize quantiles using \code{\link{meffil.normalize.quantiles}()}
#' (Default: NULL).
#' @return A "QC object" identical to \code{qc.object} but containing
#' only information that is needed for normalizing quantiles.
#'
#' @examples \dontrun{
#' random.effects = c("Slide")
#' 
#' qc.objects1 = meffil.qc(samplesheet1)
#' dat1.objects = lapply(qc.objects1, meffil.shrink.qc.object, keep.vars=random.effects)
#' 
#' qc.objects2 = meffil.qc(samplesheet2)
#' dat2.objects = lapply(qc.objects2, meffil.shrink.qc.object, keep.vars=random.effects)
#' 
#' qc.objects = c(dat1.objects, dat2.objects)
#' norm.objects = meffil.normalize.quantiles(qc.objects, number.pcs=5, random.effects=random.effects)
#' dat1.norms = norm.objects[names(dat1.objects)]
#' dat2.norms = norm.objects[names(dat2.objects)]
#' 
#' norm.objects1 = mapply(meffil.expand.norm.object, dat1.norms, qc.objects1, SIMPLIFY=F)
#' meth1 = meffil.normalize.samples(norm.objects1)
#' 
#' norm.objects2 = mapply(meffil.expand.norm.object, dat2.norms, qc.objects2, SIMPLIFY=F)
#' meth2 = meffil.normalize.samples(norm.objects2)
#' }
#' 
#' @export
meffil.shrink.qc.object <- function(qc.object,keep.vars=NULL) {
  stopifnot(is.qc.object(qc.object))
  qc.object$samplesheet = qc.object$samplesheet[,unique(keep.vars),drop=F]  
  qc.object[c(
    "class",
    "sample.name",
    "version",
    "controls",
    "featureset",
    "chip",
    "quantiles",
    "dye.intensity",
    "intensity.R",
    "intensity.G",
    "predicted.sex",
    "median.m.signal",
    "median.u.signal",
    "bad.probes.detectionp",
    "bad.probes.beadnum",
    "bad.probes.detectionp.threshold",
    "bad.probes.beadnum.threshold",
    "samplesheet")]
}

#' Restore information needed for sample normalization
#' 
#' Restore information to normalized QC objects derived from shrunk QC objects 
#'
#' @param norm.object An output from \code{\link{meffil.normalize.quantiles}()}.
#' @param qc.object The output of \code{\link{meffil.create.qc.object}()}
#' that was used to create \code{norm.object}.
#' @return A normalized "QC object" identical to \code{norm.object}
#' but containing the full complement of information found in \code{qc.object}. 
#'
#' @examples \dontrun{
#' random.effects = c("Slide")
#' 
#' qc.objects1 = meffil.qc(samplesheet1)
#' dat1.objects = lapply(qc.objects1, meffil.shrink.qc.object, keep.vars=random.effects)
#' 
#' qc.objects2 = meffil.qc(samplesheet2)
#' dat2.objects = lapply(qc.objects2, meffil.shrink.qc.object, keep.vars=random.effects)
#' 
#' qc.objects = c(dat1.objects, dat2.objects)
#' norm.objects = meffil.normalize.quantiles(qc.objects, number.pcs=5, random.effects=random.effects)
#' dat1.norms = norm.objects[names(dat1.objects)]
#' dat2.norms = norm.objects[names(dat2.objects)]
#' 
#' norm.objects1 = mapply(meffil.expand.norm.object, dat1.norms, qc.objects1, SIMPLIFY=F)
#' meth1 = meffil.normalize.samples(norm.objects1)
#' 
#' norm.objects2 = mapply(meffil.expand.norm.object, dat2.norms, qc.objects2, SIMPLIFY=F)
#' meth2 = meffil.normalize.samples(norm.objects2)
#' }
#' 
#' @export
meffil.expand.norm.object <- function(norm.object,original) {
  stopifnot(is.qc.object(norm.object))
  stopifnot(is.qc.object(original))
  stopifnot(norm.object$sample.name == original$sample.name)
  for (item in setdiff(names(original), names(norm.object)))
    norm.object[[item]] = original[[item]]
  norm.object
}
