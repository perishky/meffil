# ## post normalisation data


# use beta values to calculate PCs on 20k most variable probes

# plot each PC against each of chiprow, chip column, plate, slide IF there is a significant lm
# - see if technical variance is gone





#' Plot scree plot of control matrix
#'
#' @param  norm.objects From \code{meffil.normalize.objects}
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}

meffil.plot.control.scree <- function(norm.objects)
{

}


#' Plot each PC against each of chiprow, chip column, plate, slide IF there is a significant lm
#'
#' @param samplesheet From \code{read.450k.sheet}
#' @param  norm.objects From \code{meffil.normalize.objects}
#' @param  threshold P-value threshold for whether to plot. Default 0.05
#' @param  pcs Which PCs to plot. Default first 10
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.pc.batch <- function(samplesheet, norm.objects, threshold = 0.05, pcs=1:10)
{

}
