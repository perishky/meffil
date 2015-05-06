#' Cell count estimates
#'
#' @param qc.objects List of objects obtained from \code{\link{meffil.qc}()}
#' or \code{\link{meffil.create.qc.object}()}.
#' 
#' @export 
meffil.cell.count.estimates <- function(qc.objects) {
    stopifnot(sapply(qc.objects, is.qc.object))
    cell.counts <- NULL
    try(cell.counts <- sapply(qc.objects, function(object) object$cell.counts$counts))
    if (!is.matrix(cell.counts)) return(NULL)
    cell.counts
}
