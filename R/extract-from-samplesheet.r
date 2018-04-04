# select columns identified by \code{names}
# from the samplesheet used to create the \code{qc.objects}. 
extract.from.samplesheet <- function(qc.objects, names) {
    ret <- do.call(rbind, lapply(qc.objects, function(object) {
        stopifnot(is.qc.object(object))
        object$samplesheet[, unique(names), drop = F]
    }))
    ret
}
