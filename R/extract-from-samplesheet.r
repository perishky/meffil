extract.from.samplesheet <- function(qc.objects, names) {
    ret <- do.call(rbind, lapply(qc.objects, function(object) {
        object$samplesheet[, unique(names), drop = F]
    }))
    for (name in colnames(ret))
        stopifnot(length(unique(na.omit(ret[,name]))) > 1)
    ret
}
