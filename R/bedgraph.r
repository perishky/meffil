#' Save EWAS effect estimates to bedgraph file
#'
#' Saves EWAS effect estimates  to a bedgraph file
#' for viewing on a genome browser.
#' More file format details can be found here:
#' https://genome.ucsc.edu/goldenPath/help/bedgraph.html
#'
#' @param ewas.object Object returned by \code{\link{meffil.ewas()}}.
#' @param filename Filename for output, typically with a 'bed' extension. 
#' @param analysis The particular EWAS analysis from which to obtain summary statistics.
#' This should be one of \code{names(ewas.object$analyses)}.
#' @param name Text name to be included in the bedgraph header. 
#' @param description Text description to be included in the bedgraph header. 
#' @param header Bedgraph header. The default header uses the name and description provided. 
#' @export
meffil.ewas.bedgraph <- function(ewas.object,                               
                                 filename,
                                 analysis,
                                 name,
                                 description,
                                 header) {
    stopifnot(is.ewas.object(ewas.object))
    if (!analysis %in% names(ewas.object$analyses)) 
        stop("'analysis' should be one of: ", paste(names(ewas.object$analyses), collapse=", "))

    x <- ewas.object$analyses[[analysis]]$table
    x$start <- x$end <- x$position    
    x$chr <- as.character(x$chr)
    x <- x[which(!is.na(x$chr) & nchar(x$chr) > 0),]
    x <- x[order(x$chr, x$start, decreasing=F),]

    if (missing(name))
        name <- basename(filename)
    if (missing(description))
        description <- name    
    if (missing(header))
        header <- paste("track type=bedGraph name='", name, "' description='", description,"'", sep="")
    
    writeLines(c(header,
                 paste(x$chr,
                       " ", format(x$start, scientific=F,trim=T,digits=5),
                       " ", format(x$end, scientific=F, trim=T,digits=5),
                       if ("name" %in% colnames(x)) x$name else "",
                       " ", format(x$coefficient, scientific=F, trim=T,digits=5), sep="")),
               con=filename)
}
