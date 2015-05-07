#' Knit report using working environment and convert to HTML
#'
#' Warning: It is quite likely that this will be called within an RMD file
#' implying a recursive call to knit(). this will generate "duplicate label"
#' errors for unlabelled chunks. To avoid this, all RMD files in this
#' package should contain only named chunks.
#' Supposedly this error can also be avoided by setting the following option:
#'      options(knitr.duplicate.label='allow')
#' I tried this but it didn't seem to help.
#' 
#' @param  normalization.summary Output from \code{meffil.normalization.summary}.
#' @param  output.file Default = "meffil.normalization.report.html"
#' If specified then a html report will be generated summarising the normalization.
#' @param  ... Arguments to be passed to \code{\link{knitr::knit}}
#' @return NULL
rmd2html <- function(input.filename, output.filename, ...) {
    output.filename <- normalizePath(output.filename)

    output.dir <- dirname(output.filename)
    if (!file.exists(output.dir))
        dir.create(output.dir)

    current.dir <- getwd()
    on.exit(setwd(current.dir))
    setwd(output.dir)
    
    name <- gsub("\\.[^.]+$", "", basename(output.filename))
    md.filename <- paste(name, "md", sep=".")
    knit(input.filename, output=md.filename, envir=parent.frame(), ...)
    markdownToHTML(md.filename, basename(output.filename))
}
