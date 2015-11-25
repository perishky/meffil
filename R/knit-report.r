#' Knit report using working environment
#'
#' Warning: It is quite likely that this will be called within an RMD file
#' implying a recursive call to knit(). This will generate "duplicate label"
#' errors for unlabelled chunks. To avoid this, all code chunks
#' in \code{rmd.filename} should be named.
#' Supposedly this error can also be avoided by setting the following option:
#'      options(knitr.duplicate.label='allow')
#' I tried this but it didn't seem to help.
#'
#' @param rmd.filename RMD file.
#' @param output.filename Markdown or HTML output file.  An HTML file
#' is specified using the .htm, .html, .HTM or .HTML file extension.
#' When html is specified, a similarly named markdown file is also
#' generated.
#' All output files including cache and figures will appear in the
#' same folder as \code{output.filename}.
#' 
#' @param  ... Arguments to be passed to \code{\link{knitr::knit}}
#' @return NULL
knit.report <- function(input.filename, output.filename, ...) {
    input.filename <- normalizePath(input.filename)

    output.dir <- dirname(output.filename)
    if (!file.exists(output.dir))
        dir.create(output.dir, recursive=T)

    output.dir <- normalizePath(output.dir)
    output.filename <- file.path(output.dir, basename(output.filename))
    
    current.dir <- getwd()
    on.exit(setwd(current.dir))
    setwd(output.dir)

    name <- gsub("\\.[^.]+$", "", basename(output.filename))
    suffix <- gsub(".*\\.([^.]+)$", "\\1", output.filename)
    is.html <- tolower(suffix) %in% c("htm","html")

    if (is.html)
        md.filename <- paste(name, "md", sep=".")
    else
        md.filename <- basename(output.filename)

    ## the rmd code can reference the output directory, e.g. to save a file there.
    assign("output.dir", output.dir, envir=parent.frame())
    
    knit(input.filename, output=md.filename, envir=parent.frame(), ...)

    if (is.html)
        markdownToHTML(md.filename, basename(output.filename))
}

