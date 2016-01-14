#' Knit report using working environment
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
    
    name <- gsub("\\.[^.]+$", "", basename(output.filename))
    suffix <- gsub(".*\\.([^.]+)$", "\\1", output.filename)
    is.html <- tolower(suffix) %in% c("htm","html")

    if (is.html)
        md.filename <- file.path(output.dir, paste(name, "md", sep="."))
    else
        md.filename <- output.filename

    current.dir <- opts_knit$get("output.dir")
    on.exit(opts_knit$set(output.dir=current.dir))
    opts_knit$set(output.dir=output.dir)

    cwd <- getwd()
    on.exit(setwd(cwd), add=TRUE)

    if (is.null(options("knitr.in.progress")[[1]])) {
        setwd(output.dir)
        knit(input.filename, output=md.filename, envir=parent.frame(), ...)
    }
    else {
        opts_knit$set(progress=FALSE)
        out <- knit_child(input.filename, envir=parent.frame(), quiet=T)
        writeLines(out, con=md.filename)
    }

    lines <- readLines(md.filename)
    lines <- gsub("![plot", "\n\n![plot", lines, fixed=T)
    writeLines(lines, md.filename)

    if (is.html)
        markdownToHTML(md.filename, output.filename)
}


