library(knitr)
library(markdown)

options(markdown.HTML.options=union('toc', getOption("markdown.HTML.options")))

args <- commandArgs(TRUE)

rmd.filename <- args[1]

stopifnot(file.exists(rmd.filename))

knit(rmd.filename)

name <- sub("\\.[^.]+$", "", rmd.filename)
markdownToHTML(paste(name, "md", sep="."), paste(name, "html", sep="."))
