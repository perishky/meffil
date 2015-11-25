library(knitr)
library(markdown)

rmd2html <- function(basename) {
    cat(date(), "Knitting", basename, "-------------------\n")
    knit(paste(basename, "rmd", sep="."))
    markdownToHTML(paste(basename, "md", sep="."), paste(basename, "html", sep="."))
}

if (!file.exists("450k-demo")) source("450k-demo.r")

if (!file.exists("minfi")) rmd2html("minfi")

if (!file.exists("ecpi-demo")) source("epic-demo.r")

if (!file.exists("450k-and-epic")) source("450k-and-epic.r")

if (!file.exists("cnv")) rmd2html("cnv") ........ create cnv.rmd based on cnv.r

if (!file.exists("ewas")) source("ewas.r")

if (!file.exists("random")) rmd2html("random") ......... create  random.rmd based on random.r

clean.up.outputs <- function() {
    unlink("450k-demo", recursive = TRUE)
    unlink(c("minfi","minfi.md","minfi.html"), recursive=T)
    unlink("epic-demo", recursive=T)
    unlink("450k-and-epic", recursive=T)
    unlink(c("cnv", "cnv.md", "cnv.html"), recursive=T)
    unlink("ewas", recursive=T)
    unlink(c("random", "random.md", "random.html"), recursive=T)
    unlink(c("figure","cache"), recursive=T)
}

clean.up.data <- function() {
    unlink("data-450k-demo", recursive=T)
    unlink("data-epic-demo", recursive=T)
}

## cp -rv random random.{md,html} cnv.{md,html} 450k-and-epic 450k-demo epic-demo ewas minfi minfi.{md,html} figure OUTPUT_DIR
