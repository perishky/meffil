library(knitr)
library(markdown)

run.test <- function(filename) {
    extension <- sub("^.*(\\.[^.]+)$", "\\1", filename)
    name <- sub("\\.[^.]+$", "", filename)

    if (!file.exists(name) && !file.exists(paste(name, "html", sep="."))) {
        if (extension == ".r") {
            cat(date(), "Sourcing", filename, "------------------\n")
            source(filename, echo=T)
        }
        else {
            cat(date(), "Knitting", filename, "-------------------\n")
            knit(paste(name, "rmd", sep="."))
            markdownToHTML(paste(name, "md", sep="."), paste(name, "html", sep="."))
        }
    }
}

run.test("450k-demo.r")
run.test("epic-demo.r")
run.test("450k-and-epic.r")
run.test("ewas.r")
run.test("cnv.rmd")
run.test("random.rmd")
run.test("minfi.rmd")

clean.up.outputs <- function() {
    unlink("450k-demo", recursive=T)
    unlink("epic-demo", recursive=T)
    unlink("450k-and-epic", recursive=T)
    unlink("ewas", recursive=T)
    unlink(c("cnv", "cnv.md", "cnv.html"), recursive=T)
    unlink(c("random", "random.md", "random.html"), recursive=T)
    unlink(c("minfi","minfi.md","minfi.html"), recursive=T)
    unlink(c("figure","cache"), recursive=T)
}

clean.up.data <- function() {
    unlink("data-450k-demo", recursive=T)
    unlink("data-epic-demo", recursive=T)
}

## To check for any errors:
##  grep -e Error *.md */*.md

## To copy outputs to OUTPUT_DIR:
##  cp -rv random random.{md,html} cnv cnv.{md,html} 450k-and-epic 450k-demo epic-demo ewas minfi minfi.{md,html} figure OUTPUT_DIR
