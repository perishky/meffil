# ---- download-450k-demo-dataset ----

download.450k.demo.dataset <- function() {
    dir.create(path <- "data-450k-demo")
    
    if (length(list.files(path, "*.idat$")) == 0) {
        filename <-  file.path(path, "gse55491.tar")
        download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55491&format=file", filename)
        cat(date(), "Extracting files from GEO archive.\n")
        system(paste("cd", path, ";", "tar xvf", basename(filename)))
        unlink(filename)
        cat(date(), "Unzipping IDAT files.\n")
        system(paste("cd", path, ";", "gunzip *.idat.gz"))

        library(GEOquery)
        geo <- getGEO("GSE55491", GSEMatrix=F)
        geo <- lapply(geo@gsms, function(gsm) unlist(gsm@header))
        geo <- do.call(rbind, geo)
        geo <- as.data.frame(geo, stringAsFactors=F)
        geo$group <- geo$characteristics_ch13
        geo$sex <-   geo$characteristics_ch11
        write.csv(geo, file=file.path(path, "samples.csv"))
    }
    
    path
}
