# ---- download-450k-lead-dataset ----

download.450k.lead.dataset <- function() {
    dir.create(path <- "data-450k-lead")
    
    if (length(list.files(path, "*.idat$")) == 0) {
        filename <-  file.path(path, "gse69633.tar")
        download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE69633&format=file", filename)
        cat(date(), "Extracting files from GEO archive.\n")
        system(paste("cd", path, ";", "tar xvf", basename(filename)))
        unlink(filename)
        cat(date(), "Unzipping IDAT files.\n")
        system(paste("cd", path, ";", "gunzip *.idat.gz"))
        
        library(GEOquery)
        geo <- getGEO("GSE69633", GSEMatrix=F)
        ## this output from geo is really messed up.
        xx <- lapply(geo@gsms, function(x) x@header)
        characteristics <- t(sapply(xx, function(x) x$characteristics_ch1))
        colnames(characteristics) <- sub("([^:]+):.*", "\\1", characteristics[1,])
        rownames(characteristics) <- names(xx)
        characteristics <- apply(characteristics, 2, function(x) sub("[^:]+: (.*)", "\\1", x))
        characteristics <- as.data.frame(characteristics, stringsAsFactors=F)
        for (name in c("socioeconomic score", "gestational age", "birth weight", "pbconc (ng/dl)"))
            characteristics[[name]] <- as.numeric(characteristics[[name]])
        write.csv(characteristics, file=file.path(path, "samples.csv"))
    }
    
    path
}



