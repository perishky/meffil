# ---- download-epic-demo-dataset ----

download.epic.demo.dataset <- function() {
    dir.create(path <- "data-epic-demo")
    
    if (length(list.files(path, "*.idat$")) == 0) {
        download.zip <- function(url, path) {
            filename <- file.path(path, "data.zip")
            download.file(url, filename)
            filenames <- unzip(filename, junkpaths=T, exdir=path)
            unlink(filename)
            invisible(filenames)
        }    
        
        ftp.url <- "ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/methylationEPIC"
                
        download.zip(file.path(ftp.url, "infinium-methylationepic-demo-dataset.zip"), path)
    }
    
    path
}
