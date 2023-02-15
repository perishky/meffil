# ---- download-epic2-demo-dataset ----

download.epic2.demo.dataset <- function() {
    dir.create(path <- "data-epic2-demo")
    
    if (length(list.files(path, "*.idat$")) == 0) {
        download.zip <- function(url, path) {
            filename <- file.path(path, "data.zip")
            download.file(url, filename)
            filenames <- unzip(filename, junkpaths=T, exdir=path)
            unlink(filename)
            invisible(filenames)
        }    

        url <- "https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/DemoDataEPIC_v2.zip"
        download.zip(url, file.path(path, "data.zip"))
    }
    
    path
}
