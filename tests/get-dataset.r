## Retrieve the data from the Gene Expression Omnibus (GEO) website
## (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55491).

dir.create(path <- "data")
if (length(list.files(path, "*.idat$")) == 0) {
  filename <-  file.path(path, "gse55491.tar")
  download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55491&format=file", filename)
  cat(date(), "Extracting files from GEO archive.\n")
  system(paste("cd", path, ";", "tar xvf", basename(filename)))
  unlink(filename)
  cat(date(), "Unzipping IDAT files.\n")
  system(paste("cd", path, ";", "gunzip *.idat.gz"))
}
