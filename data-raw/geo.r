
get.substring <- function(x, pattern, type) {
    value <- sub(pattern, "\\1", as.character(x))
    if (!missing(type))
        as(value, type)
    else
        value
}

get.sentrix <- function(x) {
    get.substring(x, "[^_]+_([^_]+_[^_]+)_.+")
}

get.characteristic <- function(x, name, type="character") {
    name <- gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", name))
    get.substring(x, paste0(".*", name, ": ([^;]+).*"), type)
}

get.characteristic1 <- function(x, ...)
    get.characteristic(x$characteristics_ch1, ...)

merge.columns <- function(x, column) {
    idx <- which(colnames(x) == column)
    if (length(idx) <= 1) return(x)
    new.col <- apply(x[,idx], 1, paste, collapse=";")
    x[,idx[1]] <- as.character(new.col)
    x <- x[,-idx[-1],drop=F]
    x
}

geo.matrix.url <- function(gse) {
    file.path("ftp://ftp.ncbi.nlm.nih.gov/geo/series",
              paste(substring(gse, 1, nchar(gse)-3), "nnn", sep=""),
              gse,
              "matrix",
              paste(gse, "series", "matrix.txt.gz", sep="_"))
}

readLinesUntil <- function(pattern, con, ...) {
    n <- 10
    lines <- c()
    finished <- F
    while (!finished) {
        new <- readLines(con=con, n=n, ...)
        lines <- c(lines, new)
        finished <- (length(grep(pattern, new)) > 0
                     || length(new) < n)
        n <- 2*n
    }
    lines
}

geo.samples <- function(gse) {
    con <- gzcon(url(geo.matrix.url(gse)))
    open(con)
    on.exit(close(con))
    lines <- readLinesUntil("^!series_matrix_table_begin", con)

    idx <- grep("^!Sample_", lines)
    samples <- read.table(file=textConnection(lines[idx]), sep="\t", header=T)
    
    column.names <- sub("!Sample_", "", samples[,1])

    samples <- samples[,-1]
    samples <- t(samples)
    colnames(samples) <- column.names
    samples <- as.data.frame(samples)
    samples <- merge.columns(samples, "characteristics_ch1")
    samples <- merge.columns(samples, "supplementary_file")
    samples
}

