msg <- function(..., verbose=T) {
    if (verbose) {
        x <- paste(list(...))
        name <- sys.call(sys.parent(1))[[1]]
        cat(paste("[", name, "]", sep=""), date(), x, "\n")
    }
}







