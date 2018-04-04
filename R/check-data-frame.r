# check that data frame \code{x} contains the required
# columns and types or contents.
# \code{columns} is a list whose names correspond to column names
# and the values of each item is either a type (e.g. "character", "numeric")
# or a vector of possible values.
check.data.frame <- function(x, columns) {
    stopifnot(is.data.frame(x))
    missing.columns <- setdiff(names(columns), colnames(x))
    if (length(missing.columns) > 0)
        stop(paste("Missing columns:", paste(missing.columns, collapse=", ")))
    
    for (column in names(columns)) {
        if (length(columns[[column]]) == 1) {
            if (!is(x[[column]], columns[[column]]))
                stop(paste("Column", column, "is not of type", columns[[column]]))
        }
        else {
            if (!all(na.omit(x[[column]]) %in% columns[[column]]))
                stop(paste("Column", column, "valid values:",
                           paste(columns[[column]], collapse=",")))
        }
    }
}
