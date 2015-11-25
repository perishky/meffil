check.data.frame <- function(x, columns) {
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
