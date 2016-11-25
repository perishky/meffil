#' Convert a variable to a format that can be included in a design matrix.
#' Characters and unordered factors are converted to model matrices with binary dummy variables.
#' Ordered factors are converted to integers corresponding to ordered levels.
#' Other variables are left asis.
simplify.variable <- function(v) {
    if (is.character(v))
        v <- as.factor(v)
    if (is.factor(v)) {
        v <- droplevels(v)
        if (is.ordered(v) || length(levels(v)) <= 2)
            as.integer(v) - 1
        else
            model.matrix(~ v, model.frame(~ v, na.action=na.pass))[,-1,drop=F]
    } else
        v
}
