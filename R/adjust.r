#' Fit linear models with each column of \code{y}
#' as dependent variables and the fixed and random
#' effects as independent variables.
#' Independent variables lacking variation are omitted.
#' Returns the residual matrix.
adjust.columns <- function(y, fixed.effects=NULL, random.effects=NULL) {
    stopifnot(is.matrix(y))
    stopifnot((is.matrix(fixed.effects) || is.data.frame(fixed.effects)) && nrow(y) == nrow(fixed.effects))
    stopifnot(is.null(random.effects) || nrow(y) == nrow(random.effects))
    
    remove.invariant.columns <- function(x) {
        if (is.null(x)) return(x)
        is.variable <- apply(x, 2, function(x) length(unique(x))) > 1
        var.idx <- which(is.variable)
        if (length(var.idx) == 0) NULL
        else x[,var.idx,drop=F]
    }

    fixed.effects <- remove.invariant.columns(fixed.effects)
    random.effects <- remove.invariant.columns(random.effects)

    if (is.null(fixed.effects) && is.null(random.effects)) return(y)

    if (is.null(random.effects)) {
        if (is.data.frame(fixed.effects))
            fixed.effects <- do.call(cbind, lapply(fixed.effects, simplify.variable))        
        fit <- lm.fit(x=fixed.effects, y=y)
        return(residuals(fit))
    }
    
    data <- data.frame(random.effects, stringsAsFactors=F)
    formula <- "y ~"
    if (!is.null(fixed.effects)) {
        formula <- paste(formula, paste(colnames(fixed.effects), collapse=" + "), "+")
        data <- data.frame(data, fixed.effects, stringsAsFactors=F)
    }
    formula <- paste(formula, paste("(1 |", colnames(random.effects), ")", collapse=" + "))
    
    ret <- sapply(1:ncol(y), function(i) {
        data$y <- y[,i]
        tryCatch({
            residuals(lme4::lmer(formula, data=data))
            ## chose lme4 because it is faster than nlme
        }, error=function(e) {
            print(e)
            cat("For variable", i, "ignoring random effects.\n")
            tryCatch({
                residuals(lm(y ~ ., data=data))
            }, error=function(e) {
                print(e)
                cat("For variable", i, "setting all values to missing.\n")
                rep(NA, nrow(data))
            })
        })
    })
    dimnames(ret) <- dimnames(y)
    ret
}
