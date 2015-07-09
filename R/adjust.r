adjust.columns <- function(y, fixed.effects=NULL, random.effects=NULL) {
    stopifnot(nrow(y) == nrow(fixed.effects))
    stopifnot(is.null(random.effects) || nrow(y) == nrow(random.effects))

    if (is.null(random.effects)) {
        if (is.null(fixed.effects)) return(y)
        fit <- lm.fit(x=fixed.effects, y=y)
        ret <- residuals(fit)
    } else {
        formula <- "y ~"
        if (!is.null(fixed.effects)) 
            formula <- paste(formula, paste(colnames(fixed.effects), collapse=" + "), "+")
        formula <- paste(formula, paste("(1 |", colnames(random.effects), ")", collapse=" + "))

        if (is.null(fixed.effects))
            data <- data.frame(random.effects, stringsAsFactors=F)
        else
            data <- data.frame(fixed.effects, random.effects, stringsAsFactors=F)
        
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
    }
    ret
}
