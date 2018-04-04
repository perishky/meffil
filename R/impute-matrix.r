# to impute missing values with row means
# x <- impute.matrix(x,1)
# to impute missing values with column means
# x <- impute.matrix(x,2)
impute.matrix <- function(x, margin=1, fun=function(x) mean(x, na.rm=T)) {
    if (margin == 2) x <- t(x)
    
    idx <- which(is.na(x) | !is.finite(x), arr.ind=T)
    if (length(idx) > 0) {
        na.idx <- unique(idx[,"row"])
        v <- apply(x[na.idx,,drop=F],1,fun) ## v = summary for each row
        v[which(is.na(v))] <- fun(v)      ## if v[i] is NA, v[i] = fun(v)
        x[idx] <- v[match(idx[,"row"],na.idx)] ##
        stopifnot(all(!is.na(x)))
    }

    if (margin == 2) x <- t(x)
    x
}
