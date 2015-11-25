winsorize <- function(x, pct=0.05) {
    if (is.matrix(x)) {
        if (any(is.na(x)))
            x <- impute.matrix(x)
        low <- rowQ(x, which=floor(ncol(x)*pct)) 
        high <- rowQ(x, which=floor(ncol(x)*(1-pct))) 
        idx <- which(x < low, arr.ind=T)
        x[idx] <- low[idx[,1]]
        idx <- which(x > high, arr.ind=T)
        x[idx] <- high[idx[,1]]
    }
    else {
        low <- quantile(x, probs=pct, na.rm=T)
        high <- quantile(x, probs=1-pct, na.rm=T)
        idx <- which(x < low)
        x[idx] <- low
        idx <- which(x > high)
        x[idx] <- high
    }        
    x
}    

