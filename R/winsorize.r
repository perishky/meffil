winsorize <- function(x, pct=0.05) {
    if (is.matrix(x)) {
        quantiles <- rowQuantiles(x, probs=c(pct,1-pct), na.rm=T)
        low <- quantiles[,1]
        high <- quantiles[,2]
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

