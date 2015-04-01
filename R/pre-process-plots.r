# ## pre normalisation data

# sample missmatches concordance with genotype data
# - select snps with 99% concordance and then samples with 80% concordance

# gender
# - check predicted sex against manifest files

# methylated vs unmethylated
# - plot raw control probes and fit linear regression, remove samples that have sd(y - yhat) > mean*3

# calculate means for each sample from control probes of each
# - dyebias
# - oob ratio
# - bisulphide conversion
# - spec2red2
# - hybe3
# - remaining 37

# plot all of them - sample against mean ,  colour by plate


# plot detection p values from idat files

# plot number of beads per sample

# manhattan plot of detection pval per probe - percentage with pvalue < 0.01
# manhattan plot of number of beads by probe - percentage of probes with beads < 3 for each sample


# plot scree plot of control matrix

# plot each PC against each of chiprow, chip column, plate, slide IF there is a significant lm


# ## post normalisation data


# use beta values to calculate PCs on 20k most variable probes

# plot each PC against each of chiprow, chip column, plate, slide IF there is a significant lm
# - see if technical variance is gone
















# Sex plots


## Return gender covariate if provided in the covariates
returnGivenGender <- function(covariates){
    givenGender <- NULL
    
    if (!is.null(covariates)){
        possibilities <- c("gender","Gender","sex","Sex","GENDER","SEX")
        sum <- sum(possibilities %in% colnames(covariates))
        if (sum!=0){
            goodColumn <- possibilities[possibilities %in% colnames(covariates)][1]
            givenGender <- as.character(covariates[,goodColumn])
            givenGender <- substr(toupper(givenGender),1,1)
            givenGender <- as.character(givenGender)
        }}
    
    return(givenGender)
}

## Predict gender with X and Y chr intensities
returnPredictedGender <- function(cnQuantiles, cutoff = (-0.3)){
    ## It assumes that cnQuantiles is a list containing $X and $Y
    X <- cnQuantiles$X
    Y <- cnQuantiles$Y
    med <- floor(nrow(X)/2)
    
    x = X[med,]
    y = Y[med,]
    diff <- log2(y)-log2(x)
    n <- length(diff)
    
    predictedGender = rep("M", n)
    predictedGender[which(diff < cutoff)] <- "F"     
    names(predictedGender) <- colnames(X)   
    
    predictedGender <- as.character(predictedGender)
    return(predictedGender)
}

plotPredictedGender <- function(cnQuantiles,cutoff =(-0.3),
                                color = NULL, legend=TRUE, bty="o"){
    X <- cnQuantiles$X
    Y <- cnQuantiles$Y
    med <- floor(nrow(X)/2)
    x <- X[med,]
    y <- Y[med,]
    diff <- log2(y)-log2(x)
    n <- length(diff)
    
    if (is.null(color)){
        color <-  rep("lightskyblue", n)
        color[which(diff < cutoff)] <- "orange"
    }
    
    plot(diff, 
         jitter(rep(0,n),factor=1.5), 
         ylim = c(-1,1), 
         pch=18, 
         cex=2,
         col= color, 
         yaxt="n", 
         xlab="median CN(Y)  - median CN(X)", 
         ylab="", bty=bty
         )            
    abline(v=as.numeric(cutoff),lty=3,lwd=2)
    
    
    if (legend){
        legend("topright",
               c("Predicted male","Predicted female"),
               cex=2, 
               pch=18, 
               col=c("lightskyblue","orange"),
               bty="n"
               )
    }
    
}

plotDiscrepancyGenders <- function(cutoff = (-0.3), covariates,
                                   cnQuantiles, bty="o"){
    predictedGender <- returnPredictedGender(cutoff = cutoff,
                                             cnQuantiles = cnQuantiles)
    givenGender     <- returnGivenGender(covariates)
    ## In the case a gender information was not provided:
    if (is.null(givenGender)){
        plotPredictedGender(cnQuantiles = cnQuantiles, 
                            cutoff = cutoff, 
                            color = NULL,
                            legend = TRUE, bty = bty)
    } else {
    ## In the case a gender information WAS provided:
        X <- cnQuantiles$X
        Y <- cnQuantiles$Y
        med <- floor(nrow(X)/2)
        x = X[med,]
        y = Y[med,]
        diff <- log2(y)-log2(x)
        n <- length(diff)
        
        color <- predictedGender
        color[color == "M"] <- "lightskyblue"
        color[color == "F"] <- "orange"
        
        ## For the samples with different predicted gender:
        non.matching <- predictedGender  != givenGender
        color[non.matching] <- "black"
        
                                        # For the samples with missing values:
        na.positions <- is.na(givenGender)
        color[na.positions] <- "grey"
        
        plotPredictedGender(cnQuantiles = cnQuantiles, 
                            cutoff = cutoff, 
                            color = color,
                            legend = FALSE, bty=bty)
        
        ## To add the sample names in the plot for the discrepancy samples:
        ## if (sum(color=="black" & color!="grey")>=1){
        ## indices <- which(color=="black" & color != "grey")
        ## text(diff[indices],-0.2, names(X)[indices])
        
        ## To add the legend: 
        legend.colors <- c("lightskyblue","orange")
        legend.names  <- c("Predicted male","Predicted female")
        
        ## In the case there are unmatching samples:
        if (sum(color=="black" & color!="grey")>=1){
            legend.colors <- c(legend.colors, "black")
            legend.names  <- c(legend.names, "Unmatching samples")
        }
        
        ## In the case there are missing values in the provided sex covariate:
        if (sum(color=="grey")>=1){
            legend.colors <- c(legend.colors, "grey")
            legend.names  <- c(legend.names, "Sex not provided by user")
        }
        
        legend( "topright",
               legend = legend.names,
               col = legend.colors,
               cex = 1.5, 
               pch = 18,
               bty="n")
    }
}





# shiny methyl plots

plotPCA <- function(pca, pc1, pc2, col, covariates, selectedCov, bty="o"){
    xMin <- min(pca[,as.numeric(pc1)])
    xMax <- max(pca[,as.numeric(pc1)])
    xRange <- xMax - xMin
    xlim <- c(xMin-0.05*xRange, xMax+0.20*xRange)
    xlab <- paste("PC",as.numeric(pc1), " scores", sep="")
    ylab <- paste("PC",as.numeric(pc2), " scores", sep="")
    
    plot(pca[,as.numeric(pc1)], pca[,as.numeric(pc2)],
         col = col, pch = 18, cex = 2, xlab = xlab,
         ylab = ylab, xlim = xlim,
         main = "Principal component analysis (PCA)",
         cex.main = 1.5, cex.lab = 1.5, bty = bty)
    uColor <- unique(col)
    uCov   <- unique(covariates[,match(selectedCov, colnames(covariates))])
    
    legend("bottomright", legend = uCov, pch = 18, col = uColor,
           cex = 1.5, title = selectedCov, bty = "n")
    grid()
}


# plot qc

plotQC <- function(unmethQuantiles, methQuantiles, 
                   sampleNames, col, bty="o"){
    
    slideNames <- substr(sampleNames,1,10)
    
    med <- as.integer(nrow(methQuantiles[[3]])/2)
    x <- log2(unlist(unmethQuantiles[[3]][med,]))
    y <- log2(unlist(methQuantiles[[3]][med,]))
    
    range.x <- max(x)-min(x)
    range.y <- max(y)-min(y)
    xlim1 <- min(x)-0.2*range.x
    xlim2 <- max(x)+0.2*range.x
    ylim1 <- min(y)-0.2*range.y
    ylim2 <- max(y)+0.2*range.y
    
    plot(x, y, xlim = c(xlim1, xlim2), ylim = c(ylim1, ylim2),
         cex = 1, pch = 20, col = col, main = "QC Plot", bty=bty)
    grid()
}

addHoverQC <- function(y, selectedSamples = c(),
                       unmethQuantiles, methQuantiles){
    med <- as.integer(nrow(methQuantiles[[3]])/2)
    mediansU <- unlist(unmethQuantiles[[3]][med,])
    mediansM <- unlist(methQuantiles[[3]][med,])
    
    ## To put a circle around the last entry: 
    n <- length(selectedSamples)
    if (n>=1){
        points(log2(mediansU[selectedSamples[n]]), 
               log2(mediansM[selectedSamples[n]]),
               col = "black", cex=3, pch=1, lwd=2)
    }
    ## Make the points black:
    points(log2(mediansU[selectedSamples]), 
           log2(mediansM[selectedSamples]),
           col = "black", cex = 1, pch = 17)
}


# Some more changes