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















## Latest from Josine

#outlier detection, PCA and lm and plots controlmatrix

library(meffil)

#samples in rows, batches in columns
manifestdata<-read.table('manifestData_genotypeQC2.txt',sep='\t',header=TRUE)

#load output from "meffil.compute.normalization.object"
load("ARIES_funnorm.Robj")
length(norm.objects)

#get sample id and sort manifest file 
basename<-unlist(lapply(norm.objects,function(x) x$basename))
m<-match(basename,manifestdata$sentrix)
manifestdata<-manifestdata[m,]

#outlier detection

cols=rainbow(length(unique(manifestdata$BCD_plate)))
df2<-data.frame(unique(manifestdata$BCD_plate),cols)
m<-match(manifestdata$BCD_plate,df2[,1])
manifestdata<-data.frame(manifestdata,platecols=df2[m,2])

cols=rainbow(length(unique(manifestdata$Slide)))
df2<-data.frame(unique(manifestdata$Slide),cols)
m<-match(manifestdata$Slide,df2[,1])
manifestdata<-data.frame(manifestdata,slidecols=df2[m,2])

o<-order(manifestdata$Slide)
manifestdata<-manifestdata[o,]
o<-order(manifestdata$BCD_plate)
manifestdata<-manifestdata[o,]

###
control.matrix <- meffil.control.matrix(norm.objects)
control.matrix<-t(control.matrix)

outlier.out<-data.frame()
pdf("controlprobe_type_coloredbyplate.pdf",width=8,height=8)
par(mfrow = c(2, 2))
for (controls in 1:dim(control.matrix)[2]){

o1<-(mean(control.matrix[,controls]))+(5*sd(control.matrix[,controls]))
o2<-(mean(control.matrix[,controls]))-(5*sd(control.matrix[,controls]))
outlier<-which(control.matrix[,controls]>o1|control.matrix[,controls]<o2)
outlier
if(length(outlier)>0) {outlier.out<-rbind(outlier.out,data.frame(manifestdata[outlier,"sentrix"],controlType=rep(colnames(control.matrix)[controls],length(outlier)),control.matrix[outlier,controls]))}

plot(control.matrix[,controls],main=colnames(control.matrix)[controls],xlab="sample",ylab="signal",col=manifestdata$platecols,cex=0.8,pch=16,cex.main=0.8,cex.axis=0.8)
abline(h=o1,col="black")
abline(h=o2,col="black")
}

dev.off()

o<-order(manifestdata$Slide)
manifestdata<-manifestdata[o,]

pdf("controlprobe_type_coloredbyslide.pdf",width=8,height=8)
par(mfrow = c(2, 2))
for (controls in 1:dim(control.matrix)[2]){
plot(control.matrix[,controls],main=colnames(control.matrix)[controls],xlab="sample",ylab="signal",col=manifestdata$slidecols,cex=0.8,pch=16,cex.main=0.8,cex.axis=0.8)
abline(h=o1,col="black")
abline(h=o2,col="black")
}

dev.off()

m1<-manifestdata[manifestdata$sentrix%in%outlier.out[,1],]
m1<-merge(outlier.out[,1:2],m1,by.x=1,by.y="sentrix",all.x=T)
write.table(m1,"5sdoutlier_controlmatrix.txt",sep="\t",col.names=T,row.names=F,quote=F)


###variance explained controlmatrix
load("ARIES_funnorm.Robj")
length(norm.objects)
#get sample id
basename<-unlist(lapply(norm.objects,function(x) x$basename))
control.matrix <- meffil.control.matrix(norm.objects)
number.pcs=dim(control.matrix)[1]

#
control.matrix <- meffil.control.matrix(norm.objects)
control.components <- prcomp(t(control.matrix))
prop = control.components$sdev^2/sum(control.components$sdev^2)
#cumprop<-cumsum((control.components$sdev)^2) / sum(control.components$sdev^2)

pdf("variance_explained_controlprobes_meffil.pdf",height=6,width=6)
barplot(prop*100,main="variance explained",ylab="percentage of variance",xlab="PC")
screeplot(control.components, npcs = number.pcs, type = "lines")
dev.off()
####
#PCA and lm plots on controlmatrix

plotPCA <- function(pca, pc1, pc2, col, covariates, selectedCov, bty="o"){
    pheno.col=as.numeric(as.factor(covariates[,match(selectedCov,colnames(covariates))]))
    colors <- c("#4477AA","#CC6677","#DDCC77","#117733","#88CCEE","#AA4499","#44AA99","#999933","#882255","#661100","#6699CC","#AA4466")
    if(length(unique(pheno.col))>length(colors)){colors=rainbow(length(unique(pheno.col)))}
    
    u.pheno.col<-data.frame(unique(pheno.col),colors[1:length(unique(pheno.col))])
    m<-match(pheno.col,u.pheno.col[,1])
    col<-u.pheno.col[m,2]
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
         cex.main = 1.0, cex.lab = 1.0, bty = bty)
    
    uColor <- unique(col)
    uCov   <- unique(covariates[,match(selectedCov, colnames(covariates))])
    
    legend("topright", legend = uCov, pch = 18, col = uColor,
           cex = 0.8, title = selectedCov, bty = "n")
    grid()
}

manifestdata<-read.table('manifestData_genotypeQC2.txt',sep='\t',header=TRUE)
manifestdata$BCD_plate<-gsub("-Corrected","",manifestdata$BCD_plate)
manifestdata$BCD_plate<-gsub("_Corrected","",manifestdata$BCD_plate)
manifestdata$BCD_plate<-gsub("CORR","",manifestdata$BCD_plate)
manifestdata$BCD_plate<-gsub("_C","",manifestdata$BCD_plate)
manifestdata$BCD_plate<-gsub("Pre","",manifestdata$BCD_plate)

load("ARIES_funnorm.Robj")
length(norm.objects)
#get sample id
basename<-unlist(lapply(norm.objects,function(x) x$basename))
m<-match(basename,manifestdata$sentrix)
manifestdata<-manifestdata[m,]
control.matrix <- meffil.control.matrix(norm.objects)
control.components <- prcomp(t(control.matrix))

batch<-c("sample_type2","Slide","BCD_plate","time_point","time_code","chiprow","chipcolumn","gender")

pdf(paste("PCAplot_controlmatrix.pdf",sep=""),height=6,width=6)
for (b in 1:length(batch)){
for (pcomp in 2:dim(control.matrix)[1]){
plotPCA(pca=control.components$x,pc1=1,pc2=pcomp,bty="o",covariates=manifestdata,selectedCov=batch[b],col=c("grey"))
}}
dev.off()

pdf(paste("pca.controlmatrix_lm.pdf",sep=""),height=6,width=6)
for (b in 1:length(batch)){
m<-match(batch[b],names(manifestdata))
manifestdata[,m]<-as.factor(manifestdata[,m])
batch.out<-data.frame()    
for (pcomp in 1:dim(control.matrix)[1]){
options(contrasts=c("contr.sum","contr.poly"))
r1<-summary(lm(control.components$x[,pcomp]~1+manifestdata[,m]))
batch.out<-rbind(batch.out,data.frame(pc=pcomp,f=levels(manifestdata[,m]),r1$coefficients))
pca.pc<-control.components$x[,pcomp]
boxplot(pca.pc~manifestdata[,m],main=paste("PC",pcomp,".",batch[b],sep=""),cex.names=0.5,las=2)
write.table(batch.out,paste("controlmatrix_pca_lm.",batch[b],".txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
}
}
dev.off()



#PCA and lm plots normalized data; use of 20000 mot variable probes only
library(meffil)
library(matrixStats)
pc=2
load(paste("ARIES_funnorm.pc",pc,".Robj",sep=""))

y<-meffil.probe.info(array = "IlluminaHumanMethylation450k",annotation = "ilmn12.hg19")
m<-match(row.names(norm.data),y$name)
y<-y[m,] #485512
autMatrix<-norm.data[which(y$chr!="chrX"&y$chr!="chrY"),]

#[1] 473864
  numPositions = 20000
             
              o <- order(-rowVars(autMatrix))[1:numPositions]
              pca <- prcomp(t(autMatrix[o,]))
prop = pca$sdev^2/sum(pca$sdev^2)
cumprop<-cumsum((pca$sdev)^2) / sum(pca$sdev^2)

pdf(paste("variance_explained_funnorm_meffil_pc",pc,".pdf",sep=""),height=6,width=6)
barplot(prop[1:10]*100,main="variance explained",ylab="percentage of variance",xlab="PC")
screeplot(pca, npcs = 10, type = "lines")
dev.off()

save(pca,file=paste("pca.pc",pc,".Robj",sep=""))            

plotPCA <- function(pca, pc1, pc2, col, covariates, selectedCov, bty="o"){
    pheno.col=as.numeric(as.factor(covariates[,match(selectedCov,colnames(covariates))]))
    colors <- c("#4477AA","#CC6677","#DDCC77","#117733","#88CCEE","#AA4499","#44AA99","#999933","#882255","#661100","#6699CC","#AA4466")
    if(length(unique(pheno.col))>length(colors)){colors=rainbow(length(unique(pheno.col)))}
    
    u.pheno.col<-data.frame(unique(pheno.col),colors[1:length(unique(pheno.col))])
    m<-match(pheno.col,u.pheno.col[,1])
    col<-u.pheno.col[m,2]
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
         cex.main = 1.0, cex.lab = 1.0, bty = bty)
    
    uColor <- unique(col)
    uCov   <- unique(covariates[,match(selectedCov, colnames(covariates))])
    
    legend("topright", legend = uCov, pch = 18, col = uColor,
           cex = 0.8, title = selectedCov, bty = "n")
    grid()
}

manifestdata<-read.table('manifestData_genotypeQC2.txt',sep='\t',header=TRUE)
number.pcs=10
manifestdata$BCD_plate<-gsub("-Corrected","",manifestdata$BCD_plate)
manifestdata$BCD_plate<-gsub("_Corrected","",manifestdata$BCD_plate)
manifestdata$BCD_plate<-gsub("CORR","",manifestdata$BCD_plate)
manifestdata$BCD_plate<-gsub("_C","",manifestdata$BCD_plate)
manifestdata$BCD_plate<-gsub("Pre","",manifestdata$BCD_plate)

basename<-row.names(pca$x)
m<-match(basename,manifestdata$sentrix)
manifestdata<-manifestdata[m,]

batch<-c("sample_type2","Slide","BCD_plate","time_point","time_code","chiprow","chipcolumn","gender")

pdf(paste("PCAplot_funnorm_meffil_pc",pc,".pdf",sep=""),height=6,width=6)
for (b in 1:length(batch)){
for (pcomp in 2:number.pcs){
plotPCA(pca=pca$x,pc1=1,pc2=pcomp,bty="o",covariates=manifestdata,selectedCov=batch[b],col=c("grey"))
}}
dev.off()

pdf(paste("pca.funnorm_lm_pc",pc,".pdf",sep=""),height=6,width=6)
for (b in 1:length(batch)){
m<-match(batch[b],names(manifestdata))
manifestdata[,m]<-as.factor(manifestdata[,m])
batch.out<-data.frame()    
for (pcomp in 1:number.pcs){
options(contrasts=c("contr.sum","contr.poly"))
r1<-summary(lm(pca$x[,pcomp]~1+manifestdata[,m]))
batch.out<-rbind(batch.out,data.frame(pc=pc,f=levels(manifestdata[,m]),r1$coefficients))
pca.pc<-pca$x[,pcomp]
boxplot(pca.pc~manifestdata[,m],main=paste("PC",pcomp,".",batch[b],sep=""),cex.names=0.5,las=2)
write.table(batch.out,paste("funnorm_pc",pc,"_lm",batch[b],".txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
}
}
dev.off()










# sandpit


dir.create(path <- "data")

if (length(list.files(path, "*.idat$")) == 0) {
  filename <-  file.path(path, "gse55491.tar")
  download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55491&format=file", filename)
  cat(date(), "Extracting files from GEO archive.\n")
  system(paste("cd", path, ";", "tar xvf", basename(filename)))
  unlink(filename)
  cat(date(), "Unzipping IDAT files.\n")
  system(paste("cd", path, ";", "gunzip *.idat.gz"))
}





# Plots



