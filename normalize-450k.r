
meffil.control.matrices <- function(basenames, verbose=T) {

    probe.info <- function(type) {
        require(IlluminaHumanMethylation450kmanifest)
        getProbeInfo(IlluminaHumanMethylation450kmanifest, type =type)
    }

    probe.names <- function(type=c("I","II")) {
        if (missing(type))
            c(probe.info("I")$Name, probe.info("II")$Name)
        else {
            type <- match.arg(type)
            probe.info(type)$Name
        }   
    }

    probe.locations <- function(array="IlluminaHumanMethylation450k",annotation="ilmn12.hg19") {
        annotation <- paste(array, "anno.", annotation, sep="")
        require(annotation,character.only=T)
        data(list=annotation)
        as.data.frame(get(annotation)@data$Locations)
    }

    read.idat <- function(filename, verbose=T) {
        msg("Reading", filename, verbose=verbose)
        
        if (!file.exists(filename))
            stop("Filename does not exist:", filename)
        readIDAT(filename)$Quants[,"Mean"]
    }

    probe.types <- function() {
        type1.R <- probe.info("I-Red")
        type1.G <- probe.info("I-Green")
        type2 <- probe.info("II")
        controls <- probe.info("Control")

        ret <- rbind(data.frame(type="M", dye="R", address=type1.R$AddressB, name=type1.R$Name),
                     data.frame(type="M", dye="G", address=type1.G$AddressB, name=type1.G$Name),
                     data.frame(type="M", dye="G", address=type2$AddressA, name=type2$Name),
                     
                     data.frame(type="U", dye="R", address=type1.R$AddressA, name=type1.R$Name),
                     data.frame(type="U", dye="G", address=type1.G$AddressA, name=type1.G$Name),
                     data.frame(type="U", dye="R", address=type2$AddressA, name=type2$Name),
                     
                     data.frame(type="OOB", dye="G", address=type1.R$AddressA, name=NA),
                     data.frame(type="OOB", dye="G", address=type1.R$AddressB, name=NA),
                     data.frame(type="OOB", dye="R", address=type1.G$AddressA, name=NA),
                     data.frame(type="OOB", dye="R", address=type1.G$AddressB, name=NA),
                     data.frame(type=controls$Type,dye="R",address=controls$Address, name=NA),
                     data.frame(type=controls$Type,dye="G",address=controls$Address, name=NA))
        
        for (col in colnames(ret)) ret[,col] <- as.character(ret[,col])
        ret
    }
    
    rg.to.mu <- function(rg) {
        probes <- probe.types()
        probes.M.R <- probes[which(probes$type == "M" & probes$dye == "R"),]
        probes.M.G <- probes[which(probes$type == "M" & probes$dye == "G"),]
        probes.U.R <- probes[which(probes$type == "U" & probes$dye == "R"),]
        probes.U.G <- probes[which(probes$type == "U" & probes$dye == "G"),]

        M <- c(rg$R[probes.M.R$address], rg$G[probes.M.G$address])
        U <- c(rg$R[probes.U.R$address], rg$G[probes.U.G$address])

        names(M) <- c(probes.M.R$Name, probes.M.G$Name)
        names(U) <- c(probes.U.R$Name, probes.U.G$Name)

        U <- U[names(M)]
        list(M=M,U=U)
    }

    background.correct <- function(rg, offset=15) {
        require(MASS)
        require(limma)

        probes <- probe.types()
        
        lapply(c(R="R",G="G"), function(dye) {
            addresses <- probes$address[which(probes$type != "OOB" & probes$dye == dye)]
            xf <- rg[[dye]][addresses]
            xf[which(xf <= 0)] <- 1

            addresses <- probes$address[which(probes$type == "OOB" & probes$dye == dye)]
            oob <- rg[[dye]][addresses]
            
            ests <- MASS::huber(oob) 
            mu <- ests$mu
            sigma <- log(ests$s)
            alpha <- log(max(huber(xf)$mu - mu, 10))
            xf.bkg <- limma::normexp.signal(as.numeric(mu,sigma,alpha), xf) + offset
            names(xf.bkg) <- names(xf)
            xf.bkg
        })
    }

    compute.dye.ratio <- function(rg) {
        probes <- probes.types()
        cg.controls <- probes$address[which(probes$type %in%  c("NORM_C", "NORM_G"))]
        at.controls <- probes$address[which(probes$type %in%  c("NORM_A", "NORM_T"))]
        mean(rg$R[at.controls])/mean(rg$G[cg.controls])
    }


    require(illuminaio)

    basenames <- sub("_Grn\\.idat$", "", basenames)
    basenames <- sub("_Red\\.idat$", "", basenames)

    sapply(basenames, function(basename) {
        rg <- list(G=read.idat(paste(basename, "_Grn.idat", sep = ""), verbose=T),
                   R=read.idat(paste(basename, "_Red.idat", sep=""), verbose=T))
        rg <- background.correct(rg)
        dye.ratio <- compute.dye.ratio(rg)
        mu <- rg.to.mu(rg)

        
        ## extract control matrix from R,G
        ## compute 1,50,99 percentiles of oob R,G
        ## compute median of X chromosome, median of Y chromosome copy numbers = log2(M + U) -- probe.locations
        ## compute 500 quantiles of mu$M and mu$U for
        ##   (type1.R,type1.G,type2) times (autosomal,X,Y)
    })
}



