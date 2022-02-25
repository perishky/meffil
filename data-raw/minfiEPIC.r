extractFromRGSetEPIC <- function (rgSet) {
    rgSet <- updateObject(rgSet)
    controlType <- c(
        "BISULFITE CONVERSION I", "BISULFITE CONVERSION II", 
        "EXTENSION", "HYBRIDIZATION", "NEGATIVE", "NON-POLYMORPHIC", 
        "NORM_A", "NORM_C", "NORM_G", "NORM_T", "SPECIFICITY I",
        ## "SPECIFICITY II",
        "TARGET REMOVAL", "STAINING")
    array <- annotation(rgSet)[["array"]]
    ctrls <- getProbeInfo(rgSet, type = "Control")
    if (!all(controlType %in% ctrls$Type)) 
        stop("The `rgSet` does not contain all necessary control probes")
    ctrlsList <- split(ctrls, ctrls$Type)[controlType]
    redControls <- getRed(rgSet)[ctrls$Address, , drop = FALSE]
    redControls <- lapply(
        ctrlsList,
        function(ctl) redControls[ctl$Address, 
                                  , drop = FALSE])
    greenControls <- getGreen(rgSet)[ctrls$Address, , drop = FALSE]
    greenControls <- lapply(
        ctrlsList,
        function(ctl) greenControls[ctl$Address, 
                                    , drop = FALSE])
    oobRaw <- getOOB(rgSet)
    probs <- c(0.01, 0.5, 0.99)
    greenOOB <- t(colQuantiles(oobRaw$Grn, na.rm = TRUE, probs = probs))
    redOOB <- t(colQuantiles(oobRaw$Red, na.rm = TRUE, probs = probs))
    oob <- list(greenOOB = greenOOB, redOOB = redOOB)
    return(list(greenControls = greenControls, redControls = redControls, 
                oob = oob, ctrlsList = ctrlsList, array = array))
}
buildControlMatrixEPIC <- function (extractedData) {
    getCtrlsAddr <- function(exType, index) {
        ctrls <- ctrlsList[[index]]
        addr <- ctrls$Address
        names(addr) <- ctrls$ExtendedType
        na.omit(addr[exType])
    }
    array <- extractedData$array
    greenControls <- extractedData$greenControls
    redControls <- extractedData$redControls
    controlNames <- names(greenControls)
    ctrlsList <- extractedData$ctrlsList
    index <- match("BISULFITE CONVERSION II", controlNames)
    redControls.current <- redControls[[index]]
    bisulfite2 <- colMeans2(redControls.current, na.rm = TRUE)
    index <- match("BISULFITE CONVERSION I", controlNames)
    if (array == "IlluminaHumanMethylation450k") {
        addr <- getCtrlsAddr(exType = sprintf("BS Conversion I%sC%s", 
                                 c(" ", "-", "-"), 1:3), index = index)
    }
    else {
        addr <- getCtrlsAddr(exType = sprintf("BS Conversion I%sC%s", 
                                 c("-", "-"), 1:2), index = index)
    }
    greenControls.current <- greenControls[[index]][addr, , drop = FALSE]
    if (array == "IlluminaHumanMethylation450k") {
        addr <- getCtrlsAddr(exType = sprintf("BS Conversion I-C%s", 
                                 4:6), index = index)
    }
    else {
        addr <- getCtrlsAddr(exType = sprintf("BS Conversion I-C%s", 
                                 3:5), index = index)
    }
    redControls.current <- redControls[[index]][addr, , drop = FALSE]
    if (nrow(redControls.current) == nrow(greenControls.current)) {
        bisulfite1 <- colMeans2(redControls.current
                                + greenControls.current, 
                                na.rm = TRUE)
    }
    else {
        bisulfite1 <- colMeans2(redControls.current, na.rm = TRUE) + 
            colMeans2(greenControls.current, na.rm = TRUE)
    }
    index <- match("STAINING", controlNames)
    addr <- getCtrlsAddr(exType = "Biotin (High)", index = index)
    stain.green <- t(greenControls[[index]][addr, , drop = FALSE])
    addr <- getCtrlsAddr(exType = "DNP (High)", index = index)
    stain.red <- t(redControls[[index]][addr, , drop = FALSE])
    index <- match("EXTENSION", controlNames)
    addr <- getCtrlsAddr(exType = sprintf("Extension (%s)", c("A", 
                             "T")), index = index)
    extension.red <- t(redControls[[index]][addr, , drop = FALSE])
    colnames(extension.red) <- paste0("extRed", 1:ncol(extension.red))
    addr <- getCtrlsAddr(exType = sprintf("Extension (%s)", c("C", 
                             "G")), index = index)
    extension.green <- t(greenControls[[index]][addr, , drop = FALSE])
    colnames(extension.green) <- paste0("extGrn", 1:ncol(extension.green))
    index <- match("HYBRIDIZATION", controlNames)
    hybe <- t(greenControls[[index]])
    colnames(hybe) <- paste0("hybe", 1:ncol(hybe))
    index <- match("TARGET REMOVAL", controlNames)
    targetrem <- t(greenControls[[index]])
    colnames(targetrem) <- paste0("targetrem", 1:ncol(targetrem))
    index <- match("NON-POLYMORPHIC", controlNames)
    addr <- getCtrlsAddr(exType = sprintf("NP (%s)", c("A", "T")), 
                         index = index)
    nonpoly.red <- t(redControls[[index]][addr, , drop = FALSE])
    colnames(nonpoly.red) <- paste0("nonpolyRed", 1:ncol(nonpoly.red))
    addr <- getCtrlsAddr(exType = sprintf("NP (%s)", c("C", "G")), 
                         index = index)
    nonpoly.green <- t(greenControls[[index]][addr, , drop = FALSE])
    colnames(nonpoly.green) <- paste0("nonpolyGrn", 1:ncol(nonpoly.green))
    ##index <- match("SPECIFICITY II", controlNames)
    ##greenControls.current <- greenControls[[index]]
    ##redControls.current <- redControls[[index]]
    ##spec2.green <- t(greenControls.current)
    ##colnames(spec2.green) <- paste0("spec2Grn", 1:ncol(spec2.green))
    ##spec2.red <- t(redControls.current)
    ##colnames(spec2.red) <- paste0("spec2Red", 1:ncol(spec2.red))
    ##spec2.ratio <- (colMeans2(greenControls.current, na.rm = TRUE)
    ##                /colMeans2(redControls.current, na.rm = TRUE))
    index <- match("SPECIFICITY I", controlNames)
    addr <- getCtrlsAddr(
        exType = sprintf("GT Mismatch %s (PM)", 
            1:3), index = index)
    greenControls.current <- greenControls[[index]][addr, , drop = FALSE]
    redControls.current <- redControls[[index]][addr, , drop = FALSE]
    spec1.green <- t(greenControls.current)
    colnames(spec1.green) <- paste0("spec1Grn", 1:ncol(spec1.green))
    spec1.ratio1 <- (colMeans2(redControls.current, na.rm = TRUE)
                     /colMeans2(greenControls.current, na.rm = TRUE))
    index <- match("SPECIFICITY I", controlNames)
    addr <- getCtrlsAddr(exType = sprintf("GT Mismatch %s (PM)", 
                             4:6), index = index)
    greenControls.current <- greenControls[[index]][addr, , drop = FALSE]
    redControls.current <- redControls[[index]][addr, , drop = FALSE]
    spec1.red <- t(redControls.current)
    colnames(spec1.red) <- paste0("spec1Red", 1:ncol(spec1.red))
    spec1.ratio2 <- (colMeans2(greenControls.current, na.rm = TRUE)
                     /colMeans2(redControls.current, na.rm = TRUE))
    spec1.ratio <- (spec1.ratio1 + spec1.ratio2)/2
    index <- match(c("NORM_A"), controlNames)
    normA <- colMeans2(redControls[[index]], na.rm = TRUE)
    index <- match(c("NORM_T"), controlNames)
    normT <- colMeans2(redControls[[index]], na.rm = TRUE)
    index <- match(c("NORM_C"), controlNames)
    normC <- colMeans2(greenControls[[index]], na.rm = TRUE)
    index <- match(c("NORM_G"), controlNames)
    normG <- colMeans2(greenControls[[index]], na.rm = TRUE)
    dyebias <- (normC + normG)/(normA + normT)
    oobG <- extractedData$oob$greenOOB
    oobR <- extractedData$oob$redOOB
    oob.ratio <- oobG[2, ]/oobR[2, ]
    oobG <- t(oobG)
    colnames(oobG) <- paste0("oob", c(1, 50, 99))
    model.matrix <- cbind(
        bisulfite1, bisulfite2, extension.green, 
        extension.red, hybe, stain.green, stain.red, nonpoly.green, 
        nonpoly.red, targetrem, spec1.green, spec1.red,
        ##spec2.green, spec2.red,
        spec1.ratio1, spec1.ratio,
        ##spec2.ratio,
        spec1.ratio2, 
        normA, normC, normT, normG, dyebias, oobG, oob.ratio)
    for (colindex in 1:ncol(model.matrix)) {
        if (any(is.na(model.matrix[, colindex]))) {
            column <- model.matrix[, colindex]
            column[is.na(column)] <- mean(column, na.rm = TRUE)
            model.matrix[, colindex] <- column
        }
    }
    model.matrix <- scale(model.matrix)
    model.matrix[model.matrix > 3] <- 3
    model.matrix[model.matrix < (-3)] <- -3
    model.matrix <- scale(model.matrix)
    return(model.matrix)
}    
normalizeFunnormEPIC <- function (
    object, extractedData, nPCs, sex, verbose = TRUE){
    normalizeQuantiles <- function(matrix, indices, sex = NULL) {
        matrix <- matrix[indices, , drop = FALSE]
        oldQuantiles <- t(colQuantiles(matrix, probs = probs))
        if (is.null(sex)) {
            newQuantiles <- minfi:::.returnFit(
                controlMatrix = model.matrix, 
                quantiles = oldQuantiles, nPCs = nPCs)
        }
        else {
            newQuantiles <- minfi:::.returnFitBySex(
                controlMatrix = model.matrix, 
                quantiles = oldQuantiles, nPCs = nPCs, sex = sex)
        }
        minfi:::.normalizeMatrix(matrix, newQuantiles)
    }
    indicesList <- minfi:::.getFunnormIndices(object)
    ##model.matrix <- .buildControlMatrix450k(extractedData)
    model.matrix <- buildControlMatrixEPIC(extractedData)
    probs <- seq(from = 0, to = 1, length.out = 500)
    Meth <- getMeth(object)
    Unmeth <- getUnmeth(object)
    if (nPCs > 0) {
        for (type in c("IGrn", "IRed", "II")) {
            indices <- indicesList[[type]]
            if (length(indices) > 0) {
                if (verbose) 
                    message(sprintf("[normalizeFunnormEPIC] Normalization of the %s probes", 
                                    type))
                Unmeth[indices, ] <- normalizeQuantiles(
                    Unmeth, 
                    indices = indices, sex = NULL)
                Meth[indices, ] <- normalizeQuantiles(
                    Meth, indices = indices, 
                    sex = NULL)
            }
        }
        indices <- indicesList[["X"]]
        if (length(indices) > 0) {
            if (verbose) 
                message("[normalizeFunnormEPIC] Normalization of the X-chromosome")
            Unmeth[indices, ] <- normalizeQuantiles(
                Unmeth, indices = indices, 
                sex = sex)
            Meth[indices, ] <- normalizeQuantiles(
                Meth, indices = indices, 
                sex = sex)
        }
    }
    indices <- indicesList[["Y"]]
    if (length(indices) > 0) {
        if (verbose) 
            message("[normalizeFunnormEPIC] Normalization of the Y-chromosome")
        sex <- as.character(sex)
        levels <- unique(sex)
        nSexes <- length(levels)
        if (nSexes == 2) {
            level1 <- levels[1]
            level2 <- levels[2]
        }
        if (nSexes == 2) {
            if (sum(sex == level1) > 1) {
                Meth[indices, sex == level1] <- preprocessCore::normalize.quantiles(Meth[indices, 
                                  sex == level1, drop = FALSE])
                Unmeth[indices, sex == level1] <- preprocessCore::normalize.quantiles(Unmeth[indices, 
                                    sex == level1, drop = FALSE])
            }
            if (sum(sex == level2) > 1) {
                Meth[indices, sex == level2] <- preprocessCore::normalize.quantiles(Meth[indices, 
                                  sex == level2, drop = FALSE])
                Unmeth[indices, sex == level2] <- preprocessCore::normalize.quantiles(Unmeth[indices, 
                                    sex == level2, drop = FALSE])
            }
        }
        else {
            Meth[indices, ] <- preprocessCore::normalize.quantiles(Meth[indices, 
                                                                        ])
            Unmeth[indices, ] <- preprocessCore::normalize.quantiles(Unmeth[indices, 
                                                                            ])
        }
    }
    assay(object, "Meth") <- Meth
    assay(object, "Unmeth") <- Unmeth
    return(object)
}
