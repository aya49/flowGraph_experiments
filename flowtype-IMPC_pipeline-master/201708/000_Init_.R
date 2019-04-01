#Last date modified 20151228 by Alice Yue
#run the RUNME.R script up to the part you are starting in this script
#Must run in chronological order; user can stop between steps but must ensure required variables exist

#3iTcells Project Assumptions: all Singlets isolated already, therefore removes columns FSC & SSC & Time, just leaves marker columns
#reading and writing files (uncomment code so that ONLY if the amount of FCS files you are handling is smaller than your memory, use fclean/ft-Frames to store all FCS into memory - a lot faster :)

#Install packages before loading libraries (if not yet installed):
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("foreach", "arrayMvout", "flowCore", "flowClean", "flowDensity", "flowType", "RchyOptimyx", "sva"))
# install.packages(c("extremevalues", "foreach", "igraph", "fitdistrplus", "corrplot", "dendextend", "stringr", "sva", "MDR", "arules", "arulesViz"))

options(stringsAsFactors=FALSE)
options(device="cairo")
root <- "~/projects/IMPC"
setwd(root)
centreL <- c("SangerTCellSPLEEN","SangerTCellMLN","CIPHE","TCP","H")
centre <- centreL[3]
source("Code/_3iTcellfunctions.R")
source("Code/_funcAlice.R")
# source("Code/load_manual_results.R") # for CIPHE

#General
# libr(extremevalues)
libr(Matrix)
libr(stringr)
# libr(foreach)

libr(flowCore)

#Gating
libr(flowClean)
libr(flowDensity)
libr(flowType)

#Analysis
libr(sva)
libr(MDR)
libr(arules)
libr(arulesViz)
libr(fitdistrplus)
libr(RchyOptimyx)
libr(igraph)



#004_USE.ME_LoadMatrixCellCount; for done centres; row(phen) x col(samples) ------------------------------------------

start <- Sys.time()

for (centre in centreL[c(4)]) {
  start1 <- Sys.time()
  
  setwd(root)
  ftPath <- paste("Justin/IMPC/", centre, sep="")
  matrixCellCount <- NULL
  sampleMeta <- NULL
  if (grepl("SangerTCell", centre)) {
    ftFilePath <- sort(dir(ftPath, pattern=".Rdata", all.files=TRUE, full.names=TRUE, recursive=TRUE) )
    ftFilePathcsv <- sort(dir(ftPath, pattern=".csv", all.files=TRUE, full.names=TRUE, recursive=TRUE) )
    ftGT <- folderNames(ftFilePath)
    ftFileNames <- fileNames(ftFilePath, "Rdata")
    ft <- get(load(ftFilePath[1]))
    markers <- c("CD44","CD62L","CD25","CD8","KLRG","CD5","CD45","CD161","CD4","GITR","TCRd")
  } else {
    ftFilePath <- sort(dir(ftPath, pattern="_All.Rdata", all.files=TRUE, full.names=TRUE, recursive=TRUE) )
    ftFilePathcsv <- sort(dir(ftPath, pattern="_prop.csv", all.files=TRUE, full.names=TRUE, recursive=TRUE) )
    ftGT <- fileNames(ftFilePath, "Rdata")
    ft <- get(load(ftFilePath[1]))[[1]]$ftype
    markers <- getMarkers(names(ft@CellFreqs)[length(ft@CellFreqs)])
  }
  save(markers, file=paste("Results", centre, "/markers.Rdata", sep=""))
  
  require(flowType)
  phenotype <- unlist(lapply(ft@PhenoCodes, function(x){return( decodePhenotype(x, markers, ft@PartitionsPerMarker) )}))
  phenoMeta <- data.frame(phenotype)
  phenoMeta$phenotype <- as.character(phenoMeta$phenotype); rm(phenotype)
  phenoMeta$phenoCode <- ft@PhenoCodes
  phenoMeta$phenoLevel <- unlist(lapply(phenoMeta$phenoCode, function(x){return(length(markers) - charOccurences("0", x)) } ))
  save(phenoMeta, file=paste("Results", centre, "/phenoMeta.Rdata", sep=""))
  
  if (grepl("SangerTCell",centre)) { # Sanger TCell panel MLN & SPLEEN organ
    
    for (i in 1:length(ftFilePath)) {
      ft <- get(load(ftFilePath[i]))
      npheno <- length(ft@CellFreqs)
      cat("\n", i, ". Loading file: ", ftGT[i], "/", ftFileNames[i], "; # of phenotypes x samples: ", npheno, " x ", length(ft), "", sep ="")
      matrixCellCount <- cbind(matrixCellCount, ft@CellFreqs)
      rm(ft)
    }
    colnames(matrixCellCount) <- ftFileNames
    
    sampleMetaTemp <- data.frame(lapply(read.csv(ftFilePathcsv[1])[,c(3,2,6,7)], as.character), stringsAsFactors=FALSE)
    for (i in 1:length(ftFileNames)) { sampleMeta <- rbind(sampleMeta, sampleMetaTemp[which(sampleMetaTemp[,1]==gsub("FT_","", ftFileNames[i]) & sampleMetaTemp[,2]==ftGT[i]),]) }
    sampleMeta <- as.data.frame(sampleMeta); rm(sampleMetaTemp)
    colnames(sampleMeta) <- c("fileName", "gene", "date", "gender")
    sampleMeta$barcode <- str_extract(sampleMeta$fileName, "L[0-9]+")
    sampleMetaTemp2 <- data.frame(lapply(read.csv(ftFilePathcsv[2])[,c(1:21,246:253)], as.character))
    ind <- match(sampleMeta$barcode,sampleMetaTemp2$Label.Barcode)
    sampleMetaTemp3 <- data.frame(lapply(read.csv(ftFilePathcsv[3]), as.character))
    
    # Sanger MLN TCell panel has different date formats
    require(lubridate)
    sampleMeta$date <- dmy(sampleMeta$date)
    if (grepl("MLN",centre)) {
      require(lubridate)
      sampleMeta$date <- as.character(ymd(sampleMeta$date))
      if (grepl("SangerTCellMLN", centre)) {
        d <- as.character(sampleMeta$date)
        for (i in 1:length(d)) {
          if (nchar(sampleMeta$date[i])==8) {
            d[i] <- paste0("20",d[i],collapse="")
          }
        }
        sampleMeta$date <- d
      }
    }
    
  } else { # Justin's files; all files of a gene merged together; CIPHE has different number of markers (64 genes); TCP (61 genes)
    # 08. 1-1-1-Baiap2l2_01_All; # of phenotypes x samples: 59049 x 6; different cell pop #; used markers (CD5, CD161, CD4, TCRd, GITR, CD44, CD62L, CD25, CD24, KLRG)
    # 16. 1-1-1-Cxcr7_01_All; # of phenotypes x samples: 59049 x 6	; different cell pop #; used markers (CD5, CD161, CD4, TCRd, GITR, CD44, CD62L, CD25, CD24, KLRG)
    # 21. 1-1-1-Dnalc4_01_All; # of phenotypes x samples: 6561 x 6	; different cell pop #; used markers (CD5, CD161, CD4, TCRd, CD44, CD62L, CD25, CD24)
    # 41. 1-1-1-Otub1_01_All; # of phenotypes x samples: 19683 x 6	; different cell pop #; used markers (CD5, CD161, CD4, CD8, TCRd, CD44, CD62L, CD25, CD24)
    # 43. 1-1-1-Pax4_01_All; # of phenotypes x samples: 6561 x 6	  ; different cell pop #; used markers (CD5, CD161, CD4, TCRd, CD44, CD62L, CD25, CD24)
    # 47. 1-1-1-Ptdss1_01_All; # of phenotypes x samples: 19683 x 6	; different cell pop #; used markers (CD5, CD161, CD4, CD8, TCRd, CD44, CD62L, CD25, CD24)
    # 48. 1-1-1-Rab19_01_All; # of phenotypes x samples: 19683 x 6	; different cell pop #; used markers (CD5, CD161, CD4, CD8, TCRd, CD44, CD62L, CD25, CD24)
    # 52. 1-1-1-Setbp1_01_All; # of phenotypes x samples: 19683 x 6	; different cell pop #; used markers (CD5, CD161, CD4, CD8, TCRd, CD44, CD62L, CD25, CD24)
    # 57. 1-1-1-Snx29_01_All; # of phenotypes x samples: 59049 x 6	; different cell pop #; used markers (CD5, CD161, CD4, TCRd, GITR, CD44, CD62L, CD25, CD24, KLRG)
    
    
    fileName <- NULL
    ftGenotype <- NULL
    marker <- NULL
    excludedPC <- NULL
    excludedGT <- NULL
    
    for (i in c(1:(length(ftFilePath)))) {
      ft <- get(load(ftFilePath[i]))
      if (i==1) {
        phenotype <- names(ft[[1]]$ftype@CellFreqs)
        phenoMeta <- data.frame(phenotype)
        phenoMeta$phenotype <- as.character(phenoMeta$phenotype); rm(phenotype)
        phenoMeta$phenoCode <- ft[[1]]$ftype@PhenoCodes
        phenoMeta$phenoLevel <- unlist(lapply(phenoMeta$phenoCode, function(x){return(length(markers) - charOccurences("0", x)) } ))
        save(phenoMeta, file=paste("Results", centre, "/phenoMeta.Rdata", sep=""))
      }
      if (i!=1) {
        if (npheno!=length(ft[[1]]$ftype@CellFreqs)) {
          cat("\t; different cell pop #, skipped; used markers ("); cat(getMarkers(names(ft[[1]]$ftype@CellFreqs)[length(ft[[1]]$ftype@CellFreqs)]), sep=", "); cat(")", sep="")
          excludedGT <- c(excludedGT, ftGT[i])
          next()
        }
      }
      npheno <- length(ft[[1]]$ftype@CellFreqs)
      cat("\n", i, ". Loading file: ", ftGT[i], "; # of phenotypes x samples: ", npheno, " x ", length(ft), "", sep ="")
      fileNameTemp <- as.character(names(read.csv(ftFilePathcsv[i], nrows=1)))
      fileNameTemp <- fileNameTemp[-which(is.na(fileNameTemp) | fileNameTemp=="X")]
      for (j in 1:length(ft)) {
        ftGenotype <- c(ftGenotype, ftGT[i])
        fileName <- c(fileName, fileNameTemp[j])
        #       if (all( match(markers,getMarkers(names(ft[[j]]$ftype@CellFreqs)[npheno]))==c(1:length(markers)) )) {
        matrixCellCount <- cbind(matrixCellCount, ft[[j]]$ftype@CellFreqs)
        #       } else {
        #         if (j==1) { markeri <- markers }
        #         if (!all( match(markeri,getMarkers(names(ft[[j]]$ftype@CellFreqs)[length(ft[[j]]$ftype@CellFreqs)]))==c(1:length(markeri)) )) { # if not the same as previous markers (or markers if j==1)
        #           marker[[length(fileName)]] <- markeri <- getMarkers(names(ft[[j]]$ftype@CellFreqs)[length(ft[[j]]$ftype@CellFreqs)])
        #           cat("\t; different cell pop #; used markers ("); cat(markeri, sep=", "); cat(")", sep="")
        #           m <- phenoCC(phenoMeta$phenoCode, markers, ft[[j]]$ftype@PhenoCodes, markeri)
        #         }
        #         matrixCellCount <- cbind(matrixCellCount, ft[[j]]$ftype@CellFreqs[m[[1]]])
        #         excludedi <- which(c(1:length(ft[[j]]$ftype@PhenoCodes))%in%m[[1]])
        #         if (length(excludedi)>0) { excludedPC[[length(fileName)+1]] <- ft[[j]]$ftype@CellFreqs[excludedi] }
        #       }
      }
      rm(ft)
    }
    colnames(matrixCellCount) <- fileName <- gsub("[.]", "-", gsub(".fcs", "", gsub("^X", "", fileName)))
    
    sampleMeta <- data.frame(fileName=fileName); rm(fileName)
    sampleMeta$gene <- ftGenotype; rm(ftGenotype)
    
    if (grepl("CIPHE", centre)) {
      sampleMeta$date <- as.Date(sapply(c(1:nrow(sampleMeta)), function(x) {return(strsplit(sampleMeta$fileName[x],"[_]")[[1]][1])} ), "%y-%B-%d")
      sampleMeta$sampleNum <- as.numeric(sapply(c(1:nrow(sampleMeta)), function(x) {return(strsplit(sampleMeta$fileName[x],"[_]")[[1]][3])} ))
    } else if (grepl("TCP", centre)) {
      matrixCSV <- read.csv(paste(ftPath, "/date.code.spreadsheet.csv", sep=""))
      matrixCSV[,4] <- gsub(".labelled","-labelled", gsub(".fcs","",matrixCSV[,4]))
      ind <- match(sampleMeta$fileName, matrixCSV[,4])
      sampleMeta$date <- matrixCSV[ind,1]
      sampleMeta$code <- matrixCSV[ind,2]
      sampleMeta$mouseEar <- matrixCSV[ind,3]
    } else if (grepl("H", centre)) {
      matrixCSV <- read.csv(dir(ftPath, pattern="HsummaryTable.csv", all.files=TRUE, full.names=TRUE, recursive=TRUE))
      matrixCSV[,1] <- gsub(".fcs","", gsub("[-%]",".", matrixCSV[,1]))
      ind <- match(sampleMeta$fileName,matrixCSV[,1])
      sampleMeta$date <- matrixCSV[ind,5]
      sampleMeta$sampleNum <- matrixCSV[ind,7]
      sampleMeta$gender <- matrixCSV[ind,8]; for (gi in 1:nrow(sampleMeta)) {if (sampleMeta$gender[gi]=="F") { sampleMeta$gender[gi] <- 0 } else { sampleMeta$gender[gi] <- 1 } }
      sampleMeta$baseline <- matrixCSV[ind,9]
    }
    
    # CIPHE Manula count import
    #   CSVfile <- paste("Results", centre, "/IMPC_20142015_merged_Panel1.csv", sep="")
    #   matrixCSV <- read.csv(CSVfile)
    #   matrixManualCounts <- matrix(0, nrow=55, ncol=ncol(matrixCellCount))
    #   mcIndex <- c(19:ncol(matrixCSV))[which(c(19:ncol(matrixCSV)) %% 2 != 0)]
    #   rownames(matrixManualCounts) <- sapply(mcIndex, function(x) {return( gsub("neg", "-", gsub("pos", "+", strsplit(colnames(matrixCSV)[x], "\\.\\.")[[1]][1]) ) ) })
    #   
    #   missingMetaIndex <- NULL
    #   for (i in 1:nrow(sampleMeta)) {
    #     row <- which(matrixCSV[,1]==sampleMeta$sampleNum[i] & as.Date(matrixCSV[,3], "%d-%m-%Y")==sampleMeta$date[i])
    #     if (length(row)==0) {missingMetaIndex <- c(missingMetaIndex, i); next()}
    #     sampleMeta[i,c(5:7)] <- matrixCSV[row,c(2,10,9)]
    #     matrixManualCounts[,i] <- as.vector( as.numeric(matrixCSV[row,mcIndex] ))
    #   }
    #   save(matrixManualCounts, file=paste("Results", centre, "/matrixManualCounts.Rdata", sep=""))
    
  }
  
  #phenotypes/rows with all 0
  countLim <- 5
  allLittleCountPhen <- which(apply(matrixCellCount[,-1], 1, function(x) all(x<=countLim))==T)
  all0Phen <- which(apply(matrixCellCount[,-1], 1, function(x) all(x==0))==T); if (length(all0Phen)==0) all0Phen <- c()
  save(all0Phen, file=paste("Results", centre, "/all0Phen.Rdata", sep=""))
  if (exists("excludedGT")) { save(excludedGT, file=paste("Results", centre, "/excludedGT.Rdata", sep="")); rm(excludedGT)}
  save(sampleMeta, file=paste("Results", centre, "/sampleMeta.Rdata", sep=""))
  save(matrixCellCount, file=paste("Results", centre, "/matrixCellCount.Rdata", sep=""))
  write.table(matrixCellCount, paste("Results", centre, "/matrixCellCount.csv", sep=""), sep=",", row.names=F)
  if (exists("marker")) { if (!is.null(marker)) { save(marker, file=paste("Results", centre, "/NotSameMarker.Rdata", sep="")) }}
  if (exists("excluededPC")) { if (!is.null(excludedPC)) { save(excludedPC, file=paste("Results", centre, "/excludedPC.Rdata", sep="")) }}
  
  matrixCellProp <- matrixCellCount%*%diag(1/matrixCellCount[1,])
  colnames(matrixCellProp) <- colnames(matrixCellCount)
  save(matrixCellProp, file=paste("Results", centre, "/matrixCellProp.Rdata", sep=""))
  
  cat("Total time used to load cell count: ", TimeOutput(start), "\n", sep="")
}

cat("Total time used to load cell count: ", TimeOutput(start), "\n", sep="")


#004.1_USE.ME_Normalize cell count -------------------------------------------
## lol I thought I found something, ended up that someone had already done an article on count normalization over at RNA seq differential expression analysis
## http://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25
## Trimmed Mean of M values - reference x sqrt(factor); col / sqrt(factor)
## taksource("https://bioconductor.org/biocLite.R")
# assume most cells aren't differentially expressed, use marjority ratio x count

# require(edgeR)
# f <- calcNormFactors(matrixCellCount[-1085,], lib.size=matrixCellCount[1085,], refColumn=NULL, p=.75, logratioTrim=.3, doWeighting=TRUE, Acutoff=-1e10)

start <- Sys.time()

for (centre in centreL[c(5)]) { 
  start1 <- Sys.time()
  
  sampleMeta <- get(load(paste("Results", centre, "/sampleMeta.Rdata",sep="")))
  cat("\nCalculating ", nrow(sampleMeta), " normalization factors for ",centre,": ",sep="")
  suppressWarnings(dir.create (paste("Results", centre, "/cellCountNormFactor",sep="")))
  matrixCellCount <- get(load(paste("Results", centre, "/matrixCellCount.Rdata",sep="")))
  if (nrow(matrixCellCount)<ncol(matrixCellCount)) { matrixCellCount <- t(matrixCellCount) }
  
  
  # Taken from TMM; done manually
  x <- as.matrix(matrixCellCount)[-unique(c(1,get(load(paste("Results", centre, "/all0Phen.Rdata",sep=""))))),] #take out total cell count
  
  p=0.75
  lib.size <- matrixCellCount[1,]
  ftWTGT <- c("+_+","+_Y"); if (!grepl("Sanger",centre)) ftWTGT <- "1-1-1-WildType_01_All"; if (centre=="H") ftWTGT <- "FR-FCM-ZYCB-WildType_01_All"
  refColumn <- which( matrixCellCount[1,]==median(matrixCellCount[1,which(sampleMeta$gene%in%ftWTGT)]) )
  #   f75 <- apply(t(t(x)/lib.size),2,function(data) quantile(data,p=p))
  #   refColumn <- which.min(abs(f75-mean(f75)))
  
  f <- rep(NA,ncol(x))
  fdiff <- rep(NA,ncol(x)) #diff between density peak and value (note: logged)
  ref <- x[,refColumn]
  nR <- lib.size[refColumn]
  doWeighting <- T
  Acutoff <- -1e10
  logratioTrim <- .3
  sumTrim <- 0.05
  for(i in ncol(x):1) { cat(i," ",sep="")
    obs <- x[,i]
    nO <- lib.size[i]
    logR <- log2((obs/nO)/(ref/nR))			# log ratio of expression, accounting for libr size
    absE <- (log2(obs/nO) + log2(ref/nR))/2	# absolute expression
    v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref	 # estimated asymptotic variance
    
    #	remove infinite values, cutoff based on A
    fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
    logR <- logR[fin]
    absE <- absE[fin]
    v <- v[fin]
    
    if(max(abs(logR)) < 1e-6) { f[i] <- 1; next }
    
    #	taken from the original mean() function
    n <- length(logR)
    loL <- floor(n * logratioTrim) + 1
    hiL <- n + 1 - loL
    loS <- floor(n * sumTrim) + 1
    hiS <- n + 1 - loS
    
    #	keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
    #	a fix from leonardo ivan almonacid cardenas, since rank() can return
    #	non-integer values when there are a lot of ties
    keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)
    
    if(doWeighting) {
      f[i] <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
    } else { f[i] <- mean(logR[keep], na.rm=TRUE) }
    
    #	Results will be missing if the two libraries share no features with positive counts
    #	In this case, return unity
    if(is.na(f[i])) f[i] <- 0
    
    # check if close to peak; if not, switch to peak
    require(pracma)
    d <- density(log2((matrixCellCount[,i])/matrixCellCount[,refColumn]), na.rm=T)
    p <- as.matrix(findpeaks(d$y)); if(ncol(p)==1) p <- t(p)
    p1 <- d$x[p[which.max(p[,1]),2]]
    fdiff[i] <- p1-f[i]
    
    pngname <- paste("Results", centre, "/cellCountNormFactor/", "cellCountNormFactor_",str_pad(i, 4, pad = "0"), "_",sampleMeta$gene[i], ".png", sep = "" )
    png (file=pngname , width=700, height=1800)
    par(mfrow=c(3,1), mar=(c(5, 5, 4, 2) + 0.1))
    
    plot(d); abline(v=f[i], col="red"); abline(v=p1, col="blue"); 
    
    plot((matrixCellCount[,i]+matrixCellCount[,refColumn])/2, log2(matrixCellCount[,i]/matrixCellCount[,refColumn]), cex=.5, main=paste("mean count vs. log2 fold change: ", sampleMeta$gene[i]," over refColumn ", sampleMeta$gene[refColumn],": f=",f[i], sep=""))
    abline(h=f[i], col="red")
    
    # if f[i] too far from peak
    if (abs(f[i]-p1)>.5) {
      abline(h=p1, col="blue")
      # f[i] <- p1
    }
    
    f[i] <- 1/2^f[i]
    
    plot((matrixCellCount[,i]+matrixCellCount[,refColumn])/2, log2((matrixCellCount[,i]*f[i])/matrixCellCount[,refColumn]), cex=.5, main=paste("AFTER CHANGE: mean count vs. log2 fold change: ", sampleMeta$gene[i]," over refColumn ", sampleMeta$gene[refColumn],": f=",f[i], sep=""))
    abline(h=0, col="red")
    dev.off()
    
  }
  # multiple of 1
  rm(x)
  save(f, file=paste("Results", centre, "/normFactor.Rdata", sep=""))
  save(fdiff, file=paste("Results", centre, "/normFactorDiffDensLogged.Rdata", sep=""))
  
  matrixCellCountAdj <- sapply(c(1:ncol(matrixCellCount)), function(x) {matrixCellCount[,x]*f[x]})
  colnames(matrixCellCountAdj) <- colnames(matrixCellCount)
  save(matrixCellCountAdj, file=paste("Results", centre, "/matrixCellCountAdj.Rdata",sep=""))
  write.csv(matrixCellCountAdj, file=paste("Results", centre, "/matrixCellCountAdj.csv",sep=""))
  
  cat("\nTime to calculate norm factors: ", TimeOutput(start1), "\n", sep="")
}

cat("\nTime to calculate norm factors: ", TimeOutput(start), "\n", sep="")



#004.2_Get means of samples for each KO gene --------------
matrixCellCountMean <- matrixM(matrixCellCountAdj, matrixType="mean") #gets mean for all GT
matrixCellCountMeanDiff <- matrixMDiff(matrixCellCountMean)
rownames(matrixCellCountMeanDiff) <- rownames(matrixCellCountAdj)
write.table(matrixCellCountMean, file=paste("Results", centre, "/matrixCellCountMean.csv", sep=""), col.names=TRUE, row.names=TRUE)
save(matrixCellCountMean, file=paste("Results", centre, "/matrixCellCountMean.Rdata", sep=""))
save(matrixCellCountMeanDiff, file=paste("Results", centre, "/matrixCellCountMeanDiff.Rdata", sep=""))

#004.2_Create Mean/Median matrixCellCount differences between mean WT and every other sample; CIPHE --------------
start <- Sys.time()
# matrixCPMed <- matrixM(matrixCP, matrixType="med") #gets median for all GT
# matrixCPMedDiff <- matrixMDiff(matrixCPMed)
# write.table(matrixCPMed, file=paste("Results", centre, "/matrixCPMed.csv", col.names=TRUE, row.names=TRUE)
# save(matrixCPMed, file=paste("Results", centre, "/matrixCPMed.Rdata")
# save(matrixCPMedDiff, file=paste("Results", centre, "/matrixCPMedDiff.Rdata")

matrixCellCountDiff <- matrixMDiff(cbind( matrixM(matrixCellCount[,ftWTIndex], matrixType="mean"), matrixCellCount[,unlist(ftKOIndex)] ))
save(matrixCellCountDiff, file=paste("Results", centre, "/matrixCellCountDiff.Rdata", sep=""))
TimeOutput(start)

#004.2_Barcode: Sample means, WT KO difference, adjust difference, pvalue -----------------------------------------------------------

if(!exists("matrixCellCount")) { matrixCellCount <- get(load(paste("Results", centre, "/matrixCellCount.Rdata", sep=""))) }

#Use MatrixCellPop to get pVal Matrix
# matrixPValFULL <- matrixPValBC(matrixCellCount, append(list(ftWTIndex),ftKOIndex), cbind(rep(1,length(ftKOGT)), c(2:(length(ftKOGT)+1))), ftKOGT) #2+hrs
# save(matrixPValFULL, file="Results/matrixPValFULL.Rdata")

# #make total cell count same for each each WT, KO; same as comparing proportions, but retain total count.
# sameTotalCount <- TRUE #only do once; TRUE=done
# if (!sameTotalCount) {
#   for (j in 2:length(ftWTIndex)) {
#     m <- matrixCellCount[1,ftWTIndex[1]]/matrixCellCount[1,ftWTIndex[j]]
#     matrixCellCount[,ftWTIndex[j]] <- m*matrixCellCount[,ftWTIndex[j]]
#   }
#   for (i in 1:length(ftKOIndex)) {
#     for (j in 2:length(ftKOIndex[[i]])) {
#       m <- matrixCellCount[1,ftKOIndex[[i]][1]]/matrixCellCount[1,ftWTIndex[j]]
#       matrixCellCount[,ftKOIndex[[i]][j]] <- m*matrixCellCount[,ftKOIndex[[i]][j]]
#     }
#   }
# }
# 
# png(paste("Results", centre, "/cellCount.png", sep=""), width=6000, height=6000)
# matplot(matrixCellCount[,ftWTIndex], type="l")
# dev.off()

#005.0_Confounding factors - create residual expression matrix

#extract meta data for each sample


#004.3_USE.ME_Get pVal Matrix --------------
start <- Sys.time()
# matrixPValFULL <- matrixPValBC(matrixCellCount, append(list(ftWTIndex),ftKOIndex), cbind(rep(1,length(ftKOGT)), c(2:(length(ftKOGT)+1))), ftKOGT) #2+hrs
# save(matrixPValFULL, file=paste("Results", centre, "/matrixPValFULL.Rdata", sep=""))

adjust="BH"
for (centre in centreL[5]) {
  for (cp in 1:2) {
    start1 <- Sys.time()
    if (cp==1) {
      matrixCellCountAdj <- get(load(paste("Results", centre, "/matrixCellCountAdj.Rdata",sep="")))
    } else {
      matrixCellCountAdj <- get(load(paste("Results", centre, "/matrixCellProp.Rdata",sep="")))
    }
    if (nrow(matrixCellCountAdj)<ncol(matrixCellCountAdj)) { matrixCellCountAdj <- t(matrixCellCountAdj) }
    sampleMeta <- get(load(paste("Results", centre, "/sampleMeta.Rdata",sep="")))
    
    ftWTGT <- c("+_+","+_Y"); if (!grepl("Sanger",centre)) ftWTGT <- "1-1-1-WildType_01_All"; if (centre=="H") ftWTGT <- "FR-FCM-ZYCB-WildType_01_All"
    g <- getGTindex(sampleMeta$gene, ftWTGT)
    ftGT <- g$ftGT; ftWTIndex <- g$ftWTIndex; ftKOGT <- g$ftKOGT; ftKOIndex <- g$ftKOIndex; ftKOgoodGTIndex <- g$ftKOgoodGTIndex; rm(g)  #Index of unique KO genotypes that has 3+ samples available ftKOIndex(ftKOGTIndexIndex)
    
    colCombos <- NULL
    for (i in length(ftKOIndex):1) {
      colCombos[[i]][[1]] <- ftWTIndex
      colCombos[[i]][[2]] <- ftKOIndex[[i]]
    }
    
    matrixPValFULL <- matrixPValBC(matrixCellCountAdj, colCombos=colCombos, colLabel=ftKOGT, adjust=adjust) #2+hrs
    mpvf <- matrixPValFULL[[1]]
    mpvfa <- matrixPValFULL[[2]]
    rm(matrixPValFULL)
    
    pValThres <- .05 #delete phenotypes/rows without any significant changes from the pVal matrix
    insigPhenIndex <- which(apply(mpvf[,-1], 1, function(x) any(x<pValThres))==FALSE)
    insigPhenIndexa <- which(apply(mpvfa[,-1], 1, function(x) any(x<pValThres))==FALSE)
    
    setwd(root)
    if (cp==1) {
      save(mpvf, file=paste("Results", centre, "/matrixPValFULL.Rdata", sep=""))
      save(mpvfa, file=paste("Results", centre, "/matrixPValFULL", adjust,".Rdata", sep=""))
      save(insigPhenIndex, file=paste("Results", centre, "/insigPhenIndex.Rdata", sep=""))
      save(insigPhenIndexa, file=paste("Results", centre, "/insigPhenIndex", adjust,".Rdata", sep=""))
    } else {
      save(mpvf, file=paste("Results", centre, "/matrixPValFULLProp.Rdata", sep=""))
      save(mpvfa, file=paste("Results", centre, "/matrixPValFULL", adjust,"Prop.Rdata", sep=""))
      save(insigPhenIndex, file=paste("Results", centre, "/insigPhenIndexProp.Rdata", sep=""))
      save(insigPhenIndexa, file=paste("Results", centre, "/insigPhenIndex", adjust,"Prop.Rdata", sep=""))
    }
    
    cat("\nTime taken to calculate p values & barcode matrices is: ",TimeOutput(start1), "\n", sep="") #3iTcell ~40min
  }
}

cat("\nTime taken to calculate p values & barcode matrices is: ",TimeOutput(start), "\n", sep="") #3iTcell ~40min

#004.3B_Barcode: Trim p value and barcode matrix-----------------------------------------------------------

if(!exists("matrixCellCountAdj")) { matrixCellCountAdj <- as.matrix(get(load(paste("Results", centre, "/matrixCellCountAdj.Rdata", sep="")))) }
if(!exists("matrixPValFULL")) { matrixPValFULL <- get(load(paste("Results", centre, "/matrixPValFULL.Rdata", sep=""))) }

#TRIM: Delete phenotypes & KOgenotypes with insignificant pVal for all KO genes; Can embed into function...
pValThres <- .05 #delete phenotypes/rows without any significant changes from the pVal matrix
insigPhenIndex <- which(apply(matrixPValFULL[,-1], 1, function(x) any(x<pValThres))==FALSE)
if (length(which(insigPhenIndex==1))>0) { insigPhenIndex <- insigPhenIndex[-which(insigPhenIndex==1)] } #keep first row, all cells
# insigPhenNotInGSIndex <- insigPhenIndex[-intersect(insigPhenIndex,GSPhenIndex)]
insigKOGTIndex <- which(apply(matrixPValFULL[-1,], 2, function(x) any(x<pValThres))==FALSE)
# matrixPVal <- matrixTrim(matrixPValFULL, pValThres=pValThres) #Does same as below, returns trimmed matrix
if (length(insigKOGTIndex)>0) {
  matrixPVal <- matrixPValFULL[,-insigKOGTIndex]
  insigKOIndex <- which(colnames(matrixCellProp)%in%names(insigKOGTIndex)==TRUE)
  matrixCC <- as.matrix(matrixCellCount[,-insigKOIndex])
}
if (length(insigPhenIndex)>0) {
  matrixPVal <- matrixPValFULL[-insigPhenIndex,]
  matrixCC <- matrixCellCount[-insigPhenIndex,]
  phenoMeta <- phenoMeta[-insigPhenIndex,]
}

save(matrixPVal, file=paste("Results", centre, "/matrixPVal.Rdata",sep=""))
write.csv(matrixPVal, paste("Results", centre, "/matrixPVal.csv",sep=""), row.names=T)
save(matrixCC, file=paste("Results", centre, "/matrixCC.Rdata",sep=""))
write.csv(matrixCC, paste("Results", centre, "/matrixCC.csv",sep=""), row.names=T)
save(phenoMeta, file=paste("Results", centre, "/phenoMeta.Rdata",sep=""))

# plotRchyPDF(path=paste(pathResults[8], "test_", iGT, ".pdf", sep=""), popMetaB$WCpValBH, popMetaB$PhenoCodes, startPhenotype=popMetaB$PhenoCodes[sigPopIndex15])

#Which phenotypes are < CellCountThreshold for both KO & WT that are still in the matrix
# uniqueSigKOGT <- as.character(colnames(matrixPVal))
# cellCountThres <- 100
#list phenotypes (rows) where all WT Genotype are < cellCountThres
# lowCountWTPhenotypesIndex <- which(apply(matrixCC[,ftWTIndex], 1, function(x) !all(x<cellCountThres))==FALSE)
# lowCountKOofWTPhenotypesIndex <- NULL
# for (i in 1:length(uniqueSigKOGT)) {
#   lowCountSigKOofWTPhenotypesIndex[[i]] <- which(apply(matrixCC[lowCountWTPhenotypesIndex, which(colnames(matrixCC)==uniqueSigKOGT[i])], 1, function(x) !all(x<cellCountThres))==FALSE)
# }
# lowCountPhenotypesIndex <- lowCountWTPhenotypesIndex[] #!!!!!!!!!!!!!!!!!

#004.0_USE.ME_Get Genotype indeces -----------------------------------------------------

#3iTcell Just change the file paths and WT phenotypes
# matrixCellCount <- get(load("Results/ftALLCount.Rdata"))

# GSPhenotypes <- GatingStrategyPhenotypes() #Gating Strategy Phenotypes; 3iTCells only
# GSPhenIndex <- which(rownames(matrixCellCount)%in%GSPhenotypes ==TRUE)
# ftGT <- colnames(matrixCellCount)
# ftWTGT <- c("+_+", "+_Y") #Unique WT Genotypes
# ftWTIndex <- which(ftGT%in%ftWTGT==TRUE)
# ftKOGT <- unique(ftGT[-ftWTIndex]) #Unique KO Genotypes
# ftKOIndex <- NULL #Index of unique KO Genotypes (one genotype per vector in the list)
# for (i in 1:length(ftKOGT)) { ftKOIndex[[i]] <- which(ftGT==ftKOGT[i]) }
# ftKOgoodGTIndex <- which( unlist(lapply(ftKOIndex, length)) > 2 ) #Index of unique KO genotypes that has 3+ samples available ftKOIndex(ftKOGTIndexIndex)


#CIPHE & Sanger & TCP Meta data
# GSPhenotypes <- GatingStrategyPhenotypes() #Gating Strategy Phenotypes; 3iTCells only
# GSPhenIndex <- which(rownames(matrixCellCount)%in%GSPhenotypes ==TRUE)
ftWTGT <- c("+_+","+_Y"); if (!grepl("Sanger",centre)) ftWTGT <- "1-1-1-WildType_01_All"; if (centre=="H") ftWTGT <- "FR-FCM-ZYCB-WildType_01_All"
g <- getGTindex(sampleMeta$gene, ftWTGT)
ftGT <- g$ftGT; ftWTIndex <- g$ftWTIndex; ftKOGT <- g$ftKOGT; ftKOIndex <- g$ftKOIndex; ftKOgoodGTIndex <- g$ftKOgoodGTIndex #Index of unique KO genotypes that has 3+ samples available ftKOIndex(ftKOGTIndexIndex)
rm(g)






#006.0_GetSinglePhen ---------------------------------------------------------

#Load prop & count only for phenName

# phenName <- "KLRG+" # big wave
# phenName <- "CD25-CD8-CD5-TCRd+" # split in the middle
phenName <- "TCRd+" # split in the middle
# phenName <- "CD8+"

for (centre in centreL[c(1:5)]) { #CIPHE and TCP doesn't have gender
  setwd(root); setwd(paste("Results",centre,sep=""))
  b <- getSinglePhen(centre, phenName, get(load("sampleMeta.Rdata")), get(load("phenoMeta.Rdata")), get(load("markers.Rdata")), get(load("matrixCellCountAdj.Rdata")))
  b <- sortByDate(b)
  require(lubridate)
  # if (!exists("phenData")) {
    phenData <- b
  # } else { phenData <- rbind(phenData,b) }
    ftWTGT <- c("+_+","+_Y"); if (!grepl("Sanger",centre)) ftWTGT <- "1-1-1-WildType_01_All"; if (centre=="H") ftWTGT <- "FR-FCM-ZYCB-WildType_01_All"
    wt <- which(phenData$gene%in%ftWTGT)
    fm <- which(phenData$gender==1)
    
      require(changepoint)
      y <- phenData$count[wt]
      plot(y, main=phenName)
      
      methods <- c("BinSeg","PELT") #AMOC at most one change, PELT is poisson (quick exact), BinSeg (quick, approx), SegNeigh (slow, exact)
      
      png(filename=paste("ChangePoint_",phenName,sep=""), width=length(methods)*800, height=8*400)
      layout(matrix(c(1,2,4,6,8,9,11,13,1,3,5,7,8,10,12,14),ncol=2))
      for (i in 1:2) {
        if (i==1) {
          y <- phenData$count[wt]
          plot(y, main=phenName, pch=19,cex=.4)
        } else {
          require(FKF)
          ## Set constant parameters:
          dt <- ct <- matrix(0) 
          Zt <- Tt <- matrix(1)
          a0 <- y[1]            # Estimation of the first year flow 
          P0 <- matrix(100)     # Variance of 'a0'
          
          ## Estimate parameters 23min if TS
          fit.fkf <- optim(c(HHt = var(y, na.rm = TRUE) * .5,
                             GGt = var(y, na.rm = TRUE) * .5),
                           fn = function(par, ...) -fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik,
                           yt = rbind(y), a0 = a0, P0 = P0, dt = dt, ct = ct,
                           Zt = Zt, Tt = Tt, check.input = FALSE)
          
          ## Filter Nile data with estimated parameters:
          fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, HHt = matrix(fit.fkf$par[1]), GGt = matrix(fit.fkf$par[2]), yt = rbind(y))
          y0 <- y
          y <- fkf.obj$att[1,]
          
          ## Compare with the stats' structural time series implementation: 5min if TS
          fit.stats <- StructTS(y0, type = "level")
#           fit.fkf$par
#           fit.stats$coef
#           
          ## Plot the flow data together with fitted local levels:
          plot(y0, main = "kalman filter", pch=19,cex=.4)
          lines(fitted(fit.stats), col = "green")
          lines(y, col = "blue")
          legend("top", c("Actual datapoints", "Local level (StructTS)", "Local level (fkf)"), col = c("black", "green", "blue"), lty = 1)
        }
        for (i in 1:length(methods)) {
          mvalue <- cpt.mean(y,method=methods[i], Q=20, penalty="MBIC")
          plot(mvalue, main=paste("WT only; mean change: ",methods[i], "; Penalty MBIC; ",phenName, sep=""))
        }
        for (i in 1:length(methods)) {
          vvalue = cpt.var(diff(y), method=methods[i], penalty="MBIC")
          plot(vvalue, main=paste("variance change: ",methods[i],sep=""))
        }
        for (i in 1:length(methods)) {
          vvalue = cpt.meanvar(diff(y), method=methods[i], penalty="MBIC")
          plot(vvalue, main=paste("variance/mean change: ",methods[i],sep=""))
        }
      }
      dev.off()
      
  
    
}
# sort by time
phenData <- orderByDate(phenData)

# perturb dates so its all unique
dup <- duplicated(phenData$date)
phenData2 <- phenData
phenData2$date <- as.character(phenData$date)
phenData2$date[which(!dup)] <- paste(phenData2$date[which(!dup)], " 00:00:00", sep="")
startingDate <- which(!dup)+1; if (startingDate[length(startingDate)]>nrow(phenData)) startingDate[-length(startingDate)]
while (length(startingDate)>0) {
  phenData2$date[startingDate] <- as.character(ymd_hms(phenData2$date[startingDate-1])+seconds(1))
  startingDate <- startingDate+1; if (startingDate[length(startingDate)]>nrow(phenData2)) startingDate[-length(startingDate)]
  startingDate <- startingDate[which(dup[startingDate])]
}

# Plot
spln <- which(phenData$centre=="SangerTCellSPLEEN")
ftWTGT <- c("+_+","+_Y"); if (!grepl("Sanger",centre)) ftWTGT <- "1-1-1-WildType_01_All"; if (centre=="H") ftWTGT <- "FR-FCM-ZYCB-WildType_01_All"
wt <- which(phenData$gene%in%ftWTGT)
fm <- which(phenData$gender==1)
plot(phenData[,3], as.numeric(phenData[,6]), cex=.5, col="blue", pch=16); points(phenData[spln,3], phenData[spln,6], col="red", cex=.5, pch=16); points(phenData[wt,3], phenData[wt,6]); points(phenData[fm,3], phenData[fm,6], col="green")
x11(); plot(phenData[,3], as.numeric(phenData[,7]), cex=.5, col="blue", pch=16); points(phenData[spln,3], phenData[spln,7], col="red", cex=.5, pch=16); points(phenData[wt,3], phenData[wt,7]); points(phenData[fm,3], phenData[fm,6], col="green")

# line plot
require(plotly)
p <- plot_ly(economics, x = date, y = uempmed, name = "unemployment")
p %>% add_trace(y = fitted(loess(uempmed ~ as.numeric(date), economics)), x = date)

p <- plot_ly(phenData[intersect(wt,spln),], x = date, y = count, name = "count")
p %>% add_trace(y = fitted(loess(count ~ as.numeric(date), phenData[intersect(wt,spln),])), x = date)
htmlwidgets::saveWidget(as.widget(p), "index.html")



#006.0_GetSinglePhen20160901 ---------------------------------------------------

# phenName <- "KLRG+" # big wave
# phenName <- "CD25-CD8-CD5-TCRd+" # split in the middle
phenName <- "TCRd+" # split in the middle

for (centre in centreL[c(1:5)]) { #CIPHE and TCP doesn't have gender
  
  setwd(root); setwd(paste("Results",centre,sep=""))
  load("markers.Rdata")
  load("sampleMeta.Rdata")
  load("phenoMeta.Rdata")
  load("matrixCellCountAdj.Rdata")
  load("matrixCellProp.Rdata")
  
  png(filename="samplePhenWT_singleMarker.png", width=1000, height=length(markers)*500*2)
  par(mfrow=c(length(markers)*2,1))
  
  for (j in 1:length(markers)) {
    phenName <- paste(markers[j],"+",sep="")
    phenIndex <- which(phenoMeta$phenotype==phenName)
    if (length(phenIndex)==0) {
      phenNameAlph <- strsplit(phenName,"[-+]")[[1]]
      phenNameEtc <- strsplit(phenName,"[0-9A-z]")[[1]]
      phenNameEtc <- phenNameEtc[which(phenNameEtc!="")]
      orders <- NULL
      for (i in 1:length(phenNameAlph)) {
        orders[i] <- which(markers==phenNameAlph[i])
      }
      orders <- order(orders)
      phenName2 <- NULL
      for (i in 1:length(orders)) {
        phenName2 <- paste(phenName2, phenNameAlph[orders[i]], phenNameEtc[orders[i]], sep="")
      }
      phenIndex <- which(phenoMeta$phenotype==phenName2)
    }
    b <- sampleMeta
    b$centre <- rep(centre,nrow(sampleMeta))
    b$count <- as.numeric(matrixCellCountAdj[phenIndex,])
    b$prop <- as.numeric(matrixCellProp[phenIndex,])
    
    require(lubridate)
    b$date <- as.character(ymd(b$date))
    if (grepl("SangerTCellMLN", centre)) {
      d <- as.character(b$date)
      for (i in 1:length(d)) {
        if (nchar(b$date[i])==8) {
          d[i] <- paste0("20",d[i],collapse="")
        }
      }
      b$date <- d
    }
    g <- rep(0,nrow(b))
    g[which(b$gender=="Male")] <- 1
    b$gender <- g # Female=0; Male=1
    b$date <- ymd(b$date)
    
    ftWTGT <- c("+_+","+_Y"); if (!grepl("Sanger",centre)) ftWTGT <- "1-1-1-WildType_01_All"; if (centre=="H") ftWTGT <- "FR-FCM-ZYCB-WildType_01_All"
    wt <- which(b$gene%in%ftWTGT)
    fm <- which(b$gender==0)
    plot(b$date[wt],b$count[wt], cex=.5, pch=16, main=paste(phenName," normalized count, WT only (green=female)",sep="")); points(b$date[intersect(fm,wt)],b$count[intersect(fm,wt)], col="green")
    plot(b$date[wt],b$prop[wt], cex=.5, pch=16, main=paste(phenName," proportion, WT only (green=female)",sep="")); points(b$date[intersect(fm,wt)],b$prop[intersect(fm,wt)], col="green")
  }
  graphics.off()
  
  if (grepl("Sanger",centre)) {
    gs <- getGSsanger()
    for (i in 1:2) {
      png(filename=paste("samplePhenWT_GSpops_",i,".png",sep=""), width=1000, height=length(gs)*500)
      par(mfrow=c(length(gs)*2,1))
      for (j in 1:length(gs)) {
        phenName <- gs[j]
        phenIndex <- which(phenoMeta$phenotype==phenName)
        if (length(phenIndex)==0) {
          phenNameAlph <- strsplit(phenName,"[-+]")[[1]]
          phenNameEtc <- strsplit(phenName,"[0-9A-z]")[[1]]
          phenNameEtc <- phenNameEtc[which(phenNameEtc!="")]
          orders <- NULL
          for (i in 1:length(phenNameAlph)) {
            orders[i] <- which(markers==phenNameAlph[i])
          }
          orders <- order(orders)
          phenName2 <- NULL
          for (i in 1:length(orders)) {
            phenName2 <- paste(phenName2, phenNameAlph[orders[i]], phenNameEtc[orders[i]], sep="")
          }
          phenIndex <- which(phenoMeta$phenotype==phenName2)
        }
        if (length(phenIndex)==0) next()
        b <- sampleMeta
        b$centre <- rep(centre,nrow(sampleMeta))
        b$count <- as.numeric(matrixCellCountAdj[phenIndex,])
        b$prop <- as.numeric(matrixCellProp[phenIndex,])
        
        require(lubridate)
        b$date <- as.character(ymd(b$date))
        if (grepl("SangerTCellMLN", centre)) {
          d <- as.character(b$date)
          for (i in 1:length(d)) {
            if (nchar(b$date[i])==8) {
              d[i] <- paste0("20",d[i],collapse="")
            }
          }
          b$date <- d
        }
        g <- rep(0,nrow(b))
        g[which(b$gender=="Male")] <- 1
        b$gender <- g # Female=0; Male=1
        b$date <- ymd(b$date)
        
        ftWTGT <- c("+_+","+_Y"); if (!grepl("Sanger",centre)) ftWTGT <- "1-1-1-WildType_01_All"; if (centre=="H") ftWTGT <- "FR-FCM-ZYCB-WildType_01_All"
        wt <- which(b$gene%in%ftWTGT)
        fm <- which(b$gender==0)
        plot(b$date[wt],b$count[wt], cex=.5, pch=16, main=paste(phenName," normalized count, WT only (green=female)",sep="")); points(b$date[intersect(fm,wt)],b$count[intersect(fm,wt)], col="green")
        plot(b$date[wt],b$prop[wt], cex=.5, pch=16, main=paste(phenName," proportion, WT only (green=female)",sep="")); points(b$date[intersect(fm,wt)],b$prop[intersect(fm,wt)], col="green")
      }
      graphics.off()
    }
  }
}



#006.05_fitDistribution-------------------------------------------
libr(fitdistrplus)
libr(logspline)
fit.weibull <- fitdist(x, "weibull")
plot(fit.weibull)
fit.norm <- fitdist(y, "norm")
plot(fit.norm)
fit.poisson <- fitdist(y, "poisson")


#006.1_change point detection; dlm/filter---------------------------------------------

# create time series object
libr(its) #irregular time seires
require(zoo) ## convert to a zoo object, with order given by the `datefield`
phenDataTS <- with(data.frame(count=phenData2$count,date=ymd_hms(phenData2$date)), zoo(count, order.by=date))
start <- Sys.time(); phenDataTS <- as.ts(phenDataTS); TimeOutput(start) #SangerTCellSPLEEN; one Phen = 4min; 2000 to 51580803 data points





require(KFK)
y <- phenData$count

## Set constant parameters:
dt <- ct <- matrix(0) 
Zt <- Tt <- matrix(1)
a0 <- y[1]            # Estimation of the first year flow 
P0 <- matrix(100)     # Variance of 'a0'

## Estimate parameters 23min if TS
fit.fkf <- optim(c(HHt = var(y, na.rm = TRUE) * .5,
                   GGt = var(y, na.rm = TRUE) * .5),
                 fn = function(par, ...) -fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik,
                 yt = rbind(y), a0 = a0, P0 = P0, dt = dt, ct = ct,
                 Zt = Zt, Tt = Tt, check.input = FALSE)

## Filter Nile data with estimated parameters:
fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, HHt = matrix(fit.fkf$par[1]), GGt = matrix(fit.fkf$par[2]), yt = rbind(y))
y0 <- fkf.obj$att[1,]

## Compare with the stats' structural time series implementation: 5min if TS
fit.stats <- StructTS(y, type = "level")
fit.fkf$par
fit.stats$coef

## Plot the flow data together with fitted local levels:
plot(y, main = "Nile flow", pch=19,cex=.4)
lines(fitted(fit.stats), col = "green")
lines(y0, col = "blue")
legend("top", c("Nile flow data", "Local level (StructTS)", "Local level (fkf)"), col = c("black", "green", "blue"), lty = 1)









libr(dlm) # Bayes
#constant model: time invariant (doesn't depend on time)
buildFun <- function(x) {
  m <- dlmModPoly(1, dV = exp(x[1])) #nth order poly
#   m <- dlmModARMA() #potentially multivar ARMA process (ar autoregressive coeff list of matrices; ma moving average coeff list of matrices)
#   m <- dlmModReg() #linear regression
#   m <- dlmModSeas(frequency=NumOfSeasons, dV=exp(x[1])) #periodic/seasonal factors
#   m <- dlmModTrig(period, dV = exp(x[1])) #periodic/trigonometric form
  
#   m$JW <- matrix(1)
#   m$X <- matrix(exp(x[2]), nc = 1, nr = length(Nile))
#   j <- which(time(Nile) == 1899) #must specify changepoint
#   m$X[j,1] <- m$X[j,1] * (1 + exp(x[3]))
#   return(m)
}

start <- Sys.time(); fit <- dlmMLE(y, parm = c(10,0,0), build=buildFun); TimeOutput(start)
dlmphenDataTS <- buildFun(fit$par)
V(dlmphenDataTS)
phenDataTSfilt <- dlmFilter(y, dlmphenDataTS)
plot(y, type ='o', col = "seagreen")
lines(dropFirst(phenDataTSfilt$m), type ='o',pch = 20, col = "brown")

# fit <- dlmMLE(Nile, parm = c(10,0,0), build=buildFun)
# dlmNileJump <- buildFun(fit$par)
# V(dlmNileJump)
# dlmNileJump$X[c(1, which(time(Nile) == 1899)), 1] # take out 1899, pretty much linear
# nileJumpFilt <- dlmFilter(Nile, dlmNileJump)
# plot(Nile, type ='o', col = "seagreen")
# lines(dropFirst(nileJumpFilt$m), type ='o',pch = 20, col = "brown")

attach(nileJumpFilt)
v <- unlist(dlmSvd2var(U.C, D.C))

# upper and lower lines
pl <- dropFirst(m) + qnorm(0.05, sd = sqrt(v[-1]))
pu <- dropFirst(m) + qnorm(0.95, sd = sqrt(v[-1]))
detach()
lines(pl, lty = 2, col = "brown")
lines(pu, lty = 2, col = "brown")

#006.2_ChangePoint analysis & Align means----------------------

require(changepoint)
y <- phenData$count[wt]
plot(y, main=phenName)

methods <- c("BinSeg","PELT") #AMOC at most one change, PELT is poisson (quick exact), BinSeg (quick, approx), SegNeigh (slow, exact)

png(filename=paste("ChangePoint_",phenName,sep=""), width=length(methods)*800, height=3*400)
par(mfrow=c(3,length(methods)))
for (i in 1:length(methods)) {
  mvalue <- cpt.mean(y,method=methods[i], Q=20, penalty="MBIC")
  plot(mvalue, main=paste("mean change: ",methods[i], "; Penalty MBIC; ",phenName, sep=""))
}
for (i in 1:length(methods)) {
  vvalue = cpt.var(diff(y), method=methods[i], penalty="MBIC")
  plot(vvalue, main=paste("variance change: ",methods[i],sep=""))
}
for (i in 1:length(methods)) {
  vvalue = cpt.meanvar(diff(y), method=methods[i], penalty="MBIC")
  plot(vvalue, main=paste("variance/mean change: ",methods[i],sep=""))
}
dev.off()

x11()
d = decompose(y)
plot(d)


#006.4_MultiLevelModel ---------------------------------------------------------
require(brms)
require(multilevel) #MASS, nlme
require(brms)
require(lme4)

require(comato)
#hopkin's index; nonrandomness


# MLM
model0 <- lmer(count ~ (1 | barcode), data=phenData, REML=F)
model1 <- lmer(count ~ (1 | centre), data=phenData, REML=F)
anova(model0, model1) # centre matter?

model2 <- lmer(count ~ gender + (1 | centre), data=phenData, REML=F)
anova(model1, model2) # gender matter?

model2.5 <- lmer(count ~ gender + (1 + gender | centre), data=phenData, REML=F)
anova(model2, model2.5) # gender matters differently in different centres?
mydata <- data.frame(cohort90=phenData$gender, schoolid=phenData$centre)
predCount <- fitted(model2.5) #avg fitted regression line + uj
datapred <- unique(data.frame(cbind(predCount = predCount, cohort90 = mydata$cohort90, schoolid = mydata$schoolid)))
xyplot(predCount ~ cohort90, data = datapred, groups = schoolid, type = c("p", "l"), col = "blue")
# assume effect same for all schools therefore same slopes

model3 <- lmer(count ~ gender + (gender | gene/centre), data=phenData, REML=F)
anova(model1, model3) # gene nested in centre matter?

model3 <- lmer(count ~ gender + (1 | gene/centre), data=phenData, REML=F)
model3 <- lmer(count ~ date + (date | gender/gene/centre), data=phenData, REML=F)
model5 <- brm(formula = count ~ date * (date | gender/gene/centre), data = phenData, family = lognormal(), #family = disribution for error
              prior = c(set_prior("normal(0,5)", class = "b"), set_prior("cauchy(0,2)", class = "sd"), set_prior("lkj(2)", class = "cor")), warmup = 1000, iter = 2000, chains = 4, control = list(adapt_delta = 0.95), 
              nonlinear=)

require(forecast)
require(FKF)












# plotting

libr(ggplot2)
libr(lme4)
libr(multcomp)
dataset <- expand.grid(experiment = factor(seq_len(10)), status = factor(c("N", "D", "R"), levels = c("N", "D", "R")), reps = seq_len(10))
dataset$value <- rnorm(nrow(dataset), sd = 0.23) + with(dataset, rnorm(length(levels(experiment)), sd = 0.256)[experiment] + ifelse(status == "D", 0.205, ifelse(status == "R", 0.887, 0))) + 2.78
model <- lmer(value~status+(1|experiment), data = dataset)
tmp <- as.data.frame(confint(glht(model, mcp(status = "Tukey")))$confint)
tmp$Comparison <- rownames(tmp)
ggplot(tmp, aes(x = Comparison, y = Estimate, ymin = lwr, ymax = upr)) + geom_errorbar() + geom_point()

tmp <- as.data.frame(confint(glht(model))$confint)
tmp$Comparison <- rownames(tmp)
ggplot(tmp, aes(x = Comparison, y = Estimate, ymin = lwr, ymax = upr)) + geom_errorbar() + geom_point()

model <- lmer(value ~ 0 + status + (1|experiment), data = dataset)
tmp <- as.data.frame(confint(glht(model))$confint)
tmp$Comparison <- rownames(tmp)
ggplot(tmp, aes(x = Comparison, y = Estimate, ymin = lwr, ymax = upr)) + geom_errorbar() + geom_point()

#006.4_plotMixedModel--------------------------

set.seed(101)
dataset <- expand.grid(experiment = factor(seq_len(10)), 
                       status = factor(c("N", "D", "R"), levels = c("N", "D", "R")), 
                       reps = seq_len(10))
X <- model.matrix(~status,dataset)
dataset <- transform(dataset, 
                     value=rnorm(nrow(dataset), sd = 0.23) +   ## residual
                       rnorm(length(levels(experiment)), sd = 0.256)[experiment] +  ## block effects
                       X %*% c(2.78,0.205,0.887))  ## fixed effects

install.packages("coefplot2",repos="http://r-forge.r-project.org")
libr(coefplot2)
coefplot2(model)

#006.4_Multilevel modelling Example--------------------------------------------------
libr(merTools)
m1 <- lmer(y ~ service + lectage + studage + (1|d) + (1|s), data=InstEval)
# shinyMer(m1, simData = InstEval[1:100, ]) # just try the first 100 rows of data
libr(lme4)
libr(arm)

libr(nlme)


lmm.data <- read.table("http://researchsupport.unt.edu/class/Jon/R_SC/Module9/lmm.data.txt", header = TRUE, sep = ",", na.strings = "NA", dec = ".", strip.white = TRUE)
head(lmm.data)

# simple OLS regression; dependant variable is extroversion
OLSexamp <- lm(extro ~ open + agree + social, data = lmm.data)
display(OLSexamp)

# Generalized liniar model; model fit via max likelihood estimation
MLexamp <- glm(extro ~ open + agree + social, data = lmm.data)
display(MLexamp)
AIC(MLexamp) # poor model fit

# fit varying intercept model
MLexamp.2 <- glm(extro ~ open + agree + social + class, data = lmm.data)
display(MLexamp.2)
AIC(MLexamp.2)
MLexamp.3 <- glm(extro ~ open + agree + social + school, data = lmm.data)
display(MLexamp.3)
AIC(MLexamp.3) #School effect much better fit
MLexamp.4 <- glm(extro ~ open + agree + social + school:class, data = lmm.data) # : for interaction between school&class
display(MLexamp.4)
AIC(MLexamp.4)
MLexamp.5 <- glm(extro ~ open + agree + social + school * class - 1, data = lmm.data)
display(MLexamp.5)
AIC(MLexamp.5)

# Random slopes: fit seperate model for each school/class combo & explore par variation among them
require(plyr)
modellist <- dlply(lmm.data, .(school, class), function(x) glm(extro ~ open + agree + social, data = x))
display(modellist[[1]])
display(modellist[[2]])

# group level variable = (1|school); fit linear model with varying intercept group effect using school
MLexamp.6 <- lmer(extro ~ open + agree + social + (1 | school), data = lmm.data)
display(MLexamp.6)
MLexamp.7 <- lmer(extro ~ open + agree + social + (1 | school) + (1 | class), data = lmm.data)
display(MLexamp.7)
MLexamp.8 <- lmer(extro ~ open + agree + social + (1 | school/class), data = lmm.data) # nested group effects: fit mixed effect term for varying intercepts by schools and for classes inside schools
display(MLexamp.8)

# fit varying slope model
MLexamp.9 <- lmer(extro ~ open + agree + social + (1 + open | school/class), data = lmm.data)
display(MLexamp.9)














dendlist <- dendlist()
for(i in seq_along(dis)) {
  try( d <- vegdist(matrixCellCountAdj, method=dis[i]) ) 
  hc <- hclust(d, method = "ward.D")   
  dendlist <- dendlist(dendlist, as.dendrogram(hc))
}
names(dendlist) <- dis

dendlist_cor <- cor.dendlist(dendlist)
require(corrplot)
png (file = paste("Results", centre, "/hclust/corr_", dis[i], "_ward.D.png", sep = "" ), width=900, height=900)
corrplot::corrplot(dendlist_cor, "pie", "lower")
dev.off()


libr(cvxbiclustr)
cobra(t(m))

# X <- lung
X <- matrixCellCountAdj
m <- mean(X[which(colnames(X))=="1-1-1-CIPHE_All",])
X <- log( X / m )
## Create annotation for heatmap
types <- colnames(matrixCellCountAdj)
ty <- as.numeric(factor(types))
cols <- rainbow(length(unique(ty)))
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)
## Construct weights and edge-incidence matrices
phi <- 0.5
k <- 5
wts <- gkn_weights(X,phi=phi,k_row=k,k_col=k)
## Connected Components of Row and Column Graphs
wts$nRowComp
wts$nColComp
#### Initialize path parameters and structures
nGamma <- 5
gammaSeq <- 10**seq(0,3,length.out=nGamma)
## Generate solution path
sol <- cobra_validate(X, wts$E_row, wts$E_col, wts$w_row, wts$w_col, gammaSeq)
ix <- 4
heatmap(sol$U[[ix]],col=hmcols,labRow=NA,labCol=NA,ColSideCol=cols[ty])



#007.1_distance: hclust, tsne (005.0 normalized cell count matrix) -------------------------

dodist <- F
doHC <- T
doTsne <- T
require(colorspace)
dis <- c("canberra","bray","kulczynski", "binomial","cao", "canberra2") # "jaccard","mahalanobis", "gower","altGower","morisita","horn", "euclidean", ""manhattan", "maximum", "binary", "minkowski", "raup", "mountford" (all at same level),)
link <- c("ward.D", "ward.D2", "mcquitty") #, "single", "complete") # , "median", "average", "centroid"

start <- Sys.time()

# Distance Objects
for (centre in centreL[c(1:5)]) { cat("\nCentre: ",centre,sep="")
  start0 <- Sys.time()
  
  setwd(root); setwd(paste("Results",centre,sep=""))
  suppressWarnings(dir.create ("dist"))
  load("sampleMeta.Rdata")
  load("phenoMeta.Rdata")
  load("all0Phen.Rdata")
  if (exists("insigPhenIndex")) {
    if (length(insigPhenIndex)>0 | length(all0Phen)>0) phenoMeta <- phenoMeta[-unique(c(all0Phen,insigPhenIndex)),]
  } else if (length(all0Phen)>0) { phenoMeta <- phenoMeta[-all0Phen,] }
  if (file.exists("insigPhenIndex.Rdata")) load("insigPhenIndex.Rdata")
  # if only do whole matrix i.e. k=max(phenolevel) only
  k <- max(phenoMeta$phenoLevel)+1
  for (mcp in 1:2) {
    # load & fix matrices
    if (mcp==1) {dname <- paste("dist/", dis, "_layer-", str_pad(k, 2, pad = "0"), ".Rdata", sep = "" )
    } else {dname <- paste("dist/", dis, "_layer-", str_pad(k, 2, pad = "0"), "_prop.Rdata", sep = "" )}
    
    if (!all(dname%in%dir("dist/")) | dodist) {
      if (mcp==1) {matrixCellCountAdj <- get(load(paste("matrixCellCountAdj.Rdata",sep="")))
      } else {matrixCellCountAdj <- get(load(paste("matrixCellProp.Rdata",sep="")))}
      if (nrow(matrixCellCountAdj)>ncol(matrixCellCountAdj)) matrixCellCountAdj <- t(matrixCellCountAdj) #transpose; row = sample
      rownames(matrixCellCountAdj) <- sampleMeta$gene
      if (exists("insigPhenIndex")) {
        if (length(insigPhenIndex)>0 | length(all0Phen)>0) matrixCellCountAdj <- matrixCellCountAdj[,-unique(c(all0Phen,insigPhenIndex))]
      } else if (length(all0Phen)>0) { matrixCellCountAdj <- matrixCellCountAdj[,-all0Phen] }
    }
    
    for (i in 1:length(dis)) { 
      start1 <- Sys.time()
      for (k in max(phenoMeta$phenoLevel)+1) {
        
        # get distance object i
        if (mcp==1) {dname <- paste("dist/", dis[i], "_layer-", str_pad(k, 2, pad = "0"), ".Rdata", sep = "" )
        } else {dname <- paste("dist/", dis[i], "_layer-", str_pad(k, 2, pad = "0"), "_prop.Rdata", sep = "" )}
        if (file.exists(dname) & !dodist) { cat("\nLoading Distance Object #", length(dis)-i+1, " ", dis[i]," ",sep="")
          load(dname)
        } else { cat("\nCalculating Distance Object #", length(dis)-i+1, " ", dis[i], sep="")
          # m <- matrixCellCountAdj[,which(phenoMeta$phenoLevel<=k)]
          m <- matrixCellCountAdj
          start2 <- Sys.time()
          #Dist matrix
          require(vegan) # libr(proxy)
          if (dis[i]=="canberra2") {
            m[which(m<0)] <- 1
            try( d <- vegdist(m, method="canberra") )
          } else {
            try( d <- vegdist(m, method=dis[i]) )
          }
          save(d,file=dname)
          cat(" ",TimeOutput(start2), sep="")
        }
        
        start2 <- Sys.time()
        
        # create hclust for distance object i with all links
        if (doHC) {
          require(fastcluster)
          require(dendextend)
          cat("; HClust: ")
          for (j in 1:length(link)) { cat(length(link)-j+1," ",sep="") 
            #Hclust
            try({
              hc <- hclust(d, method=link[j])
              dend <- as.dendrogram(hc)
              dend <- rotate(dend, 1:nrow(sampleMeta))
            })
            
            numplots <- 2; if (!grepl("CIPHE",centre) & !grepl("TCP",centre)) numplots <- numplots+1
            if (mcp==1) {pngname <- paste("dist/", dis[i], "_", link[j], "_layer-", str_pad(k, 2, pad = "0"), "_hclust.png", sep = "" )
            } else {pngname <- paste("dist/", dis[i], "_", link[j], "_layer-", str_pad(k, 2, pad = "0"), "_hclust_prop.png", sep = "" )}
            png (file=pngname , width=1000*numplots, height=8000)
            par(mar=(c(5,5,5,40) + 0.1), mfrow=c(1,numplots))
            
            dend0 <- dend
            gene <- factor(sampleMeta$gene)
            suppressWarnings({ labels_colors(dend) <- rainbow_hcl(length(rev(levels(gene))))[ sort_levels_values( as.numeric(gene)[order.dendrogram(dend)] ) ]  }) 
            labels(dend) <- paste(as.character(gene)[order.dendrogram(dend)], sep = "")
            dendh <- hang.dendrogram(dend,hang_height=0.1)
            plot(dendh, main = paste(centre," samples (labels = true gene); \ndist=",dis[i],"; link=",link[j],sep=""), horiz=TRUE,  nodePar=list(cex = .007))
            # dend <- color_branches(dend, k=length(gene))
            # legend("topleft", legend = gene, fill = rainbow_hcl(length(gene)))
            
            dend <- dend0
            date <- factor(sampleMeta$date)
            suppressWarnings({ labels_colors(dend) <- rainbow_hcl(length(rev(levels(date))))[ sort_levels_values( as.numeric(date)[order.dendrogram(dend)] ) ]  }) 
            labels(dend) <- paste(as.character(date)[order.dendrogram(dend)], sep = "")
            dendh <- hang.dendrogram(dend,hang_height=0.1)
            plot(dendh, main = paste(centre," samples (labels = date); \ndist=",dis[i],"; link=",link[j],sep=""), horiz=TRUE,  nodePar=list(cex = .007))
            
            if (!grepl("CIPHE",centre) & !grepl("TCP",centre)) {
              dend <- dend0
              gender <- factor(sampleMeta$gender)
              suppressWarnings({ labels_colors(dend) <- rainbow_hcl(length(rev(levels(gender))))[ sort_levels_values( as.numeric(gender)[order.dendrogram(dend)] ) ]  }) 
              labels(dend) <- paste(as.character(gender)[order.dendrogram(dend)], sep = "")
              dendh <- hang.dendrogram(dend,hang_height=0.1)
              plot(dendh, main = paste(centre," samples (labels = gender); \ndist=",dis[i],"; link=",link[j],sep=""), horiz=TRUE,  nodePar=list(cex = .007))
            }
            
          }
          graphics.off()
        }
        
        
        # Rtsne plot for distance object i
        if (doTsne) {
          require(Rtsne)
          cat("; Tsne theta: ")
          for (theta in c(0,.5)) { cat(theta," ", sep="")
            tryCatch({
              tsne <- Rtsne(d, is_distance=T, theta=theta)
              colnames(tsne$Y) <- c("x","y")
              # tsnem <- Rtsne(m)
              # rownames(tsnem$Y) <- sampleMeta$gene
              # palette <- choose_palette(pal=rainbow_hcl, n=length(unique(sampleMeta$gene)))
              
              plotTsne <- function(x, continuous=F, main, colBarWidth=.08, colBarTitle="", leg=T) {
                x1 <- rownames(x)
                require(lubridate)
                if (is.Date(x1) | is.timepoint(x1)) {
                  x1 <- as.numeric(x1-min(x1))
                }
                x <- x[order(x1),]
                x1 <- sort(x1)
                if (continuous) {
                  ts <- as.numeric(factor(x1))
                  colour <- heat.colors(length(unique(ts)))
                  plot(x, col=colour[ts], pch=20, main=main)
                  legend.col(col = colour, lev = ts)
                } else {
                  c1 <- 7
                  if (length(unique(x1))/7>25) {
                    c1 <- ceiling(length(unique(x1))/25)
                  }
                  colour <- rainbow(c1)
                  plot(x, t='n', main=main)
                  g <- 1
                  for (d in c(20,1:19,21:25)) {
                    for (c in 1:c1) {
                      points(x[which(x1%in%unique(x1)[g]),], col=colour[c], pch=d)
                      if (g==length(unique(x1))) break
                      g <- g+1
                    }
                    if (g==length(unique(x1))) break
                  }
                  if (leg) legend("topleft", legend=unique(x1), fill=colour[rep(c(1:c1),g)], pch=rep(c(20,1:19,21:25),each=c1)[c(1:g)])
                }
              }
              
              if (mcp==1) {pngname <- paste("dist/", dis[i], "_layer-", str_pad(k, 2, pad = "0"), "_tsne",theta,".png", sep = "" )
              } else {pngname <- paste("dist/", dis[i], "_layer-", str_pad(k, 2, pad = "0"), "_tsne",theta,"_prop.png", sep = "" )}
              width <- 700; height <- 700 # of each plot in png
              perplot <- height/12
              
              
              numplots <- ceiling( length(unique(sampleMeta$gene))/perplot ) # each legend row takes up 18 pixels; prevent legend from overspilling
              col <- numplots+2; if (numplots>1) col <- col+1 # all genes + WT only, >3 samples
              row <- 3; if(!grepl("TCP",centre) & !grepl("CIPHE",centre) & col<6) row <- row+1; 
              
              png(pngname, width=width*col, height=height*row)
              par(mfrow=c(row,col))
              
              ftWTGT <- c("+_+","+_Y"); if (!grepl("Sanger",centre)) ftWTGT <- "1-1-1-WildType_01_All"; if (centre=="H") ftWTGT <- "FR-FCM-ZYCB-WildType_01_All"
              sampleNo <- 3
              g <- getGTindex(sampleMeta$gene, ftWTGT)
              ftGT <- g$ftGT; ftWTIndex <- g$ftWTIndex; ftKOGT <- g$ftKOGT; ftKOIndex <- g$ftKOIndex; ftKOgoodGTIndex <- g$ftKOgoodGTIndex; rm(g)  #Index of unique KO genotypes that has 3+ samples available ftKOIndex(ftKOGTIndexIndex)
              ftKOIndexL <- NULL; for (l in 1:length(ftKOIndex)) { ftKOIndexL[l] <- length(ftKOIndex[[l]]) }
              
              rownames(tsne$Y) <- sampleMeta$gene
              for (leg in c(T,F)) {
                plotTsne(tsne$Y[which(sampleMeta$gene%in%ftWTGT),], leg=leg, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta,"; WT only", sep=""))
                plotTsne(tsne$Y[sampleMeta$gene%in%c(ftWTGT,as.character(ftKOGT)[which(ftKOIndexL>3)]),], leg=leg, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta,"; gene with ",sampleNo,"< samples", sep=""))
                if (numplots>1) plotTsne(tsne$Y, leg=leg, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta,"; all gene", sep=""))
                for (j  in 1:numplots) {
                  e <- min(j*perplot, length(unique(sampleMeta$gene)))
                  gt <- unique(sampleMeta$gene)[c(((j-1)*perplot+1):e)]
                  plotTsne(tsne$Y[sampleMeta$gene%in%gt,], leg=leg, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta, sep=""))
                }
              }
              
              require(lubridate)
              rownames(tsne$Y) <- ymd(sampleMeta$date)
              plotTsne(tsne$Y[which(sampleMeta$gene%in%ftWTGT),], continuous=T, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta,"; \ndays since ", min(sampleMeta$date),"; WT only", sep=""))
              plotTsne(tsne$Y[sampleMeta$gene%in%c(ftWTGT,as.character(ftKOGT)[which(ftKOIndexL>3)]),], continuous=T, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta,"; \ndays since ", min(sampleMeta$date),"; gene with ",sampleNo,"< samples", sep=""))
              plotTsne(tsne$Y, continuous=T, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta, "\n days since ", min(sampleMeta$date), sep=""))
              
              if (!grepl("TCP",centre) & !grepl("CIPHE",centre) ) {
                rownames(tsne$Y) <- as.character(sampleMeta$gender)
                plotTsne(tsne$Y[which(sampleMeta$gene%in%ftWTGT),], main=paste("TSNE ",dis[i]," distance matrix; gender; WT only", sep=""))
                plotTsne(tsne$Y[sampleMeta$gene%in%c(ftWTGT,as.character(ftKOGT)[which(ftKOIndexL>3)]),], main=paste("TSNE ",dis[i]," distance matrix; gender; gene with ",sampleNo,"< samples", sep=""))
                plotTsne(tsne$Y, main=paste("TSNE ",dis[i]," distance matrix; gender", sep=""))
              }
              
              graphics.off()
              
              # +_Y at 353; +_+ at 2
              require(cluster)
              require(geometry)
              
              #           cNo <- 7
              #           km <- kmeans(tsne$Y,cNo)
              #           for (i in 1:cNo) {
              #             line(convhulln(tsne$Y[which(km$cluster==i),]))
              #           }
              #           convhulln
              #           
              #           require(RDRToolbox)
              #           isomap <- Isomap(data=matrixCellCountAdj, dims=2, k=10)
              #           ecb(isomap$dim2, unique(sampleMeta$gene), paste("Isomap ",dis[i]," distance matrix; theta=0", sep="") )
              #           ecb(isomap$dim2, ftWTGT, paste("Isomap ",dis[i]," distance matrix; theta=0; WildType only", sep=""))
              
            }, error = function(err) { graphics.off()})
          }
        }
        
        cat(" viz ",TimeOutput(start2)," ",sep="")
      }
      cat(" dis ",TimeOutput(start1)," ",sep="")
    }
  }
  cat(" centre ",TimeOutput(start0)," ",sep="")
}
cat("\ntotal ",TimeOutput(start)," ",sep="")




require(RDRToolbox)
libr(rgl)
## Isomap ##



# libr(dbscan)
# oc <- optics(d, eps=10, minPts=3) # no results (gravitate towards each other)



#007.2_distance manual: FPM -----------------------------------------

# get significant phenotypes from PVal matrix for comparison
pval <- .05
sigPheno <- NULL
for (i in 1:ncol(matrixPVal)) {
  sigIndex <- which(matrixPVal[,i]<pval)
  sigIndex <- sigIndex[order(matrixPVal[sigIndex,i])]
  if (length(sigIndex)>0) {
    sigPheno[[i]] <- sigIndex
  }
}

#Manually make distance matrix (distance per single phen) See what is most different
dis <- "canberra" # (0,a)==>1; (0,0)==>NA

ftWTGT <- c("+_+", "+_Y")
ftWTGTi <- which(sampleMeta$gene%in%ftWTGT)
ftKOGT <- unique(sampleMeta$gene[-ftWTGTi])
ftKOGTi <- NULL
for (i in length(ftKOGT):1) { ftKOGTi[[i]] <- c(which(sampleMeta$gene==ftKOGT[i])) }
for (i in ftKOGT) {
  d <- featDistList(matrixCellCountAdj[ftWTGTi,], matrixCellCountAdj[which(sampleMeta$gene==i),], dis)
  save(d,file=paste("featDistList_",ftWTGT,".Rdata",sep=""))
  do <- NULL
  for (j in length(d)) {
    for(k in length(d[[j]])) {
      do <- order
    }
  }
  #sort all lists
  #find most representative list
  #find outliers
  #plot each list onto a single plot
}

# #Manually make distance matrix (hardcoded per single phen) and get list for each pairwise comparison
start <- Sys.time()
canberra2 <- 2 #Make 0 into 1; so distance isn't just 1
for (centre in centreL[c(5)]) {
  setwd(root); setwd(paste("Results",centre,sep=""))
  suppressWarnings(dir.create("dist/phenoTypes"))
  
  m <- get(load("matrixCellCountAdj.Rdata")); if(nrow(m)>ncol(m)) m <- t(m); if (canberra2==2) m[which(m<1)] <- 1
  for (i in nrow(m):2) { cat(i," ",sep="") #compare each row; e.g. 6561 pheno, 100 samples = 9seconds; 200 samples = 39 seconds
    t <- NULL
    for (j in 1:(i-1)) {
      t[[j]] <- abs(m[i,]-m[j,])/(m[i,]+m[j,])
      t[[j]][is.nan(t[[j]])] <- 0
      # plot(density(t))
      # plot(sort(m[j,which(m[j,]>0 & t==1)]))
      # points(sort(m[i,which(m[i,]>0 & t==1)]))
    }
    save(t, file=paste("dist/phenoTypes/", str_pad(i, 4, pad = "0"), ".Rdata", sep="") )
  }
  
  d <- matrix(0,nrow=nrow(m), ncol=nrow(m))
  tall <- NULL
  for (i in nrow(m):2) {
    load(paste("dist/phenoTypes/", str_pad(i, 4, pad = "0"), ".Rdata", sep=""))
    for (j in 1:(i-1)) {
      d[i,j] <- d[j,i] <- sum(t[[j]])
      tall <- cbind(tall, t[[j]])
    }
  }
  d <- as.dist(d)
  save(d,file="dist/phenoTypes/d.Rdata")
  save(tall, file="dist/phenoTypes/tall.R")
}
TimeOutput(start)

load("markers.Rdata")
load("phenoMeta.Rdata")
require(arules)
tr <- as(createTransactionMatrix(phenoMeta$phenoCode, markers), "transactions")
for (i in nrow(m):2) {
  load(paste("dist/phenoTypes/", str_pad(i, 4, pad = "0"), ".Rdata", sep=""))
  for (j in 1:(i-1)) {
    transactionInfo(t) <- data.frame(phenoCode=phenoMeta$phenoCode, weight=t[[j]] )
    s <- weclat(t, parameter=list(support=0, maxLen=length(markers), type="absolute"), control=list(verbose=T))
    inspect(sort(s))
  }
}




















#SampleMeta stats -------------------------------------------

#male female distribution in each gene
fPerc <- NULL
for (i in 1:length(unique(sampleMeta$gene))) {
  rows <- which(sampleMeta$gene==unique(sampleMeta$gene[i]))
  fPerc[i] <- length(which(grepl("F",sampleMeta$gender[rows]))) / length(rows)
}
plot(unique(sampleMeta$gene),fPerc)
noFGT <- unique(sampleMeta$gene)[which(fPerc==0)]
noMGT <- unique(sampleMeta$gene)[which(fPerc==1)]

#date distribtution in each gene
plot(sampleMeta$gene[order(sampleMeta$gene)],sampleMeta$date[order(sampleMeta$gene)])

#number of samples of each gene
ftKOIndexL <- NULL
less3samp <- NULL
minSamp <- 3
for (i in 1:length(ftKOIndex)) {
  ftKOIndexL[i] <- length(ftKOIndex[[i]])
  if (ftKOIndexL[i]<minSamp) less3samp <- c(less3samp,ftKOIndex[[i]])
}
names(ftKOIndexL) <- ftKOGT


phenoNo <- 654
#distribution of 1 cell population count across all genes
plot(factor(sampleMeta$gene),matrixCellCountAdj[,phenoNo],ylim=c(0,quantile(matrixCellCountAdj[,phenoNo],probs=c(.99))))
#distribution of 1 cell population count across all genes on dates
m <- data.frame(gene=factor(sampleMeta$gene), date=as.Date(sampleMeta$date, "%d-%b-%y"), phenCount=matrixCellCountAdj[phenoNo,])
colour <- rainbow(length(unique(m$gene)))
plot(m$date, m$phenCount, cex=.3)
axis(1, m$date, format(m$date, "%b %d"), cex.axis = .7)
for (i in 1:length(unique(m$gene))[-which(unique(m$gene)%in%ftWTGT)]) { # +_+ = 156, +_Y = 356; Sanger SPLEEN
  geneIndex <- which(m$gene==unique(m$gene)[i])
  geneIndex <- geneIndex[order(m$date[geneIndex])]
  points(m[geneIndex,c(2,3)], col=colour[i], cex=.3)
}
points(m[which(m$gene%in%ftWTGT),c(2,3)], cex=1)

refColumn <- which( matrixCellCountAdj[1,]==median(matrixCellCountAdj[1,which(sampleMeta$gene%in%ftWTGT)]) )
m <- matrixCellCountAdj[order(matrixCellCountAdj[,refColumn]),]





require(arrayMvout)
GoodFakeData <- NULL
for ( x in 1:4 ) { GoodFakeData <- cbind(GoodFakeData, runif(nrow(m),5,5.5))  }
outInd <- as.numeric(rownames( ArrayOutliers(data.frame(cbind(m1,GoodFakeData[c(1:nrow(m)),])), alpha=0.5)$outl ))

