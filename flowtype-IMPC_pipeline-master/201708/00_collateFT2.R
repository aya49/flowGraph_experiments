# aya43@sfu.ca 20151228 -- still need to standardize, don't use for now...
# Reads in IMPC flowtype files and collates/save cell counts into matrix

root = "~/projects/IMPC"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
panelL = c("P2")
centreL = c("Sanger_SPLEEN")#,"Sanger_MLN","CIPHE","TCP","H")
setwd(root)

options(stringsAsFactors=F)
options(device="cairo")


#Input
data_dir = "/mnt/f/Brinkman group/current/Alice/IMPC/data" #Sanger saves flowtype differently from all other centres
ft_dir = paste0(data_dir, "/", panelL, "/", centreL, "/ft")
csv_dir = paste0(data_dir, "/", panelL, "/", centreL, "/sampleMeta.csv") #differenct for each centre, columns aren't standardized yet either, fix after standardization
countcutoff = 0 #delete columns/cellPops with all values equal or below countcutoff
levelcutoff = 8 #Inf if delete no layers; >levelcutoff are deleted

#Output
markers_dir = paste0(result_dir, "/", panelL, "/", centreL, "/markers.Rdata")
phenoMeta_dir = paste0(result_dir, "/", panelL, "/", centreL, "/phenoMeta")
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta", sep="")

matrixCount_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixCount", sep="")
matrixProp_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixProp", sep="")
excludedGT_dir = paste(result_dir, "/", panelL, "/", centreL, "/excludedGT.Rdata", sep="")

libr(Matrix)
libr(stringr)
libr(lubridate)
libr(flowCore)
libr(flowType)
libr(foreach)
libr(doMC)
source("~/projects/IMPC/code/_funcAlice.R")
# source("code/load_manual_results.R") # for CIPHE



start = Sys.time()

no_cores = detectCores() - 1
registerDoMC(no_cores)

for (ci in 1:length(paste0(panelL,centreL))) {
  centre = paste0(panelL," ",centreL)[ci]
  start1 = Sys.time()
  cat("\n",centre)
  #make file paths, markers, hardcode markers, fix after standardization --------------------------------
  ft_dirall = sort(dir(ft_dir[ci], pattern=".Rdata", all.files=TRUE, full.names=TRUE, recursive=TRUE) )
  if (grepl("Sanger", centre)) {
    ftGT = folderNames(ft_dirall)
    ftFileNames = fileNames(ft_dirall, "Rdata")
    ft = get(load(ft_dirall[1]))
    markers = c("CD44","CD62L","CD25","CD8","KLRG","CD5","CD45","CD161","CD4","GITR","TCRd")
  } else {
    ftGT = NULL
    ftGenotype = fileNames(ft_dirall, "Rdata")
    ftFileNames = NULL
    ft = get(load(ft_dirall[1]))[[1]]$ftype
    markers = getMarkers(names(ft@CellFreqs)[length(ft@CellFreqs)])
    csv_dirall = sort(dir(ft_dir[ci], pattern=".csv", all.files=TRUE, full.names=TRUE, recursive=TRUE) )
  }
  save(markers, file=markers_dir[ci])
  
  #loading & compiling flowTypes --------------------------------
  excludedGT = NULL
  if (grepl("Sanger",centre)) { #Sanger TCell panel MLN & SPLEEN organ
    matrixCountList = foreach(i=1:length(ft_dirall)) %dopar% {
      ft = get(load(ft_dirall[i]))
      npheno = length(ft@CellFreqs)
      cat("\n", i, ". Loaded file: ", ftGT[i], "/", ftFileNames[i], "; # of phenotypes x samples: ", npheno, " x ", length(ft), "", sep ="")
      return(ft@CellFreqs)
    }
    ft = get(load(ft_dirall[1]))
  } else {
    # all files of a gene merged together; CIPHE has different number of markers (64 genes); TCP (61 genes)
    # 08. 1-1-1-Baiap2l2_01_All; # of phenotypes x samples: 59049 x 6; different cell pop #; used markers (CD5, CD161, CD4, TCRd, GITR, CD44, CD62L, CD25, CD24, KLRG)
    # 16. 1-1-1-Cxcr7_01_All; # of phenotypes x samples: 59049 x 6	; different cell pop #; used markers (CD5, CD161, CD4, TCRd, GITR, CD44, CD62L, CD25, CD24, KLRG)
    # 21. 1-1-1-Dnalc4_01_All; # of phenotypes x samples: 6561 x 6	; different cell pop #; used markers (CD5, CD161, CD4, TCRd, CD44, CD62L, CD25, CD24)
    # 41. 1-1-1-Otub1_01_All; # of phenotypes x samples: 19683 x 6	; different cell pop #; used markers (CD5, CD161, CD4, CD8, TCRd, CD44, CD62L, CD25, CD24)
    # 43. 1-1-1-Pax4_01_All; # of phenotypes x samples: 6561 x 6	  ; different cell pop #; used markers (CD5, CD161, CD4, TCRd, CD44, CD62L, CD25, CD24)
    # 47. 1-1-1-Ptdss1_01_All; # of phenotypes x samples: 19683 x 6	; different cell pop #; used markers (CD5, CD161, CD4, CD8, TCRd, CD44, CD62L, CD25, CD24)
    # 48. 1-1-1-Rab19_01_All; # of phenotypes x samples: 19683 x 6	; different cell pop #; used markers (CD5, CD161, CD4, CD8, TCRd, CD44, CD62L, CD25, CD24)
    # 52. 1-1-1-Setbp1_01_All; # of phenotypes x samples: 19683 x 6	; different cell pop #; used markers (CD5, CD161, CD4, CD8, TCRd, CD44, CD62L, CD25, CD24)
    # 57. 1-1-1-Snx29_01_All; # of phenotypes x samples: 59049 x 6	; different cell pop #; used markers (CD5, CD161, CD4, TCRd, GITR, CD44, CD62L, CD25, CD24, KLRG)
    
    # row = 1
    result = foreach (i = 1:length(ft_dirall)) %dopar% {
      excludedGT = NULL
      ftFileNames = NULL
      matrixCountList = NULL
      ftt = get(load(ft_dirall[i]))
      cat("\n", i, ". Loaded file: ", ftGenotype[i], "; # of phenotypes x samples: ", length(ftt[[1]]$ftype@CellFreqs), " x ", length(ftt), "", sep ="")
      if (i!=1) { #check if same number of markers
        if (npheno!=length(ftt[[1]]$ftype@CellFreqs)) {
          cat("\t; different cell pop #, skipped; used markers ("); cat(getMarkers(names(ftt[[1]]$ftype@CellFreqs)[length(ftt[[1]]$ftype@CellFreqs)]), sep=", "); cat(")", sep="")
          excludedGT <- ftGT[i]
        }
      }
      if(is.null(excludedGT)) {
        npheno = length(ftt[[1]]$ftype@CellFreqs)
        fn = as.character(names(read.csv(csv_dirall[i], nrows=1)))
        ftFileNames = gsub("[.]", "-", gsub(".fcs", "", gsub("^X", "", fn[-which(is.na(fn) | fn=="X")])))
        for (j in 1:length(ftt)) {
          ft = ftt[[j]]$ftype
          ftGT = c(ftGT, ftGenotype[i])
          matrixCountList[[j]] = ft@CellFreqs
          # row = row+1
        }
        matrixCountList = do.call(rbind, matrixCountList)
        
        return(list(matrixCountList=matrixCountList, ftFileNames=ftFileNames, excludedGT=excludedGT))
      }
      ft = get(load(ft_dirall[1]))[[1]]$ftype
    }
    
    excludedGT = list()
    ftFileNames = list()
    matrixCountList = list()
    for (i in 1:length(result)) {
      if (!is.null(result[[i]]$excludedGT)) excludedGT[[i]] = result[[i]]$excludedGT
      if (!is.null(result[[i]]$ftFileNames)) ftFileNames[[i]] = result[[i]]$ftFileNames
      if (!is.null(result[[i]]$matrixCountList)) matrixCountList[[i]] = result[[i]]$matrixCountList
    }
    excludedGT <- unlist(excludedGT[which(vapply(excludedGT, Negate(is.null), NA))])
    ftFileNames <- unlist(ftFileNames[which(vapply(ftFileNames, Negate(is.null), NA))])
    matrixCountList <- matrixCountList[which(vapply(matrixCountList, Negate(is.null), NA))]
  }
  matrixCount = do.call(rbind, matrixCountList)
  rownames(matrixCount) = ftFileNames
  matrixProp = matrixCount/matrixCount[,1]
  colnames(matrixCount) <- colnames(matrixProp) <- phenotype <- unlist(lapply(ft@PhenoCodes, function(x){return( decodePhenotype(x, markers, ft@PartitionsPerMarker) )}))
  
  #make phenoMeta --------------------------------
  phenoMeta = data.frame(phenotype,stringsAsFactors=F)
  phenoMeta$phenotype = as.character(phenoMeta$phenotype); rm(phenotype)
  phenoMeta$phenocode = ft@PhenoCodes
  phenoMeta$phenolevel = unlist(lapply(phenoMeta$phenocode, function(x){return(length(markers) - charOccurences("0", x)) } ))
  
  #trim matrix/phenoMeta & save ------------------------------
  trimColIndex1 <- NULL
  trimColIndex2 <- NULL
  tr1 = which(apply(matrixCount[-1,], 2, function(x) all(x<=countcutoff))==T) #delete cols of all 0
  tr2 = which(phenoMeta$phenolevel>levelcutoff) #delete cols of too high level
  if (length(tr1)>0) trimColIndex1 = tr1
  if (length(tr2)>0) trimColIndex2 = tr2
  trimColIndex = union(trimColIndex1,trimColIndex2)
  if (length(trimColIndex)>0) {
    matrixCount <- matrixCount[,-trimColIndex]
    matrixProp <- matrixProp[,-trimColIndex]
    phenoMeta <- phenoMeta[-trimColIndex,]
  }
  
  save(phenoMeta, file=paste0(phenoMeta_dir[ci],".Rdata"))
  
  save(matrixCount, file=paste0(matrixCount_dir[ci],".Rdata"))
  write.csv(matrixCount, file=paste0(matrixCount_dir[ci],".csv"))
  save(matrixProp, file=paste0(matrixProp_dir[ci],".Rdata"))
  write.csv(matrixProp, file=paste0(matrixProp_dir[ci],".csv"))
  if (!is.null(excludedGT)) save(excludedGT, file=excludedGT_dir[ci])
  
  
  
  
  #make sampleMeta, hardcode column names of wanted columns, fix after standardization --------------------------------
  sampleMeta = NULL
  sampleMetaList = NULL
  if (grepl("Sanger",centre)) { #Sanger TCell panel MLN & SPLEEN organ
    sampleMetaTemp = data.frame(lapply(read.csv(csv_dir[ci])[,c("Genotype","Assay.Date","Gender","Label.Barcode", "Colony.Prefix","Birth.Date","Coat.Colour")], as.character), stringsAsFactors=F)
    sampleMetaTemp[,"Genotype"] = gsub("/","_",sampleMetaTemp[,"Genotype"])
    barcodes <- str_extract(ftFileNames, "L[0-9]+")
    for (i in 1:length(barcodes)) {
      row = intersect(grep(barcodes[i],sampleMetaTemp[,"Label.Barcode"]), which(sampleMetaTemp[,"Genotype"]==ftGT[i]))[1]
      sampleMetaList[[i]] = c(ftFileNames[i], unlist(sampleMetaTemp[row,]))
    }
    sampleMeta = rbind(sampleMeta, do.call(rbind, sampleMetaList))
    sampleMeta = as.data.frame(sampleMeta, stringsAsFactors=F)
    colnames(sampleMeta) = c("fileName", "gene", "date", "gender", "barcode", "colony", "birth_date", "fur")
    sampleMeta$barcode = barcodes
    sampleMeta$date = dmy(sampleMeta$date)
    sampleMeta$birth_date = dmy(sampleMeta$birth_date)
    
    #check duplicate and delete
    duplicate.index <- which(duplicated(sampleMeta$barcode)==TRUE)
    if (length(duplicate.index)>0) {
      sampleMeta = sampleMeta[-duplicate.index,]
      matrixCount = matrixCount[-duplicate.index,]
    }
    
    # Sanger MLN TCell panel has different date formats
    # sampleMeta$date = dmy(sampleMeta$date)
    # if (grepl("MLN",centre)) {
    #   sampleMeta$date = as.character(ymd(sampleMeta$date))
    #   if (grepl("SangerTCellMLN", centre)) {
    #     d = as.character(sampleMeta$date)
    #     for (i in 1:length(d)) {
    #       if (nchar(sampleMeta$date[i])==8) {
    #         d[i] = paste0("20",d[i],collapse="")
    #       }
    #     }
    #     sampleMeta$date = d
    #   }
    # }
    
  } else if (grepl("CIPHE", centre)) {
    sampleMetaTemp = data.frame(lapply(read.csv(csv_dir[ci])[,c("Date","SEX..M.F.","Sample.N.")], as.character))
    sampleMetaTempDate = dmy(sampleMetaTemp[,"Date"])
    ftFileNamesSplit = lapply(c(1:length(ftFileNames)), function(x) {return(strsplit(ftFileNames[x],"[_]")[[1]])} )
    date = ymd(sapply(c(1:length(ftFileNamesSplit)), function(x) {return(ftFileNamesSplit[[x]][1])} ))
    sample = as.numeric(sapply(c(1:length(ftFileNamesSplit)), function(x) {return(ftFileNamesSplit[[x]][3])} ))
    for (i in 1:length(ftFileNames)) {
      row = which(sampleMetaTempDate==date[i] & sampleMetaTemp[,"Sample.N."]==sample[i])
      sampleMetaList[[i]] = c(ftFileNames[i], ftGT[i], unlist(sampleMetaTemp[row,]))
    }
    sampleMeta = rbind(sampleMeta, do.call(rbind, sampleMetaList))
    sampleMeta = as.data.frame(sampleMeta, stringsAsFactors=F); rm(sampleMetaTemp)
    colnames(sampleMeta) = c("filename", "gene", "date", "gender", "sample")
    sampleMeta$date = dmy(sampleMeta$date)
    
  } else if (grepl("TCP", centre)) {
    sampleMetaTemp = data.frame(lapply(read.csv(csv_dir[ci])[,c("expDate", "sex", "SampleID", "strainCode")], as.character))
    sampleMetaTemp2 = data.frame(lapply(read.csv(csv_dirtcp)[,c("dates", "list.date.short")], as.character))
    sampleMetaTemp[,"expDate"] = mdy(sampleMetaTemp[,"expDate"])
    
    sampleMetaTemp2[,"list.date.short"] = gsub(".labelled","-labelled", gsub(".fcs","",sampleMetaTemp2[,"list.date.short"]))
    rows = match(ftFileNames, sampleMetaTemp2[,"list.date.short"])
    
    ftFileNamesSplit = lapply(c(1:length(ftFileNames)), function(x) {return(strsplit(ftFileNames[x],"[_]")[[1]])} )
    strain = sapply(c(1:length(ftFileNamesSplit)), function(x) {return(ftFileNamesSplit[[x]][3])} )
    sample = sapply(c(1:length(ftFileNamesSplit)), function(x) {return(ftFileNamesSplit[[x]][4])} )
    date = ymd(sampleMetaTemp2[rows,"dates"])
    gender = NULL
    for (i in 1:length(ftFileNames)) {
      row = which(sampleMetaTemp[,"expDate"]==date[i] & sampleMetaTemp[,"SampleID"]==sample[i] & sampleMetaTemp[,"strainCode"]==strain[i])[1]
      gender = append(gender, sampleMetaTemp[row,"sex"])
    }
    sampleMeta = as.data.frame(cbind(ftFileNames, ftGT, date, gender, sample, strain), stringsAsFactors=F)
    colnames(sampleMeta) = c("filename", "gene", "date", "gender", "sample", "strain")
    sampleMeta$date = date
    
  } else if (grepl("H", centre)) {
    sampleMetaTemp = data.frame(lapply(read.csv(csv_dir[ci])[,c("Filename", "Date","Gender","IMPCSpecimenCode")], as.character))
    sampleMetaTemp[,"Filename"] = gsub("[.]", "-", gsub(".fcs","", gsub("[-%]",".", sampleMetaTemp[,"Filename"])))
    rows = match(ftFileNames, sampleMetaTemp[,"Filename"])
    sampleMeta = cbind(sampleMetaTemp[rows,"Filename"], ftGT, sampleMetaTemp[rows,c("Date","Gender","IMPCSpecimenCode")])
    sampleMeta = as.data.frame(sampleMeta, stringsAsFactors=F); rm(sampleMetaTemp)
    colnames(sampleMeta) = c("filename", "gene", "date", "gender", "specimen")
    sampleMeta$date = ymd(sampleMeta$date)
    
    # CIPHE Manula count import
    #   CSVfile = paste("Results", centre, "/IMPC_20142015_merged_Panel1.csv", sep="")
    #   matrixCSV = read.csv(CSVfile)
    #   matrixManualCounts = matrix(0, nrow=55, ncol=ncol(matrixCount))
    #   mcIndex = c(19:ncol(matrixCSV))[which(c(19:ncol(matrixCSV)) %% 2 != 0)]
    #   rownames(matrixManualCounts) = sapply(mcIndex, function(x) {return( gsub("neg", "-", gsub("pos", "+", strsplit(colnames(matrixCSV)[x], "\\.\\.")[[1]][1]) ) ) })
    #   
    #   missingMetaIndex = NULL
    #   for (i in 1:nrow(sampleMeta)) {
    #     row = which(matrixCSV[,1]==sampleMeta$sampleNum[i] & as.Date(matrixCSV[,3], "%d-%m-%Y")==sampleMeta$date[i])
    #     if (length(row)==0) {missingMetaIndex = c(missingMetaIndex, i); next()}
    #     sampleMeta[i,c(5:7)] = matrixCSV[row,c(2,10,9)]
    #     matrixManualCounts[,i] = as.vector( as.numeric(matrixCSV[row,mcIndex] ))
    #   }
    #   save(matrixManualCounts, file=paste("Results", centre, "/matrixManualCounts.Rdata", sep=""))
    
  }
  g = rep(1,nrow(sampleMeta)); g[grep("F",sampleMeta$gender)] = 0; g[which(is.na(sampleMeta$gender))] = NA
  sampleMeta$gender = g # Female=0; Male=1
  
  save(sampleMeta, file=paste0(sampleMeta_dir[ci],".Rdata"))
  write.csv(sampleMeta, file=paste0(sampleMeta_dir[ci],".csv"))
  
  cat("\nTotal time used to load cell count for ", centre, ": ", TimeOutput(start), "\n", sep="")
}

cat("\nTotal time used to load cell count: ", TimeOutput(start), "\n", sep="")

