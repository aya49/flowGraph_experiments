## input: flowtype file,  meta data paths
## output: feat_file_cell_count, meta_file, meta_cell
## process: 
## - takes flowtype output files and compiles them together to create cell count matrix
## - reformats meta file (meta info for fcm files)
## - creates cell meta file (meta info for fcm cell populations)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metrics"
setwd(root)

panel = c("panel1") # panel2 is sanger-spleen only
centre = c("sanger-spleen")#,"sanger-mln","ciphe","tcp","harwell")

result_dir = paste0(root, "/result/impc_",panel,"_",centre); suppressWarnings(dir.create (result_dir, recursive=T))
# result_dir = paste0("result/", panelL, "/", centreL); suppressWarnings(dir.create (result_dir))
data_dir = paste0("/mnt/f/Brinkman group/current/Alice/IMPC/data/", panel, "/", centre) #Sanger saves flowtype differently from all other centres


## input directory
ft_dir = paste0(data_dir, "/ft")
csv_dir = paste0(data_dir, "/sampleMeta.csv") #different for each centre, columns aren't standardized yet either, fix after standardization
if (centre=="tcp") csv_dirtcp = paste0(data_dir, "/date.code.spreadsheet.csv") #second meta spreadsheet for tcp...


## output directory (see last section, output split by panel)
meta_dir = paste0(result_dir,"/meta"); dir.create(meta_dir, showWarnings=F)
meta_cell_dir = paste(meta_dir, "/cell", sep="")
meta_file_dir = paste(meta_dir, "/file", sep="")
markers_dir = paste(meta_dir, "/cell_markers", sep="")

feat_dir = paste(result_dir, "/feat", sep=""); dir.create(feat_dir, showWarnings=F)
feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")
feat_file_cell_prop_dir = paste(feat_dir, "/file-cell-prop", sep="")

excludedGT_dir = paste(result_dir, "/excludedGT.Rdata", sep="") #IMPC ONLY


## libraries
source("source/_funcAlice.R")
# source("code/load_manual_results.R") # for CIPHE
libr(c("Matrix", "stringr", "lubridate",
     "flowCore", "flowType",
     "foreach", "doMC"))


## cores
no_cores = detectCores() - 10
registerDoMC(no_cores)

## options
writecsv = F
options(stringsAsFactors=F)
options(device="cairo")
countThres = 0 #delete columns/rows where all values equal or below countThres
levelThres = 8 #Inf if delete no layers; >levelcutoff are deleted

IMPC_control = "[+]_[+]|[+]_Y"



start = Sys.time()

## make file paths, markers, hardcode markers, fix after standardization --------------------------------
ft_dirall = sort(list.files(ft_dir, pattern=".Rdata", all.files=TRUE, full.names=TRUE, recursive=TRUE) )
if (grepl("sanger", centre)) {
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
  csv_dirall = sort(dir(ft_dir, pattern=".csv", all.files=TRUE, full.names=TRUE, recursive=TRUE) )
}
save(markers, file=paste0(markers_dir,".Rdata"))


## loading & compiling flowTypes --------------------------------
excludedGT = NULL
if (grepl("sanger",centre)) { #Sanger TCell panel MLN & SPLEEN organ
  feat_file_cell_countList = foreach(i=1:length(ft_dirall)) %dopar% {
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
    feat_file_cell_countList = NULL
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
        feat_file_cell_countList[[j]] = ft@CellFreqs
        # row = row+1
      }
      feat_file_cell_countList = do.call(rbind, feat_file_cell_countList)
      
      return(list(feat_file_cell_countList=feat_file_cell_countList, ftFileNames=ftFileNames, excludedGT=excludedGT))
    }
    ft = get(load(ft_dirall[1]))[[1]]$ftype
  }
  
  excludedGT = list()
  ftFileNames = list()
  feat_file_cell_countList = list()
  for (i in 1:length(result)) {
    if (!is.null(result[[i]]$excludedGT)) excludedGT[[i]] = result[[i]]$excludedGT
    if (!is.null(result[[i]]$ftFileNames)) ftFileNames[[i]] = result[[i]]$ftFileNames
    if (!is.null(result[[i]]$feat_file_cell_countList)) feat_file_cell_countList[[i]] = result[[i]]$feat_file_cell_countList
  }
  excludedGT <- unlist(excludedGT[which(vapply(excludedGT, Negate(is.null), NA))])
  ftFileNames <- unlist(ftFileNames[which(vapply(ftFileNames, Negate(is.null), NA))])
  feat_file_cell_countList <- feat_file_cell_countList[which(vapply(feat_file_cell_countList, Negate(is.null), NA))]
}
feat_file_cell_count = do.call(rbind, feat_file_cell_countList)
rownames(feat_file_cell_count) = ftFileNames
feat_file_cell_prop = feat_file_cell_count/feat_file_cell_count[,1]
colnames(feat_file_cell_count) <- colnames(feat_file_cell_prop) <- phenotype <- unlist(lapply(ft@PhenoCodes, function(x){return( decodePhenotype(x, markers, ft@PartitionsPerMarker) )}))



## make meta_cell --------------------------------
meta_cell = getPhen(phenotype)


## make meta_file, hardcode column names of wanted columns, fix after standardization --------------------------------
meta_file = NULL
meta_fileList = NULL
if (grepl("sanger",centre)) { #Sanger TCell panel MLN & SPLEEN organ
  meta_fileTemp = data.frame(lapply(read.csv(csv_dir)[,c("Genotype","Assay.Date","Gender","Label.Barcode", "Colony.Prefix","Birth.Date","Coat.Colour")], as.character), stringsAsFactors=F)
  meta_fileTemp[,"Genotype"] = gsub("/","_",meta_fileTemp[,"Genotype"])
  barcodes <- str_extract(ftFileNames, "L[0-9]+")
  for (i in 1:length(barcodes)) {
    row = intersect(grep(barcodes[i],meta_fileTemp[,"Label.Barcode"]), which(meta_fileTemp[,"Genotype"]==ftGT[i]))[1]
    meta_fileList[[i]] = c(ftFileNames[i], unlist(meta_fileTemp[row,]))
  }
  meta_file = rbind(meta_file, do.call(rbind, meta_fileList))
  meta_file = as.data.frame(meta_file, stringsAsFactors=F)
  colnames(meta_file) = c("id", "gene", "date", "gender", "barcode", "colony", "birth_date", "fur")
  meta_file$barcode = barcodes
  meta_file$date = dmy(meta_file$date)
  meta_file$birth_date = dmy(meta_file$birth_date)
  
  #check duplicate and delete
  duplicate.index <- which(duplicated(meta_file$barcode)==TRUE)
  if (length(duplicate.index)>0) {
    meta_file = meta_file[-duplicate.index,]
    feat_file_cell_count = feat_file_cell_count[-duplicate.index,]
  }
  
  # Sanger MLN TCell panel has different date formats
  # meta_file$date = dmy(meta_file$date)
  # if (grepl("MLN",centre)) {
  #   meta_file$date = as.character(ymd(meta_file$date))
  #   if (grepl("SangerTCellMLN", centre)) {
  #     d = as.character(meta_file$date)
  #     for (i in 1:length(d)) {
  #       if (nchar(meta_file$date[i])==8) {
  #         d[i] = paste0("20",d[i],collapse="")
  #       }
  #     }
  #     meta_file$date = d
  #   }
  # }
  
  gene0 = meta_file$gene_ori = meta_file$gene
  # gene0 = gene0[!grepl(" KO$",gene0)]
  gene0[grepl(" KO$",gene0)] = sapply(str_split(gene0[grepl(" KO$",gene0)]," "), function(x) x[1])
  # gene0 = gene0[!grepl("failed", gene0)]
  gene0[grepl("failed", gene0)] = "failed"
  controlind = grepl(IMPC_control, gene0)
  gene = str_split(gene0,"_")
  gene1 = sapply(gene, function(x) x[1])
  gene2 = sapply(gene, function(x) x[2])
  # identical(gene1,gene2)
  # sort(unique(gene1[gene2=="+"]))
  # for (g in unique(gene1[gene2=="+"])) cat("\n(",sum(gene1==g & gene2=="+"),"/",sum(gene1==g),") ", g, sep="")
  # sort(unique(gene1[gene2=="Y"]))
  # for (g in unique(gene1[gene2=="Y"])) cat("\n(",sum(gene1==g & gene2=="Y"),"/",sum(gene1==g),") ", g, sep="")
  # gene = unique(gene1[gene2!="+" & gene2!="Y"]) #delete 118 genes because mostly 1/2 or 1/1 -/+
  # gene = gene1
  gene = str_split(gene1,"[(]") #take out "(b)" at the end of gene name, number of genes dont decrease
  gene = sapply(gene, function(x) x[1])
  gene = toupper(gene)
  # write.table(gene,file=gene_dir,sep="\t",row.names=F,col.names=F,quote=F)
  gene[controlind] = meta_file$gene[controlind]
  meta_file$gene = gene
  
} else if (grepl("ciphe", centre)) {
  meta_fileTemp = data.frame(lapply(read.csv(csv_dir)[,c("Date","SEX..M.F.","Sample.N.")], as.character))
  meta_fileTempDate = dmy(meta_fileTemp[,"Date"])
  ftFileNamesSplit = lapply(c(1:length(ftFileNames)), function(x) {return(strsplit(ftFileNames[x],"[_]")[[1]])} )
  date = ymd(sapply(c(1:length(ftFileNamesSplit)), function(x) {return(ftFileNamesSplit[[x]][1])} ))
  sample = as.numeric(sapply(c(1:length(ftFileNamesSplit)), function(x) {return(ftFileNamesSplit[[x]][3])} ))
  for (i in 1:length(ftFileNames)) {
    row = which(meta_fileTempDate==date[i] & meta_fileTemp[,"Sample.N."]==sample[i])
    meta_fileList[[i]] = c(ftFileNames[i], ftGT[i], unlist(meta_fileTemp[row,]))
  }
  meta_file = rbind(meta_file, do.call(rbind, meta_fileList))
  meta_file = as.data.frame(meta_file, stringsAsFactors=F); rm(meta_fileTemp)
  colnames(meta_file) = c("id", "gene", "date", "gender", "sample")
  meta_file$date = dmy(meta_file$date)
  
} else if (grepl("TCP", centre)) {
  meta_fileTemp = data.frame(lapply(read.csv(csv_dir)[,c("expDate", "sex", "SampleID", "strainCode")], as.character))
  meta_fileTemp2 = data.frame(lapply(read.csv(csv_dirtcp)[,c("dates", "list.date.short")], as.character))
  meta_fileTemp[,"expDate"] = mdy(meta_fileTemp[,"expDate"])
  
  meta_fileTemp2[,"list.date.short"] = gsub(".labelled","-labelled", gsub(".fcs","",meta_fileTemp2[,"list.date.short"]))
  rows = match(ftFileNames, meta_fileTemp2[,"list.date.short"])
  
  ftFileNamesSplit = lapply(c(1:length(ftFileNames)), function(x) {return(strsplit(ftFileNames[x],"[_]")[[1]])} )
  strain = sapply(c(1:length(ftFileNamesSplit)), function(x) {return(ftFileNamesSplit[[x]][3])} )
  sample = sapply(c(1:length(ftFileNamesSplit)), function(x) {return(ftFileNamesSplit[[x]][4])} )
  date = ymd(meta_fileTemp2[rows,"dates"])
  gender = NULL
  for (i in 1:length(ftFileNames)) {
    row = which(meta_fileTemp[,"expDate"]==date[i] & meta_fileTemp[,"SampleID"]==sample[i] & meta_fileTemp[,"strainCode"]==strain[i])[1]
    gender = append(gender, meta_fileTemp[row,"sex"])
  }
  meta_file = as.data.frame(cbind(ftFileNames, ftGT, date, gender, sample, strain), stringsAsFactors=F)
  colnames(meta_file) = c("id", "gene", "date", "gender", "sample", "strain")
  meta_file$date = date
  
} else if (grepl("harwell", centre)) {
  meta_fileTemp = data.frame(lapply(read.csv(csv_dir)[,c("Filename", "Date","Gender","IMPCSpecimenCode")], as.character))
  meta_fileTemp[,"Filename"] = gsub("[.]", "-", gsub(".fcs","", gsub("[-%]",".", meta_fileTemp[,"Filename"])))
  rows = match(ftFileNames, meta_fileTemp[,"Filename"])
  meta_file = cbind(meta_fileTemp[rows,"Filename"], ftGT, meta_fileTemp[rows,c("Date","Gender","IMPCSpecimenCode")])
  meta_file = as.data.frame(meta_file, stringsAsFactors=F); rm(meta_fileTemp)
  colnames(meta_file) = c("id", "gene", "date", "gender", "specimen")
  meta_file$date = ymd(meta_file$date)
  
  # CIPHE Manula count import
  #   CSVfile = paste("Results", centre, "/IMPC_20142015_merged_Panel1.csv", sep="")
  #   matrixCSV = read.csv(CSVfile)
  #   matrixManualCounts = matrix(0, nrow=55, ncol=ncol(feat_file_cell_count))
  #   mcIndex = c(19:ncol(matrixCSV))[which(c(19:ncol(matrixCSV)) %% 2 != 0)]
  #   rownames(matrixManualCounts) = sapply(mcIndex, function(x) {return( gsub("neg", "-", gsub("pos", "+", strsplit(colnames(matrixCSV)[x], "\\.\\.")[[1]][1]) ) ) })
  #   
  #   missingMetaIndex = NULL
  #   for (i in 1:nrow(meta_file)) {
  #     row = which(matrixCSV[,1]==meta_file$sampleNum[i] & as.Date(matrixCSV[,3], "%d-%m-%Y")==meta_file$date[i])
  #     if (length(row)==0) {missingMetaIndex = c(missingMetaIndex, i); next()}
  #     meta_file[i,c(5:7)] = matrixCSV[row,c(2,10,9)]
  #     matrixManualCounts[,i] = as.vector( as.numeric(matrixCSV[row,mcIndex] ))
  #   }
  #   save(matrixManualCounts, file=paste("Results", centre, "/matrixManualCounts.Rdata", sep=""))
}
g = rep(1,nrow(meta_file)); g[grep("F",meta_file$gender)] = 0; g[which(is.na(meta_file$gender))] = NA
meta_file$gender = g # Female=0; Male=1


## trim matrix/meta_cell ------------------------------
rowIndex = apply(feat_file_cell_count, 1, function(x) any(x > countThres)) #delete rows of all 0 or too little count
colIndex1 = apply(feat_file_cell_count, 2, function(x) any(x > countThres)) #delete cols of all 0 or too little count
colIndex2 = meta_cell$phenolevel <= levelThres #delete cols of too high level
colIndex = colIndex1 & colIndex2

feat_file_cell_count <- feat_file_cell_count[rowIndex,colIndex]
feat_file_cell_prop <- feat_file_cell_prop[rowIndex,colIndex]
meta_cell <- meta_cell[colIndex,]
meta_file <- meta_file[rowIndex,]

#rename classes, so there is a control group
meta_file$gene[grepl(IMPC_control,meta_file$gene)] = "control"
meta_file$gene[grepl("wildtype",meta_file$gene, ignore.case=T)] = "control"
colnames(meta_file)[colnames(meta_file)=="gene"] = "class"
if ("gender" %in% colnames(meta_file)) { 
  meta_file$gender = as.character(meta_file$gender)
  meta_file$gender[meta_file$gender=="0"] = "control"
  meta_file$gender[meta_file$gender=="1"] = "male"
} 


## save ----------------------------------------------------
save(meta_cell, file=paste0(meta_cell_dir,".Rdata"))
if (writecsv) write.csv(meta_cell, file=paste0(meta_cell_dir,".csv"))

save(meta_file, file=paste0(meta_file_dir,".Rdata"))
if (writecsv) write.csv(meta_file, file=paste0(meta_file_dir,".csv"))

feat_file_cell_count = as.matrix(feat_file_cell_count)
save(feat_file_cell_count, file=paste0(feat_file_cell_count_dir,".Rdata"))
if (writecsv) write.csv(feat_file_cell_count, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)

feat_file_cell_prop = as.matrix(feat_file_cell_prop)
save(feat_file_cell_prop, file=paste0(feat_file_cell_prop_dir,".Rdata"))
if (writecsv) write.csv(feat_file_cell_prop, file=paste0(feat_file_cell_prop_dir,".csv"), row.names=T)
if (!is.null(excludedGT)) save(excludedGT, file=excludedGT_dir)


time_output(start)

