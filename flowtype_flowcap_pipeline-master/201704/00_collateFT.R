# aya43@sfu.ca 20161220
# Reads in flowCAP-II AML flowtype files and collates/save cell counts into matrix

root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
data_dir = "/mnt/f/Brinkman group/current/Alice/flowCAP-II/data"
ft_dir = paste0(data_dir,"/FT")
csv_dir = paste0(data_dir,"/AML.csv")

#Output
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
phenoMeta_dircsv = paste(result_dir, "/phenoMeta.csv", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
sampleMeta_dircsv = paste(result_dir, "/sampleMeta.csv", sep="")
markers_dir = paste(result_dir, "/markers.Rdata", sep="")
matrixCount_dir = paste(result_dir, "/matrixCount.Rdata", sep="")
matrixCount_dircsv = paste(result_dir, "/matrixCount.csv", sep="")
matrixProp_dir = paste(result_dir, "/matrixProp.Rdata", sep="")
matrixProp_dircsv = paste(result_dir, "/matrixProp.csv", sep="")

libr(flowCore)
libr(flowType)
libr(foreach)
libr(doMC)
source("~/projects/IMPC/code/_funcAlice.R")

no_cores = 6#detectCores()-3
registerDoMC(no_cores)


start = Sys.time()

ftFile_dir = sort( dir(ft_dir, pattern<-".Rda", all.files<-TRUE, full.names<-TRUE, recursive<-TRUE) )
ftGT = folderNames(ftFile_dir)
ftFileNames = fileNames(ftFile_dir, "Rda")

sampleMetatemp = data.frame(read.csv(csv_dir))
colnames(sampleMetatemp) = c("fileName", "tube", "specimen", "aml")

ft = get(load(ftFile_dir[1]))
markers = ft@MarkerNames
save(markers, file=markers_dir)

phenotype = rownames(ft@MFIs)
phenocode = unlist(lapply(phenotype, function(x){return( encodePhenotype(x, markers) )}))
phenolevel = unlist(lapply(phenocode, function(x){return(length(markers) - charOccurences("0", x)) } ))
phenoMeta = data.frame(phenotype, phenocode, phenolevel, stringsAsFactors=F)

matrixCount = NULL
sampleMeta = NULL

result = foreach (i = 1:length(ftFile_dir), .combine="rbind") %dopar% {
  ft = get(load(ftFile_dir[i]))
  print(ft@MarkerNames)
  npheno = length(ft@CellFreqs)
  #matrixCount = rbind(matrixCount, ft@CellFreqs)
  cat("\n", i, "Loaded file: ", ftGT[i], "/", ftFileNames[i], "; # of phenotypes x samples: ", npheno, " x ", length(ft), "", sep <-"")
  
  fn = ftFileNames[i]
  fn = strsplit(fn, "\\.")[[1]][1]
  fn = strsplit(fn, "S")[[1]]
  tube = substr(fn[1],2,nchar(fn[1]))
  specimen = gsub("FT","",fn[2])
  #sampleMeta = rbind(sampleMeta, sampleMetatemp[which(as.numeric(sampleMetatemp$specimen)==as.numeric(specimen) & as.numeric(sampleMetatemp$tube)==as.numeric(tube)),])
  return(c(sampleMetatemp[which(as.numeric(sampleMetatemp$specimen)==as.numeric(specimen) & as.numeric(sampleMetatemp$tube)==as.numeric(tube)),], ft@CellFreqs))
  
  rm(ft)
}

matrixCount = as.matrix(apply(result[,(ncol(sampleMetatemp)+1):ncol(result)], 2, as.numeric))
sm = result[,1:(ncol(sampleMetatemp))]

#delete phenotypes with no cells
all0Phen = which(apply(matrixCount[,-1], 1, function(x) all(x==0))==T)
if (length(all0Phen)>0) {
  matrixCount = matrixCount[,-all0Phen]
  phenoMeta[-all0Phen,]
}

#order matrix phenotypes
phenoorder = order(phenoMeta$phenocode)
phenoMeta = phenoMeta[phenoorder,]
matrixCount = matrixCount[,phenoorder]

rownames(matrixCount) = ftFileNames
sampleMeta = as.data.frame(ftFileNames)
sampleMeta$tube = unlist(sm[,2])
sampleMeta$specimen = unlist(sm[,3])
sampleMeta$aml = unlist(sm[,4])
colnames(sampleMeta) = colnames(sampleMetatemp)

rownames(matrixCount) = sampleMeta[,1] = ftFileNames
colnames(matrixCount) = phenoMeta$phenotype
matrixProp = matrixCount/matrixCount[,1]
dimnames(matrixProp) = dimnames(matrixCount)

save(phenoMeta, file=phenoMeta_dir)
write.csv(phenoMeta, file=phenoMeta_dircsv, row.names=F)
save(sampleMeta, file=sampleMeta_dir)
write.csv(sampleMeta, file=sampleMeta_dircsv, row.names=F)
save(matrixCount, file=matrixCount_dir)
write.csv(matrixCount, file=matrixCount_dircsv, row.names=F)
save(matrixProp, file=matrixProp_dir)
write.csv(matrixProp, file=matrixProp_dircsv, row.names=F)

TimeOutput(start)






