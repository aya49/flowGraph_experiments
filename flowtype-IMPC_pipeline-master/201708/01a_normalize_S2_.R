## Input: original count matrix --> Output: normalized count matrix
#aya43@sfu.ca 20151228

#Directory
root = "~/projects/IMPC/SangerP2"
result_dir = "Results"; suppressWarnings(dir.create (result_dir))
setwd(root)

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

#Input
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
matrixCount_dir = paste(result_dir, "/matrixCount.Rdata", sep="")

#Output
matrixCountAdj_dir = paste(result_dir, "/matrixCountAdj",sep="")
matrixCountAdjlog_dir = paste(result_dir, "/matrixCountAdjlog",sep="")
norm_dir = paste(result_dir, "/cellCountNormFactor",sep=""); suppressWarnings(dir.create(norm_dir))
normFactor_dir = paste(result_dir, "/normFactor.Rdata", sep="")
normFactorDiffLog_dir = paste(result_dir, "/normFactorDiffDensLogged.Rdata", sep="")
norm_fdiffplot_dir = paste0(norm_dir,"/cellCountNormFactor_f_diff_from_peak_abs__tube")

#Libraries/Functions
libr(stringr)
libr(pracma)
libr(foreach)
libr(doMC)
libr(flowDensity)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)











#Options for script
cutoff = c(Inf) #c(.6) #if TMM-peak>cutoff, then apply peak instead of TMM; run this script and look at norm_fdiffplot plot to determine this number
layer = 0 #calculate TMM using only phenotypes in this layer; set to 0 if do for all layers
cellCountThres = 200 #don't use phenotypes with cell count lower than cellCountThres
splitby = NULL #column to split by -- FlowCAP-II normalizes only between files of same tube/panel; set to NULL if don't split
target_col = "gene" #column with control/experiment
control = "+_+|+_Y" #control value in target_col column



#Prepare data
sampleMeta0 = get(load(sampleMeta_dir))
# phenoMeta = get(load(phenoMeta_dir))
matrixCount0 = get(load(matrixCount_dir))










start = Sys.time()

#calculate TMM
fdiff0 = f0 = rep(0,nrow(sampleMeta0))

for (ti in ifelse(is.null(splitby), NULL, unique(sampleMeta0[,splitby]))) {
  start1 = Sys.time()
  
  #split by ti
  tubeind = ifelse(is.null(ti), rep(T,nrow(sampleMeta0)), sampleMeta0[,splitby]==ti)
  matrixCount = matrixCount0[tubeind,]
  sampleMeta = sampleMeta0[tubeind,]
  
  #prepare matrixCounts
  x = x0 = as.matrix(matrixCount)[,-1] #take out total cell count
  if (layer>0) x = as.matrix(x0[,colnames(x0)%in%phenoMeta$phenotype[phenoMeta$phenolevel==layer] & sapply(1:ncol(x0), function(y) any(x0[,y]>cellCountThres))])
  lib.size = matrixCount[,1]
  refColumn = which.min(abs( lib.size-median(lib.size[grep(control,sampleMeta[,target_col])]) )) #reference column: median total count out of all control files
  
  #prepare plot paths/titles
  loop.ind = 1:nrow(x)
  pngnames = sapply(loop.ind, function(i) paste0(norm_dir, "/cellCountNormFactor_",str_pad(i, 4, pad = "0"), "_",sampleMeta[i,target_col], "_", sampleMeta$fileName[i], ".png")) 
  mains = sapply(loop.ind, function(i) paste0("mean count vs. ln fold change:\n", sampleMeta[i,target_col]," over refColumn ", sampleMeta[refColumn,target_col], "___layer-",layer))
  
  ## calculate absolute count TMM, mostly taken from TMM
  fresult = tmm(x,x0,lib.size,refColumn,cutoff=Inf,plotimg=T,pngnames=pngnames,mains=mains,no_cores=no_cores,samplesOnCol=F)
  f0[tubeind] = fresult$f
  fdiff0[tubeind] = fresult$fdiff
  
  #plot difference between TMM and peak for all files
  pngname = ifelse(is.null(ti), paste0(norm_fdiffplot_dir,".png"), paste0(norm_fdiffplot_dir,"-",ti,".png"))
  png(file=pngname , width=700, height=700)
  plot(sort(abs(fdiff0[tubeind])), cex=.4, ylim=c(0,3))
  lines(sort(abs(fdiff0[tubeind])), col="blue")
  graphics.off()
  
  TimeOutput(start1)
}

matrixCountAdj = sapply(c(1:nrow(matrixCount0)), function(x) {matrixCount0[x,]*f0[x]})
matrixCountAdj = t(matrixCountAdj)
colnames(matrixCountAdj) = colnames(matrixCount0)
rownames(matrixCountAdj) = rownames(matrixCount0)
#phenotype on cols
matrixCountAdjlog = log(matrixCountAdj)
matrixCountAdjlog[which(matrixCountAdjlog<0)] = 0

#save
save(f0, file=normFactor_dir)
save(fdiff0, file=normFactorDiffLog_dir)
save(matrixCountAdj, file=paste0(matrixCountAdj_dir,".Rdata"))
write.csv(matrixCountAdj, file=paste0(matrixCountAdj_dir,".csv"), row.names=F)
save(matrixCountAdjlog, file=paste0(matrixCountAdjlog_dir,".Rdata"))
write.csv(matrixCountAdjlog, file=paste0(matrixCountAdjlog_dir,".csv"), row.names=F)

TimeOutput(start)






