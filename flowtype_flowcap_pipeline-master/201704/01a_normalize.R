# Normalizes cell count matrix
# aya43@sfu.ca 20151228

root <- "~/projects/flowCAP-II"
result_dir <- "result"; suppressWarnings(dir.create (result_dir))
setwd(root)


options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir <- paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir <- paste(result_dir, "/sampleMeta.Rdata", sep="")
matrixCount_dir <- paste(result_dir, "/matrixCount.Rdata", sep="")
cutoff = c(Inf)#c(.6) #run this script and look at plot to determine this number
layer = 4
cellCountThres = 200

#Output
matrixCountAdj_dir <- paste(result_dir, "/matrixCountAdj",sep="")
matrixCountAdjlog_dir <- paste(result_dir, "/matrixCountAdjlog",sep="")
norm_dir <- paste(result_dir, "/cellCountNormFactor",sep="")
normFactor_dir <- paste(result_dir, "/normFactor.Rdata", sep="")
normFactorDiffLog_dir <- paste(result_dir, "/normFactorDiffDensLogged.Rdata", sep="")

libr(stringr)
libr(pracma)
libr(foreach)
libr(doMC)
libr(flowDensity)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")



start <- Sys.time()

no_cores <- detectCores() - 1
registerDoMC(no_cores)

control <- "normal"

suppressWarnings(dir.create(norm_dir))
sampleMeta0 <- get(load(sampleMeta_dir))
phenoMeta <- get(load(phenoMeta_dir))
matrixCount0 <- get(load(matrixCount_dir))
#phenotype on rows
# if (ncol(matrixCount0)!=nrow(sampleMeta0)) { matrixCount0 <- t(matrixCount0) }

fdiff0 = f0 = rep(0,nrow(sampleMeta0))

for (ti in unique(sampleMeta0$tube)) {
  start1 = Sys.time()
  
  tubeind = which(sampleMeta0$tube==ti)
  matrixCount = matrixCount0[tubeind,]
  sampleMeta = sampleMeta0[tubeind,]
  
  ## Taken from TMM ------------------------------------------------------
  x <- x0 <- as.matrix(matrixCount)[,-1] #take out total cell count
  if (layer>0) x = as.matrix(x0[,colnames(x0)%in%phenoMeta$phenotype[phenoMeta$phenolevel==layer] & sapply(1:ncol(x0), function(y) any(x0[,y]>cellCountThres))])
  lib.size <- matrixCount[,1]
  refColumn <- which.min(abs( lib.size-median(lib.size[grep(control,sampleMeta$aml)]) )) #Get reference column: median normal patient
  #p=0.75
  #  f75 <- apply(t(t(x)/lib.size),2,function(data) quantile(data,p=p))
  #  refColumn <- which.min(abs(f75-mean(f75)))
  
  loop.ind <- 1:nrow(x)
  pngnames <- sapply(loop.ind,function(i) paste(norm_dir, "/cellCountNormFactor_",str_pad(i, 4, pad = "0"), "_",sampleMeta$aml[i], "_", sampleMeta$fileName[i], ".png", sep = "" )) 
  mains = sapply(loop.ind,function(i) paste0("mean count vs. ln fold change:\n", sampleMeta$aml[i]," over refColumn ", sampleMeta$aml[refColumn], "___layer-",layer,sep=""))
  
  fresult = tmm(x,x0,lib.size,refColumn,cutoff=Inf,plotimg=T,pngnames=pngnames,mains=mains,no_cores=no_cores,samplesOnCol=F)
  
  f0[tubeind] = f = fresult$f
  fdiff0[tubeind] = fdiff = fresult$fdiff
  
  pngname <- paste(norm_dir, "/cellCountNormFactor_f_diff_from_peak_abs__tube-",ti,".png", sep = "" )
  png (file=pngname , width=700, height=700)
  plot(sort(abs(fdiff)), cex=.4, ylim=c(0,3))
  lines(sort(abs(fdiff)), col="blue")
  dev.off()
  
  TimeOutput(start1)
}

matrixCountAdj <- sapply(c(1:nrow(matrixCount0)), function(x) {matrixCount0[x,]*f0[x]})
matrixCountAdj <- t(matrixCountAdj)
colnames(matrixCountAdj) <- colnames(matrixCount0)
rownames(matrixCountAdj) <- rownames(matrixCount0)
#phenotype on cols
matrixCountAdjlog = log(matrixCountAdj)
matrixCountAdjlog[which(matrixCountAdjlog<0)] = 0

save(f0, file=normFactor_dir)
save(fdiff0, file=normFactorDiffLog_dir)
save(matrixCountAdj, file=paste0(matrixCountAdj_dir,".Rdata"))
write.csv(matrixCountAdj, file=paste0(matrixCountAdj_dir,".csv"), row.names=F)
save(matrixCountAdjlog, file=paste0(matrixCountAdjlog_dir,".Rdata"))
write.csv(matrixCountAdjlog, file=paste0(matrixCountAdjlog_dir,".Rdata"), row.names=F)


# #max change in slope
# cutoff = sort(abs(fdiff))[which.max(diff(sort(abs(fdiff))))]

# #Find inflection point of fdiff curve to define cut-off for when to switch to highest density method
# x = seq(1,length(fdiff))
# y = sort(abs(fdiff))
# #The exact inflection point is ip=5.0
# #Because of the total symmetry we expect inflection point near the middle of x-range:
# A<-findiplist(x,y,0);A;#Our expectation came true
# #Let's make some ESE and EDE iterations and plot them:
# a<-findipiterplot(x,y,0,TRUE,TRUE,TRUE);
# a$first;#Show first solution
# a$BESE;#Show ESE iterations
# a$CRESE;#Show cross ESE iterations from EDE iterations
# a$esm;#Show all ESE iterations
# a$BEDE;#Show EDE iterations
# a$CREDE;#Show cross EDE iterations from ESE iterations
# a$edm;#Show all EDE iterations
# a$aesmout;#Statistics and 95%c c.i. for ESE
# a$besmout;#Statistics and 95%c c.i. for cross EsE
# a$aedmout;#Statistics and 95%c c.i. for EDE
# a$bedmout;#Statistics and 95%c c.i. for cross EDE
# a$esmout;#Statistics and 95%c c.i. for all ESE results
# a$edmout;#Statistics and 95%c c.i. for all EDE results
# a$ipall;#Statistics and 95%c c.i. for results
# #Close the 2 previously opened devises








TimeOutput(start) #21min parallel, 5 centres






