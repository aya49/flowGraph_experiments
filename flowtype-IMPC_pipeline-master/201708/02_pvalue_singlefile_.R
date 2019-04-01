# aya43@sfu.ca 20151228
# Calculate PValues

#root = "~/projects/IMPC"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN","Sanger_MLN","CIPHE","TCP","H")

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrix", sep="")
matrix_type = c("CountAdj")#, "Prop")
#adjust = "bonferroni"
# test = "wilcox" #pvalue test
# adjust = "BH" #pvalue adjustment
cellCountThres = 10 #insignificant if count under
pvalThres = .025 #delete phenotypes/rows without any significant changes from the pVal matrix
# sampleCountThres = 3 #only compare if >=3 samples available

#Output
# matrixPval_dir = paste(result_dir, "/matrixPval1",sep="")
# matrixPvalAdj_dir = paste(result_dir, "/matrixPvalAdj",sep="")
matrixPvalFULL_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixPvalFULL1",sep="")
# matrixPvalAdjFULL_dir = paste(result_dir, "/matrixPvalAdjFULL",sep="")
matrixLogRatio_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixLogRatio1",sep="")

libr(stringr)
libr(foreach)
libr(doMC)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")


start = Sys.time()

no_cores = 6#detectCores() - 2
registerDoMC(no_cores)

#split every patient into 2 groups of samples
# half = sampleMeta$specimen
# specimen_max = max(sampleMeta$specimen)
# for (i in unique(sampleMeta$specimen)) {
#   iind = which(sampleMeta$specimen==i)
#   half[iind[1:floor(length(iind)/2)]] = specimen_max+i
# }
# 

for (ci in length(paste0(panelL,centreL)):1) {
  centre = paste0(panelL,centreL)[ci]

  control = c("+_+","+_Y"); if (!grepl("Sanger",centre)) control = "WildType"
  sampleMeta = get(load(sampleMeta_dir[ci]))
  colnames(sampleMeta)[1] = "filename"
  save(sampleMeta,file=sampleMeta_dir[ci])
  
  g = getGTindex(sampleMeta$gene, control, cellCountThres, sampleMeta$filename)
  ftGT = g$attlist; ftWTIndex = g$controlIndex; ftKOGT = g$exp; ftKOIndex = g$expIndex; ftKOgoodGTIndex = g$goodexpIndex; rm(g)  #Index of unique KO genotypes that has 3+ samples available ftKOIndex(ftKOGTIndexIndex)
  
  rowcombos = NULL
  for (i in length(ftKOIndex):1) {
    rowcombos[[i]][[1]] = ftWTIndex
    rowcombos[[i]][[2]] = ftKOIndex[[i]]
  }
  
  for (mcp in matrix_type) {
    start1 = Sys.time()
    
    #load matrices
    m = get(load(paste0(matrix_dir[ci], mcp,".Rdata")))
    
    
    cat(paste(centre, " getting pValues of ", ncol(m), " possible cell populations for ", length(rowcombos), " experiments ", sep="")) #3iTCell specific
    
    loop.ind = 1:ncol(m)
    result = foreach(k = loop.ind) %dopar% { #for each phenotype
      #for (k in 1:ncol(m)){ cat(paste(" ", j, sep="")) {
      pvalcol = ratiocol = rep(1,length(rowcombos))
      for (j in 1:length(rowcombos)) { #for each KO gene
        compare1 = as.numeric(m[ rowcombos[[j]][[1]],k ])
        compare2 = as.numeric(m[ rowcombos[[j]][[2]],k ])
        if (!(median(compare1)<cellCountThres & median(compare2)<cellCountThres) | !grepl("Count",mcp)) { #TRIM: if both WT and KO medians (to avoid outler influence) is < cell count threshold then change is not significant (set as 1)
          # if (test=="wilcox") {
          pvalcol[j] = t.test.single(compare1, compare2)
          # } else if (test=="ttest") {
          #   try({ pvalcol[j] = t.test(compare1, compare2)$p.value })
          # }
          ratiocol[j] = log( compare2/exp(mean(log(compare1))) ) #over geo mean, exp it to get non-log ratio
        }
      }
      return(list(p=pvalcol,r=ratiocol))
    }
    
    matrixPvalFULL = matrixLogRatio = matrix(0,nrow=length(rowcombos),ncol=length(result))
    for (k in 1:length(result)) {
      matrixPvalFULL[,k] = result[[k]]$p
      matrixLogRatio[,k] = result[[k]]$r
    }
    rm(result)
    matrixPvalFULL[is.nan(matrixPvalFULL)] = 1
    colnames(matrixPvalFULL) = colnames(matrixLogRatio) = colnames(m)
    rownames(matrixPvalFULL) = rownames(matrixLogRatio) = ftKOGT
    
    # matrixPvalAdj = NULL
    # for (i in 1:nrow(matrixPval)) { matrixPvalAdj = rbind(matrixPvalAdj, p.adjust(matrixPval[i,], method=adjust)) }
    
    # #trim pvalues
    # matrixPvalFULL = matrixPval
    
    # matrixPval = matrixTrim(matrixPvalFULL, thres=pvalThres, trimGreater=T)
    # matrixPvalAdj = matrixTrim(matrixPvalAdjFULL, thres=pvalThres, trimGreater=T)
    
    try({ save(matrixPvalFULL, file=paste0(matrixPvalFULL_dir[ci], ".Rdata")) })
    write.csv(matrixPvalFULL, file=paste0(matrixPvalFULL_dir[ci], ".csv"))
    # save(matrixPvalFULL, file=paste0(matrixPvalFULL_dir, mcp, "control-", control, "_", col, ".Rdata"))
    # save(matrixPvalAdj, file=paste0(matrixPvalAdj_dir, mcp, "control-", control, "_", col, ".Rdata"))
    # save(matrixPvalAdjFULL, file=paste0(matrixPvalAdjFULL_dir, mcp, "control-", control, "_", col, ".Rdata"))
    try({ save(matrixLogRatio, file=paste0(matrixLogRatio_dir[ci], ".Rdata")) })
    write.csv(matrixLogRatio, file=paste0(matrixLogRatio_dir[ci], ".csv"))
    
    TimeOutput(start1) #3iTcell ~40min
  }
}
TimeOutput(start)

