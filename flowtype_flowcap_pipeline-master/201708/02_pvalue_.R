# aya43@sfu.ca 20151228
# Calculate PValues

root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrix", sep="")
matrix_type = c("CountAdj", "Prop")
#adjust = "bonferroni"
test = "wilcox" #pvalue test
adjust = "BH" #pvalue adjustment
cellCountThres = 10 #insignificant if count under
pvalThres = .05 #delete phenotypes/rows without any significant changes from the pVal matrix
sampleCountThres = 3 #only compare if >=3 samples available

#Output
matrixPval_dir = paste(result_dir, "/matrixPval",sep="")
matrixPvalAdj_dir = paste(result_dir, "/matrixPvalAdj",sep="")
matrixPvalFULL_dir = paste(result_dir, "/matrixPvalFULL",sep="")
matrixPvalAdjFULL_dir = paste(result_dir, "/matrixPvalAdjFULL",sep="")


libr(stringr)
libr(foreach)
libr(doMC)
source("code/_funcAlice.R")



start = Sys.time()

no_cores = detectCores() - 2
registerDoMC(no_cores)
sampleMeta = get(load(sampleMeta_dir))

#split every patient into 2 groups of samples
half = sampleMeta$specimen
specimen_max = max(sampleMeta$specimen)
for (i in unique(sampleMeta$specimen)) {
  iind = which(sampleMeta$specimen==i)
  half[iind[1:floor(length(iind)/2)]] = specimen_max+i
}

for (ci in 1:2) {
  if (ci==1) { control = "normal"
  } else { control = "aml" }
  
  for (col in c("all","half")) {
    if (col=="all") {
      g = getGTindex(sampleMeta$aml, control, sampleCountThres, sampleMeta$specimen)
    } else if (col=="half") {
      g = getGTindex(sampleMeta$aml, control, sampleCountThres, half)
    }
    ftGT = g$attlist; ftWTIndex = g$controlIndex; ftKOGT = g$exp; ftKOIndex = g$expIndex; ftKOgoodGTIndex = g$goodexpIndex; rm(g)  #Index of unique KO genotypes that has 3+ samples available ftKOIndex(ftKOGTIndexIndex)
    
    if (col=="half") {
      newhalf = which(ftKOGT>specimen_max)
      ftKOGT[newhalf] = ftKOGT[newhalf]-specimen_max
    }
    
    for (mcp in matrix_type) {
      start1 = Sys.time()
      
      #load matrices
      m = get(load(paste0(matrix_dir[ci], mcp,".Rdata")))
      
      rowcombos = NULL
      for (i in length(ftKOIndex):1) {
        rowcombos[[i]][[1]] = ftWTIndex
        rowcombos[[i]][[2]] = ftKOIndex[[i]]
      }
      
      cat(paste("getting pValues of ", ncol(m), " possible cell populations for ", length(rowcombos), " experiments ", sep="")) #3iTCell specific
      
      loop.ind = 1:ncol(m)
      matrixPval = foreach(k = loop.ind, .combine = 'cbind', .multicombine = T) %dopar% { #for each phenotype
        #for (k in 1:ncol(m)){ cat(paste(" ", j, sep="")) {
        pvalcol = rep(1,length(rowcombos))
        for (j in 1:length(rowcombos)) { #for each KO gene
          compare1 = as.numeric(m[ rowcombos[[j]][[1]],k ])
          compare2 = as.numeric(m[ rowcombos[[j]][[2]],k ])
          if (!(median(compare1)<cellCountThres & median(compare2)<cellCountThres) | !grepl("Count",mcp)) { #TRIM: if both WT and KO medians (to avoid outler influence) is < cell count threshold then change is not significant (set as 1)
            if (test=="wilcox") {
              pvalcol[j] = wilcox.test(compare1, compare2)$p.value
            } else if (test=="ttest") {
              try({ pvalcol[j] = t.test(compare1, compare2)$p.value })
            }
          }
        }
        return(pvalcol)
      }
      matrixPval[is.nan(matrixPval)] = 1
      colnames(matrixPval) = colnames(m)
      
      matrixPvalAdj = NULL
      for (i in 1:nrow(matrixPval)) { matrixPvalAdj = rbind(matrixPvalAdj, p.adjust(matrixPval[i,], method=adjust)) }
      
      #trim pvalues
      matrixPvalFULL = matrixPval
      matrixPvalAdjFULL = matrixPvalAdj
      
      rownames(matrixPvalFULL) = rownames(matrixPvalAdjFULL) = ftKOGT
      
      matrixPval = matrixTrim(matrixPvalFULL, thres=pvalThres, trimGreater=T)
      matrixPvalAdj = matrixTrim(matrixPvalAdjFULL, thres=pvalThres, trimGreater=T)
      
      save(matrixPval, file=paste0(matrixPval_dir, mcp, "control-", control, "_", col, ".Rdata"))
      save(matrixPvalFULL, file=paste0(matrixPvalFULL_dir, mcp, "control-", control, "_", col, ".Rdata"))
      save(matrixPvalAdj, file=paste0(matrixPvalAdj_dir, mcp, "control-", control, "_", col, ".Rdata"))
      save(matrixPvalAdjFULL, file=paste0(matrixPvalAdjFULL_dir, mcp, "control-", control, "_", col, ".Rdata"))
      
      cat("\nTime taken to calculate p values & barcode matrices is: ",TimeOutput(start1), "\n", sep="") #3iTcell ~40min
    }
    
  }
}
TimeOutput(start)

