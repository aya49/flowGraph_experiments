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
matrix_type = c("CountAdj", "Prop")
#adjust = "bonferroni"
test = "wilcox" #pvalue test
adjust = "BH" #pvalue adjustment
cellCountThres = 10 #insignificant if count under
pvalThres = .05 #delete phenotypes/rows without any significant changes from the pVal matrix
sampleCountThres = 3 #only compare if >=3 samples available

#Output
matrixPval_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixPval",sep="")
matrixPvalAdj_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixPvalAdj",sep="")
matrixPvalFULL_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixPvalFULL",sep="")
matrixPvalAdjFULL_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixPvalAdjFULL",sep="")


libr(stringr)
libr(foreach)
libr(doMC)
source("~/projects/IMPC/code/_funcAlice.R")



start = Sys.time()

no_cores = detectCores() - 2
registerDoMC(no_cores)

for (ci in 1:length(paste0(panelL,centreL))) {
  centre = paste0(panelL,centreL)[ci]
  
  for (mcp in matrix_type) {
    start1 = Sys.time()
    
    #load matrices
      m = get(load(paste0(matrix_dir[ci], mcp,".Rdata")))
    sampleMeta = get(load(sampleMeta_dir[ci]))
    
    #wildtypes
    ftWTGT = c("+_+","+_Y"); if (!grepl("Sanger",centre)) ftWTGT = "WildType"
    g = getGTindex(sampleMeta$gene, ftWTGT, sampleCountThres)
    ftGT = g$attlist; ftWTIndex = g$controlIndex; ftKOGT = g$exp; ftKOIndex = g$expIndex; ftKOgoodGTIndex = g$goodexpIndex; rm(g)  #Index of unique KO genotypes that has 3+ samples available ftKOIndex(ftKOGTIndexIndex)
    
    rowcombos = NULL
    for (i in length(ftKOIndex):1) {
      rowcombos[[i]][[1]] = ftWTIndex
      rowcombos[[i]][[2]] = ftKOIndex[[i]] 
    }
    
    cat(paste("getting pValues of ", ncol(m), " possible cell populations for ", length(rowcombos), " genotypes ", sep="")) #3iTCell specific
    
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
          } 
          else if (test=="ttest") {
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
    matrixPval = matrixTrim(matrixPvalFULL, thres=pvalThres, trimGreater=T)
    matrixPvalAdj = matrixTrim(matrixPvalAdjFULL, thres=pvalThres, trimGreater=T)
    
    rownames(matrixPvalFULL) = rownames(matrixPvalAdjFULL) = rownames(matrixPval) = rownames(matrixPvalAdj) = rownames(matrixPval) = rownames(matrixPvalAdj) = ftKOGT
    
    save(matrixPval, file=paste0(matrixPval_dir[ci], "_", mcp,".Rdata"))
    save(matrixPvalFULL, file=paste0(matrixPvalFULL_dir[ci], "_", mcp,".Rdata"))
    save(matrixPvalAdj, file=paste0(matrixPvalAdj_dir[ci], "_", mcp,".Rdata"))
    save(matrixPvalAdjFULL, file=paste0(matrixPvalAdjFULL_dir[ci], "_", mcp,".Rdata"))
    
    cat("\nTime taken to calculate p values & barcode matrices is: ",TimeOutput(start1), "\n", sep="") #3iTcell ~40min
  }
}

cat("\nTime taken to calculate p values & barcode matrices is: ",TimeOutput(start), "\n", sep="") #3iTcell ~40min
