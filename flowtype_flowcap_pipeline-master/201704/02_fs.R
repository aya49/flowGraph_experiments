# aya43@sfu.ca 20161220
# Uses different distance measures to calculate distance & plot samples (for pvalue matrices)

root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
matrixPval_dir = list.files(result_dir, pattern="Pval", full.names=T)
matrixPval_dir = matrixPval_dir[-grep("FULL",matrixPval_dir)]
matrixPval_names = gsub("matrix","",fileNames(matrixPval_dir, ext="Rdata")) 


#Output
dist_dir = paste(result_dir, "/dist", sep=""); suppressWarnings(dir.create (dist_dir))

libr(stringr)
libr(colorspace)
libr(vegan) # libr(proxy)
libr(prabclus)
libr(fastcluster)
libr(dendextend)
libr(circlize)
libr(Rtsne)
libr(foreach)
libr(doMC)
libr(lubridate) #if there are date variables
source("code/_funcAlice.R")



dodist = F
doHC = F
doTsne = T
pvalThres = .05

dis = c("canberra","bray","kulczynski", "binomial","cao", "jaccardInd") # "jaccard","mahalanobis", "gower","altGower","morisita","horn", "euclidean", ""manhattan", "maximum", "binary", "minkowski", "raup", "mountford" (all at same level),)

start = Sys.time()

no_cores = detectCores() - 2
registerDoMC(no_cores)

load(sampleMeta_dir)
load(phenoMeta_dir)

k0 = c(1,max(phenoMeta[,3])+1) # how many markers to consider i.e. k=max(phenolevel) only






























# remove redundant features -----------------------------------

# ensure the results are repeatable
set.seed(7)
# load the libr
libr(mlbench)
libr(caret)
# load the data
data(PimaIndiansDiabetes)
# calculate correlation matrix
correlationMatrix <- cor(PimaIndiansDiabetes[,1:8])
# summarize the correlation matrix
print(correlationMatrix)
# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.5)
# print indexes of highly correlated attributes
print(highlyCorrelated)









