# aya43@sfu.ca 20161220
# Uses different distance measures to calculate distance & plot samples

root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
matrixCount_dir = paste(result_dir, "/matrixCount.Rdata", sep="")
matrixCountAdj_dir = paste(result_dir, "/matrixCountAdj.Rdata", sep="")
matrixProp_dir = paste(result_dir, "/matrixProp.Rdata", sep="")

#Output
py_dir = paste(result_dir, "_py", sep=""); suppressWarnings(dir.create (py_dir))
matrixCount_dirpy = paste(py_dir, "/matrixCount.csv", sep="")
matrixProp_dirpy = paste(py_dir, "/matrixProp.csv", sep="")
matrixCountAdj_dirpy = paste(py_dir, "/matrixCountAdj.csv", sep="")
matrixCount_cortrim_dirpy = paste(py_dir, "/matrixCount_cortrim.csv", sep="")
matrixProp_cortrim_dirpy = paste(py_dir, "/matrixProp_cortrim.csv", sep="")
matrixCountAdj_cortrim_dirpy = paste(py_dir, "/matrixCountAdj_cortrim.csv", sep="")
aml_dirpy = paste(py_dir, "/aml.csv", sep="")


libr(stringr)
source("code/_funcAlice.R")



dodist = F
doHC = F
doTsne = T


start = Sys.time()

load(sampleMeta_dir)
write.csv(sampleMeta$aml, file=aml_dirpy, row.names=F)


for (mcp in 1:3) { # Load & fix cell count/countAdj/proportion matrix; 1:3 are regular matrices, 4... are pvalues
  if (mcp==1) {
    m = get(load(matrixCountAdj_dir))
    write.csv(m, file=matrixCountAdj_dirpy, row.names=F)
    
    cormatrix = cor(m)
    #cormatrix_n = cor(m[which(sampleMeta$aml=="normal"),])
    #cormatrix_a = cor(m[which(sampleMeta$aml=="aml"),])
    cor = findCorrelation(cormatrix)
    #cor_n = findCorrelation(cormatrix_n)
    #cor_a = findCorrelation(cormatrix_a)
    
    write.csv(m[,-cor], file=matrixCountAdj_cortrim_dirpy, row.names=F)
  } else if (mcp==2) {
    m = get(load(matrixCount_dir))
    write.csv(m, file=matrixCount_dirpy, row.names=F)
  } else if (mcp==3) {
    m = get(load(matrixProp_dir))
    write.csv(m, file=matrixProp_dirpy, row.names=F)
  }
}


TimeOutput(start)





