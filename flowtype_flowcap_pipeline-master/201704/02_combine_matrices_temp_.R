# aya43@sfu.ca 20170213
# combine all centre matrices

#root = "~/projects/IMPC"
result_dir = "result"
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN","Sanger_MLN","CIPHE","TCP","H")

options(stringsAsFactors=FALSE)
#options(device="cairo")
options(na.rm=T)

#Input
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
phenoMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoMeta.Rdata", sep="")
markers_dir = paste(result_dir, "/", panelL, "/", centreL, "/markers.Rdata", sep="")
matrixAdj_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixCountAdj.Rdata", sep="")

#Output
sampleMeta0_dir = paste(result_dir, "/", panelL, "/sampleMeta.Rdata", sep="")
sampleMeta0_dircsv = paste(result_dir, "/", panelL, "/sampleMeta.csv", sep="")
phenoMeta0_dir = paste(result_dir, "/", panelL, "/phenoMeta.Rdata", sep="")
markers0_dir = paste(result_dir, "/", panelL, "/markers.Rdata", sep="")
matrixAdj0_dir = paste(result_dir, "/", panelL, "/matrixCountAdj.Rdata", sep="")
matrixAdj0_dircsv = paste(result_dir, "/", panelL, "/matrixCountAdj.csv", sep="")

libr(stringr)
libr(foreach)
libr(doMC)
libr(lubridate) #if there are date variables
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")

start = Sys.time()

no_cores = detectCores()-4
registerDoMC(no_cores)

sampleMeta0 = NULL
for (ci in 5:1) { #5th centre has least markers
  sampleMeta = get(load(sampleMeta_dir[ci]))
  sampleMeta_temp = cbind(rep(centreL[ci],nrow(sampleMeta)), sampleMeta[,c("fileName","gene","date","gender")])
  colnames(sampleMeta_temp)[1] = "centre"
  sampleMeta0 = rbind(sampleMeta0, sampleMeta_temp)
  phenoMeta = get(load(phenoMeta_dir[ci]))
  m = get(load(matrixAdj_dir[ci]))
  if (ci==5) {
    m_all = m
    phenoMeta0 = phenoMeta
    markers0 = get(load(markers_dir[ci]))
    colnames(sampleMeta0)[1] = "centre"
  } else {
    cols = match(phenoMeta0$phenotype,phenoMeta$phenotype)
    m_all = rbind(m_all, m[,cols])
  }
  
}
save(markers0,file=markers0_dir)
save(phenoMeta0,file=phenoMeta0_dir)
save(sampleMeta0,file=sampleMeta0_dir)
write.csv(sampleMeta0,file=sampleMeta0_dircsv)
save(m_all,file=matrixAdj0_dir)
write.table(m_all,file=matrixAdj0_dircsv, sep=",", col.names=F,row.names=F)






