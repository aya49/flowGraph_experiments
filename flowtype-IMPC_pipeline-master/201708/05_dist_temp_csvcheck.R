# aya43@sfu.ca 20170116
# Evaluate all distance metrics

#root = "~/projects/IMPC"
result_dir = "result"
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#"Sanger_MLN","CIPHE","TCP","H")

options(na.rm=T)

#Input
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
dist_dir = paste(result_dir, "/", panelL, "/", centreL, "/dist", sep="")

libr(stringr)
libr(foreach)
libr(doMC)
source("~/projects/IMPC/code/_funcAlice.R")


start = Sys.time()

no_cores = 4#detectCores()-1
registerDoMC(no_cores)

for (ci in 1:length(paste0(panelL,centreL))) {
  start2 = Sys.time()
  cat("\n", paste0(panelL," ",centreL)[ci], sep="")
  centre = paste0(panelL," ",centreL)[ci]

  distmfile = list.files(dist_dir[ci], recursive=F, full.names=T, pattern=".Rdata")
  result = foreach(i=1:length(distmfile)) %dopar% {
    start1 = Sys.time()
    cat("\n", distmfile[i], ", loading dist matrix", sep="")
    d = get(load(distmfile[i]))
    d = as.matrix(d)
    try ({
      if (sum(is.na(d))>0) distmfile[i] = paste0(distmfile[i],"_NA")
      if (sum(is.nan(d))>0) distmfile[i] = paste0(distmfile[i],"_NaN")
      if (sum(d==Inf)>0) distmfile[i] = paste0(distmfile[i],"_Inf")
    })
    write.csv(d, file=paste0(distmfile[i],".csv"))
  }
}

TimeOutput(start)









