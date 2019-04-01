# aya43@sfu.ca 20170116
# Evaluate all distance metrics

root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

options(na.rm=T)

#Input
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
dist_dir = paste(result_dir, "/dist", sep="")

libr(stringr)
libr(foreach)
libr(doMC)
source("~/projects/IMPC/code/_funcAlice.R")


start = Sys.time()

no_cores = detectCores()-1
registerDoMC(no_cores)

distmfile = list.files(dist_dir, recursive=F, full.names=T, pattern=".Rdata")
result = foreach(i=1:length(distmfile)) %dopar% {
  start1 = Sys.time()
  cat("\n", distmfile[i], ", loading dist matrix", sep="")
  d = get(load(distmfile[i]))
  d = as.matrix(d)
  if (sum(is.na(d))>0) distmfile[i] = paste0(distmfile[i],"_NA")
  if (sum(is.nan(d))>0) distmfile[i] = paste0(distmfile[i],"_NaN")
  if (!is.na(sum(d==Inf))) distmfile[i] = paste0(distmfile[i],"_Inf")
  write.csv(d, file=paste0(distmfile[i],".csv"))
}

TimeOutput(start)









