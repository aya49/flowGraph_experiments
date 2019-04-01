# CountStats
# aya43@sfu.ca 20161220

root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

#Input
matrix_dir = paste(result_dir, "/matrix", sep="")
mcp = c("CountAdj") #countAdj comes first, to set which columns are deleted based on cellcountthres

#Output
png_dir = paste(result_dir, "/count_stats.png", sep="")

libr(foreach)
libr(doMC)
libr(colorspace)
source("~/projects/IMPC/code/_funcAlice.R")

no_cores = 6#detectCores() - 4
registerDoMC(no_cores)

start = Sys.time()

mm = get(load(paste0(matrix_dir, mcp, ".Rdata")))
markers = unlist(strsplit(colnames(mm)[which.max(nchar(colnames(mm)))],"[+-]"))
phenolevel = getPhen(colnames(mm), phenotype=F, markers=markers)$phenoLevel

start1 = Sys.time()
countThres = seq(20,2000,20)
underCountlist = foreach(c = 1:length(countThres)) %dopar% {
  return(which(apply(mm,2,function(x) all(x<=countThres[c]))))
}
TimeOutput(start1)

start1 = Sys.time()
k0 = seq(1,max(phenolevel))
underklist = foreach(k = 1:length(k0)) %dopar% {
  return(which(phenolevel<=k0[k]))
}
TimeOutput(start1)

start1 = Sys.time()
underBoth = foreach(c = 1:length(countThres)) %dopar% {
  ub = rep(0,length(k0))
  for (k in 1:length(k0)) {
    l = length(intersect(underCountlist[[c]],underklist[[k]]))
    if (!(length(l)>0)) l = 0
    ub[k] = l
  }
  return(ub)
}
underBoth = do.call(rbind, underBoth)
TimeOutput(start1)

png(png_dir, width=400, height=400)
colour = rainbow_hcl(ncol(underBoth))
plot(countThres, underBoth[,ncol(underBoth)], type="l", col=colour[ncol(underBoth)], xlab="Count threshold", ylab="# of cell populations with count <= Count threshold, for all samples")
for(i in (ncol(underBoth)-1):1) {
  lines(countThres, underBoth[,i], col=colour[i])
}
legend("topleft",legend=paste0("<= phenolevel ", c(1:ncol(underBoth))), fill=colour)
graphics.off()

TimeOutput(start)





