## input: countadj
## output: some useful plots; change/add if want


## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metrics"
setwd(root)
for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  # result_dir = paste0(root, "/result/impc_panel1_sanger-spleen") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  
## input directories
feat_dir = paste(result_dir, "/feat", sep="")

## output directories
stat_dir = paste(result_dir, "/stats", sep=""); dir.create(stat_dir, showWarnings=F)
png_dir = paste(stat_dir, "/stats_count.png", sep="")

## libraries
source("source/_funcAlice.R")
libr(c("stringr",
       "foreach","doMC",
       "colorspace"))

## cores
no_cores = 10# detectCores()-1
registerDoMC(no_cores)











## options
options(stringsAsFactors=F)
# options(device="cairo")
options(na.rm=T)

countThres = seq(20,2000,20)
feat_count = c("file-cell-countAdj") #specify countmatrix to use



#Prepare data
mm = get(load(paste0(feat_dir, "/", feat_count, ".Rdata")))
markers = unlist(strsplit(colnames(mm)[which.max(nchar(colnames(mm)))],"[+-]"))
phenolevel = str_count(colnames(mm), "[+-]")

k0 = seq(1,max(phenolevel)) #layers to plot











start = Sys.time()


# #trim low count phenotypes
# underCountlist = foreach(c = 1:length(countThres)) %dopar% { return(sum(apply(mm,2,function(x) all(x<=countThres[c])))) }
# 
# #trim layers
# underklist = foreach(k = 1:length(k0)) %dopar% { return(which(phenolevel<=k0[k])) }

# TimeOutput(start)
# start = Sys.time()

## collate data according to count and layer thresholds
underBoth = foreach(c = 1:length(countThres)) %dopar% {
  ub = rep(0,length(k0))
  for (k in 1:length(k0)) {
    ub[k] = sum( apply(mm,2,function(x) all(x<=countThres[c])) & phenolevel<=k0[k] )
  }
  return(ub)
}
underBoth = Reduce("rbind", underBoth)

TimeOutput(start)

## plot
png(png_dir, width=400, height=400)
colour = rainbow_hcl(ncol(underBoth))
plot(countThres, underBoth[,ncol(underBoth)], type="l", col=colour[ncol(underBoth)], 
     xlab="Count threshold", ylab="# of cell populations with count <= Count threshold, for all samples")
for(i in (ncol(underBoth)-1):1) { lines(countThres, underBoth[,i], col=colour[i]) }
legend("topleft",legend=paste0("<= phenolevel ", c(1:ncol(underBoth))), fill=colour)
graphics.off()

TimeOutput(start)

}



