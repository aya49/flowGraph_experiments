## Input: original features --> Output: prints whether there are irregular values in features
# aya43@sfu.ca 20170316

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr("Matrix")

start = Sys.time()

tab = NULL
result_dirs = list.dirs(paste0(root, "/result"), full.names=T, recursive=F)
for (result_dir in result_dirs) {
  ## input directories
  feat_dir = paste(result_dir, "/feat", sep="") # feature files
  data = fileNames(result_dir)
  
  start1 = Sys.time()
  
  # feature paths
  feat_types = list.files(path=feat_dir,pattern=".Rdata")
  feat_types = gsub(".Rdata","",feat_types)

  result = ldply(feat_types, function(feat_type) {
    ## upload and prep feature matrix + meta
    m = as.matrix(Matrix(get(load(paste0(feat_dir,"/", feat_type,".Rdata")))))
    return(data.frame(data=data, feat=feat_type, nrow=nrow(m), ncol=ncol(m), inf=sum(is.infinite(m)), na=sum(is.na(m)), nan= sum(is.nan(m)), neg=sum(m<0), pos=sum(m>0), zero=sum(m==0), max=max(m[is.finite(m)])))
  })
  tab = rbind(tab,result)
}

write.csv(tab, file=paste0(root,"/featstats.csv"))
time_output(start)
