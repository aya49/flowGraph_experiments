## input: features + meta_file
## output: features normalized for patient; 
## process: for the countAdj feature (prop made afterwards),
## - note: each patient has 4 classes or is sampled at 4 time points throughout pregnancy, we use the first time point as our control to calculate p values in the previous p value script
## - take each patient and calculate the mean feature vector of her 4 fcm files
## - get the difference between her 4 fcm files and the mean; save this as her new fcm file features
## - output features marked with "-paired"


## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)


## libraries
source("source/_func.R")
libr(c("foreach","doMC",
       "stringr","plyr","Matrix"))


## cores
no_cores = detectCores()-1
registerDoMC(no_cores)


## options for script
options(stringsAsFactors=FALSE)
options(na.rm=T)

writecsv = F

count_feature = "file-cell-countAdj"



result_dirs = list.dirs(paste0(root, "/result"), full.names=T, recursive=F)
for (result_dir in result_dirs) {
  if (!grepl("pregnant|bodenmiller",result_dir)) next()
  print(result_dir)
  
  start = Sys.time()
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta") # meta files directory
  meta_file_dir = paste(meta_dir, "/file", sep="") #meta for rows (sample)
  feat_dir = paste(result_dir, "/feat", sep="") #feature files directory
  
  
  ## ouput directories
  feat_dir_ = paste0(feat_dir,"_unpaired"); dir.create(feat_dir_, showWarnings=F)
  
  
  ## start: take difference from mean of samples from one subject
  
  # load meta
  meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
  
  # feat paths
  feat_types = gsub(".Rdata","",list.files(path=feat_dir, full.names=F, pattern=".Rdata"))

  a = llply (feat_types, function(feat_type) {
    cat("\n", feat_type, " ",sep="")
    start2 = Sys.time()
    
    if (file.exists(paste0(feat_dir_,"/", feat_type,".Rdata"))) {
      m = Matrix(as.matrix(get(load(paste0(
        feat_dir_,"/", feat_type,".Rdata")))))
    } else {
      m = Matrix(as.matrix(get(load(paste0(
        feat_dir,"/", feat_type,".Rdata")))))
      save(m, file=paste0(feat_dir_,"/", feat_type,".Rdata"))
    }
    meta_file = meta_file0[match(rownames(m),meta_file0$id),]
    
    for (pi in meta_file$subject) {
      pii = meta_file$patient==pi
      mp = m[pii,]
      mpm = colMeans(as.matrix(mp))
      m[pii,] = as.matrix(ldply(1:nrow(mp), function(i) mp[i,]-mpm ))
    }
    save(m, file=paste0(feat_dir,"/", feat_type,".Rdata"))
    if (writecsv) write.csv(m, file=paste0(feat_dir,"/", feat_type,"-paired.csv"))
  }, .parallel=T)
  
  time_output(start)
}


