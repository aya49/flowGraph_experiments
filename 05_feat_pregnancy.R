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

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)

## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)


overwrite = T
writecsv = F

count_feature = "file-cell-countAdj"

for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  if (!grepl("pregnant|bodenmiller",result_dir)) next()
  # result_dir = paste0(root, "/result/impc_panel1_sanger-spleen") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  print(result_dir)
  
  start = Sys.time()
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta") # meta files directory
  meta_file_dir = paste(meta_dir, "/file", sep="") #meta for rows (sample)
  feat_dir = paste(result_dir, "/feat", sep="") #feature files directory
  
  feat_dir_ = paste0(feat_dir,"_unpaired"); dir.create(feat_dir_, showWarnings=F)
  
  # meta
  meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
  
  #data paths
  feat_types = gsub(".Rdata","",list.files(path=feat_dir, full.names=F, pattern=".Rdata"))

  
  a= llply (feat_types, function(feat_type) {
    cat("\n", feat_type, " ",sep="")
    start2 = Sys.time()
    
    if (file.exists(paste0(feat_dir_,"/", feat_type,".Rdata"))) {
      m = m0 = Matrix(as.matrix(get(load(paste0(feat_dir_,"/", feat_type,".Rdata")))))
    } else {
      m = m0 = Matrix(as.matrix(get(load(paste0(feat_dir,"/", feat_type,".Rdata")))))
      save(m0, file=paste0(feat_dir_,"/", feat_type,".Rdata"))
    }
    meta_file = meta_file0[match(rownames(m0),meta_file0$id),]
    
    
    for (pi in meta_file$patient) {
      pii = meta_file$patient==pi
      mp = m[pii,]
      mpm = colMeans(as.matrix(mp))
      m[pii,] = as.matrix(ldply(1:nrow(mp), function(i) mp[i,]-mpm ))
    }
    save(m, file=paste0(feat_dir,"/", feat_type,".Rdata"))
    if (writecsv) write.csv(m, file=paste0(feat_dir,"/", feat_type,"-paired.csv"))
  }, .parallel=T)
  time_output(start, result_dir)
  
}


