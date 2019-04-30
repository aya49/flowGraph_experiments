## Input: original features --> Output: random walk features (edge matrix + paths -- paths need to be saved one at a time, memory will overload otherwise)
# aya43@sfu.ca 20161220

## root directory
root = "~/projects/flowtype_metrics"
setwd(root)

## libraries
source("source/_funcAlice.R")
libr(c("foreach","doMC",
       "stringr","plyr","Matrix"))

#Setup Cores
no_cores = 3 #detectCores()-6
registerDoMC(no_cores)

## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)


overwrite = T
writecsv = F

for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  if (!grepl("pregnancy",result_dir)) next()
  # result_dir = paste0(root, "/result/impc_panel1_sanger-spleen") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta") # meta files directory
  meta_file_dir = paste(meta_dir, "/file", sep="") #meta for rows (sample)
  feat_dir = paste(result_dir, "/feat", sep="") #feature files directory

  ## output directories
  

  #data paths
  feat_types = gsub(".Rdata","",list.files(path=feat_dir,pattern=glob2rx("*.Rdata")))
  feat_types = feat_types[!grepl("paired",feat_types)]
  
  

  start = Sys.time()
  
  meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))

  for (feat_type in feat_types) {
      cat("\n", feat_type, " ",sep="")
      start2 = Sys.time()

      m = m0 = Matrix(as.matrix(get(load(paste0(feat_dir,"/", feat_type,".Rdata")))))
      meta_file = meta_file0[match(rownames(m0),meta_file0$id),]
      
      for (pi in meta_file$patient) {
        pii = meta_file$patient==pi
        mp = m[pii,]ø
        mpm = colMeans(mp)
        m[pii,] = foreach(i=1:nrow(mp),.combine="rbind") %do% { return(mp[i,]-mpm) }
      }
      save(m, file=paste0(feat_dir,"/", feat_type,"-paired.Rdata"))
      if (writecsv) write.csv(m, file=paste0(feat_dir,"/", feat_type,"-paired.csv"))
      
  }
  
  time_output(start)
  
}


