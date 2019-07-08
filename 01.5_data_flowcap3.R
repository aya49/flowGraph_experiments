## input: count adjusted matrix (flowcap data set only) (must run 00_data_flowcap.R before running this, every time)
## output: count adjusted matrix with only control samples

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)


## libraries
source("source/_func.R")
libr(c("stringr",
       "foreach","doMC"))
# libr(flowDensity)

## cores
no_cores = 10#detectCores()-1
registerDoMC(no_cores)



## options
writecsv = F
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

id_col = "id"
target_col = "class" #column with control/experiment
control = "control" #control value in target_col column

cutoff = c(Inf) #c(.6) #if TMM-peak>cutoff, then apply peak instead of TMM; run this script and look at norm_fdiffplot plot to determine this number
layer_norm = 4 #0 #calculate TMM using only phenotypes in this layer; set to 0 if do for all layers
cellCountThres = 1000 #don't use phenotypes with cell count lower than cellCountThres

amlprop = 1/2 # proportion of control sample to make as "aml"


randomind = NULL
for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  if (!grepl("flowcap",result_dir) | grepl("artificial",result_dir) | grepl("ctrl",result_dir)) next
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta")
  meta_file_dir = paste(meta_dir, "/file", sep="")
  meta_cell_dir = paste(meta_dir, "/cell", sep="")
  feat_dir = paste(result_dir, "/feat", sep="")
  feat_file_cell_countAdj_dir = paste(feat_dir, "/file-cell-countAdj", sep="")
  
  ## output directories
  meta_dir_ = paste0(result_dir,".ctrl/meta"); dir.create(meta_dir_, recursive=T)
  meta_file_dir_ = paste(meta_dir_, "/file", sep="")
  meta_cell_dir_ = paste(meta_dir_, "/cell", sep="")
  feat_dir_ = paste(result_dir, ".ctrl/feat", sep=""); dir.create(feat_dir_, recursive=T)
  feat_file_cell_countAdj_dir_ = paste(feat_dir_, "/file-cell-countAdj", sep="")
  feat_file_cell_prop_dir_ = paste(feat_dir_, "/file-cell-prop", sep="")
  
  
  #Prepare data
  meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
  meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
  feat_file_cell_countAdj0 = get(load(paste0(feat_file_cell_countAdj_dir,".Rdata")))
  
  
  feat_file_cell_countAdj = feat_file_cell_countAdj0[meta_file0$class=="control",]
  meta_file = meta_file0[meta_file0$class=="control",]
  
  start = Sys.time()
  
  #choose matching samples based on similar total cell count; do once
  if (is.null(randomind)) 
    randomind = sample(1:nrow(meta_file), floor(amlprop*nrow(meta_file)))
  meta_file$class[randomind] = "exp"
  

  feat_file_cell_prop = feat_file_cell_countAdj/feat_file_cell_countAdj[,1]
  
  save(feat_file_cell_prop, file=paste0(feat_file_cell_prop_dir_,".Rdata"))
  save(feat_file_cell_countAdj, file=paste0(feat_file_cell_countAdj_dir_,".Rdata"))
  save(meta_file, file=paste0(meta_file_dir_,".Rdata"))

  time_output(start)
}