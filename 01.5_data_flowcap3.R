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
no_cores = detectCores()-1
registerDoMC(no_cores)


## options
writecsv = F
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

amlprop = 1/2 # proportion of control sample to make as "aml"



start = Sys.time()

randomind = NULL
result_dirs =list.dirs(paste0(root,"/result"),full.names=T,recursive=F)
for (result_dir in result_dirs) {
  if (!grepl("flowcap",result_dir) | grepl("pos",result_dir) | grepl("ctrl",result_dir)) next
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta")
  meta_file_dir = paste(meta_dir, "/file", sep="")
  meta_cell_dir = paste(meta_dir, "/cell", sep="")
  feat_dir = paste(result_dir, "/feat", sep="")
  feat_file_cell_countAdj_dir = paste(feat_dir, "/file-cell-countAdj", sep="")
  
  
  ## output directories
  result_dir_ = paste0(result_dir,"_ctrl")
  meta_dir_ = paste0(result_dir_,"/meta"); dir.create(meta_dir_, recursive=T)
  meta_file_dir_ = paste(meta_dir_, "/file", sep="")
  meta_cell_dir_ = paste(meta_dir_, "/cell", sep="")
  feat_dir_ = paste(result_dir_, "/feat", sep=""); dir.create(feat_dir_)
  feat_file_cell_countAdj_dir_ = paste(feat_dir_, "/file-cell-countAdj", sep="")
  feat_file_cell_prop_dir_ = paste(feat_dir_, "/file-cell-prop", sep="")
  
  
  ## prepare data
  file.copy(file.path(meta_dir,list.files(meta_dir)), meta_dir_)
  meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
  feat_file_cell_countAdj0 = get(load(paste0(feat_file_cell_countAdj_dir,".Rdata")))
  
  
  ## extract controls
  m0 = feat_file_cell_countAdj0[meta_file0$class=="control",]
  p0 = m0/m0[,1]
  
  ## meta/file
  meta_file = meta_file0[meta_file0$class=="control",]
  if (is.null(randomind)) 
    randomind = sample(1:nrow(meta_file), floor(amlprop*nrow(meta_file)))
  meta_file$class[randomind] = "exp"
  

  ## save
  save(p0, file=paste0(feat_file_cell_prop_dir_,".Rdata"))
  save(m0, file=paste0(feat_file_cell_countAdj_dir_,".Rdata"))
  save(meta_file, file=paste0(meta_file_dir_,".Rdata"))
}

time_output(start)
