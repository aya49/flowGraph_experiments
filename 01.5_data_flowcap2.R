## input: count adjusted matrix (flowcap data set only) (must run 00_data_flowcap.R before running this, every time)
## output: count adjusted matrix + frankenstein samples
## process: 
## - add a 3rd category of samples as a linear combination (weight = 1/(2*matchsamples)) of randomly sampled (5) control & (5) aml samples -- sampling done once, combinations same for all flowcap panels
## - one frankenstein sample for every control sample is made

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

id_col = "id"
target_col = "class" #column with control/experiment
control = "control" #control value in target_col column

cutoff = c(Inf) #c(.6) #if TMM-peak>cutoff, then apply peak instead of TMM; run this script and look at norm_fdiffplot plot to determine this number
layer_norm = 4 #0 #calculate TMM using only phenotypes in this layer; set to 0 if do for all layers
cellCountThres = 1000 #don't use phenotypes with cell count lower than cellCountThres

matchsamples = 5 # number of samples from each normal and aml to mix


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
  meta_dir_ = paste0(result_dir,".artificial/meta"); dir.create(meta_dir_, recursive=T)
  meta_file_dir_ = paste(meta_dir_, "/file", sep="")
  meta_cell_dir_ = paste(meta_dir_, "/cell", sep="")
  feat_dir_ = paste(result_dir, ".artificial/feat", sep=""); dir.create(feat_dir_, recursive=T)
  feat_file_cell_countAdj_dir_ = paste(feat_dir_, "/file-cell-countAdj", sep="")
  feat_file_cell_prop_dir_ = paste(feat_dir_, "/file-cell-prop", sep="")
  
  
  #Prepare data
  meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
  meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
  feat_file_cell_countAdj0 = get(load(paste0(feat_file_cell_countAdj_dir,".Rdata")))
  
  
  
  start = Sys.time()
  
  #choose matching samples based on similar total cell count; do once
  normali = which(meta_file0$class=="control")
  amli = which(meta_file0$class=="aml")
  if (is.null(randomind)) {
    for (i in 1:length(amli)) {
      randomind[[i]] = list()
      randomind[[i]]$normal = sample(normali, matchsamples)
      randomind[[i]]$aml = sample(amli, matchsamples)
    }
  }
  
  weight = 1/(2*matchsamples)
  
  feat_file_cell_countAdj2 = foreach (randomi = randomind, .combine="rbind") %dopar% {
    frankenstein = weight*colSums(feat_file_cell_countAdj0[append(randomi$normal,randomi$aml),])
    return(frankenstein)
  }
  rownames(feat_file_cell_countAdj2) = c((1+nrow(feat_file_cell_countAdj0)):(nrow(feat_file_cell_countAdj2)+nrow(feat_file_cell_countAdj0)))
  feat_file_cell_countAdj = rbind(feat_file_cell_countAdj0,feat_file_cell_countAdj2)
  meta_file = rbind(meta_file0, data.frame(
    id=max(meta_file0$id)+c(1:length(randomind)), 
    tube=rep(as.numeric(gsub("p","",str_extract(result_dir,"p[0-9]$"))), length(randomind)), 
    specimen=rep(0, length(randomind)), 
    class=rep("frankenstein", length(randomind)), 
    type=append(rep("test",floor(length(randomind)/2)), rep("train",ceiling(length(randomind)/2)))))

  feat_file_cell_prop = feat_file_cell_countAdj/feat_file_cell_countAdj[,1]
  
  save(feat_file_cell_prop, file=paste0(feat_file_cell_prop_dir_,".Rdata"))
  save(feat_file_cell_countAdj, file=paste0(feat_file_cell_countAdj_dir_,".Rdata"))
  save(meta_file, file=paste0(meta_file_dir_,".Rdata"))
  save(meta_cell, file=paste0(meta_cell_dir_,".Rdata"))
  
  time_output(start)
}