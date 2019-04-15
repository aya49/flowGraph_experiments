## Input: count adjusted matrix --> Output: add a third category of samples such that it is a mix of normals and aml (must run 00_data_flowcap.R before running this, every time)
#aya43@sfu.ca 20151228

## root directory
root = "~/projects/flowtype_metrics"
setwd(root)

randomind = NULL
for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  if (!grepl("flowcap",result_dir) | grepl("artificial",result_dir)) next
  
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
  
  
  ## libraries
  source("source/_funcAlice.R")
  source("source/_funcdist.R")
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
  
  
  #Prepare data
  meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
  meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
  feat_file_cell_countAdj0 = get(load(paste0(feat_file_cell_countAdj_dir,".Rdata")))
  
  
  
  start = Sys.time()
  
  #choose matching samples based on similar total cell count; do once
  normali = which(meta_file0$class=="normal")
  amli = which(meta_file0$class=="aml")
  if (is.null(randomind)) {
    for (i in 1:length(normali)) {
      randomind[[i]] = list()
      randomind[[i]]$normal = sample(length(normali), matchsamples)
      randomind[[i]]$aml = sample(length(normali), matchsamples)
    }
  }
  
  weight = 1/(2*matchsamples)
  
  feat_file_cell_countAdj2 = foreach (randomi = randomind, .combine="rbind") %dopar% {
    ni = normali[randomi$normal]
    ai = amli[randomi$aml]
    frankenstein = weight*colSums(feat_file_cell_countAdj0[append(ni,ai),])
    return(frankenstein)
  }
  feat_file_cell_countAdj = rbind(feat_file_cell_countAdj0,feat_file_cell_countAdj2)
  meta_file = rbind(meta_file0, data.frame(id=max(meta_file0$id)+c(1:length(randomind)), class=rep("frankenstein", length(randomind))))
  
  feat_file_cell_prop = feat_file_cell_countAdj/feat_file_cell_countAdj[,1]
  
  save(feat_file_cell_prop, file=paste0(feat_file_cell_prop_dir_,".Rdata"))
  save(feat_file_cell_countAdj, file=paste0(feat_file_cell_countAdj_dir_,".Rdata"))
  save(meta_file, file=paste0(meta_file_dir_,".Rdata"))
  save(meta_cell, file=paste0(meta_cell_dir_,".Rdata"))
  
  time_output(start)
}