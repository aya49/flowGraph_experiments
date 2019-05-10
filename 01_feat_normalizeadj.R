## Input: original count matrix --> Output: normalized count matrix
#aya43@sfu.ca 20151228

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metrics"
setwd(root)

## libraries
source("source/_funcAlice.R")
source("source/_funcdist.R")
libr(c("stringr", "pracma",
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

for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  # result_dir = paste0(root, "/result/flowcap_panel6") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  if (grepl("artificial",result_dir)) next
  
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta")
  meta_cell_dir = paste(meta_dir, "/cell", sep="")
  meta_file_dir = paste(meta_dir, "/file", sep="")
  feat_dir = paste(result_dir, "/feat", sep="")
  feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")
  
  ## output directories
  feat_file_cell_countAdj_dir = paste(feat_dir, "/file-cell-countAdj", sep="")
  norm_dir = paste(result_dir, "/cell_count_norm",sep=""); dir.create(norm_dir,showWarnings=F)
  norm_factor_dir = paste(norm_dir, "/norm_factor", sep=""); dir.create(norm_factor_dir,showWarnings=F) #plot of norm factor for each file
  norm_factor_diff_dir = paste(norm_dir, "/norm_factor_diff", sep="")
  
  
  #Prepare data
  meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
  meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
  feat_file_cell_count0 = get(load(paste0(feat_file_cell_count_dir,".Rdata")))
  
  
  
  start = Sys.time()
  
  #calculate TMM
  
  #save original objects for testing
  feat_file_cell_count = feat_file_cell_count0
  meta_file = meta_file0
  
  #prepare feat_file_cell_counts
  x = x0 = as.matrix(feat_file_cell_count)[,-1] #take out total cell count
  if (layer_norm>0) x = as.matrix(x0[,colnames(x0)%in%meta_cell$phenotype[meta_cell$phenolevel==layer_norm] & sapply(1:ncol(x0), function(y) any(x0[,y]>cellCountThres))])
  lib.size = feat_file_cell_count[,1]
  refColumn = which.min(abs( lib.size - median(lib.size[grepl(control,meta_file[,target_col])]) )) #reference column: median total count out of all control files
  
  #prepare plot paths/titles
  pngnames = sapply(1:nrow(x), function(i) {
    paste0(norm_factor_dir,"/", meta_file[i,target_col], "_", meta_file[i,id_col],".png")
  })
  mains = sapply(1:nrow(x), function(i) paste0("mean count vs. ln fold change:\n", meta_file[i,target_col]," over refColumn ", meta_file[refColumn,target_col], "___layer-",layer_norm))
  
  ## calculate absolute count TMM, mostly taken from TMM
  fresult = tmm(x,x0,lib.size,refColumn,cutoff=Inf,plotimg=T,pngnames=pngnames,mains=mains,no_cores=no_cores,samplesOnCol=F)
  f0 = fresult$f
  fdiff0 = fresult$fdiff
  
  #plot difference between TMM and peak for all files
  pngname = paste0(norm_factor_dir,"/all.png")
  png(file=pngname , width=700, height=700)
  plot(sort(abs(fdiff0)), cex=.4, ylim=c(0,3), main="cell-count-norm-factor_f_diff_from_peak_abs")
  lines(sort(abs(fdiff0)), col="blue")
  graphics.off()
  
  
  
  feat_file_cell_countAdj = t(sapply(c(1:nrow(feat_file_cell_count0)), function(x) {feat_file_cell_count0[x,]*f0[x]}))
  colnames(feat_file_cell_countAdj) = colnames(feat_file_cell_count0)
  rownames(feat_file_cell_countAdj) = rownames(feat_file_cell_count0)
  
  #save
  save(feat_file_cell_countAdj, file=paste0(feat_file_cell_countAdj_dir,".Rdata"))
  if (writecsv) write.csv(feat_file_cell_countAdj, file=paste0(feat_file_cell_countAdj_dir,".csv"), row.names=T)
  save(f0, file=paste0(norm_factor_dir,".Rdata"))
  if (writecsv) write.csv(f0, file=paste0(norm_factor_dir,".csv"), row.names=T)
  save(fdiff0, file=paste0(norm_factor_diff_dir,".Rdata"))
  if (writecsv) write.csv(fdiff0, file=paste0(norm_factor_diff_dir,".csv"), row.names=T)
  
  time_output(start)
}





