## input: original count feature matrix 
## output: normalized countadj feature matrix
## process:
## - normalize count using TMM trimmed mean normalization
## - reference fcm file is the one with median total count amongst the control files
## - plots for each fcm file, the kernel density estimate of its counts ratio over the reference file's counts
## - if the normalization factor is too far from peak ratio, directly use peak ratio, see output plots for which is used
## - remember to trim matrices!

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr(c("stringr", "pracma", "fitdistrplus",
       "plyr",
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
  print(result_dir)
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta")
  meta_cell_dir = paste(meta_dir, "/cell", sep="")
  meta_file_dir = paste(meta_dir, "/file", sep="")
  feat_dir = paste(result_dir, "/feat", sep="")
  feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")
  feat_file_cell_prop_dir = paste(feat_dir, "/file-cell-prop", sep="")
  
  
  ## output directories
  feat_file_cell_countAdj_dir = paste(feat_dir, "/file-cell-countAdj", sep="")
  norm_dir = paste(result_dir, "/cell_count_norm",sep=""); dir.create(norm_dir,showWarnings=F)
  norm_factor_dir = paste(norm_dir, "/norm_factor", sep=""); dir.create(norm_factor_dir,showWarnings=F) #plot of norm factor for each file
  norm_factor_diff_dir = paste(norm_dir, "/norm_factor_diff", sep="")
  
  
  #Prepare data
  meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
  meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
  feat_file_cell_count0 = get(load(paste0(feat_file_cell_count_dir,".Rdata")))
  feat_file_cell_prop0 = get(load(paste0(feat_file_cell_prop_dir,".Rdata")))
  
  
  
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
  
  
  
  ## trim ---------------------------------------------
  
  #trim columns with too many 0's or low values
  minclassn = min(table(meta_file[,target_col]))
  col_min0 = apply(feat_file_cell_countAdj, 2, function(x) sum(x>0)>(minclassn*.5))
  col_cnts = apply(feat_file_cell_countAdj, 2, function(x) any(x>cellCountThres))
  
  
  #trim columns with too small of a mean/sd
  try ({
    
    small = ldply(1:ncol(feat_file_cell_countAdj), function(i) {
      fit = NULL
      try({
        fit = fitdist(unlist(feat_file_cell_countAdj[,i]), "nbinom")
      })
      if (is.null(fit)) return(rep(0,4))
      # plot(density(feat_file_cell_countAdj[,i]))
      return( c(fit$estimate[1], fit$sd[1],fit$estimate[2], fit$sd[2]))
    }, .parallel=T)
    errsum = sum(apply(small,1,function(x) all(x==0)))
    muuu = small[,3:4]; size = small[,1:2]; rm(small)
    colnames(size)[2] = colnames(muuu)[2] = "sd"
    musd = muuu[,1]/muuu[,2] # mean/sd for every column
    
    png(file=paste0(feat_dir,"_countadj_negbin.png"), width=300, height=800)
    par(mfrow=c(3,1))
    plot(density(musd, na.rm=T), main="negative binomial per cell pop, mu/sd")
    plot(muuu, main="negative binomial per cell pop, mu vs sd")
    plot(size, main="negative binomial per cell pop, size vs sd")
    graphics.off()
    
    layers = max(meta_cell$phenolevel):min(meta_cell$phenolevel)
    png(file=paste0(feat_dir,"_countadj_negbin_bylayer.png"), width=300, height=300*length(layers))
    par(mfrow=c(length(layers),1))
    for (l in layers) { if (l>0)
      lind = meta_cell$phenolevel==l
      plot(density(musd[lind], na.rm=T), main=paste0("negative binomial per cell pop, mu/sd for layer ", l))
    }
    graphics.off()
    
  })
  
  
  
  #save
  finalinds = col_min0 & col_cnts
  
  feat_file_cell_countAdj_ = feat_file_cell_countAdj[,finalinds]
  save(feat_file_cell_countAdj_, file=paste0(feat_file_cell_countAdj_dir,".Rdata"))
  if (writecsv) write.csv(feat_file_cell_countAdj_, file=paste0(feat_file_cell_countAdj_dir,".csv"), row.names=T)
  
  feat_file_cell_count_ = feat_file_cell_count0[,finalinds]
  save(feat_file_cell_count_, file=paste0(feat_file_cell_countAdj_dir,".Rdata"))
  if (writecsv) write.csv(feat_file_cell_count_, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)
  
  feat_file_cell_prop_ = feat_file_cell_prop0[,finalinds]
  save(feat_file_cell_prop_, file=paste0(feat_file_cell_countAdj_dir,".Rdata"))
  if (writecsv) write.csv(feat_file_cell_prop, file=paste0(feat_file_cell_prop_dir,".csv"), row.names=T)
  
  meta_cell_ = meta_cell[finalinds,]
  save(meta_cell_, file=paste0(meta_cell_dir,".Rdata"))
  
  
  save(f0, file=paste0(norm_factor_dir,".Rdata"))
  if (writecsv) write.csv(f0, file=paste0(norm_factor_dir,".csv"), row.names=T)
  save(fdiff0, file=paste0(norm_factor_diff_dir,".Rdata"))
  if (writecsv) write.csv(fdiff0, file=paste0(norm_factor_diff_dir,".csv"), row.names=T)
  
  time_output(start)
}





