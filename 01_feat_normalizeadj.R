## input: original count feature matrix 
## output: normalized countadj feature matrix, meta_cell, meta_cell_childpn/parent, and these converted into graphs (edge list /vertices) convertable into igraph with graph_from_data_frame function
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
libr(c("stringr", "plyr", 
       "pracma", "fitdistrplus", "flowType",
       "foreach","doMC"))


## cores
no_cores = detectCores()-1
registerDoMC(no_cores)


## options
writecsv = F
options(stringsAsFactors=FALSE)
options(na.rm=T)

cutoff = c(Inf) #c(.6) #if TMM-peak>cutoff, then apply peak instead of TMM; run this script and look at norm_fdiffplot plot to determine this number
layer_norm = c(2,4) #0 #calculate TMM using only phenotypes in this layer; set to 0 if do for all layers
cellCountThres = .01 # don't use phenotypes with cell count all lower than cellCountThres



result_dirs =list.dirs(paste0(root,"/result"),full.names=T,recursive=F)
for (result_dir in result_dirs) {
  if (grepl("ctrl|pos",result_dir)) next
  print(fileNames(result_dir))
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta")
  meta_file_dir = paste(meta_dir, "/file", sep="")
  meta_cell_dir = paste(meta_dir, "/cell", sep="")
  feat_dir = paste(result_dir, "/feat", sep="")
  feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count",sep="")
  
  
  ## output directories
  feat_file_cell_countAdj_dir = paste(feat_dir, "/file-cell-countAdj", sep="")
  norm_dir = paste(result_dir, "/cell_count_norm",sep=""); dir.create(norm_dir,showWarnings=F)
  norm_factor_dir = paste(norm_dir, "/norm_factor", sep=""); dir.create(norm_factor_dir,showWarnings=F) #plot of norm factor for each file
  norm_factor_diff_dir = paste(norm_dir, "/norm_factor_diff", sep="")
  
  
  ## load data
  meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
  feat_file_cell_count0 = as.matrix(get(load(paste0(feat_file_cell_count_dir,".Rdata"))))
  meta_cell0 = getPhen(colnames(feat_file_cell_count0))

  
  
  start = Sys.time()
  
  
  ## prepare data
  
  #save original objects for testing
  feat_file_cell_count = feat_file_cell_count0
  meta_file = meta_file0
  
  # prepare feat_file_cell_counts
  x = x0 = as.matrix(feat_file_cell_count)[,-1] # take out total cell count
  maxx = max(x0[is.finite(x0)])
  rootc = feat_file_cell_count[,1]
  refsample = which.min(abs( rootc-median(rootc[meta_file$class=="control"]) )) #reference column: median total count out of all control files
  
  # extract cell populations that would define TMM (layer/count)
  if (layer_norm[1]>0) 
    x = x[, colnames(x0) %in% meta_cell0$phenotype[
      meta_cell0$phenolevel>=layer_norm[1] & 
        meta_cell0$phenolevel<=layer_norm[2]] ]
  x = x[,sapply(1:ncol(x), function(y) any(x[,y]>cellCountThres*maxx))]
  
  # prepare paths to plot to
  pngnames = paste0(norm_factor_dir,"/", meta_file$class, "_", meta_file$id,".png")
  mains = paste0("mean count vs. ln fold change over ref sample ", 
                 meta_file$class[refsample], " on layer ", 
                 paste(layer_norm,collapse="-"))
  
  
  ## calculate absolute count TMM
  fresult = tmm(x,x0,rootc,refsample,cutoff=Inf,plotimg=T,pngnames=pngnames,mains=mains,no_cores=no_cores,samplesOnCol=F)
  f0 = fresult$f
  fdiff0 = fresult$fdiff
  
  m00 = as.matrix(ldply(c(1:nrow(feat_file_cell_count0)), function(x) feat_file_cell_count0[x,]*f0[x]))
  dimnames(m00) = dimnames(feat_file_cell_count0)

  # plot difference between TMM and peak for all files
  png(paste0(norm_factor_dir,"/all.png") , width=700, height=700)
  plot(sort(abs(fdiff0)), cex=.4, ylim=c(0,3), main="cell-count-norm-factor_f_diff_from_peak_abs")
  lines(sort(abs(fdiff0)), col="blue")
  graphics.off()
  

  
  #trim columns with too small of a mean/sd -- not used
  # try ({
  #   small = ldply(1:ncol(feat_file_cell_countAdj), function(i) {
  #     fit = NULL
  #     try({
  #       fit = fitdist(unlist(feat_file_cell_countAdj[,i]), "nbinom")
  #     })
  #     if (is.null(fit)) return(rep(0,4))
  #     # plot(density(feat_file_cell_countAdj[,i]))
  #     return( c(fit$estimate[1], fit$sd[1],fit$estimate[2], fit$sd[2]))
  #   }, .parallel=T)
  #   errsum = sum(apply(small,1,function(x) all(x==0)))
  #   muuu = small[,3:4]; size = small[,1:2]; rm(small)
  #   colnames(size)[2] = colnames(muuu)[2] = "sd"
  #   musd = muuu[,1]/muuu[,2] # mean/sd for every column
  #   
  #   png(file=paste0(feat_dir,"_countadj_negbin.png"), width=300, height=800)
  #   par(mfrow=c(3,1))
  #   plot(density(musd, na.rm=T), main="negative binomial per cell pop, mu/sd")
  #   plot(muuu, main="negative binomial per cell pop, mu vs sd")
  #   plot(size, main="negative binomial per cell pop, size vs sd")
  #   graphics.off()
  #   
  #   layers = max(meta_cell$phenolevel):min(meta_cell$phenolevel)
  #   png(file=paste0(feat_dir,"_countadj_negbin_bylayer.png"), width=300, height=300*length(layers))
  #   par(mfrow=c(length(layers),1))
  #   for (l in layers) { if (l>0)
  #     lind = meta_cell$phenolevel==l
  #     plot(density(musd[lind], na.rm=T), main=paste0("negative binomial per cell pop, mu/sd for layer ", l))
  #   }
  #   graphics.off()
  #   
  # })
  
  
  
  ## trim & save
  
  #trim columns with too many 0's or low values
  minclassn = min(table(meta_file$class))
  col_min0 = apply(m00, 2, function(x) sum(x>0)>(minclassn*.5))
  maxx = max(x[is.finite(x)])
  col_cnts = apply(m00, 2, function(x) any(x>cellCountThres*maxx))
  finalinds = col_min0 & col_cnts
  
  
  m0 = m00[,finalinds]
  save(m0, file=paste0(feat_file_cell_countAdj_dir,".Rdata"))
  if (writecsv) write.csv(m0, file=paste0(feat_file_cell_countAdj_dir,".csv"), row.names=T)
  
  c0 = feat_file_cell_count0[,finalinds]
  save(c0, file=paste0(feat_file_cell_count_dir,".Rdata"))
  if (writecsv) write.csv(c0, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)
  
  
  save(f0, file=paste0(norm_factor_dir,".Rdata"))
  if (writecsv) write.csv(f0, file=paste0(norm_factor_dir,".csv"), row.names=T)
  save(fdiff0, file=paste0(norm_factor_diff_dir,".Rdata"))
  if (writecsv) write.csv(fdiff0, file=paste0(norm_factor_diff_dir,".csv"), row.names=T)
  
  time_output(start)
}



