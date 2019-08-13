## input: features 
## output: p values and their corelation between train/test sets; p value qq plots, recall/precision, number of significant features

set.seed(10)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)


## libraries
source("source/_func.R")
libr(c("stringr","Matrix", "plyr",
       "foreach","doMC"))


## cores
no_cores = detectCores()-1
registerDoMC(no_cores)



## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

# cvn-fold cross validation

overwrite = F #overwrite?
writecsv = F

pthress = c(.05,.025,.01) # p value sig threshold for t test

## calcuate p values!

# format of p values:
# tri = foldsip[[uc]]$indices$train
# tei = foldsip[[uc]]$indices$test
# pv_tr_ = foldsip[[uc]]$p[[ptype]][[adj]]$train
# pv_tr2_ = foldsip[[uc]]$p[[ptype]][[adj]]$p_train2
# pv_te_ = foldsip[[uc]]$p[[ptype]][[adj]]$p_test
# pv_all_ = foldsip[[uc]]$p[[ptype]][[adj]]$p_all


start = Sys.time()


## make table for stats
table = NULL
result_dirs = list.dirs(paste0(root, "/result"), full.names=T, recursive=F)
for (result_dir in result_dirs) {
  print(result_dir)
  data = fileNames(result_dir)
  
  ## input directories
  sum_dir = paste0(result_dir,"/feat_summaries")
  sum_p_dir = paste0(result_dir,"/feat_pval")
  sum_m_dir = paste0(result_dir,"/feat_mean")
  
  p_paths = list.files(sum_p_dir,full.names=T)
  m_paths = list.files(sum_m_dir,full.names=T)
  
  pp = str_split(gsub(".Rdata","",fileNames(p_paths)),"_")
  ppta = str_split(sapply(pp, function(x) x[3]),"-")
  mm = str_split(gsub(".Rdata","",fileNames(m_paths)),"_")
  
  fe = sapply(pp,function(x) x[1])
  table = rbind(table, data.frame(
    pathp=p_paths, pathm=m_paths[match(fe,mm)], data=data,
    feat=fe,
    class=sapply(pp,function(x) x[2]),
    test=sapply(ppta,function(x) x[1]),
    test_adj=sapply(ppta,function(x) x[2])))
} # result



## calculate corr
for (i in 1:nrow(table)) {
  pvs = get(load(table$pathp[i]))
  
  ## calculate correlations between p values test & train/2
  pv_trl = -log(pvs$train)
  pv_tel = -log(pvs$test)
  # pv_alll
  
  table$pcorr[i] = cor(pv_trl,pv_tel, method="spearman")
  table$pcorrp[i] = cor.test(pv_trl,pv_tel, method="spearman")$p.value
  if (is.na(table$pcorr[i])) { 
    table$pcorr[i] = 1
    table$pcorrp[i] = 0
  }
  # pcorr2 = cor(pv_trl2,pv_tel, method="spearman")
  # pcorrp2 = cor.test(pv_trl2,pv_tel, method="spearman")$p.value
}


## calculate p value sig stats
tbl = NULL
for (i in 1:nrow(table)) {
  for (pthres in pthress) {
    # get sig
    table$pthres[i] = pthres
    
    pvs = get(load(table$pathp[i]))
    pv_tr_ = pvs$train
    pv_te_ = pvs$test
    pv_all_ = pvs$all
    
    tr_sig = pv_tr_<pthres
    # tr2_sig = pv_tr2_<pthres
    te_sig = pv_te_<pthres
    all_sig = pv_all_<pthres
    
    table$p_thres[i]=pthres
    table$m_all_sig[i]=sum(all_sig)
    table$m_all[i]=length(all_sig)
    table$m_all_perc[i]=sum(all_sig)/length(all_sig)
    table$m_train_sig[i]=sum(tr_sig)
    table$m_train[i]=length(tr_sig)
    table$m_train_perc[i]=sum(tr_sig)/length(tr_sig)
    table$m_test_sig[i]=sum(te_sig)
    table$m_test[i]=length(te_sig)
    table$m_test_perc[i]=sum(te_sig)/length(te_sig)
    table$m_overlap_sig[i] = overlap = sum(tr_sig & te_sig)
    if (overlap==0) {
      table$rec[i] = table$prec[i] = table$f[i] = 0
    } else {
      table$rec[i] = rec = overlap/sum(te_sig)
      table$prec[i] = prec = overlap/sum(tr_sig)
      table$f[i] = 2*((prec*rec)/(prec+rec))
    }
    
    tbl = rbind(tbl,table[i,])
  }
}
save(tbl, file=paste0(root,"/pvals.Rdata")) # table

time_output(start)


