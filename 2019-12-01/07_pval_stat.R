## input: pval etc. directories 
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

pthress = c(.05,.01) # p value sig threshold for t test

## calcuate p values!

# format of p values:
# tri = foldsip[[uc]]$indices$train
# tei = foldsip[[uc]]$indices$test
# pv_tr_ = foldsip[[uc]]$p[[ptype]][[adj]]$train
# pv_tr2_ = foldsip[[uc]]$p[[ptype]][[adj]]$p_train2
# pv_te_ = foldsip[[uc]]$p[[ptype]][[adj]]$p_test
# pv_all_ = foldsip[[uc]]$p[[ptype]][[adj]]$p_all


start = Sys.time()


## make table for p values paths
table = NULL
result_dirs = list.dirs(paste0(root, "/result"), full.names=T, recursive=F)
for (result_dir in result_dirs) {
  if (grepl("bodenmiller$|pregnancy$",result_dir)) next
  data = fileNames(result_dir)
  
  ## input directories
  sum_dir = paste0(result_dir,"/feat_summaries")
  sum_p_dir = paste0(result_dir,"/feat_pval")

  p_paths = list.files(sum_p_dir,full.names=T)

  pp = str_split(gsub(".Rdata","",fileNames(p_paths)),"_")
  ppta = str_split(sapply(pp, function(x) x[3]),"-")

  fe = sapply(pp,function(x) x[1])
  table = rbind(table, data.frame(
    pathp=p_paths,
    data=data,
    feat=fe,
    class=sapply(pp,function(x) x[2]),
    test=sapply(ppta,function(x) x[1]),
    test_adj=sapply(ppta,function(x) x[2])))
} # result

# table = table[table$test=="t" & table$test_adj%in%c("BY","none"),]


## add overlap paths
for (i in 1:nrow(table)) 
  table$patho[i] = paste0(root,"/result/",table$data[i],"/feat_overlap/",table$feat[i],"_",table$class[i],".Rdata")


## calculate corr
start1 = Sys.time()
# l_ply(loopInd(1:nrow(table),no_cores), function(ii) {
for (i in 1:nrow(table)) { 
  # try ({
  
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
# })
}
# },.parallel=T)
time_output(start1)


## calculate p value sig stats
start1 = Sys.time()
tbl = llply(loopInd(1:nrow(table),no_cores), function(ii) {
  Reduce(rbind,
         llply(ii, function(i) { try ({
           # for (i in ii) {
           pvs = get(load(table$pathp[i]))
           Reduce(rbind,
                  llply(pthress, function(pthres) {
                    # get sig
                    tt = table[i,]
                    
                    tt$pthres = pthres
                    
                    pv_tr_ = pvs$train
                    pv_te_ = pvs$test
                    pv_all_ = pvs$all
                    
                    tr_sig = pv_tr_<pthres
                    # tr2_sig = pv_tr2_<pthres
                    te_sig = pv_te_<pthres
                    all_sig = pv_all_<pthres
                    
                    tt$p_thres=pthres
                    tt$m_all_sig=sum(all_sig)
                    tt$m_all=length(all_sig)
                    tt$m_all_perc=sum(all_sig)/length(all_sig)
                    tt$m_train_sig=sum(tr_sig)
                    tt$m_train=length(tr_sig)
                    tt$m_train_perc=sum(tr_sig)/length(tr_sig)
                    tt$m_test_sig=sum(te_sig)
                    tt$m_test=length(te_sig)
                    tt$m_test_perc=sum(te_sig)/length(te_sig)
                    tt$m_overlap_sig = overlap = sum(tr_sig & te_sig)
                    if (overlap==0) {
                      tt$rec = tt$prec = tt$f = 0
                    } else {
                      tt$rec = rec = overlap/sum(te_sig)
                      tt$prec = prec = overlap/sum(tr_sig)
                      tt$f = 2*((prec*rec)/(prec+rec))
                    }
                    
                    # calculate corr again
                    pv_trl = -log(pv_tr_[tr_sig | te_sig])
                    pv_tel = -log(pv_te_[tr_sig | te_sig])
                    # pv_alll
                    
                    if (length(pv_trl)>1) {
                      tt$pcorr_sigonly = cor(pv_trl,pv_tel, method="spearman")
                      tt$pcorrp_sigonly = cor.test(pv_trl,pv_tel, method="spearman")$p.value
                    } else {
                      tt$pcorr_sigonly = 1
                      tt$pcorrp_sigonly = 0
                    }
                    # pcorr2 = cor(pv_trl2,pv_tel, method="spearman")
                    # pcorrp2 = cor.test(pv_trl2,pv_tel, method="spearman")$p.value
                    
                    
                    return(tt)
                  })
           )
         }) })
  )
},.parallel=T)
tbl = Reduce(rbind,tbl)
save(tbl, file=paste0(root,"/pvals.Rdata")) # table
time_output(start1)

time_output(start)


