## input: features 
## output: p values and their corelation between train/test

set.seed(10)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr(c("stringr","colorspace", "Matrix", "plyr",
       "vegan", # libr(proxy)
       "foreach","doMC",
       "kernlab"))



#Setup Cores
no_cores = 6 # detectCores()-5
registerDoMC(no_cores)



## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

# cvn-fold cross validation

overwrite = T #overwrite distances?
writecsv = F

readcsv = F

adjust = c("none","BY","BH","bonferroni") #pvalue adjustment
pthres = .05 # p value sig threshold for t test
good_count = 5
# minfold = 5 # minimum number of samples in each fold/class
good_sample = minfold = 5 # each class must have more thatn good_sample amount of samples or else prune
cvn0 = 10 # cvn-fold cross validation
testn = 1/5 #proportion of samples to make into test sample if none specified
# countThres = 1000 #insignificant if count under

id_col = "id"
target_col = "class"
order_cols = NULL

control = "control"


for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  print(result_dir)
  # result_dir = paste0(root, "/result/impc_panel1_sanger-spleen") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta") # meta files directory
  meta_file_dir = paste(meta_dir, "/file", sep="") #meta for rows (sample)
  meta_cell_dir = paste(meta_dir, "/cell", sep="") #meta for rows (sample)
  feat_dir = paste(result_dir, "/feat", sep="") #feature files directory
  
  
  ## output directories
  pval_dir = paste0(result_dir,"/pval"); suppressWarnings(dir.create (pval_dir))
  
  
  
  #data paths
  feat_count = "file-cell-countAdj"
  feat_types = list.files(path=feat_dir,pattern=".Rdata")
  # feat_types = feat_types[!feat_types=="file-cell-count.Rdata"]
  feat_types = gsub(".Rdata","",feat_types)
  feat_types = feat_types[!grepl("KO|Max",feat_types)]
  # feat_types = feat_types[!grepl("Freqp_orig",feat_types)]
  
  meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
  
  
  start = Sys.time()
  
  mc0 = Matrix(get(load(paste0(feat_dir,"/", feat_count,".Rdata"))))
  
  #load different features matrix and calculate distance matrix
  # for (feat_type in feat_types_) {
  foldres0 = ldply(feat_types, function(feat_type) {
    tryCatch({
      cat("\n", feat_type, " ",sep="")
      start2 = Sys.time()
      
      ## upload and prep feature matrix
      m0 = as.matrix(Matrix(get(load(paste0(feat_dir,"/", feat_type,".Rdata")))))
      # iscell = grepl("_",colnames(m0))[10]
      
      sm = meta_file[match(rownames(m0),meta_file$id),]
      # good_sample
      tcl = table(sm$class); 
      delrow = rownames(m0)%in%sm$id[sm$class%in%names(tcl)[tcl<=minfold]]
      sm = sm[!delrow,]
      m = m0[!delrow,]
      if (all(sm$class==sm$class[1])) return(F)
      # if (feat_type=="file-cell-lnpropexpect") 
      m[is.infinite(m)] = min(m[!is.infinite(m)])
      m[is.nan(m) | is.na(m)] = 0
      controli = sm$class=="control"
      controln = sum(controli)
      
      foldsip = list() # save original
      foldres = NULL # return as table
      for (uc in unique(sm$class[!controli])) {
        uci = sm$class==uc
        if (sum(uci)<minfold) next
        
        # test/train indices
        if ("type"%in%colnames(sm)) {
          foldsip[[uc]]$test = sm$id[sm$type=="test" & uci]
          foldsip[[uc]]$train = sm$id[sm$type=="train" & uci]
        } else {
          testni = sample(which(uci), max(minfold,testn*sum(uci)))
          foldsip[[uc]]$test = sm$id[testni]
          foldsip[[uc]]$train = sm$id[!c(1:nrow(sm))%in%testni & uci]
        }
        
        # n-fold cv train indices; if only two folds? then still train/test; number of folds rounded
        trind = sample(foldsip[[uc]]$train)
        if (length(trind)>(minfold*3)) {
          cvn = floor(min(cvn0,length(trind)/minfold))
          foldsip[[uc]]$train = split(trind, cut(seq_along(trind), cvn, labels=F))
        }
      }
      
      # calculate p value and plot; plot needs clean up, right now plots same number of plots per row, will want to fix so that plot only up until a nmber of folds per row so it's a unique class per row
      
      nh = 2*2*length(adjust) # number of p value methods
      nw = length(foldsip)
      png(paste0(pval_dir,"/",feat_type,".png"), 
          width=nw*500,height=nh*500)
      par(mfcol=c(nh,nw), mar=c(3,2,6,1))
      
      for (uc in names(foldsip)) {
        tri = foldsip[[uc]]$train
        tei = foldsip[[uc]]$test
        
        pvs = llply(1:ncol(m), function(j) {
          cm = m[controli,j]; 
          trm = m[unlist(tri),j]; 
          trm2 = NULL; if (is.list(tri)) trm2 = lapply(tri, function(trii) m[trii,j])
          tem = m[tei,j]
          
          t = wilcox = list()
          t$test = t$train = t$train2 = 
            wilcox$test = wilcox$train = wilcox$train2 = 1
          if (!all(tem==tem[1])) {
            t$test = t.test(cm, tem)$p.value
            wilcox$test = wilcox.test(cm, tem)$p.value
          }
          if (!all(trm==trm[1])) {
            t$train = t.test(cm, trm)$p.value
            wilcox$train = wilcox.test(cm, trm)$p.value
            if (!is.null(trm2)) {
              t$train2 = max(sapply(trm2, function(trm2i) 
                t.test(trm2i)$p.value))
              wilcox$train2 = max(sapply(trm2, function(trm2i) 
                wilcox.test(trm2i)$p.value))
            }
          }
          return(list(t=t, wilcox=wilcox))
        })
        
        for (ptype in names(pvs[[1]])) {
          pv_tr = sapply(pvs, function(x) x[[ptype]]$train)
          pv_tr2 = sapply(pvs, function(x) x[[ptype]]$train2)
          pv_te = sapply(pvs, function(x) x[[ptype]]$test)
          names(pv_tr) = names(pv_tr2) = names(pv_te) = colnames(m)
          # adjust p values
          for (adj in adjust) {
            pv_tr_ = p.adjust(pv_tr, method=adj)
            pv_tr2_ = p.adjust(pv_tr2, method=adj)
            pv_te_ = p.adjust(pv_te, method=adj)
            
            foldsip[[uc]][[ptype]][[adj]]$p_train = pv_tr_
            foldsip[[uc]][[ptype]][[adj]]$p_train2 = pv_tr2_
            foldsip[[uc]][[ptype]][[adj]]$p_test = pv_te_
            
            # calculate correlations between p values test & train/2
            pv_trl = -log(pv_tr_)
            pv_trl2 = -log(pv_tr2_)
            pv_tel = -log(pv_te_)
            
            foldsip[[uc]][[ptype]][[adj]]$pcorr = pcorr = 
              cor(pv_trl,pv_tel, method="pearson")
            foldsip[[uc]][[ptype]][[adj]]$pcorrp = pcorrp =
              cor.test(pv_trl,pv_tel, method="pearson")$p.value
            foldsip[[uc]][[ptype]][[adj]]$pcorr2 = pcorr2 = 
              cor(pv_trl2,pv_tel, method="pearson")
            foldsip[[uc]][[ptype]][[adj]]$pcorrp2 = pcorrp2 = 
              cor.test(pv_trl2,pv_tel, method="pearson")$p.value
            
            # plot
            phens = names(pv_trl)
            if(grepl("_",phens[10])) phens = sapply(str_split(names(pv_trl),"_"), function(x) x[1])
            
            mcm = colMeans(mc0[match( as.character(sm$id[sm$class=="control"]), rownames(mc0)), match(phens,colnames(mc0))])
            mcm_ = log(mcm); mcm_ = mcm_/max(mcm_) # for plot dot size
            
            # dcol = rep("black",length(pv_trl))
            tr_sig = pv_tr_<pthres
            tr2_sig = pv_tr2_<pthres
            te_sig = pv_te_<pthres
            # dcol[tr_sig] = "blue"
            # dcol[te_sig] = "red"
            # dcol[tr_sig & te_sig] = "purple"
            
            title = paste0("class ",uc,", ", adj," adj ",ptype, " pvalues",
                           "\nsize=ln(meancout)/max(ln(meancount)); sigs=",pthres, "; (ntrain / test / ctrl=",length(unlist(tri))," / ",length(tei)," / ",sum(controln),")")
            plot(pv_trl, pv_tel, pch=16, cex=mcm_, col=rgb(0,0,0,.5), #col=dcol,
                 xlab="-ln(train)", ylab="-ln(test)",
                 main=paste0(title, "\n pearson corr / p=",round(pcorr,3)," / ",round(pcorrp,3),
                             ";   sig train / test / both=",sum(tr_sig)," / ",sum(te_sig)," / ",sum(tr_sig & te_sig)))
            abline(h=-log(pthres), v=-log(pthres))
            
            plot(pv_trl2, pv_tel, pch=16, cex=mcm_, col=rgb(0,0,0,.5), #col=dcol,
                 xlab="-ln(train) max(train folds)", ylab="-ln(test)",
                 main=paste0(title, "\n pearson corr / p=",round(pcorr2,3)," / ",round(pcorrp2,3),
                             ";   sig train / test / both=",sum(tr2_sig)," / ",sum(te_sig)," / ",sum(tr2_sig & te_sig)))
            
            abline(h=-log(pthres), v=-log(pthres))
            
            
            foldres = 
              rbind(foldres, 
                    data.frame(feature=feat_type, class=uc, sig_test=ptype, adjustment=adj, p_thres=pthres, folds_in_train=ifelse(is.list(tri),1,length(tri)), samples_train=length(unlist(tri)), samples_test=length(tei), samples_control=sum(controln),
                    nsig_test=sum(te_sig), 
                    nsig_train=sum(tr_sig), nsig_overlap=sum(tr_sig & te_sig), pcorr=pcorr, pcorr_pvalue=pcorrp, nsig_recall=sum(tr_sig & te_sig)/sum(te_sig), nsig_precision=sum(tr_sig & te_sig)/sum(tr_sig),
                    nsig_train_fold=sum(tr2_sig), nsig_overlap_fold=sum(tr2_sig & te_sig), pcorr_fold=pcorr2, pcorr_pvalue_fold=pcorrp2, nsig_recall_fold=sum(tr2_sig & te_sig)/sum(te_sig), nsig_precision_fold=sum(tr2_sig & te_sig)/sum(tr2_sig)
              ))
          } # adjust
        } # ptype
      } # uc
      graphics.off()
      
      save(foldsip, file=paste0(pval_dir,"/",feat_type,".Rdata"))
      return(foldres)
      time_output(start2, feat_type)
    }, error = function(err) { cat(paste("ERROR:  ",err)); return(T) })
  }, .parallel=T)
  write.csv(foldres0, file=paste0(pval_dir,"/result.csv"), row.names=F)
  
  time_output(start)
}


