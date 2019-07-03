## input: features 
## output: p values and their corelation between train/test

set.seed(10)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr(c("stringr","colorspace", "Matrix", "plyr",
       "lattice", # libr(proxy)
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
          cmtf = !all(tem==tem[1])
          
          t = wilcox = list()
          t$test = t$train = t$train2 = 
            wilcox$test = wilcox$train = wilcox$train2 = 1
          if (!all(tem==tem[1]) & cmtf) {
            t$test = t.test(cm, tem)$p.value
            wilcox$test = wilcox.test(cm, tem)$p.value
          }
          if (!all(trm==trm[1]) & cmtf) {
            t$train = t.test(cm, trm)$p.value
            wilcox$train = wilcox.test(cm, trm)$p.value
            if (!is.null(trm2)) {
              t$train2 = median(sapply(trm2, function(trm2i) {
                if (!all(trm2i==trm2i[1])) return(1)
                return(t.test(cm, trm2i)$p.value)
              }))
              wilcox$train2 = median(sapply(trm2, function(trm2i) {
                if (!all(trm2i==trm2i[1])) return(1)
                return(wilcox.test(cm, trm2i)$p.value)
              }))
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
            
            title = paste0("-ln(train) vs -ln(test); class ",uc,", ", adj," adj ",ptype, " pvalues",
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
            
            overlap = sum(tr_sig & te_sig)
            rec = overlap/sum(te_sig)
            prec = overlap/sum(tr_sig)
            overlap2 = sum(tr2_sig & te_sig)
            rec2 = overlap2/sum(te_sig)
            prec2 = overlap2/sum(tr2_sig)
            foldres = 
              rbind(foldres, 
                    data.frame(feature=feat_type, class=uc, sig_test=ptype, adjustment=adj, p_thres=pthres, folds_in_train=ifelse(is.list(tri),1,length(tri)), samples_train=length(unlist(tri)), samples_test=length(tei), samples_control=sum(controln),
                    nsig_test=sum(te_sig), 
                    nsig_train=sum(tr_sig), nsig_overlap=overlap, pcorr=pcorr, pcorr_pvalue=pcorrp, nsig_recall=rec, nsig_precision=prec, f=2*((prec*rec)/(prec+rec2)),
                    nsig_train_fold=sum(tr2_sig), nsig_overlap_fold=overlap2, pcorr_fold=pcorr2, pcorr_pvalue_fold=pcorrp2, nsig_recall_fold=rec2, nsig_precision_fold=prec2, f_fold=2*((prec2*rec2)/(prec2+rec2))
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
  save(foldres0, file=paste0(pval_dir,"/result.Rdata"), row.names=F)
  
  time_output(start)
}


for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  load(paste0(pval_dir,"/result.Rdata"))
  foldres0$pmethod_adjust = paste(foldres0$sig_test, foldres0$adjustment)
  extras = ifelse(any(grepl("paired",foldres0$feature)),2,1)
  # png(paste0(pval_dir,"/result.png"), 
  #     width=length(unique(foldres0$class))*500,
  #     height=extras*3*2*500)
  # par(mfcol=c(extras*3*2,length(unique(foldres0$class))))

  for (uc in unique(foldres0$class)) {
    fr = foldres0[!grepl("paired",foldres0$feature),]
    for (i in 1:2) {
      if (i==2 & !extras) next
      if (i==2 & extras) 
        fr = foldres0[grepl("paired",foldres0$feature),]
      if (nrow(fr)==0) next
      try ({
        pl_r = barchart(nsig_recall~feature, data=fr ,groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)))
        png(paste0(pval_dir,"/result_class-",uc,"_recall",ifelse(i==2, "_paired",""),".png"))
        print(pl_r)
        graphics.off()
      })
      try ({
        pl_r_ = barchart(nsig_recall_fold~feature, data=fr ,groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)))
        png(paste0(pval_dir,"/result_class-",uc,"_recall",ifelse(i==2, "_paired",""),"_.png"))
        print(pl_r_)
        graphics.off()
      })
      try ({
        pl_p = barchart(nsig_precision~feature,data=fr, groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)))
        png(paste0(pval_dir,"/result_class-",uc,"_precision",ifelse(i==2, "_paired",""),".png"))
        print(pl_p)
        graphics.off()
      })
      try ({
        pl_p_ = barchart(nsig_precision_fold~feature,data=fr, groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)))
        png(paste0(pval_dir,"/result_class-",uc,"_precision",ifelse(i==2, "_paired",""),"_.png"))
        print(pl_p_)
        graphics.off()
      })
      try ({
        pl_f = barchart(f~feature,data=fr, groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)))
        png(paste0(pval_dir,"/result_class-",uc,"_f",ifelse(i==2, "_paired",""),".png"))
        print(pl_f)
        graphics.off()
      })
      try ({
        pl_f_ = barchart(f_fold~feature,data=fr, groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)))
        png(paste0(pval_dir,"/result_class-",uc,"_f",ifelse(i==2, "_paired",""),"_.png"))
        print(pl_f_)
        graphics.off()
      })
      try ({
        pl_c = barchart(pcorr~feature,data=fr, groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)))
        png(paste0(pval_dir,"/result_class-",uc,"_corr",ifelse(i==2, "_paired",""),".png"))
        print(pl_c)
        graphics.off()
      })
      try ({
        pl_c_ = barchart(pcorr_fold~feature,data=fr, groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)))
        png(paste0(pval_dir,"/result_class-",uc,"_corr",ifelse(i==2, "_paired",""),"_.png"))
        print(pl_c_)
        graphics.off()
      })
    }
  }
  
}

