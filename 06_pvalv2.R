## input: features 
## output: p values and their corelation between train/test sets; p value qq plots, recall/precision, number of significant features

set.seed(10)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr(c("stringr","colorspace", "Matrix", "plyr",
       "lattice", "gridExtra",# libr(proxy)
       "metap",
       "foreach","doMC",
       "kernlab"))



#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)



## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

# cvn-fold cross validation

overwrite = F #overwrite?
writecsv = F

readcsv = F

adjust = c("BY")#,"fisher","none","BH","bonferroni") #pvalue adjustment; "lanc",
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
feat_count = "file-cell-countAdj"

## calcuate p values!

# format of p values:
# tri = foldsip[[uc]]$indices$train
# tei = foldsip[[uc]]$indices$test
# pv_tr_ = foldsip[[uc]]$p[[ptype]][[adj]]$train
# pv_tr2_ = foldsip[[uc]]$p[[ptype]][[adj]]$p_train2
# pv_te_ = foldsip[[uc]]$p[[ptype]][[adj]]$p_test
# pv_all_ = foldsip[[uc]]$p[[ptype]][[adj]]$p_all

start = Sys.time()
table = pvals = NULL
for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  print(result_dir)
  data = fileNames(result_dir)
  
  ## input directories
  meta_file_dir = paste(result_dir, "/meta/file", sep="") #meta for rows (sample)
  feat_dir = paste(result_dir, "/feat", sep="") #feature files directory
  
  ## output directories
  unlink(paste0(result_dir,"/pval"))
  pvalsource_dir = paste0(result_dir,"/pval/source"); suppressWarnings(dir.create (pvalsource_dir, recursive=T))

  #data paths
  feat_types = list.files(path=feat_dir,pattern=".Rdata")
  # feat_types = feat_types[!feat_types=="file-cell-count.Rdata"]
  feat_types = gsub(".Rdata","",feat_types)
  # feat_types = feat_types[!grepl("KO|Max",feat_types)]

  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
  
  
  start1 = Sys.time()
  
  #load different features matrix and calculate distance matrix
  # for (feat_type in feat_types_) {
  result = ldply(feat_types, function(feat_type) #for (feat_type in feat_types) 
  {
    start2 = Sys.time()
    
    # tryCatch({
    if (!file.exists(paste0(pvalsource_dir,"/",feat_type,".Rdata")) | overwrite) {
      cat("\n", feat_type, " ",sep="")
      
      ## upload and prep feature matrix + meta
      m0 = as.matrix(Matrix(get(load(paste0(feat_dir,"/", feat_type,".Rdata")))))
      sm = meta_file[match(rownames(m0),meta_file$id),]
      # good_sample
      tcl = table(sm$class); 
      delrow =rownames(m0)%in%sm$id[sm$class%in%names(tcl)[tcl<=minfold]]
      sm = sm[!delrow,]
      if (all(sm$class==sm$class[1])) return(F)
      m = m0[!delrow,]
      
      m[is.infinite(m)] = min(m[!is.infinite(m)])
      m[is.nan(m) | is.na(m)] = 0
      
      controli = sm$class=="control"
      controln = sum(controli)
      
      
      foldsip = foldres = NULL # save original
      for (uc in unique(sm$class[!controli])) {
        
        ## make test/train(x10) indices
        uci = sm$class==uc
        if (sum(uci)<minfold) next
        
        if ("type"%in%colnames(sm)) {
          foldsip[[uc]]$indices$test = as.character(sm$id[sm$type=="test" & uci])
          foldsip[[uc]]$indices$train = as.character(sm$id[sm$type=="train" & uci])
        } else {
          testni = sample(which(uci), max(minfold,testn*sum(uci)))
          foldsip[[uc]]$indices$test = as.character(sm$id[testni])
          foldsip[[uc]]$indices$train = as.character(sm$id[!c(1:nrow(sm))%in%testni & uci])
        }
        
        # n-fold cv train indices; if only two folds? then still train/test; number of folds rounded
        trind = sample(foldsip[[uc]]$indices$train)
        if (length(trind)>(minfold*3)) {
          cvn = floor(min(cvn0,length(trind)/minfold))
          foldsip[[uc]]$indices$train = split(trind, cut(seq_along(trind), cvn, labels=F))
        }
        
        tri = foldsip[[uc]]$indices$train
        tei = foldsip[[uc]]$indices$test
        
        ## calculate p value
        pvs = llply(1:ncol(m), function(j) {
          # for (j in 1:ncol(m)) {
          cm = m[controli,j]; 
          cmm = m[!controli,j]; 
          trm = m[unlist(tri),j]; 
          trm2 = NULL; if (is.list(tri)) trm2 = lapply(tri, function(trii) m[trii,j])
          tem = m[tei,j]
          cmtf = !all(cm==cm[1])
          
          t = wilcox = list()
          t$test = t$train = t$train2 = t$all = 
            wilcox$test = wilcox$train = wilcox$train2 = wilcox$all = 1
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
          if (!all(cmm==cmm[1]) & cmtf) {
            t$all = t.test(cm, cmm)$p.value
            wilcox$all = wilcox.test(cm, cmm)$p.value
          }
          return(list(t=t, wilcox=wilcox))
        })
        
        for (ptype in names(pvs[[1]])) {
          for (adj in adjust) {
            for (tretype in names(pvs[[1]][[ptype]])) {
              pv = sapply(pvs, function(x) x[[ptype]]$train); 
              pv[is.nan(pv) | is.na(pv)] = 1
              names(pv) = colnames(m)
              
              ## adjust or combo p values
              if (adj%in%c("lanc","fisher")) {
                # if (grepl("group",feat_type)) next()
                
                phens = names(pv)
                if(grepl("_",phens[10])) phens = sapply(str_split(names(pv),"_"), function(x) x[1])
                nonp = gsub("[-]|[+]","",phens)
                allplus = which(grepl("[-]|[+]",phens) & !duplicated(nonp))
                if (length(allplus)==length(phens)) next
                # nonpu = unique(nonp)
                groupi = match(nonp, nonp)
                groups = llply(allplus, function(i) {
                  ii = which(groupi==groupi[i])
                  if (length(ii)>1) return(ii)
                  return(NULL)
                })
                names(groups) = gsub("-","+",phens[allplus])
                groups = plyr::compact(groups)
                
                pv_ = a = rep(1,length(groups)); names(a) = names(groups)
                ap = sapply(groups, function(x) any(pv[x]<1))
                
                # don't use, not sure degrees of freedom...
                if (adj=="lanc") {
                  if (any(ap)) pv_[ap] = 
                      laply(groups[ap],function(x) invchisq(pv[x],2)$p)
                }
                if (adj=="fisher") {
                  if (any(ap)) pv_[ap] = 
                      laply(groups[ap],function(x) sumlog(pv_[x])$p)
                }
                names(pv_) = names(groups)
              } else {
                pv_ = p.adjust(pv, method=adj)
              }
              
              foldsip[[uc]]$p[[ptype]][[adj]][[tretype]] = pv_
            } # tretype
            
            pv_tr_ = foldsip[[uc]]$p[[ptype]][[adj]]$train
            pv_tr2_ = foldsip[[uc]]$p[[ptype]][[adj]]$p_train2
            pv_te_ = foldsip[[uc]]$p[[ptype]][[adj]]$p_test
            pv_all_ = foldsip[[uc]]$p[[ptype]][[adj]]$p_all
            
            ## calculate correlations between p values test & train/2
            pv_trl = -log(pv_tr_)
            pv_trl2 = -log(pv_tr2_)
            pv_tel = -log(pv_te_)
            # pv_alll
            pcorr = cor(pv_trl,pv_tel, method="spearman")
            pcorrp = cor.test(pv_trl,pv_tel, method="spearman")$p.value
            pcorr2 = cor(pv_trl2,pv_tel, method="spearman")
            pcorrp2 = cor.test(pv_trl2,pv_tel, method="spearman")$p.value
            
            
            # save table
            tr_sig = pv_tr_<pthres
            tr2_sig = pv_tr2_<pthres
            te_sig = pv_te_<pthres
            all_sig = pv_all_<pthres
            
            overlap = sum(tr_sig & te_sig)
            rec = overlap/sum(te_sig)
            prec = overlap/sum(tr_sig)
            overlap2 = sum(tr2_sig & te_sig)
            rec2 = overlap2/sum(te_sig)
            prec2 = overlap2/sum(tr2_sig)
            foldres = 
              rbind(
                foldres, data.frame(
                  data=data, feature=feat_type, class=uc, sig_test=ptype, adjust.combine=adj, p_thres=pthres, 
                  n_train_folds=ifelse(is.list(tri),1,length(tri)), n_samples_train=length(unlist(tri)), n_samples_test=length(tei), n_samples_control=sum(controln),
                  
                  m_test_sig=sum(te_sig), m_test=length(te_sig), 
                  m_test_perc=sum(te_sig)/length(te_sig),
                  m_train_sig=sum(tr_sig), m_train_sig2=sum(tr2_sig), m_train=length(tr_sig), 
                  m_train_perc=sum(tr_sig)/length(tr_sig), m_train_perc2=sum(tr2_sig)/length(tr2_sig),
                  m_overlap_sig=overlap, m_overlap_sig2=overlap2, 
                  
                  corr_spear=pcorr, corr_spear_p=pcorrp,
                  corr_spear2=pcorr2, corr_spear_p2=pcorrp2,
                  recall=rec, precision=prec, f=2*((prec*rec)/(prec+rec)),
                  recall2=rec2, precision2=prec2, f2=2*((prec2*rec2)/(prec2+rec2))
              ))
          } # adj
        } # ptype
      }
      save(foldsip, file=paste0(pvalsource_dir,"/",feat_type,".Rdata"))
      save(foldres, file=paste0(pvalsource_dir,"/",feat_type,"_table.Rdata"))
    } # if
    foldsip = get(load(paste0(pvalsource_dir,"/",feat_type,".Rdata")))
    foldres = get(load(paste0(pvalsource_dir,"/",feat_type,"_table.Rdata")))
    a = list(foldsip=foldsip, foldres=foldres)
    
    time_output(start2)
    return(a)
  }, .parallel=T)
  pvals[[data]] = llply(result, function(x) x$foldsip)
  names(pvals[[data]]) = feat_types
  table = rbind(table, ldply(result, function(x) x$foldres))
  
  time_output(start1)
} # result

mc0 = Matrix(get(load(paste0(feat_dir,"/", feat_count,".Rdata"))))


    ## build table
    foldres = NULL
    for (uc in names(foldsip)) {
      tri = pvals[[data]][[feat_type]][[uc]]$indices$train
      tei = pvals[[data]][[feat_type]][[uc]]$indices$test
      
      pan = length(foldsip[[uc]]$p) + length(foldsip[[uc]]$p[[]])
      png(paste0(pval_dir,"/",feat_type,"_",uc,".png"), 
          width=pan*500,height=2*500)
      par(mfcol=c(2,pan), mar=c(3,3,6,1))
      
      for (ptype in names(foldsip[[uc]]$p)) {
        for (adj in names(foldsip[[uc]][[ptype]])) {
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
          
          
          if (!all(pv_tel==pv_tel[1])) {
            
            title = paste0("-ln(train) vs -ln(test); class ",uc,", ", adj," adj ",ptype, " pvalues",
                           "\nsize=ln(meancout)/max(ln(meancount)); sigs=",pthres, "
                           \n(ntrain/test/ctrl = ",length(unlist(tri)),"/",length(tei),"/",sum(controln),")")
            if (!all(pv_trl==pv_trl[1])) {
              plot(pv_trl, pv_tel, pch=16, cex=mcm_, col=rgb(0,0,0,.5), #col=dcol,
                   xlab="-ln(train)", ylab="-ln(test)",
                   main=paste0("spearnman corr/p = ",round(pcorr,3),"/",round(pcorrp,3),
                               ";   sig train/test/both = ",sum(tr_sig),"/",sum(te_sig),"/",sum(tr_sig & te_sig)))
              abline(h=-log(pthres), v=-log(pthres))
            } else {
              next
            }
            
            if (!all(pv_trl==pv_trl[1])) {
              plot(pv_trl2, pv_tel, pch=16, cex=mcm_, col=rgb(0,0,0,.5), #col=dcol,
                   xlab="-ln(train) max(train folds)", ylab="-ln(test)",
                   main=paste0(title, "\n spearnman corr / p=",round(pcorr2,3)," / ",round(pcorrp2,3),
                               ";   sig train / test / both=",sum(tr2_sig)," / ",sum(te_sig)," / ",sum(tr2_sig & te_sig)))
              abline(h=-log(pthres), v=-log(pthres))
            } else {
              next
            }
          } else {
            next
          }
          
        } # adjust
      } # ptype
      graphics.off()
    } # uc
    # graphics.off()
    
    time_output(start2, feat_type)
    return(foldres)
    #}, error = function(err) { cat(paste("ERROR:  ",err)); return(T) })
  }, .parallel=T)
  save(foldres0, file=paste0(pval_dir,"/table.Rdata"))
  
  # combine results
  foldps = llply(feat_types, function(x) {
    fna = paste0(pval_dir,"/",x,".Rdata")
    a = NULL; if (file.exists(fna)) {
      a = get(load(fna))
    }
    return(a)
  })
  names(foldps) = feat_types
  foldps = plyr::compact(foldps)
  save(foldps, file=paste0(pval_dir,"/result.Rdata"))
  
  time_output(start1)
}
time_output(start)


## write table of results
table = ldply (list.dirs(paste0(root, "/result"), full.names=T, recursive=F), function(result_dir) {
  foldres0 = get(load(paste0(result_dir, "/pval/table.Rdata")))
  cbind(data=fileNames(result_dir), foldres0)
})
write.csv(table, paste0(root,"/pval.csv"))


## plot result table stuff
start = Sys.time()
tbl0 = table
pval_dir = paste0(result_dir,"/pval")

width = 700
tbl0$pmethod_adjust = paste(tbl0$sig_test, tbl0$adjust.combine)
tbl0 = tbl0[tbl0$feature!="file-cell-count",]
extras = ifelse(any(grepl("paired",tbl0$feature)),2,1)
# png(paste0(pval_dir,"/result.png"), 
#     width=length(unique(tbl$class))*500,
#     height=extras*3*2*500)
# par(mfcol=c(extras*3*2,length(unique(tbl$class))))
pl_r = pl_r_ = pl_p = pl_p_ = pl_f = pl_f_ = pl_c = pl_c_ = pl_n = pl_n_ = list()
dss = unique(tbl0$data)
dss = c(dss[grepl("ctrl",dss)],dss[grepl("pos",dss)],dss[!grepl("pos|ctrl",dss)])
for (ds in dss) {
  dsp = get(load(paste0(root,"/result/",ds,"/pval/result.Rdata")))
  tbl = tbl0[tbl0$data==ds,]
  for (uc in unique(tbl$class)) {
    fr = tbl[!grepl("paired",tbl$feature),]
    for (i in 1:2) {
      if (i==2 & !extras) next
      if (i==2 & extras) 
        fr = tbl[grepl("paired",tbl$feature),]
      if (nrow(fr)==0) next
      ploti = paste0(ds,"_",uc,ifelse(i==2,"_paired",""))
      plotnom = paste0("data: ",ds,ifelse(i==2,"_paired",""),"; class: ",uc)
      
      try ({
        pl_r[[ploti]] = barchart(nsig_recall~feature, data=fr ,groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)), main=plotnom)
        png(paste0(root,"/pval_",ds,"_class-",uc,"_recall",ifelse(i==2, "_paired",""),".png"), width=width)
        print(pl_r[[ploti]])
        graphics.off()
      })
      try ({
        pl_r_[[ploti]] = barchart(nsig_recall_fold~feature, data=fr ,groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)), main=plotnom)
        png(paste0(root,"/pval_",ds,"_class-",uc,"_recall",ifelse(i==2, "_paired",""),"_.png"), width=width)
        print(pl_r_[[ploti]])
        graphics.off()
      })
      try ({
        pl_p[[ploti]] = barchart(nsig_precision~feature,data=fr, groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)), main=plotnom)
        png(paste0(root,"/pval_",ds,"_class-",uc,"_precision",ifelse(i==2, "_paired",""),".png"), width=width)
        print(pl_p[[ploti]])
        graphics.off()
      })
      try ({
        pl_p_[[ploti]] = barchart(nsig_precision_fold~feature,data=fr, groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)), main=plotnom)
        png(paste0(root,"/pval_",ds,"_class-",uc,"_precision",ifelse(i==2, "_paired",""),"_.png"), width=width)
        print(pl_p_[[ploti]])
        graphics.off()
      })
      try ({
        pl_f[[ploti]] = barchart(f~feature,data=fr, groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)), main=plotnom)
        png(paste0(root,"/pval_",ds,"_class-",uc,"_f",ifelse(i==2, "_paired",""),".png"), width=width)
        print(pl_f[[ploti]])
        graphics.off()
      })
      try ({
        pl_f_[[ploti]] = barchart(f_fold~feature,data=fr, groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)), main=plotnom)
        png(paste0(root,"/pval_",ds,"_class-",uc,"_f",ifelse(i==2, "_paired",""),"_.png"), width=width)
        print(pl_f_[[ploti]])
        graphics.off()
      })
      try ({
        pl_c[[ploti]] = barchart(pcorr~feature,data=fr, groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)), main=plotnom)
        png(paste0(root,"/pval_",ds,"_class-",uc,"_corr",ifelse(i==2, "_paired",""),".png"), width=width)
        print(pl_c[[ploti]])
        graphics.off()
      })
      try ({
        pl_c_[[ploti]] = barchart(pcorr_fold~feature,data=fr, groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)), main=plotnom)
        png(paste0(root,"/pval_",ds,"_class-",uc,"_corr",ifelse(i==2, "_paired",""),"_.png"), width=width)
        print(pl_c_[[ploti]])
        graphics.off()
      })
      try ({
        pl_n[[ploti]] = barchart(nsig_train~feature,data=fr, groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)), main=paste0(plotnom, " ::: total feats: ", fr$nsigall_train))
        png(paste0(root,"/pval_",ds,"_class-",uc,"_nsig",ifelse(i==2, "_paired",""),".png"), width=width)
        print(pl_c[[ploti]])
        graphics.off()
      })
      try ({
        pl_n_[[ploti]] = barchart(nsig_train_fold~feature,data=fr, groups=pmethod_adjust, auto.key = list(columns=2), cex.axis=3, las=2, scales=list(x=list(rot=90,cex=0.8)), main=paste0(plotnom, " ::: total feats: ", fr$nsigall_train_fold))
        png(paste0(root,"/pval_",ds,"_class-",uc,"_nsig_",ifelse(i==2, "_paired",""),"_.png"), width=width)
        print(pl_c_[[ploti]])
        graphics.off()
      })
      
    }
  }
}

g = do.call(grid.arrange, pl_r)
ggsave(file=paste0(root,"/pval_recall.png"), g) #saves g
g = do.call(grid.arrange, pl_r_)
ggsave(file=paste0(root,"/pval_recall_.png"), g) #saves g
g = do.call(grid.arrange, pl_p)
ggsave(file=paste0(root,"/pval_prec.png"), g) #saves g
g = do.call(grid.arrange, pl_p_)
ggsave(file=paste0(root,"/pval_prec_.png"), g) #saves g
g = do.call(grid.arrange, pl_f)
ggsave(file=paste0(root,"/pval_f.png"), g) #saves g
g = do.call(grid.arrange, pl_f_)
ggsave(file=paste0(root,"/pval_f_.png"), g) #saves g
g = do.call(grid.arrange, pl_c)
ggsave(file=paste0(root,"/pval_corr.png"), g) #saves g
g = do.call(grid.arrange, pl_c_)
ggsave(file=paste0(root,"/pval_corr_.png"), g) #saves g
g = do.call(grid.arrange, pl_n)
ggsave(file=paste0(root,"/pval_nsig.png"), g) #saves g
g = do.call(grid.arrange, pl_n_)
ggsave(file=paste0(root,"/pval_nsig_.png"), g) #saves g

time_output(start)


## qq plots
start = Sys.time()
for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  pval_dir = paste0(result_dir,"/pval")
  load(paste0(pval_dir,"/result.Rdata"))
  features = names(foldps)
  ucs = names(foldps[[1]])
  ptypes = names(foldps[[1]][[1]]); ptypes = ptypes[!ptypes%in%c("train","test")]
  adjs = names(foldps[[1]][[1]][[3]])
  for (uc in ucs) {
    
    png(paste0(pval_dir,"/qq_",uc,".png"),width=length(adjs)*530, height=length(ptypes)*400)
    par(mfrow=c(length(ptypes),length(adjs)), mar=c(5.1, 4.1, 4.1, 13), xpd=TRUE) # Add extra space to right of plot area
    for (ptype in ptypes) {
      for (adj in adjs) {
        
        for (fi in 1:length(features)) {
          if (is.null(foldps[[fi]][[uc]][[ptype]][[adj]])) next
          p_all = foldps[[fi]][[uc]][[ptype]][[adj]]$p_all
          qqnorm(p_all, col=fi, ylim=c(0,1), xlim=c(-4,3), main=ifelse(fi==1,paste0("p value quantiles; ", ptype, " ", adj),""), pch=16,cex=.3)
          if (fi<length(features)) par(new=T)
        }
        qqline(c(0,1))
        legend("topright", inset=c(-.6,0), legend=features, col=1:length(features))
      }
    }
    graphics.off()
  }
}
time_output(start)
