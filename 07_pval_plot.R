## input: features 
## output: p values and their corelation between train/test sets; p value qq plots, recall/precision, number of significant features

set.seed(10)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr(c("stringr","colorspace", "Matrix", "plyr",
       "lattice", "gridExtra", "plotly",# libr(proxy)
       "metap",
       "foreach","doMC",
       "kernlab"))

pvalr_dir = paste0(root,"/pval"); dir.create(pvalr_dir, showWarnings=F, recursive=T)

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


table = get(load(paste0(root,"/pval_table.Rdata")))
pvals = get(load(paste0(root,"/pval.Rdata")))

# data and uc together, ptype and adj together, 
c_datas = names(pvals)
c_feats = unique(table$feature)
# uc, p
c_ptypes = unique(table$sig_test)
c_adjs = unique(table$adjust.combine)

tbl = table
tbl$pmethod_adj = paste(table$sig_test, table$adjust.combine)
ptl = -log(table$p_thres[1])

plot_int = function(dat, pch='.', ...) {
  colPalette = colorRampPalette(c("blue", "green", "yellow", "red"))
  col = densCols(dat, colramp=colPalette)
  graphics::plot(dat, col=col, pch=pch, ...)
}


## histograms + scatterplots
for (data in c_datas) {
  png(paste0(pvalr_dir,"/hist_",data,".png"), height=length(pvals[[data]])*300, width=length(pvals[[data]][[1]])*length(c_ptypes)*2*300)
  mv1 = matrix(c(1,1,2,3),nrow=2,byrow=F)
  for (pi in 2:length(c_ptypes))
    mv1 = cbind(mv1, matrix(c(1,1,2,3)+max(mv1),nrow=2,byrow=F))
  mv = mv1
  for (fi in 1:length(pvals[[data]]))
    mv = rbind(mv, mv1+max(mv))
  layout(mv) # scatterplot + histograms
  
  mc0 = Matrix(get(load(paste0(root,"/result/",data,"/feat/", feat_count,".Rdata"))))
  mcm = colMeans(mc0)
  mcm = -log(mcm/max(mcm))
  mcm = 2*mcm/max(mcm)
  mcm[mcm<.5] = .5

  for (uc in names(pvals[[data]][[1]])) {
    for (feat in names(pvals[[data]])) {
      for (ptype in c_ptypes) {
        for (adj in c_adjs) {
          pvs = pvals[[data]][[feat]][[uc]]$p[[ptype]][[adj]]
          
          set1 = -log(pvs$train)
          set2 = -log(pvs$test)
          set1[set1>10] = set2[set2>10] = 10
          
          plot_int(cbind(set1,set2), xlim=c(0,10), ylim=c(0,10), xlab="set 1 -ln(p values)", ylab="set 2 -ln(p values)", main=paste0("set2 class: ",uc, "; feature: ", feat, "; method: ", ptype,"/",adj), cex=mcm[match(names(set1),names(mcm))])
          abline(v=ptl,h=ptl)
          
          hist(set1, freq=F, main=paste0("set 1"), xlab="-ln(p values)", ylab="histogram", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, xlim=c(0,10), ylim=c(0,1), col=rgb(0,0,0,.5), breaks=10)
          abline(v=ptl)
          
          hist(set2, freq=F, main="set 2", xlab="-ln(p values)", ylab="histogram", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, xlim=c(0,10), ylim=c(0,1), col=rgb(0,0,0,.5), breaks=10)
          abline(v=ptl)
        }
      }
    }
  }
  graphics.off()
}

## qq plot; log both axis
for (data in c_datas) {
  png(paste0(pvalr_dir,"/hist_",data,".png"), height=length(pvals[[data]])*300, width=length(pvals[[data]][[1]])*length(c_ptypes)*2*300)
  mv1 = matrix(c(1,1,2,3),nrow=2,byrow=F)
  for (pi in 2:length(c_ptypes))
    mv1 = cbind(mv1, matrix(c(1,1,2,3)+max(mv1),nrow=2,byrow=F))
  mv = mv1
  for (fi in 1:length(pvals[[data]]))
    mv = rbind(mv, mv1+max(mv))
  layout(mv) # scatterplot + histograms
  
  mc0 = Matrix(get(load(paste0(root,"/result/",data,"/feat/", feat_count,".Rdata"))))
  mcm = colMeans(mc0)
  mcm = -log(mcm/max(mcm))
  mcm = 2*mcm/max(mcm)
  mcm[mcm<.5] = .5
  
  for (uc in names(pvals[[data]][[1]])) {
    for (feat in names(pvals[[data]])) {
      for (ptype in c_ptypes) {
        for (adj in c_adjs) {
          pvs = pvals[[data]][[feat]][[uc]]$p[[ptype]][[adj]]
          
          set1 = -log(pvs$train)
          set2 = -log(pvs$test)
          set1[set1>10] = set2[set2>10] = 10
          
          plot_int(cbind(set1,set2), xlim=c(0,10), ylim=c(0,10), xlab="set 1 -ln(p values)", ylab="set 2 -ln(p values)", main=paste0("set2 class: ",uc, "; feature: ", feat, "; method: ", ptype,"/",adj), cex=mcm[match(names(set1),names(mcm))])
          abline(v=ptl,h=ptl)
          
          hist(set1, freq=F, main=paste0("set 1"), xlab="-ln(p values)", ylab="histogram", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, xlim=c(0,10), ylim=c(0,1), col=rgb(0,0,0,.5), breaks=10)
          abline(v=ptl)
          
          hist(set2, freq=F, main="set 2", xlab="-ln(p values)", ylab="histogram", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, xlim=c(0,10), ylim=c(0,1), col=rgb(0,0,0,.5), breaks=10)
          abline(v=ptl)
        }
      }
    }
  }
  graphics.off()
}



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


## number of significant features

## correlation scatterplot; histogram of p values

## other metrics

mc0 = Matrix(get(load(paste0(feat_dir,"/", feat_count,".Rdata"))))

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


