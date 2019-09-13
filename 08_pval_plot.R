## input: features 
## output: p values and their corelation between train/test sets; p value qq plots, recall/precision, number of significant features

set.seed(10)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr(c("stringr","Matrix", "plyr",
       "lattice", "gridExtra", "plotly", "RColorBrewer", "plotrix", "ggrepel", "ggplot2", "colorspace", # libr(proxy)
       "metap", "igraph",
       "arules", "arulesViz",
       "foreach","doMC"))

pvalr_dir = paste0(root,"/pval"); dir.create(pvalr_dir, showWarnings=F, recursive=T)
graph_dir = paste0(root,"/pval_graphs.Rdata")


#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)



## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

breaks = 30# colour breaks

overwrite = F #overwrite?
writecsv = F

feat_count = "file-cell-countAdj"
overlapmin = .5 # how much <overlap there must be before a cell population is considered significant (for graphs)

start = Sys.time()


## load table
table = get(load(paste0(root,"/pvals.Rdata")))

# use tbl (trim)
# tbl = table[grepl("pos7|pos8|pos9",tbl$data),]
tbl = tbl_ = table[!grepl("raw|edge|group|entropy|short",table$feat),]
tbl_$pmethod_adj = tbl$pmethod_adj = paste(tbl$test, tbl$test_adj, sep="-")
tblc = tbl[grepl("^ctrl",tbl$data),]
tbl = tbl[!grepl("ctrl[1-9]",tbl$data),]

tbl_ = tbl_[!grepl("ctrl[1-9]",tbl_$data),]


# pars
c_datas = unique(tbl$data)
c_feats = unique(tbl$feat)
c_pas = unique(tbl$pmethod_adj)
c_pts = unique(tbl$pthres)



## default counts/props
mcms = llply(c_datas, function(data) {
  m = get(load(paste0(root,"/result/",data,"/feat_mean/file-cell-prop.Rdata")))
  llply(m, function(mc) {
    mcm = mc - min(mc)
    mcm/max(mcm)
  })
  
  # sm =   get(load(paste0(root,"/result/",data,"/meta/", "file",".Rdata")))
  # sm = sm[match(rownames(m),sm$id),]
  # mc0 = llply(unique(sm$class), function(uc) {
  #   a = apply(m[sm$class=="control",],2,mean_geo)
  #   names(a) = colnames(m)
  #   # return(a)
  #   mc = a
  #   mc = mc0$control
  #   mcm = mc - min(mc)
  #   mcm/max(mcm)
  # })
  # names(mc0) = unique(sm$class)
  # return(mc0)
}, .parallel=T)
names(mcms) = c_datas


## default all feats
mfms = llply(c_datas, function(data) {
  fm_dir = paste0(root,"/result/",data,"/feat_mean")
  fm_dirs = gsub(".Rdata","",list.files(fm_dir, full.names=F, pattern=".Rdata"))
  mf = llply(fm_dirs, function(feat) get(load(paste0(fm_dir,"/",feat,".Rdata"))) )
  names(mf) = fm_dirs
  mf
}, .parallel=T)
names(mfms) = c_datas



## default graphs
gr0s = llply(c_datas, function(data) {
  get(load(paste0(root,"/result/",data,"/meta/cell_graph.Rdata")))
})
names(gr0s) = c_datas
# for (i in length(gr0s):1) if (nrow(gr0s[[i]]$v)>2500) gr0s[[i]] = NULL
grp0s = llply(names(gr0s), function(data) {
  get(load(paste0(root,"/result/",data,"/meta/cell_graphpos.Rdata"))) 
})
names(grp0s) = names(gr0s)


# density plot
plot_int = function(dat, pch='.', ...) {
  colPalette = colorRampPalette(c("blue", "green", "yellow", "red"))
  col = densCols(dat, colramp=colPalette)
  graphics::plot(dat, col=col, pch=pch, ...)
}


# make folders
for (i in 1:nrow(tbl)) 
  dir.create(paste0(root,"/result/",tbl$data[i],"/plot_pval/",tbl$pmethod_adj[i],"/pthres-",tbl$pthres[i]),showWarnings=F,recursive=T)





## histograms p ----------------------------------
start1 = Sys.time()
l_ply(loopInd(which(tbl$pthres==.05), no_cores), function(ii) {
  for (i in ii) {
    pt = tbl$pthres[i]; ptl = -log(pt)
    classn = ifelse(tbl$class[i]=="exp", "", paste0("_",tbl$class[i]))
    pathn = paste0(root,"/result/",tbl$data[i],"/plot_pval/",tbl$pmethod_adj[i])
    pv = get(load(tbl$pathp[i]))
    
    
    png(paste0(pathn,"/hist_p",classn,"_",tbl$feat[i],".png"), height=1000, width=1000)
    mv1 = matrix(c(1,1,2,3,4,4,5,6),ncol=2,byrow=F)
    layout(mv1) # scatterplot + histograms
    par(mar=c(7,7,3,3))
    
    mcm = mcms[[tbl$data[i]]]$control*4
    
    set1 = pv$train
    set2 = pv$test
    
    if (grepl("-cell-",tbl$feat[i])) 
      mcmind = match(names(set1),names(mcm))
    if (!grepl("-cell-",tbl$feat[i])) 
      mcmind = match(sapply(strsplit(names(set1),"_"), function(x) x[1]),names(mcm))
    mcm1 = mcm[mcmind]
    
    plot_int(cbind(set1,set2), pch=1, xlim=c(0,1), ylim=c(0,1), 
             xlab="set 1 p values", ylab="set 2 p values", main="", cex=mcm1, cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5)
    abline(v=pt,h=pt, col="red")
    
    hist(set1, freq=F, main=paste0("set 1"), xlim=c(0,1), breaks=seq(0,1,.05), xlab="p values", ylab="histogram", cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5, col=rgb(0,0,0,.5))
    abline(v=pt, col="red")
    
    hist(set2, freq=F, main="set 2", xlim=c(0,1), breaks=seq(0,1,.05), xlab="p values", ylab="histogram", cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5, col=rgb(0,0,0,.5))
    abline(v=pt, col="red")
    
    
    set1 = -log(set1); 
    set1[set1==Inf] = 10
    set2 = -log(set2); 
    set2[set2==Inf] = 10
    plot_int(cbind(set1,set2), pch=1, xlim=c(0,10), ylim=c(0,10), 
             xlab="set 1 -ln(p values)", ylab="set 2 -ln(p values)", main="", cex=mcm1, cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5)
    abline(v=ptl,h=ptl, col="red")
    
    hist(set1, freq=F, main=paste0("set 1"), xlim=c(0,max(4,max(set1))), xlab="-ln(p values)", ylab="histogram", cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5, col=rgb(0,0,0,.5))
    abline(v=ptl, col="red")
    
    hist(set2, freq=F, main="set 2", xlim=c(0,max(4,max(set2))), xlab="-ln(p values)", ylab="histogram", cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5, col=rgb(0,0,0,.5))
    abline(v=ptl, col="red")
    
    graphics.off()
  }
},.parallel=T)
time_output(start1)



## histograms o ----------------------------------
start1 = Sys.time()
l_ply(loopInd(which(tbl$pthres==c_pts[1] & tbl$pmethod_adj==c_pas[1]), no_cores), function(ii) {
  for (i in ii) {
    pt = tbl$pthres[i]; ptl = -log(pt)
    classn = ifelse(tbl$class[i]=="exp", "", paste0("_",tbl$class[i]))
    pathn = paste0(root,"/result/",tbl$data[i],"/plot_pval")
    dir.create(pathn,showWarnings=F,recursive=T)
    mo = get(load(tbl$patho[i]))
    
    
    png(paste0(pathn,"/hist_o",classn,"_",tbl$feat[i],".png"), height=500, width=500)
    mv1 = matrix(c(1,1,2,3),ncol=1,byrow=F)
    layout(mv1) # scatterplot + histograms
    par(mar=c(7,7,3,3))
    
    mcm = mcms[[tbl$data[i]]]$control*4
    
    set1 = mo$train
    set2 = mo$test
    
    if (grepl("-cell-",tbl$feat[i])) 
      mcmind = match(names(set1),names(mcm))
    if (!grepl("-cell-",tbl$feat[i])) 
      mcmind = match(sapply(strsplit(names(set1),"_"), function(x) x[1]),names(mcm))
    mcm1 = mcm[mcmind]
    
    plot_int(cbind(set1,set2), pch=1, xlim=c(0,1), ylim=c(0,1), 
             xlab="set 1  overlaps/control", ylab="set 2 overlaps/control", main="", cex=mcm1, cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5)
    
    hist(set1, freq=F, main=paste0("set 1"), xlim=c(0,1), breaks=seq(0,1,.05), xlab="p values", ylab="histogram", cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5, col=rgb(0,0,0,.5))
    
    hist(set2, freq=F, main="set 2", xlim=c(0,1), breaks=seq(0,1,.05), xlab="p values", ylab="histogram", cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5, col=rgb(0,0,0,.5))
    
    graphics.off()
  }
},.parallel=T)
time_output(start1)




## control qq plot -------------------------------
classn = ""
mcm = mcms[["ctrl0"]]$control
for (pa in c_pas) {
  tb = tblc[tblc$pmethod_adj==pa & tblc$pthres==.05,]
  
  for (trtype in c("train","test","all")) {
    
  pvs = llply(unique(tb$feat), function(x) {
    unlist(llply(tb$pathp[tb$feat==x], function(y) get(load(y))[[trtype]]))
  })
  names(pvs) = unique(tb$feat)
  
  pathn = paste0(root,"/result/ctrl0/plot_pval/",pa,"/qqctrl.png")
  
  
  
  ## qq plot; log both axis -----------------------
  nqs = llply(sapply(pvs,length), function(n) 1:n/(n+1))
  cs = rainbow(length(pvs))
  
  mcm1l = llply(pvs, function(x) {
    if (grepl("_",names(x)[1])) mcmind = match(names(x),names(mcm))
    if (!grepl("_",names(x)[1]))
      mcmind = match(sapply(strsplit(names(x),"_"), function(x) x[1]),names(mcm))
    mcm[mcmind]
  })
  
  names(pvs) = names(nqs) = names(cs) = unique(tb$feat) # feats
  
  png(pathn, height=500, width=1000)
  par(mfrow=c(1,2),mar=c(7,7,3,3))
  
  plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab=paste0(pa, " p values"), xlab="theoretical quantile", main="", cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5)
  for (i in 1:length(pvs)) {
    yo = order(pvs[[i]])
    points(nqs[[i]], pvs[[i]][yo], pch=1, cex=mcm1l[[i]][yo]*4, col=adjustcolor(cs[i], alpha.f=.75))
  }
  abline(h=.05)
  # abline(h=.025)
  # abline(h=.01)
  lines(x=c(-100,100), y=c(-100,100))
  legend("topleft", legend=names(pvs), fill=cs, bg="transparent",pt.cex=1.5)
  
  nqls = llply(nqs, log)
  pvls = llply(pvs, log)
  plot(NULL, xlim=c(min(sapply(nqls,function(x)x[1])),0), ylim=c(min(sapply(nqls,function(x)x[1])),0), ylab=paste0(pa, " ln(p values)"), xlab="theoretical quantile", main="", cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5)
  for (i in 1:length(pvs)) {
    yo = order(pvls[[i]])
    points(nqls[[i]], pvls[[i]][yo], pch=1, cex=mcm1l[[i]][yo]*4, col=adjustcolor(cs[i], alpha.f=.75))
  }
  abline(h=log(.05))
  # abline(h=log(.025))
  # abline(h=log(.01))
  lines(x=c(-100,100), y=c(-100,100))
  legend("topleft", legend=names(pvs), fill=cs, bg="transparent",pt.cex=1.5)
  
  graphics.off()
  }
}




## feature comparison plots ------------------------------
start1 = Sys.time()
for (data in c_datas) {
  mcm0 = mcms[[data]]
  ucs = unique(tbl$class[tbl$data==data])
  for (uc in ucs) {
    classn = ifelse(length(ucs)==1, "", paste0("_",uc))
    mcm = mcm0[[uc]]
    for (pa in c_pas) {
      tbinds = tbl$data==data & tbl$class==uc & tbl$pmethod_adj==pa
      tb = tbl[tbinds & tbl$pthres==.05,]
      pt = tbl$pthres
      
      for (trtype in c("train","test","all")) {
        
      # load p values
      pvs = llply(tb$pathp, function(x) get(load(x))[[trtype]])
      names(pvs) = tb$feat
      
      pathn = paste0(root,"/result/",data,"/plot_pval/",pa)
      ## qq plot; log both axis -----------------------
      nqs = llply(sapply(pvs,length), function(n) 1:n/(n+1))
      cs = rainbow(length(pvs))
      
      mcm1l = llply(pvs, function(x) {
        if (grepl("_",names(x)[1])) mcmind = match(names(x),names(mcm))
        if (!grepl("_",names(x)[1]))
          mcmind = match(sapply(strsplit(names(x),"_"), function(x) x[1]),names(mcm))
        mcm[mcmind]
      })
      names(nqs) = names(cs) = names(mcm1l) = tb$feat # feats
      
      png(paste0(pathn,"/qqpl",classn,"_",trtype,".png"), height=500, width=1000)
      par(mfrow=c(1,2),mar=c(7,7,3,3))
      
      plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab=paste0(pa, " p values"), xlab="theoretical quantile", main="", cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5)
      for (i in 1:length(pvs)) {
        yo = order(pvs[[i]])
        points(nqs[[i]], pvs[[i]][yo], pch=1, cex=mcm1l[[i]][yo]*4, col=adjustcolor(cs[i], alpha.f=.75))
      }
      abline(h=.05)
      lines(x = c(-100,100), y = c(-100,100))
      legend("topleft", legend=names(pvs), fill=cs, bg="transparent")
      
      nqls = llply(nqs, log)
      pvls = llply(pvs, log)
      plot(NULL, xlim=c(min(sapply(nqls,function(x)x[1])),0), ylim=c(min(sapply(nqls,function(x)x[1])),0), ylab=paste0(pa, " ln(p values)"), xlab="theoretical quantile", main="", cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5)
      for (i in 1:length(pvs)) {
        yo = order(pvls[[i]])
        points(nqls[[i]], pvls[[i]][yo], pch=1, cex=mcm1l[[i]][yo]*4, col=adjustcolor(cs[i], alpha.f=.75))
      }
      abline(h=log(.05))
      lines(x = c(-100,100), y = c(-100,100))
      legend("topleft", legend=names(pvs), fill=cs, bg="transparent")
      
      graphics.off()
      
      
      for (pt in c_pts) { 
        ptl = -log(pt)
        tb = tbl[tbinds & tbl$pthres==pt,]
        
        pathn = paste0(root,"/result/",data,"/plot_pval/",pa,"/pthres-",pt)
        
        ## number of significant features -----------------------
        plens = laply(pvs, length)
        npts = laply(pvs, function(pv) {
          a = sum(pv<pt); ifelse (a==0,0,a/length(pv))
        })
        names(pvs) = names(npts) = names(plens) = tb$feat
        
        png(paste0(pathn,"/num",classn,"_",trtype,".png"), height=500, width=500)
        par(mar=c(10,7,5,3))
        
        maint = paste0("% of sig/features (bar/line) \n data: ", data, "/", uc, ", method: ", pa)
        
        barplot(npts, ylim=c(0,1), xlab="", las=2, ylab="% significant features", col=rgb(0,0,0,.5), main=maint, cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5)
        
        abline(h=pt)
        par(new=T)
        
        plot(plens, pch=15, xlab="",ylab="",axes=F,type="b", col="red", ylim=c(0, max(plens)+(.15*max(plens))), cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5)
        mtext("number of features",side=4,col="red",line=4)
        axis(4, col="red", col.axis="red", las=1)
        
        graphics.off()
        
        
        
        ## robustness stats -----------------------
        fr = adply(tb,1, function(x) {
          ldply(c("rec","prec","f","pcorr"), function(y)
            data.frame(x,metric=y,score=as.numeric(x[[y]])) )
        })
        
        png(paste0(pathn,"/rpfc",classn,"_",trtype,".png"), height=500, width=500)
        par(mar=c(10,4,5,4))
        
        pl = barchart(score~feat, fr, groups=metric, las=1,# cex.axis=3,
                      ylim=c(0,1),
                      auto.key = list(columns=2),
                      scales=list(x=list(rot=45,cex=0.8)),
                      main=paste0("data: ", data, "/", uc, "\ntests: f measure, precision, recall, pearson correlation"), cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5)
        print(pl)
        
        graphics.off()
        
        
        
        
        ## pos stats -----------------------
        if(grepl("pos[0-9]",data)) {
          asig = llply(pvs, function(x) x<pt)
          rsig = llply(pvs, function(x) {
            if (data=="pos1") a = grepl("^A[+|-]$",names(x))
            if (data=="pos2") a = grepl("^[A|B][+|-]$",names(x))
            if (data%in%c("pos4","pos7")) a = grepl("^A[+|-]B[+|-]$",names(x))
            if (data=="pos5") a = grepl("^A[+|-]B[+|-]$",names(x)) | grepl("^D[+|-]E[+|-]$",names(x))
            if (data%in%c("pos3","pos6","pos8")) a = grepl("^A[+|-]B[+|-]C[+|-]$",names(x))

            return(a)
          })
          hsig = llply(pvs, function(x) {
            if (data=="pos1") a = grepl("A",names(x))
            if (data=="pos2") a = grepl("A|B",names(x))
            if (data%in%c("pos4","pos7")) a = grepl("A[+|-]B[+|-][C-Z|+|-]+",names(x))
            if (data=="pos5") a = grepl("A[+|-]B[+|-][C-Z|+|-]+",names(x)) | grepl("[A-C|+|-]+D[+|-]E[+|-]",names(x))
            
            if (data%in%c("pos3","pos6","pos8")) a = grepl("A[+|-]B[+|-]C[+|-][D-Z|+|-]+",names(x))
            
            return(a)
          })
          hsig = llply(1:length(hsig), function(j) hsig[[j]] & !rsig[[j]])
          names(pvs) = names(asig) = names(rsig) = names(hsig) = tb$feat
          
          
          for (t in c("real","influenced")) {
            if (t=="real") sig = rsig
            if (t=="influenced") sig = hsig
            # sig = switch(t, real = rsig, influenced = hsig)
            sbs = ldply(tb$feat, function(i) {
              overlap = sum(asig[[i]] & sig[[i]])
              if (overlap==0) {
                rec = prec = f = 0
              } else {
                rec = overlap/sum(asig[[i]])
                prec = overlap/sum(sig[[i]])
                f = 2*((prec*rec)/(prec+rec))
              }
              return(data.frame(feature=paste0(i," (",sum(asig[[i]]),",",overlap,",",sum(sig[[i]]),")"), type=t,
                                sig=sum(sig[[i]]), actual=sum(asig[[i]]),
                                metric=c("recall","precision","f"),
                                score=c(rec,prec,f)))
            })
            
            png(paste0(pathn,"/pos-",t,classn,"_",trtype,".png"), height=500, width=500)
            pl = barchart(score~feature, sbs, groups=metric, las=1, ylim=c(0,1),
                          main="feature (actual,overlap,hypothetical) \noverlap / actual (p) vs hypothetical (r)",
                          auto.key = list(space = "top"),
                          scales = list(x = list(rot=45)), cex.lab=2, cex.main=2, cex.sub=2, cex.axis=1.5)
            print(pl)
            graphics.off()
          }
        }
      }
      }
    }
  }    
}
time_output(start1)



## sig cell population histograms --------------------
start1 = Sys.time()
inds = which(tbl$m_all_sig>0 & tbl$test_adj!="none")
l_ply(loopInd(sample(inds,length(inds)),no_cores), function(ii) {
  i_d = i_f = ""
  for (i in ii) {
    pathn = paste0(root,"/result/",tbl$data[i],"/plot_pval/",pa,"/overlap/",tbl$class[i],"/",tbl$feat[i])
    dir.create(pathn, showWarnings=F, recursive=T)
    
    # load p values
    pvs = get(load(tbl$pathp[i]))$all
    mo = get(load(tbl$patho[i]))$all
    
    i_d_ = i_d==tbl$data[i]
    i_f_ = i_f==tbl$feat[i]
    if (!i_d_) 
      sm0 = get(load(paste0(root,"/result/",tbl$data[i],"/meta/file.Rdata")))
    
    if (!i_f_ | !i_d_) {
      m0 = get(load(paste0(root,"/result/",tbl$data[i],"/feat/",tbl$feat[i],".Rdata")))
      sm = sm0[match(rownames(m0),sm0$id),]
    }
    m = m0
    mc = m[sm$class=="control",]
    me = m[sm$class==tbl$class[i],]
    
    for (j in names(pvs)[pvs<tbl$pthres[i]]) {
      moi = names(mo)==j
      if (mo[moi]==1) next()
      if (mo[moi]>overlapmin) next() # & grepl("^ctrl|^pos",tbl$data[i])) next()
      mi = colnames(m)==j
      png(paste0(pathn,"/",round(mo[moi],3),"_",j,".png"))
      mr = range(c(mc[,mi],me[,mi]))
      h = hist(mc[,mi], plot=F)
      h$counts <- h$counts / sum(h$counts)
      plot(h, col=rgb(1,0,0,0.5),xlim=c(mr[1],mr[2]), ylim=c(0,1), main="red=control; blue=experiment")
      par(new=T)
      h_ = hist(me[,mi], plot=F)
      h_$counts <- h_$counts / sum(h_$counts)
      plot(h_, col=rgb(0,0,1,0.5),xlim=c(mr[1],mr[2]), ylim=c(0,1), add=T)
      graphics.off()
    }
    
    i_d = tbl$data[i]
    i_f = tbl$feat[i]
  }
}, .parallel=T)
time_output(start1)




## graph plots ----------------------------------
tbl = table[!grepl("raw|edge|group|entropy",table$feat),]
tbl$pmethod_adj = paste(tbl$test, tbl$test_adj, sep="-")
tbl = tbl[!grepl("ctrl[1-9]",tbl$data),]

start1 = Sys.time()
# l_ply(loopInd(sample(which(tbl$m_all_sig>0 & tbl$data%in%names(grp0s) & tbl$test_adj!="none"),1050), no_cores), function(ii) {
  l_ply(loopInd(which(tbl$m_all_sig>0 & tbl$data%in%names(grp0s) & tbl$test_adj!="none" & tbl$pthres==.01 & tbl$pmethod_adj=="t-BY" & !grepl("^pos",tbl$data) & grepl("prop|countAdj|lnpropexpect",tbl$feat) & grepl("pregnancy|flowcap6$|genentech",tbl$data)), no_cores), function(ii) {
    #   for (i in ii) {
  # for (i in which(tbl$m_all_sig>0 & tbl$data%in%names(grp0s) & tbl$test_adj!="none")) {
  for (i in ii) { 
    
    # do only for cell feature types and nodes<1000 see below... nvm
    
    pt = tbl$pthres[i]; ptl = -log(pt)
    classn = ifelse(tbl$class[i]=="exp", "", paste0("_",tbl$class[i]))
    pathn = paste0(root,"/result/",tbl$data[i],"/plot_pval/",tbl$pmethod_adj[i],"/pthres-",pt)
    
    main0 = paste0("data: ",tbl$data[i],"; feat: ",tbl$feat[i], "; class: ",tbl$class[i], "; pthres: ",round(-log(tbl$pthres[i]),3))
    mfm = mfms[[tbl$data[i]]][[tbl$feat[i]]]
    mcm = mcms[[tbl$data[i]]]
    mfm_tf = grepl("lnpropexpect",tbl$feat[i])
    if (mfm_tf) mfm_ = mfms[[tbl$data[i]]][[paste0(tbl$feat[i],"-raw")]]
    mo0 = get(load(tbl$patho[i]))$all
    pv = get(load(tbl$pathp[i]))
    # ptr = pv$train
    # pte = pv$test
    p = pv$all
    if (tbl$data[i]%in%c("pos1","pos2")) p = p[!grepl("D|E|F|G|H",names(p))]
    if (tbl$data[i]%in%c("pos6")) p = p[!grepl("E|F|G|H",names(p))]
    
    p_ = p<pt
    # ptr_ = ptr<pt & pte>pt
    # pte_ = pte<pt & ptr>pt
    mfmc = mfm$control
    mfmu = mfm[[tbl$class[i]]]
    mfmc = mfmc[match(names(p),names(mfmc))]
    mfmu = mfmu[match(names(p),names(mfmu))]
    if (mfm_tf) {
      mfmc_ = mfm_$control
      mfmu_ = mfm_[[tbl$class[i]]]
      mfmc_ = mfmc_[match(names(p),names(mfmc_))]
      mfmu_ = mfmu_[match(names(p),names(mfmu_))]
    }
    mfmdiff = abs(mfmu-mfmc)/sd(mfmc)
    mcmc = mcm$control
    mcmu = mcm[[tbl$class[i]]]
    mcmc = mcmc[match(names(p),names(mcmc))]
    mcmu = mcmu[match(names(p),names(mcmu))]
    gr0 = gr0s[[tbl$data[i]]]
    gr0$v = gr0$v[match(names(p),gr0$v$name),,drop=F]
    gr0$v = data.frame(name=unlist(gr0$v))
    mo = mo0[match(names(p),names(mo0))]
    gr0$e = gr0$e[gr0$e$from%in%gr0$v$name & gr0$e$to%in%gr0$v$name,]
    gr = layout_gr(gr0$e,gr0$v)
    gr = gpdf(gr)
    gr$v$color = ifelse(mfmu>mfmc,"_increase","decrease") # node colour
    # label_ind = rep(F,length(p))
    # label_ind_ = p+mo/2
    # label_ind[which(p_)[head(order(label_ind_[p_]),5)]] = T
    if (grepl("^ctrl|^pos",tbl$data[i])) {
      label_ind = p_ & mo==0 & !grepl("[-]",names(p_))
      gr$v$label = paste0(gr$v$name,":",round(mcmu,3))#,"/",round(mfmc,3))
      gr$v$size = mfmdiff
      main = paste0(main0,"\nsize=# of sd apart; label=prop(if pos/ctrl/prop)")
      
    } else {
      label_ind = rep(F,length(p))
      # label_ind_ = p+mo/2
      label_ind[which(p_)[tail(order(mfmdiff[p_]),20)]] = T
      gr$v$label = paste0(gr$v$name,":",round(mfmu,3))#,"/",round(mfmc,3))
      gr$v$size = -log(p)
      gr$v$size[is.infinite(gr$v$size)] = max(gr$v$size[!is.infinite(gr$v$size)])
      main = paste0(main0,"\nsize=-ln(p value); label=feature value (prop if lnpropexpect)")
      
    }
    
    gpp = gggraph(
      gr, v_ind=p_, vb_ind=rep(F,nrow(gr$v)), main=main,
      e_ind=gr$e[,1]%in%gr$v$name[p_] & gr$e[,2]%in%gr$v$name[p_],
      label_ind=rep(F,nrow(gr$v)))
    # gpp = gpp + geom_point(data=gr$v[p_,],aes(x=x,y=y),size=2,color="grey")
    gp = gpp + geom_label_repel(
      data=gr$v[label_ind,],
      aes(x=x,y=y,label=label, color=color),
      nudge_x=-.1, direction="y", hjust=1, segment.size=0.2)
    
    ggsave(paste0(pathn,"/gr",classn,"_",tbl$feat[i],".png"), plot=gp, scale=1, width=9, height=9, units="in", dpi=500, limitsize=T)
    
    if (mfm_tf) {
      # if (grepl("^pos",tbl$data[i])) {
      #   gr$v$color = (mfmu_-mcmu)/mcmu # node colour
      # }
      gr$v$label = paste0(gr$v$name,":",round(mfmu_,3),"/",round(mfmc_,3))
      gp = gpp + geom_label_repel(
        data=gr$v[label_ind,],
        aes(x=x,y=y,color=color,label=label),
        nudge_x=-.1, direction="y", hjust=1, segment.size=0.2)
      
      ggsave(paste0(pathn,"/gr",classn,"_",tbl$feat[i],"-raw.png"), plot=gp, scale=1, width=9, height=9, units="in", dpi=500, limitsize=T)
    }
    
    
  
    
    
    try ({
      pn = names(p)[p_ & mo<overlapmin]
      pnse = str_extract_all(pn,"[A-Za-z0-9]+[+|-]")
      pnseu = sort(unique(unlist(pnse)))
      templ = rep(0,length(pnseu))
      pndf = Reduce(rbind,llply(loopInd(1:length(pnse),no_cores), function(ii) {
        Reduce(rbind,llply(ii, function(i) {
          x = pnse[[i]]
          templ[match(x,pnseu)] = 1
          templ
        }))
      }, .parallel=T))
      colnames(pndf) = pnseu
      rules = apriori(pndf, parameter=list(support=1/nrow(pndf), confidence=.1))#, confidence=0.5))
      # ivl = str_count(iv,"[+|-]")
      # igr0 = list(v=data.frame(name=iv), e=as.data.frame(t(apply(which(ia, arr.ind=T),1,function(x) {
      #   if (ivl[x[2]]-ivl[x[1]]==1) return(c(iv[x[1]],iv[x[2]]))
      #   return(c(NA,NA))
      # } ))))
      # colnames(igr0$e) = c("from","to")  
      # igr0$e = igr0$e[!is.na(igr0$e[,1]),]
      # igr = layout_gr(igr0$e,igr0$v)
      # igr = gpdf(igr)
      # igr$v$label = paste0(igr$v$name,":",ic)
      # gpp = gggraph(
      #   igr, v_ind=rep(T,nrow(igr$v)), vb_ind=rep(F,nrow(igr$v)), main="",
      #   e_ind=rep(T,nrow(igr$e)),
      #   label_ind=rep(F,nrow(gr$v)))
      
      
      # rules plot
      dir.create(paste0(pathn,"/arules/",tbl$class[i]), showWarnings=F, recursive=T)
      png(paste0(pathn,"/arules/",tbl$class[i],"/",tbl$feat[i],"_paracoord.png"))
      plot(rules, method="paracoord") # lines arrow
      graphics.off()
      png(paste0(pathn,"/arules/",tbl$class[i],"/",tbl$feat[i],"_graph.png"))
      plot(rules, method="graph") # dots with arrows point to them by lift
      graphics.off()
      png(paste0(pathn,"/arules/",tbl$class[i],"/",tbl$feat[i],"_scatterplot.png"))
      plot(rules, method="scatterplot")
      graphics.off()
      png(paste0(pathn,"/arules/",tbl$class[i],"/",tbl$feat[i],"_grouped.png"))
      plot(rules, method="grouped") # item in LHs by RHS
      graphics.off()
      
      
      # itemset plot
      itemsets = eclat(pndf, parameter=list(support=1/nrow(pndf)))#, confidence=.1))#, confidence=0.5))
      ic = itemsets@quality$count
      idu = itemsets@items@itemInfo$labels
      id = itemsets@items@data
      # ia = is.subset(itemsets)
      # colnames(ia) = rownames(ia) = gsub("[{|}|,]","",colnames(ia))
      iv = sapply(1:ncol(id), function(j) paste(idu[id[,j]],collapse=""))
      names(ic) = iv
      
      ivm = match(gr$v$name,names(ic))
      gr$v$size = ic[ivm]
      gr$v$label = paste0(gr$v$name,":",ic[ivm])
      lbnotna = !grepl("NA$",gr$v$label)
      label_ind[!lbnotna] = F

    main = paste0(main0,"\ncolours=p value (if pos_, %change expect/prop); size=frequent itemset count; label=see size")
    gpp = gggraph(
      gr, v_ind=lbnotna, vb_ind=rep(F,nrow(gr$v)), main=main,
      e_ind=gr$e[,1]%in%gr$v$name[lbnotna] & gr$e[,2]%in%gr$v$name[lbnotna],
      label_ind=rep(F,nrow(gr$v)))
    gp = gpp + geom_label_repel(
      data=gr$v[label_ind,],
      aes(x=x,y=y,label=label, color=color),
      nudge_x=-.1, direction="y", hjust=1, segment.size=0.2)
    
    ggsave(paste0(pathn,"/gr",classn,"_",tbl$feat[i],"_itemset.png"), plot=gp, scale=1, width=9, height=9, units="in", dpi=500, limitsize=T)
    
    gp = gpp + geom_label_repel(
      data=gr$v[lbnotna,],
      aes(x=x,y=y,label=label, color=color),
      nudge_x=-.1, direction="y", hjust=1, segment.size=0.2)
    
    ggsave(paste0(pathn,"/gr",classn,"_",tbl$feat[i],"_itemset_alllabellled.png"), plot=gp, scale=1, width=9, height=9, units="in", dpi=500, limitsize=T)
    
    
    })
    
    
    
    
    
    
    
    if (!grepl("short",tbl$feat[i])) {
      posind = !grepl("[-]",names(p))
      pp = p[posind] # p
      pp_ = pp<pt # p sig
      # ptrp_ = ptr_[posind] # p sig train only
      # ptep_ = pte_[posind] # p sig test only
      mfmcp = mfmc[match(names(pp),names(mfmc))] # prop control
      mfmup = mfmu[match(names(pp),names(mfmu))] # prop exp
      if (mfm_tf) {
        mfmcp_ = mfmc_[match(names(pp),names(mfmc_))] # prop control
        mfmup_ = mfmu_[match(names(pp),names(mfmu_))] # prop exp
      }
      mfmdiffp = abs(mfmup-mfmcp)/sd(mfmcp)
      mcmcp = mcmc[match(names(pp),names(mcmc))] # prop control
      mcmup = mcmu[match(names(pp),names(mcmu))] # prop exp
      mop = mo0[match(names(pp),names(mo0))]
      grp0 = grp0s[[tbl$data[i]]] # load graph
      grp0$v = grp0$v[match(names(pp),grp0$v$name),,drop=F]
      grp0$v = data.frame(name=unlist(grp0$v))
      grp0$e = grp0$e[grp0$e$from%in%grp0$v$name & grp0$e$to%in%grp0$v$name,]
      grp = layout_gr(grp0$e,grp0$v) # graph layout
      grp = gpdf(grp) # add graph attributes
      grp$v$label = paste0(grp$v$name,":",round(mfmup,3))
      grp$v$color = ifelse(mfmup>mfmcp,"_increase","decrease") # node colour
      # grp$v$colorb = ifelse(ptrp_,"set1","none") # node border colour
      # grp$v$colorb[ptep_] = "set2"
      # grp$v$colorb[(!ptep_ | !ptrp_) & pp_] = "all"
      # grp$v$colorb = as.factor(grp$v$colorb)
      # grp$v$fill = ifelse(ptep_ | ptrp_,F,T)
      if (grepl("^ctrl|^pos",tbl$data[i])) {
        label_ind = pp_ & mop==0 & !grepl("[-]",names(pp_))
        grp$v$label = paste0(grp$v$name,":",round(mcmup,3))#,"/",round(mfmc,3))
        grp$v$size = mfmdiffp
        main = paste0(main0,"\npositive nodes only\nsize=# of sd apart; label=prop(if pos/ctrl/prop)")
        
      } else {
        label_ind = rep(F,length(pp))
        # label_ind_ = p+mo/2
        label_ind[which(pp_)[tail(order(mfmdiffp[p_]),20)]] = T
        grp$v$label = paste0(grp$v$name,":",round(mfmup,3))#,"/",round(mfmc,3))
        grp$v$size = -log(p)
        grp$v$size[is.infinite(grp$v$size)] = max(grp$v$size[!is.infinite(grp$v$size)])
        main = paste0(main0,"\npositive nodes only\nsize=-ln(p value); label=feature value (prop if lnpropexpect)")
        
      }
      
      main = paste0(main0,"\ncolours=p value (if ctrl/pos, %change); size=1-overlap(control); label=prop(exp)\n")
      gpp = gggraph(
        grp, v_ind=pp_, vb_ind=rep(F,nrow(grp$v)), main=main,
        e_ind=grp$e[,1]%in%grp$v$name[pp_] & grp$e[,2]%in%grp$v$name[pp_],
        label_ind=rep(F,nrow(grp$v)))
      gp = gpp + geom_label_repel(
        data=grp$v[label_ind,],
        aes(x=x,y=y,label=label, color=color),
        nudge_x=-.1, direction="y", hjust=1, segment.size=0.2) 
      
      ggsave(paste0(pathn,"/grpos",classn,"_",tbl$feat[i],".png"), plot=gp, scale=1, width=9, height=9, units="in", dpi=500, limitsize=T)
      
      if (mfm_tf) {
        if (grepl("^pos",tbl$data[i])) {
          grp$v$color = (mfmup_-mcmup)/mcmup # node colour
        }
        grp$v$label = paste0(grp$v$name,":",round(mfmup_,3),"/",round(mfmcp_,3))
        gp = gpp + geom_label_repel(
          data=grp$v[label_ind,],
          aes(x=x,y=y,color=color,label=label), 
          nudge_x=-.1, direction="y", hjust=1, segment.size=.2) 
        ggsave(paste0(pathn,"/grpos",classn,"_",tbl$feat[i],"-raw.png"), plot=gp, scale=1, width=9, height=9, units="in", dpi=500, limitsize=T)
      }
      
      
      try ({
        pn = names(pp)[pp_ & mop<overlapmin]
        pnse = str_extract_all(pn,"[A-Za-z0-9]+[+|-]")
        pnseu = sort(unique(unlist(pnse)))
        templ = rep(0,length(pnseu))
        pndf = Reduce(rbind,llply(loopInd(1:length(pnse),no_cores), function(ii) {
          Reduce(rbind,llply(ii, function(i) {
            x = pnse[[i]]
            templ[match(x,pnseu)] = 1
            templ
          }))
        }, .parallel=T))
        colnames(pndf) = pnseu
        rules = eclat(pndf, parameter=list(support=1/nrow(pndf)))#, confidence=0.5))
        dir.create(paste0(pathn,"/arules/",tbl$class[i]), showWarnings=F, recursive=T)
        png(paste0(pathn,"/arules/",tbl$class[i],"/",tbl$feat[i],"_paracoord_pos.png"))
        plot(rules, method="paracoord") # lines arrow
        graphics.off()
        png(paste0(pathn,"/arules/",tbl$class[i],"/",tbl$feat[i],"_graph_pos.png"))
        plot(rules, method="graph") # dots with arrows point to them by lift
        graphics.off()
        png(paste0(pathn,"/arules/",tbl$class[i],"/",tbl$feat[i],"_scatterplot_pos.png"))
        plot(rules, method="scatterplot")
        graphics.off()
        png(paste0(pathn,"/arules/",tbl$class[i],"/",tbl$feat[i],"_grouped_pos.png"))
        plot(rules, method="grouped") # item in LHs by RHS
        graphics.off()
        
        
        # itemset plot
        itemsets = eclat(pndf, parameter=list(support=1/nrow(pndf)))#, confidence=.1))#, confidence=0.5))
        ic = itemsets@quality$count
        idu = itemsets@items@itemInfo$labels
        id = itemsets@items@data
        # ia = is.subset(itemsets)
        # colnames(ia) = rownames(ia) = gsub("[{|}|,]","",colnames(ia))
        iv = sapply(1:ncol(id), function(j) paste(idu[id[,j]],collapse=""))
        names(ic) = iv
        
        ivm = match(grp$v$name,names(ic))
        grp$v$size = ic[ivm]
        grp$v$label = paste0(grp$v$name,":",ic[ivm])
        lbnotna = !grepl(":NA",grp$v$label)
        # label_ind = rep(F,nrow(grp$v))
        label_ind[!lbnotna] = F
        
        main = paste0(main0,"\ncolours=p value (if pos_, %change expect/prop); size=frequent itemset count; label=see size")
        gpp = gggraph(
          grp, v_ind=lbnotna, vb_ind=rep(F,nrow(grp$v)), main=main,
          e_ind=grp$e[,1]%in%grp$v$name[lbnotna] & grp$e[,2]%in%grp$v$name[lbnotna],
          label_ind=rep(F,nrow(grp$v)))
        gp = gpp + geom_label_repel(
          data=grp$v[label_ind,],
          aes(x=x,y=y,label=label, color=color),
          nudge_x=-.1, direction="y", hjust=1, segment.size=0.2)
        
        ggsave(paste0(pathn,"/grpos",classn,"_",tbl$feat[i],"_itemset.png"), plot=gp, scale=1, width=9, height=9, units="in", dpi=500, limitsize=T)
        
        gp = gpp + geom_label_repel(
          data=grp$v[lbnotna,],
          aes(x=x,y=y,label=label, color=color),
          nudge_x=-.1, direction="y", hjust=1, segment.size=0.2)
        
        ggsave(paste0(pathn,"/grpos",classn,"_",tbl$feat[i],"_itemset_alllabelled.png"), plot=gp, scale=1, width=9, height=9, units="in", dpi=500, limitsize=T)
        
        
      })
      
      
    }
  }
},.parallel=T)
time_output(start1)







