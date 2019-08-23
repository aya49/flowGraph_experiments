()## input: features 
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


start = Sys.time()


## load table
table = get(load(paste0(root,"/pvals.Rdata")))

# use tbl (trim)
tbl = table
tbl = tbl[!grepl("raw|edge|group|entropy",tbl$feat),]
tbl$pmethod_adj = paste(tbl$test, tbl$test_adj, sep="-")
tblc = tbl[grepl("^ctrl",tbl$data),]
tbl = tbl[!grepl("ctrl[1-9]",tbl$data),]

pt = table$p_thres[1]
ptl = -log(table$p_thres[1])

# pars
c_datas = unique(tbl$data)
c_feats = unique(tbl$feat)
c_pas = unique(tbl$pmethod_adj)
c_pts = unique(tbl$pthres)



## default counts
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
names(mcms) = names(mcms) = c_datas


## default graphs
gr0s = llply(c_datas, function(data) {
  get(load(paste0(root,"/result/",data,"/meta/cell_graph.Rdata")))
})
names(gr0s) = c_datas
for (i in length(gr0s):1) if (nrow(gr0s[[i]]$v)>2500) gr0s[[i]] = NULL
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
    
    
    png(paste0(pathn,"/hist_p",classn,"_",tbl$feat[i],".png"), height=800, width=800)
    mv1 = matrix(c(1,1,2,3,4,4,5,6),ncol=2,byrow=F)
    layout(mv1) # scatterplot + histograms
    
    mcm = mcms[[tbl$data[i]]]$control*4
    
    set1 = pv$train
    set2 = pv$test
    
    if (grepl("-cell-",tbl$feat[i])) 
      mcmind = match(names(set1),names(mcm))
    if (!grepl("-cell-",tbl$feat[i])) 
      mcmind = match(sapply(strsplit(names(set1),"_"), function(x) x[1]),names(mcm))
    mcm1 = mcm[mcmind]
    
    plot_int(cbind(set1,set2), pch=1, xlim=c(0,1), ylim=c(0,1), 
             xlab="set 1 p values", ylab="set 2 p values", main=paste0("data:",tbl$data[i],"; class: ",tbl$class[i], "; feature: ", tbl$feat[i], "; method: ", tbl$pmethod_adj[i]), cex=mcm1)
    abline(v=pt,h=pt, col="red")
    
    hist(set1, freq=F, main=paste0("set 1"), xlim=c(0,1), breaks=seq(0,1,.05), xlab="p values", ylab="histogram", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col=rgb(0,0,0,.5))
    abline(v=pt, col="red")
    
    hist(set2, freq=F, main="set 2", xlim=c(0,1), breaks=seq(0,1,.05), xlab="p values", ylab="histogram", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col=rgb(0,0,0,.5))
    abline(v=pt, col="red")
    
    
    set1 = -log(set1); 
    set1[set1==Inf] = 10
    set2 = -log(set2); 
    set2[set2==Inf] = 10
    plot_int(cbind(set1,set2), pch=1, xlim=c(0,10), ylim=c(0,10), 
             xlab="set 1 -ln(p values)", ylab="set 2 -ln(p values)", main=paste0("set2 class: ",tbl$class[i], "; feature: ", tbl$feat[i], "; method: ", tbl$pmethod_adj[i]), cex=mcm1)
    abline(v=ptl,h=ptl, col="red")
    
    hist(set1, freq=F, main=paste0("set 1"), xlim=c(0,max(4,max(set1))), xlab="-ln(p values)", ylab="histogram", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col=rgb(0,0,0,.5))
    abline(v=ptl, col="red")
    
    hist(set2, freq=F, main="set 2", xlim=c(0,max(4,max(set2))), xlab="-ln(p values)", ylab="histogram", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col=rgb(0,0,0,.5))
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
    
    
    png(paste0(pathn,"/hist_o",classn,"_",tbl$feat[i],".png"), height=800, width=400)
    mv1 = matrix(c(1,1,2,3),ncol=1,byrow=F)
    layout(mv1) # scatterplot + histograms
    
    mcm = mcms[[tbl$data[i]]]$control*4
    
    set1 = mo$train
    set2 = mo$test
    
    if (grepl("-cell-",tbl$feat[i])) 
      mcmind = match(names(set1),names(mcm))
    if (!grepl("-cell-",tbl$feat[i])) 
      mcmind = match(sapply(strsplit(names(set1),"_"), function(x) x[1]),names(mcm))
    mcm1 = mcm[mcmind]
    
    plot_int(cbind(set1,set2), pch=1, xlim=c(0,1), ylim=c(0,1), 
             xlab="set 1  overlaps/control", ylab="set 2 overlaps/control", main=paste0("data:",tbl$data[i],"; class: ",tbl$class[i], "; feature: ", tbl$feat[i]), cex=mcm1)
    
    hist(set1, freq=F, main=paste0("set 1"), xlim=c(0,1), breaks=seq(0,1,.05), xlab="p values", ylab="histogram", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col=rgb(0,0,0,.5))
    
    hist(set2, freq=F, main="set 2", xlim=c(0,1), breaks=seq(0,1,.05), xlab="p values", ylab="histogram", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col=rgb(0,0,0,.5))
    
    graphics.off()
  }
},.parallel=T)
time_output(start1)




## control qq plot -------------------------------
classn = ""
mcm = mcms[["ctrl0"]]$control
for (pa in c_pas) {
  tb = tblc[tblc$pmethod_adj==pa & tblc$pthres==.05,]
  
  pvs = llply(unique(tb$feat), function(x) {
    unlist(llply(tb$pathp[tb$feat==x], function(y) get(load(y))$all))
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
  par(mfrow=c(1,2))
  
  plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab=paste0(pa, " p values"), xlab="theoretical quantile", main=paste0("qq plot for all contrls (lines=.01/.025/.05); test: ", pa))
  for (i in 1:length(pvs)) {
    yo = order(pvs[[i]])
    points(nqs[[i]], pvs[[i]][yo], pch=1, cex=mcm1l[[i]][yo], col=adjustcolor(cs[i], alpha.f=.5))
  }
  abline(h=.05)
  abline(h=.025)
  abline(h=.01)
  lines(x = c(-100,100), y = c(-100,100))
  legend("topleft", legend=names(pvs), fill=cs, bg="transparent")
  
  nqls = llply(nqs, log)
  pvls = llply(pvs, log)
  plot(NULL, xlim=c(min(sapply(nqls,function(x)x[1])),0), ylim=c(min(sapply(nqls,function(x)x[1])),0), ylab=paste0(pa, " p values"), xlab="theoretical quantile", main="ln(qq plot)")
  for (i in 1:length(pvs)) {
    yo = order(pvls[[i]])
    points(nqls[[i]], pvls[[i]][yo], pch=1, cex=mcm1l[[i]][yo], col=adjustcolor(cs[i], alpha.f=.5))
  }
  abline(h=log(.05))
  abline(h=log(.025))
  abline(h=log(.01))
  lines(x = c(-100,100), y = c(-100,100))
  legend("topleft", legend=names(pvs), fill=cs, bg="transparent")
  
  graphics.off()
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
      
      # load p values
      pvs = llply(tb$pathp, function(x) get(load(x))$all)
      names(pvs) = tb$feat
      
      pathn = paste0(root,"/result/",data,"/plot_pval/",pa,"/pthres-",pt)
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
      
      png(paste0(pathn,"/qqpl",classn,".png"), height=400, width=800)
      par(mfrow=c(1,2))
      
      plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab=paste0(pa, " p values"), xlab="theoretical quantile", main=paste0("qq plot | data: ", data, ", class: ", uc, ", method: ", pa))
      for (i in 1:length(pvs)) {
        yo = order(pvs[[i]])
        points(nqs[[i]], pvs[[i]][yo], pch=1, cex=mcm1l[[i]][yo], col=adjustcolor(cs[i], alpha.f=.5))
      }
      abline(h=.05)
      lines(x = c(-100,100), y = c(-100,100))
      legend("topleft", legend=names(pvs), fill=cs, bg="transparent")
      
      nqls = llply(nqs, log)
      pvls = llply(pvs, log)
      plot(NULL, xlim=c(min(sapply(nqls,function(x)x[1])),0), ylim=c(min(sapply(nqls,function(x)x[1])),0), ylab=paste0(pa, " p values"), xlab="theoretical quantile", main="ln(qq plot)")
      for (i in 1:length(pvs)) {
        yo = order(pvls[[i]])
        points(nqls[[i]], pvls[[i]][yo], pch=1, cex=mcm1l[[i]][yo], col=adjustcolor(cs[i], alpha.f=.5))
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
        
        png(paste0(pathn,"/num",classn,".png"), height=400, width=400)
        par(mar=c(10,4,5,4))
        
        maint = paste0("%/# of sig/features (bar/line) \n data: ", data, "/", uc, ", method: ", pa)
        
        barplot(npts, ylim=c(0,1), xlab="", las=2, ylab="% significant features", col=rgb(0,0,0,.5), main=maint)
        
        abline(h=pt)
        par(new=T)
        
        plot(plens, pch=15, xlab="",ylab="",axes=F,type="b", col="red", ylim=c(0, max(plens)+(.15*max(plens))))
        mtext("number of features",side=4,col="red",line=4)
        axis(4, col="red", col.axis="red", las=1)
        
        graphics.off()
        
        
        
        ## robustness stats -----------------------
        fr = adply(tb,1, function(x) {
          ldply(c("rec","prec","f","pcorr"), function(y)
            data.frame(x,metric=y,score=as.numeric(x[[y]])) )
        })
        
        png(paste0(pathn,"/rpfc",classn,".png"), height=400, width=400)
        par(mar=c(10,4,5,4))
        
        pl = barchart(score~feat, fr, groups=metric, las=1,# cex.axis=3,
                      ylim=c(0,1),
                      auto.key = list(columns=2),
                      scales=list(x=list(rot=45,cex=0.8)),
                      main=paste0("data: ", data, "/", uc, "\ntests: f measure, precision, recall, pearson correlation"))
        print(pl)
        
        graphics.off()
        
        
        
        
        ## pos stats -----------------------
        if(grepl("pos[0-9]",data)) {
          asig = llply(pvs, function(x) x<pt)
          rsig = llply(pvs, function(x) {
            if (data=="pos1") a = grepl("^A[+|-]$",names(x))
            if (data=="pos2") a = grepl("^[A|B][+|-]$",names(x))
            if (data%in%c("pos3","pos5")) a = grepl("^A[+|-]B[+|-]C[+|-]$",names(x))
            if (data=="pos4") a = grepl("^A[+|-]B[+|-]C[+|-]D[+|-]$",names(x))
            
            return(a)
          })
          hsig = llply(pvs, function(x) {
            if (data=="pos1") a = grepl("A",names(x))
            if (data=="pos2") a = grepl("A|B",names(x))
            if (data=="pos3") a = grepl("A[+|-]B[+|-]C[+|-][D-Z|+|-]+",names(x))
            if (data=="pos4") a = grepl("A[+|-]B[+|-]C[+|-]D[+|-][E-Z|+|-]+",names(x))
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
            
            png(paste0(pathn,"/pos-",t,classn,".png"), height=400, width=400)
            pl = barchart(score~feature, sbs, groups=metric, las=1, ylim=c(0,1),
                          main="feature (actual,overlap,hypothetical) \noverlap / actual (p) vs hypothetical (r)",
                          auto.key = list(space = "top"),
                          scales = list(x = list(rot = 45)))
            print(pl)
            graphics.off()
          }
        }
        
      }
    }
  }    
}
time_output(start1)





## graph stats ----------------------------------
start1 = Sys.time()
# a=llply(loopInd(which(tbl$m_all_sig>0 & tbl$data%in%names(grp0s)), no_cores), function(ii) {
#   for (i in ii) {
    for (i in which(tbl$m_all_sig>0 & tbl$data%in%names(grp0s) & tbl$test_adj!="none")) {
      
    # do only for cell feature types and nodes<1000 see below
    
    pt = tbl$pthres[i]; ptl = -log(pt)
    classn = ifelse(tbl$class[i]=="exp", "", paste0("_",tbl$class[i]))
    pathn = paste0(root,"/result/",tbl$data[i],"/plot_pval/",tbl$pmethod_adj[i],"/pthres-",pt)
    
    main0 = paste0("data: ",tbl$data[i],"; feat: ",tbl$feat[i], "; class: ",tbl$class[i], "; pthres: ",tbl$pthres[i])
    
    mcm = mcms[[tbl$data[i]]]
    mo0 = get(load(tbl$patho[i]))$all
    pv = get(load(tbl$pathp[i]))
    ptr = pv$train
    pte = pv$test
    p = pv$all
    
    p_ = p<pt
    ptr_ = ptr<pt & pte>pt
    pte_ = pte<pt & ptr>pt
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
    gr$v$label = paste0(gr$v$name,":",round(mcmu,3))
    gr$v$size = 1-mo
    gr$v$color = p # node colour
    label_ind = rep(F,length(p))
    label_ind_ = p+mo/2
    label_ind[which(p_)[head(order(label_ind_[p_]),5)]] = T
    
    main = paste0(main0,"\ncolours=p value; size=1-overlap(control); label=prop(exp)\n")
    gpp = gggraph(
      gr, v_ind=p_, vb_ind=rep(F,nrow(gr$v)), main=main,
      e_ind=gr$e[,1]%in%gr$v$name[p_] & gr$e[,2]%in%gr$v$name[p_],
      label_ind=rep(F,nrow(gr$v)))
    gp = gpp + geom_label_repel(
      data=gr$v[label_ind,],
      aes(x=x,y=y,label=label, color=color),
      nudge_x = -.2, direction = "y", hjust = 1, segment.size = 0.2) +
      geom_point(data=gr$v[ptr_|pte_,],aes(x=x,y=y),size=2,color="grey")
    
    ggsave(paste0(pathn,"/gr",classn,"_",tbl$feat[i],".png"), plot=gp, scale = 1, width =11, height =15, units = "in", dpi = 300, limitsize = TRUE)
    
    
    if (!grepl("short",tbl$feat[i])) {
      posind = !grepl("[-]",names(p))
      pp = p[posind] # p
      pp_ = pp<pt # p sig
      ptrp_ = ptr_[posind] # p sig train only
      ptep_ = pte_[posind] # p sig test only
      mcmcp = mcmc[match(names(pp),names(mcmc))] # prop control
      mcmup = mcmu[match(names(pp),names(mcmu))] # prop exp
      mop = mo0[match(names(pp),names(mo0))]
      grp0 = grp0s[[tbl$data[i]]] # load graph
      grp0$v = grp0$v[match(names(pp),grp0$v$name),,drop=F]
      grp0$v = data.frame(name=unlist(grp0$v))
      grp0$e = grp0$e[grp0$e$from%in%grp0$v$name & grp0$e$to%in%grp0$v$name,]
      grp = layout_gr(grp0$e,grp0$v) # graph layout
      grp = gpdf(grp) # add graph attributes
      grp$v$label = paste0(grp$v$name,":",round(mcmup,3))
      grp$v$size = 1-mop
      grp$v$color = pp # node colour
      # grp$v$colorb = ifelse(ptrp_,"set1","none") # node border colour
      # grp$v$colorb[ptep_] = "set2"
      # grp$v$colorb[(!ptep_ | !ptrp_) & pp_] = "all"
      # grp$v$colorb = as.factor(grp$v$colorb)
      # grp$v$fill = ifelse(ptep_ | ptrp_,F,T)
      label_ind = rep(F,length(pp))
      label_ind_ = pp+mop/2
      label_ind[which(pp_)[head(order(label_ind_[pp_]),5)]] = T
      
      main = paste0(main0,"\ncolours=p value; size=1-overlap(control); label=prop(exp)\n")
      gpp = gggraph(
        grp, v_ind=pp_, vb_ind=rep(F,nrow(grp$v)), main=main,
        e_ind=grp$e[,1]%in%grp$v$name[pp_] & grp$e[,2]%in%grp$v$name[pp_],
        label_ind=rep(F,nrow(grp$v)))
      gp = gpp + geom_label_repel(
        data=grp$v[label_ind,],
        aes(x=x,y=y,label=label, color=color),
        nudge_x = -.2, direction = "y", hjust = 1, segment.size = 0.2) +
        geom_point(data=grp$v[ptrp_|ptep_,],aes(x=x,y=y),size=2,color="grey") 
      
      ggsave(paste0(pathn,"/grpos",classn,"_",tbl$feat[i],".png"), plot=gp, scale = 1, width =11, height =15, units = "in", dpi = 300, limitsize = TRUE)
      
      
    }
  }
# },.parallel=T)
time_output(start1)







