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
tbl = tbl[!grepl("raw|edge|group",tbl$feat),]
tbl$pmethod_adj = paste(tbl$test, tbl$test_adj, sep="-")

pt = table$p_thres[1]
ptl = -log(table$p_thres[1])

# pars
c_datas = unique(tbl$data)
c_feats = unique(tbl$feat)
c_pas = unique(tbl$pmethod_adj)
c_pts = unique(tbl$pthres)



## default counts
mcms = llply(c_datas, function(data) {
  mc0 = get(load(paste0(root,"/result/",data,"/feat_mean/", feat_count,".Rdata")))
  mc = mc0$control
  mcm = mc - min(mc)
  2*mcm/max(mcm)
}, .parallel=T)
mcms_ = llply(c_datas, function(data) {
  mc0 = get(load(paste0(root,"/result/",data,"/feat_mean/", feat_count,".Rdata")))
  llply(mc0, function(mc) {
    mcm = mc - min(mc)
    mcm/max(mcm)
  })
}, .parallel=T)
names(mcms) = names(mcms_) = c_datas

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

## histograms + scatterplots

# density plot
plot_int = function(dat, pch='.', ...) {
  colPalette = colorRampPalette(c("blue", "green", "yellow", "red"))
  col = densCols(dat, colramp=colPalette)
  graphics::plot(dat, col=col, pch=pch, ...)
}

# graph base
gpbase = ggblank()

for (i in 1:nrow(tbl)) {
  pt = tbl$pthres[i]
  ptl = -log(pt)
  classn = ifelse(tbl$class[i]=="exp", "", paste0("_",tbl$class[i]))
  pathn = paste0(root,"/result/",tbl$data[i],"/plot_pval/pthres-",pt,"/",tbl$pmethod_adj[i])
  pv = get(load(tbl$pathp[i]))
  
  
  
  ## histograms
  png(paste0(pathn,"/hist",classn,"_",tbl$feat[i],".png"), height=800, width=800)
  mv1 = matrix(c(1,1,2,3,4,4,5,6),ncol=2,byrow=F)
  layout(mv1) # scatterplot + histograms
  
  mcm = mcms[[tbl$feat[i]]]
  
  set1 = pv$train
  set2 = pv$test
  
  if (grepl("-cell-",tbl$feat[i])) 
    mcmind = match(names(set1),names(mcm))
  if (!grepl("-cell-",tbl$feat[i])) 
    mcmind = match(sapply(strsplit(names(set1),"_"), function(x) x[1]),names(mcm))
  mcm1 = mcm[mcmind]
  
  plot_int(cbind(set1,set2), pch=1, xlim=c(0,1), ylim=c(0,1), 
           xlab="set 1 p values", ylab="set 2 p values", main=paste0("set2 class: ",tbl$class[i], "; feature: ", tbl$feat[i], "; method: ", tbl$pmethod_adj[i]), cex=mcm1)
  abline(v=pt,h=pt, col="red")
  
  hist(set1, freq=F, main=paste0("set 1"), xlim=c(0,1), breaks=seq(0,1,.05), xlab="p values", ylab="histogram", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col=rgb(0,0,0,.5))
  abline(v=pt, col="red")
  
  hist(set2, freq=F, main="set 2", xlim=c(0,1), breaks=seq(0,1,.05), xlab="p values", ylab="histogram", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col=rgb(0,0,0,.5))
  abline(v=pt, col="red")
  
  
  set1 = -log(set1); 
  set1[set1==Inf] = max(max(set1[set1!=Inf]),ptl+1)
  set2 = -log(set2); 
  set2[set2==Inf] = max(max(set2[set2!=Inf]),ptl+1)
  plot_int(cbind(set1,set2), pch=1, xlim=c(0,10), ylim=c(0,10), 
           xlab="set 1 -ln(p values)", ylab="set 2 -ln(p values)", main=paste0("set2 class: ",tbl$class[i], "; feature: ", tbl$feat[i], "; method: ", tbl$pmethod_adj[i]), cex=mcm1)
  abline(v=ptl,h=ptl, col="red")
  
  hist(set1, freq=F, main=paste0("set 1"), xlim=c(0,max(4,max(set1))), xlab="-ln(p values)", ylab="histogram", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col=rgb(0,0,0,.5))
  abline(v=ptl, col="red")
  
  hist(set2, freq=F, main="set 2", xlim=c(0,max(4,max(set2))), xlab="-ln(p values)", ylab="histogram", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col=rgb(0,0,0,.5))
  abline(v=ptl, col="red")
  graphics.off()
  
  
  
  
  
  ## graph stats
  # do only for cell feature types and nodes<1000 see below
  if(!tbl$data[i]%in%names(grp0s)) next
  gr0 = gr0s[[tbl$data[i]]]
  grp0 = grp0s[[tbl$data[i]]]
  
  mcm = mcms_[[tbl$data[i]]]
  
  p = pv$all
  p_ = p<pt
  pp = p[!grepl("[-]",names(p))]
  pp_ = pp<pt
  
  mcmc = mcm$control
  mcmc = mcmc[match(names(p),names(mcmc))]
  mcmcp = mcmc[match(names(pp),names(mcmc))]
  mcmu = mcm[[tbl$class[i]]]
  mcmu = mcmu[match(names(p),names(mcmu))]
  mcmup = mcmu[match(names(pp),names(mcmu))]

  gr0$v = gr0$v[match(names(p),gr0$v$name),,drop=F]
  gr0$v = data.frame(name=unlist(gr0$v))
  gr0$e = gr0$e[gr0$e$from%in%gr0$v$name & gr0$e$to%in%gr0$v$name,]
  grp0$v = grp0$v[match(names(pp),grp0$v$name),,drop=F]
  grp0$v = data.frame(name=unlist(grp0$v))
  grp0$e = grp0$e[grp0$e$from%in%grp0$v$name & grp0$e$to%in%grp0$v$name,]
  
  gr = layout_gr(gr0$e,gr0$v)
  gr = gpdf(gr)
  grp = layout_gr(grp0$e,grp0$v)
  grp = gpdf(grp)
  
  gr$v$label = paste0(gr$v$name,":",round(mcmc,3))
  gr$v$color = mcmc
  gr$v$colorb = mcmu
  grp$v$label = paste0(grp$v$name,":",round(mcmcp,3))
  grp$v$color = mcmcp
  grp$v$colorb = mcmup
  
  gr$v$sizeb = grp$v$sizeb = 5
  gr$v$size = grp$v$size = 3
  
  main = paste0("data: ",tbl$data[i],"; feat: ",tbl$feat[i], "; class: ",tbl$class[i], "; pthres: ",tbl$pthres[i],"\n sig & >.1 ratio change cell pops labelled (border/centre colours = class/control proportion)")
  if (!grepl("short",tbl$feat[i])) {
    gp = gggraph(gr, v_ind=rep(F,nrow(gr$v)), e_ind=rep(F,nrow(gr$e)), label_ind=rep(T,nrow(gr$v)))
  gp = gp + geom_label_repel(
    data=gr$v[(mcmu-mcmc)/mcmc >.1 & p_,],
    aes(x=x,y=y,label=label, color=color), nudge_y = .3) +
    ggtitle(paste0("(positive only) ",main))

  ggsave(paste0(pathn,"/gr",classn,"_",tbl$feat[i],".png"), plot=gp, scale = 1, width =11, height =7.5, units = "in", dpi = 300, limitsize = TRUE)
  }
  
  gpp = gggraph(grp, v_ind=rep(F,nrow(grp$v)), e_ind=rep(F,nrow(grp$e)), label_ind=rep(T,nrow(grp$v)))
  gp = gp + geom_label_repel(
    data=grp$v[(mcmup-mcmcp)/mcmcp >.1 & pp_,],
    aes(x=x,y=y,label=label, color=color), nudge_y = .3) +
    ggtitle(main)

  ggsave(paste0(pathn,"/gr",classn,"_",tbl$feat[i],".png"), plot=gpp, scale = 1, width =11, height =7.5, units = "in", dpi = 300, limitsize = TRUE)
  
  
  
#   c_feats_ = c_feats
#   c_feats_ = c_feats[!grepl("edge|entropy",c_feats_)]
#   
#   
#   for (feat in c_feats_) {
#     feat_ = NULL # contains original values
#     if (grepl("lnpropexpect",feat)) {
#       lno = str_extract(feat,"lnpropexpect")
#       feat_ = gsub(lno,paste0(lno,"_"),feat)
#       
#       grname = paste0(feat_,"_",uc)
#       if (!grname%in%names(grs)) {
#         feat_ = NULL
#       } else {
#         gr_e_ = grs[[grname]]$e
#         gr_v_ = grs[[grname]]$v
#       }
#     }
#     
#     grname = paste0(feat,"_",uc)
#     if (!grname%in%names(grs)) next
#     gr_e = grs[[grname]]$e
#     gr_v = grs[[grname]]$v
#     if (nrow(gr_v)>1000) next
#     
#     all_sig = pvals[[data]][[feat]][[uc]]$p[[ptype]][[adj]]$all<pt
#     all_sig_ = names(all_sig)[all_sig]
#     
#     # colour palette
#     # palblue = colorRampPalette(brewer.pal(3, "Blues"))
#     
#     al = layout_gr(gr_e,gr_v,layout.reingold.tilford)
#     gr_e = al$e
#     gr_v = al$v
#     main = paste0(data," data; class=",uc, ", feature=",feat,", method=",ptype,".",adj)
#     
#     # gr = delete.vertices(gr, V(gr)[degree(gr)<3]) # exclude low degree from graph
#     # V(gr)$color = ifelse(V(gr)$name=='CA', 'blue', 'red')
#     # V(gr)$size = degree(gr)/10
#     # cols = as.numeric(cut(append(gr_v$mean_uc,gr_v$mean_ctrl),breaks=breaks))
#     
#     # sig indices
#     gr_vsig = gr_v$name %in% all_sig_
#     gr_vhigh = str_count(gr_v[,1],"[+-]")<5
#     gr_esig = gr_e[,1] %in% all_sig_ & gr_e[,1] %in% all_sig_
#     
#     gp0 = gpbase + 
#       geom_segment(
#         data=gr_e, 
#         aes(x=from.x,xend = to.x, y=from.y,yend = to.y),
#         colour=ifelse(gr_e[,1]%in%all_sig_ & gr_e[,2]%in%all_sig_, "grey40", "grey")) +
#       geom_point(
#         data=gr_v,aes(x=x,y=y),colour='grey',size=1) +
#       geom_point(
#         data=gr_v[gr_vsig,],   # adds node border
#         aes(x=x,y=y, colour=mean_uc), size=2)
#     # geom_point(
#     #   data=gr_v[gr_vsig,], 
#     #   aes(x=x,y=y, colour=mean_ct), size=2)
#     
#     gp_sigmarker = gp0 + ggtitle(paste0(main,"\nlabel=sig; colour=mean exp value")) +
#       geom_label_repel(
#         data=gr_v[gr_vsig & gr_vhigh,],
#         aes(x=x,y=y,label=name, color=mean_uc), nudge_y = .3)
#     
#     ggsave(paste0(root,"/result/",data,"/plot_pval/",pa,"/grpl",classn,"_",feat,".png"), plot=gp_sigmarker, scale = 1, width =11, height =7.5, units = "in", dpi = 300, limitsize = TRUE)
#     if (!is.null(feat_)) { # more plots, only change label
#       gp_allmean = gp0 + ggtitle(paste0(main,"\nlabel=mean exp value"))
#       if (grepl("short",feat)) {
#         gr_v$label = paste0(gr_v_$name,":",round(gr_v_$mean_uc,3))
#         gp_allmean = gp0 + geom_label_repel(
#           data=gr_v[!gr_vsig,],
#           aes(x=x,y=y,label=label), nudge_y = .3,
#           color="grey") +
#           geom_label_repel(
#             data=gr_v[gr_vsig,], nudge_y = .3,
#             aes(x=x,y=y,label=label, color=mean_uc))
#       } else {
#         gr_v$label = round(gr_v_$mean_uc,3)
#         gp_allmean = gp0 + 
#           geom_label_repel(
#             data=gr_v[gr_vsig & gr_vhigh,], nudge_y = .3,
#             aes(x=x,y=y,label=label, color=label))
#       }
#       ggsave(paste0(root,"/result/",data,"/plot_pval/",pa,"/grpl",classn,"_",feat,"_.png"), plot=gp_allmean, scale = 1, width = 11, height = 7.5, units = c("in"), dpi = 300, limitsize = TRUE)
#       
#     }
#     
# }
# 
# 
# 
# ## feature comparison plots
# for (data in c_datas) {
#   mcm = mcms[[data]]
#   ucs = unique(tbl$class[tbl$data==data])
#   for (uc in ucs) {
#     classn = ifelse(length(ucs)==1, "", paste0("_",uc))
#     for (pa in c_pas) {
#       for (pt in c_pts) { ptl = -log(pt)
#       tb = tbl[tbl$data==data & tbl$class==uc & tbl$pmethod_adj==pa & tbl$pthres==pt,]
#       
#       pvs = llply(tb$pathp, function(x) get(load(x))$all)
#       
#       pathn = paste0(root,"/result/",data,"/plot_pval/pthres-",pt,"/",pa)
#       dir.create(pathn, showWarnings=F, recursive=T)
#       
#       
#       
#       ## qq plot; log both axis
#       nqs = llply(sapply(pvs,length), function(n) 1:n/(n+1))
#       cs = rainbow(length(pvs))
#       
#       mcm1l = llply(pvs, function(x) {
#         if (grepl("_",names(x)[1])) mcmind = match(names(x),names(mcm))
#         if (!grepl("_",names(x)[1])) 
#           mcmind = match(sapply(strsplit(names(x),"_"), function(x) x[1]),names(mcm))
#         mcm[mcmind]
#       })
#       
#       names(pvs) = names(nqs) = names(cs) = tb$feat # feats
#       
#       png(paste0(pathn,"/qqpl",classn,".png"), height=400, width=800)
#       par(mfrow=c(1,2))
#       
#       plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab=paste0(pa, " p values"), xlab="theoretical quantile", main=paste0("qq plot | data: ", data, ", class: ", uc, ", method: ", pa))
#       for (i in 1:length(pvs)) {
#         yo = order(pvs[[i]])
#         points(nqs[[i]], pvs[[i]][yo], pch=1, cex=mcm1l[[i]][yo], col=adjustcolor(cs[i], alpha.f=.5))
#       }
#       abline(h=pt)
#       lines(x = c(-100,100), y = c(-100,100))
#       legend("topleft", legend=names(pvs), fill=cs, bg="transparent")
#       
#       nqls = llply(nqs, log)
#       pvls = llply(pvs, log)
#       plot(NULL, xlim=c(min(sapply(nqls,function(x)x[1])),0), ylim=c(min(sapply(nqls,function(x)x[1])),0), ylab=paste0(pa, " p values"), xlab="theoretical quantile", main="ln(qq plot)")
#       for (i in 1:length(pvs)) {
#         yo = order(pvls[[i]])
#         points(nqls[[i]], pvls[[i]][yo], pch=1, cex=mcm1l[[i]][yo], col=adjustcolor(cs[i], alpha.f=.5))
#       }
#       abline(h=-ptl)
#       lines(x = c(-100,100), y = c(-100,100))
#       legend("topleft", legend=names(pvs), fill=cs, bg="transparent")
#       
#       graphics.off()
#       
#       
#       
#       ## number of significant features
#       plens = laply(pvs, length)
#       npts = laply(pvs, function(pv) {
#         a = sum(pv<pt); ifelse (a==0,0,a/length(pv))
#       })
#       names(pvs) = names(npts) = names(plens) = tb$feat
#       
#       png(paste0(pathn,"/num",classn,".png"), height=400, width=400)
#       par(mar=c(10,4,5,4))
#       
#       maint = paste0("%/# of sig/features (bar/line) \n data: ", data, "/", uc, ", method: ", pa)
#       
#       barplot(npts, ylim=c(0,1), xlab="", las=2, ylab="% significant features", col=rgb(0,0,0,.5), main=maint)
#       
#       abline(h=pt)
#       par(new=T)
#       
#       plot(plens, pch=15, xlab="",ylab="",axes=F,type="b", col="red", ylim=c(0, max(plens)+(.15*max(plens))))
#       mtext("number of features",side=4,col="red",line=4)
#       axis(4, col="red", col.axis="red", las=1)
#       
#       graphics.off()
#       
#       
#       
#       ## robustness stats
#       fr = adply(tb,1, function(x) {
#         ldply(c("rec","prec","f","pcorr"), function(y) 
#           data.frame(x,metric=y,score=as.numeric(x[[y]])) )
#       })
#       
#       png(paste0(pathn,"/rpfc",classn,".png"), height=400, width=400)
#       par(mar=c(10,4,5,4))
#       
#       pl = barchart(score~feat, fr, groups=metric, las=1,# cex.axis=3, 
#                     ylim=c(0,1),
#                     auto.key = list(columns=2),
#                     scales=list(x=list(rot=90,cex=0.8)),
#                     main=paste0("data: ", data, "/", uc, "\ntests: f measure, precision, recall, pearson correlation"))
#       print(pl)
#       
#       graphics.off()
#       
#       
#       
#       
#       ## pos stats
#       if(grepl("pos[0-9]",data)) {
#         asig = llply(pvs, function(x) x<pt)
#         rsig = llply(pvs, function(x) {
#           if (data=="pos1") a = grepl("A[+|-]",names(x))
#           if (data=="pos2") a = grepl("[A|B][+|-]",names(x))
#           if (data=="pos3") a = grepl("A[+|-]B[+|-]C[+|-]",names(x))
#           if (data=="pos4") a = grepl("[A|B|C|D][+|-]",names(x))
#           return(a)
#         })
#         hsig = llply(pvs, function(x) {
#           if (data=="pos1") a = grepl("A",names(x))
#           if (data=="pos2") a = grepl("A|B",names(x))
#           if (data=="pos3") a = grepl("A[+|-]B[+|-]C[+|-][D-Z|+|-]+",names(x))
#           if (data=="pos4") a = grepl("A|B|C|D",names(x))
#           return(a)
#         })
#         hsig = llply(1:length(hsig), function(i) hsig[[i]] & !rsig[[i]])
#         names(pvs) = names(asig) = names(rsig) = names(hsig) = tb$feat
#         
#         
#         for (t in c("real","influenced")) {
#           sig = switch(t, real = rsig, influenced = hsig)
#           sbs = ldply(tb$feat, function(i) {
#             overlap = sum(asig[[i]] & sig[[i]])
#             if (overlap==0) {
#               rec = prec = f = 0
#             } else {
#               rec = overlap/sum(asig[[i]])
#               prec = overlap/sum(sig[[i]])
#               f = 2*((prec*rec)/(prec+rec))
#             }
#             return(data.frame(feature=paste0(i," (",sum(asig[[i]]),",",overlap,",",sum(sig[[i]]),")"), type=t, 
#                               sig=sum(sig[[i]]), actual=sum(asig[[i]]),
#                               metric=c("recall","precision","f"),
#                               score=c(rec,prec,f)))
#           })
#           
#           png(paste0(pathn,"/pos-",t,classn,".png"), height=600, width=400)
#           pl = barchart(score~feature, sbs, groups=metric, las=1, ylim=c(0,1),
#                         main="feature (actual,overlap,hypothetical) \noverlap / actual (p) vs hypothetical (r)",
#                         auto.key = list(space = "top"),
#                         scales = list(x = list(rot = 45)))
#           print(pl)
#           graphics.off()
#         }
#         
#       }
#     }
#   }
}





