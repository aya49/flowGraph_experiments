## input: features 
## output: p values and their corelation between train/test sets; p value qq plots, recall/precision, number of significant features

set.seed(10)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr(c("stringr","colorspace", "Matrix", "plyr",
       "lattice", "gridExtra", "plotly", "RColorBrewer", "plotrix",# libr(proxy)
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

pthres = .05 # p value sig threshold for t test
good_count = 5
# minfold = 5 # minimum number of samples in each fold/class
good_sample = minfold = 5 # each class must have more thatn good_sample amount of samples or else prune
# countThres = 1000 #insignificant if count under

id_col = "id"
target_col = "class"
order_cols = NULL

control = "control"
feat_count = "file-cell-countAdj"


table = get(load(paste0(root,"/pval_table.Rdata")))
pvals = get(load(paste0(root,"/pval.Rdata")))

tbl = table
tbl$pmethod_adj = paste(table$sig_test, table$adjust.combine)
tbl = tbl[!grepl("group",tbl$feature) & tbl$feature!="file-cell-count",]
pt = table$p_thres[1]
ptl = -log(table$p_thres[1])

# data and uc together, ptype and adj together, 
c_datas = names(pvals)
c_feats = unique(tbl$feature)
# uc, p
c_ptypes = unique(tbl$sig_test)
c_adjs = unique(tbl$adjust.combine)

c_featnes = c("cell","edge")

plot_int = function(dat, pch='.', ...) {
  colPalette = colorRampPalette(c("blue", "green", "yellow", "red"))
  col = densCols(dat, colramp=colPalette)
  graphics::plot(dat, col=col, pch=pch, ...)
}


## histograms + scatterplots
for (data in c_datas) {
  # png(paste0(pvalr_dir,"/hist_",data,".png"), height=length(pvals[[data]])*300, width=length(pvals[[data]][[1]])*length(c_ptypes)*2*300)
  # mv1 = matrix(c(1,1,2,3),nrow=2,byrow=F)
  # for (pi in 2:length(c_ptypes))
  #   mv1 = cbind(mv1, matrix(c(1,1,2,3)+max(mv1),nrow=2,byrow=F))
  # mv = mv1
  # for (fi in 1:length(pvals[[data]]))
  #   mv = rbind(mv, mv1+max(mv))
  # layout(mv) # scatterplot + histograms
  
  mc0 = Matrix(get(load(paste0(root,"/result/",data,"/feat/", feat_count,".Rdata"))))
  mc0 = mc0 - min(mc0)
  mcm = colMeans(as.matrix(mc0))
  mcm = mcm/max(mcm)
  mcm = mcm*2

  dir.create(paste0(root,"/pval/",data), showWarnings=F, recursive=T)
  # dir.create(paste0(root,"/pval/",data), showWarnings=F, recursive=T)
  for (uc in names(pvals[[data]][[1]])) {
    classn = paste0(uc,"_")
    if (length(pvals[[data]][[1]])==1) classn = ""
    for (feat in c_feats) {
      for (ptype in c_ptypes) {
        for (adj in c_adjs) {
          png(paste0(root,"/pval/",data,"/hist_",classn,ptype,"-",adj,"_",feat,".png"), height=400, width=800)
          mv1 = matrix(c(1,1,2,3),nrow=2,byrow=F)
          layout(mv1) # scatterplot + histograms
          
          pvs = pvals[[data]][[feat]][[uc]]$p[[ptype]][[adj]]
          
          set1 = -log(pvs$train)
          set2 = -log(pvs$test)
          # set1[set1> 10] = set2[set2> 10] = 10
          
          if (grepl("-cell-",feat)) mcmind = match(names(set1),names(mcm))
          if (!grepl("-cell-",feat)) mcmind = match(sapply(strsplit(names(set1),"_"), function(x) x[1]),names(mcm))
          mcm1 = mcm[mcmind]
          plot_int(cbind(set1,set2), pch=1, xlim=c(0,10), ylim=c(0,10), 
                   xlab="set 1 -ln(p values)", ylab="set 2 -ln(p values)", main=paste0("set2 class: ",uc, "; feature: ", feat, "; method: ", ptype,"/",adj), cex=mcm)
          abline(v=ptl,h=ptl, col="red")
          
          hist(set1, freq=F, main=paste0("set 1"), xlim=c(0,max(4,max(set1))), xlab="-ln(p values)", ylab="histogram", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col=rgb(0,0,0,.5))
          abline(v=ptl, col="red")
          
          hist(set2, freq=F, main="set 2", xlim=c(0,max(4,max(set2))), xlab="-ln(p values)", ylab="histogram", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col=rgb(0,0,0,.5))
          abline(v=ptl, col="red")
          graphics.off()
          
        }
      }
    }
  }
  # graphics.off()
}


## qq plot; log both axis
for (data in c_datas) {
  for (uc in names(pvals[[data]][[1]])) {
    classn = paste0(uc,"_")
    if (length(pvals[[data]][[1]])==1) classn = ""
    for (ptype in c_ptypes) {
      for (adj in c_adjs) {
        # n =100
        # y = runif(n)
        ys = llply(c_feats, function(feat)
          pvals[[data]][[feat]][[uc]]$p[[ptype]][[adj]]$all)
        nqs = llply(sapply(ys,length), function(n) 1:n/(n+1))
        cs = rainbow(length(ys))
        names(ys) = names(nqs) = names(cs) = c_feats # feats
        
        png(paste0(root,"/pval/",data,"/qqpl_",classn,ptype,"-",adj,".png"), height=400, width=800)
        par(mfrow=c(1,2))
        
        plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab=paste0(ptype," test + ", adj, " p values"), xlab="theoretical quantile", main=paste0("qq plot | data: ", data, ", class: ", uc, ", method: ", ptype,"-",adj))
        for (i in 1:length(ys)) {
          points(nqs[[i]], sort(ys[[i]]), pch=16, cex=.1, col=cs[i])
        }
        abline(h=pt, col="red")
        legend("topleft", legend=names(ys), fill=cs, bg="transparent")
        
        
        yls = llply(ys, log)
        nqls = llply(nqs, log)
        plot(NULL, xlim=c(min(unlist(nqls)),0), ylim=c(min(unlist(yls)),0), ylab=paste0(ptype," test + ", adj, " p values"), xlab="theoretical quantile", main="ln(qq plot)")
        for (i in 1:length(ys)) {
          # seq(0.05,1,length(y))
          # for(i in 1:1000)
          #   points(sort(runif(length(y))),sort(y),pch=".")
          set1 = -log(nqs[[i]])
          set2 = -log(ys[[i]])
          # set1[set1< -10] = set2[set2< -10] = -10
          points(set1, sort(set2), pch=16, cex=.5, col=cs[i])#,col="red",pch=16)
        }
        abline(h=pt, col="red")
        abline(h=ptl, col="red")
        legend("topleft", legend=names(ys), fill=cs, bg="transparent")
        
        graphics.off()
        
      }
    }
  }
}


## number of significant features
## other metrics
for (data in c_datas) {
  for (uc in names(pvals[[data]][[1]])) {
    classn = paste0(uc,"_")
    if (length(pvals[[data]][[1]])==1) classn = ""
    for (ptype in c_ptypes) {
      for (adj in c_adjs) {
        ## number of significant features
        pvs = llply(c_feats, function(feat)
          pvals[[data]][[feat]][[uc]]$p[[ptype]][[adj]]$all)
        npts = laply(pvs, function(pv) {
          a = sum(pv<pt)
          if (a==0) return(0)
          return(a/length(pv))
        })
        plens = laply(pvs, length)
        names(pvs) = names(npts) = names(plens) = c_feats

        png(paste0(root,"/pval/",data,"/num_",classn,ptype,"-",adj,".png"), height=400, width=400)
        par(mar=c(10,4,5,4))
        
        maint = paste0("# of sig/features | data: ", data, "/", uc, ", method: ", ptype,"-",adj)
        
        barplot(npts, ylim=c(0,max(.5,max(npts))), xlab="", las=2, ylab="% significant features", col=rgb(0,0,0,.5), main=maint)
        
        abline(h=pt, col="green")
        par(new=T)
        
        plot(plens, pch=15, xlab="",ylab="",axes=F,type="b", col="red", ylim=c(0, max(plens)))
        mtext("number of features",side=4,col="red",line=4)
        axis(4, col="red", col.axis="red", las=1)
        
        graphics.off()
        
        
        if (data=="pos") {
          markers = unique(unlist(str_split(names(pvs[[1]]),"[+-]")))[-1]
          for (i in 1:(length(markers)-1)) {
            for (j in (i+1):length(markers)) {
              m1 = markers[i]
              m2 = markers[j]
              sbs = ldply(c_feats, function(pvi) {
                pv = pvs[[pvi]]
                ssig = grepl(m1,names(pv)) | grepl(m2,names(pv))
                # a = unlist(strsplit(names(pv),"[-]|[+]"))
                asig = pv<pt
                # data.frame(metric=c("recall","precision"), score=c(sum(ssig==asig)/sum(ssig), sum(ssig==asig)/sum(asig)), feature=pvi)
                data.frame(metric=c("recall","precision"), score=c(sum(ssig & asig)/sum(ssig), sum(ssig & asig)/sum(asig)), feature=pvi)
              })
              png(paste0(root,"/pval/",data,"/numpos_",m1,"-",m2,"_",classn,ptype,"-",adj,".png"), height=400, width=400)
              pl = barchart(score~feature, data=sbs, groups=metric, las=1,
                            main="overlap / hypothetical (r) vs actual (p) actual sig",
                            auto.key = list(space = "top"),
                            scales = list(x = list(rot = 45)))
              print(pl)
              graphics.off()
              
            }
          }
        }
        

        
        ## other stats
        fr = tbl[tbl$data==data & tbl$class==uc & tbl$sig_test==ptype & tbl$adjust.combine==adj,]
        
        pls = fr$recall
        names(pls) = fr$feature
        png(paste0(root,"/pval/",data,"/r_",classn,ptype,"-",adj,".png"))
        par(mar=c(10,4,5,4))
        barplot(pls, las=2, main=paste0("recall | data: ", data, "/", uc, ", method: ", ptype,"-",adj), ylim=c(0,1), ylab="recall")
        graphics.off()
        
        pls = fr$precision
        names(pls) = fr$feature
        png(paste0(root,"/pval/",data,"/p_",classn,ptype,"-",adj,".png"))
        par(mar=c(10,4,5,4))
        barplot(pls, las=2, main=paste0("precision | data: ", data, "/", uc, ", method: ", ptype,"-",adj), ylim=c(0,1), ylab="precision")
        graphics.off()
        
        pls = fr$f
        names(pls) = fr$feature
        png(paste0(root,"/pval/",data,"/f_",classn,ptype,"-",adj,".png"))
        par(mar=c(10,4,5,4))
        barplot(pls, las=2, main=paste0("f | data: ", data, "/", uc, ", method: ", ptype,"-",adj), ylim=c(0,1), ylab="f measure")
        graphics.off()
        
        pls = fr$corr_spear
        names(pls) = fr$feature
        png(paste0(root,"/pval/",data,"/c_",classn,ptype,"-",adj,".png"))
        par(mar=c(10,4,5,4))
        barplot(pls, las=2, main=paste0("spearman corr | data: ", data, "/", uc, ", method: ", ptype,"-",adj), ylim=c(0,1), ylab="spearman correlation")
        graphics.off()
      }
    }
  }
}


## graph stats
start1 = Sys.time()
for (data in c_datas) {
  start2 = Sys.time()
  result_dir = paste0(root,"/result/",data)
  
  # make base graph
  meta_dir = paste0(result_dir,"/meta")
  meta_cell_childpn_names_dir = paste(meta_dir, "/cell_childpn_names",sep="") #specifies a phenotypes children and splits them into +/- (only for when both -/+ exists)
  meta_cell_childpn_names = get(load(paste0(meta_cell_childpn_names_dir, ".Rdata")))
  meta_cell_childpn_names_ = ldply(names(meta_cell_childpn_names), function(x) {
    to = unlist(meta_cell_childpn_names[[x]])
    data.frame(from=rep(x,length(to)), to=to) 
  })
  gr0 = graph_from_edgelist(as.matrix(meta_cell_childpn_names_))
  gr0_v = names(V(gr0)[[]])
  # meta_cell_parent_names = get(load(paste0(meta_cell_parent_names_dir, ".Rdata")))
  
  for (uc in names(pvals[[data]][[1]])) {
    classn = paste0(uc,"_")
    if (length(pvals[[data]][[1]])==1) classn = ""
    for (ptype in c_ptypes) {
      for (adj in c_adjs) {
        png(paste0(root,"/pval/",data,"/gr_",classn,ptype,"-",adj,".png"), height=10*250, width=2*500)
        mv = matrix(c(1,1,2,2,3,4,5,6,7,8,8,9,9,10,11,12,13,14),ncol=2,byrow=F)
        layout(mv) # scatterplot + histograms
        par(cex=1.2)
        for (featne in c_featnes) {
          grt = grclusthist = grhubhist = NULL
          for (feat in c_feats[grepl(featne,c_feats)]) {
            
            ## make graph
            pv_all_ = pvals[[data]][[feat]][[uc]]$p[[ptype]][[adj]]$all
            all_sig = pv_all_<pthres
            all_sig_ = names(all_sig)[all_sig]
            gr = NULL
            if (length(all_sig_)>2) {
              # make edge list & graph
              etf = grepl("_",all_sig_[1])
              if (etf) {
                elist = as.data.frame(Reduce("rbind",str_split(all_sig_,"_")))
                colnames(elist) = c("from","to")
                gr = graph_from_data_frame(elist)
              } else {
                gr = gr0 - setdiff(gr0_v, all_sig_)
              }
            }
            if (is.null(gr)) next
            
            con = components(gr) # membership (cluster id/feat), csize (cluster sizes), no (of clusers)
            grs = decompose.graph(gr)
            vertex_connectivity(grs[[1]])
            
            
            # page_rank
            # authority_score: principal eigenvector of t(A)*A
            # eigen_centrality: first eigenvector of A (vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on))
            grhubhist[[feat]] = hub_score(gr)$vector # A*t(A); vector (Scores), value (eigenvalue of principle eigenvector)
            grclusthist[[feat]] = component_distribution(gr) # histogram for the maximal connected component sizes
            grt[feat] = con$no
          }
          if (is.null(grt)) next
          
          ## plot
          par(mar=c(10,5,3,3))
          barplot(grt,las=2, main="number of connected components")
          clusthistmax = max(sapply(grclusthist, function(x) max(which(x>0))))
          
          par(mar=c(5,5,3,3))
          plot(NULL,ylim=c(0,1.5), xlim=c(0,clusthistmax), ylab="% of cc", xlab="size of cc", main="size of connected comp cc (0=#sig/nodes)")
          cs = rainbow(length(grclusthist))
          for (i in 1:length(grclusthist)) {
            lines(grclusthist[[i]], col=cs[i])
            points(length(grhubhist[[i]]),1,col=cs[i])
          }
          legend("topleft", legend=names(grclusthist), fill=cs, bg="transparent")
          
          # get percentage histogram
          for (i in 1:length(grhubhist)) {
            h = hist(grhubhist[[i]],plot=FALSE)
            h$density = h$counts/sum(h$counts)
            plot(h,freq=F,ylim=c(0,1), ylab="% frequency of features", xlab="kleinberg's hub centrality score", main=paste0(names(grhubhist)[i], "\nhistogram of hubness scores"))
          }
        }
        graphics.off()
      }
    }
  }
  time_output(start2,data)
}
time_output(start1,"graph plots")


