## input: features 
## output: p values and their corelation between train/test sets; p value qq plots, recall/precision, number of significant features

set.seed(10)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr(c("stringr","colorspace", "Matrix", "plyr",
       "lattice", "gridExtra", "plotly", "RColorBrewer", "plotrix", "ggrepel", "ggplot2",# libr(proxy)
       "metap",
       "foreach","doMC",
       "kernlab"))

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
          
          set1 = -log(pvs$train); set1[set1==Inf] = max(max(set1[set1!=Inf]),ptl+1)
          set2 = -log(pvs$test); set2[set2==Inf] = max(max(set2[set2!=Inf]),ptl+1)
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
        for (i in 1:length(ys)) 
          points(nqs[[i]], sort(ys[[i]]), pch=16, cex=.5, col=cs[i])
        abline(h=pt)
        lines(x = c(-100,100), y = c(-100,100))
        legend("topleft", legend=names(ys), fill=cs, bg="transparent")
        
        nqls = llply(nqs, log)
        yls = llply(1:length(ys), function(xi) {
          a = log(ys[[xi]]); 
          a[a<nqls[[xi]][1]] = nqls[[xi]][1]; 
          return(a) })
        plot(NULL, xlim=c(min(sapply(nqls,function(x)x[1])),0), ylim=c(min(sapply(nqls,function(x)x[1])),0), ylab=paste0(ptype," test + ", adj, " p values"), xlab="theoretical quantile", main="-ln(qq plot)")
        for (i in 1:length(ys)) 
          points(nqls[[i]], sort(yls[[i]]), pch=16, cex=.5, col=cs[i])
        abline(h=-ptl)
        lines(x = c(-100,100), y = c(-100,100))
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
        
        maint = paste0("%/# of sig/features (bar/line) \n data: ", data, "/", uc, ", method: ", ptype,"-",adj)
        
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
  pvalgr_dir = paste0(root,"/result/",data,"/pval/graph")
  gr_dirs = list.files(pvalgr_dir, full.names=T)
  grs = llply(gr_dirs, function(x) get(load(x)))
  names(grs) = gsub(".Rdata","",fileNames(gr_dirs))
  for (uc in names(pvals[[data]][[1]])) {
    classn = paste0(uc,"_")
    if (length(pvals[[data]][[1]])==1) classn = ""
    for (ptype in c_ptypes) {
      for (adj in c_adjs) {
        mv = mvv = c(rep(1,2),rep(2,2),c(3:(2+max(sapply(c_featnes, function(x) sum(grepl(x,c_feats)))))))
        if (length(c_featnes)>0) {
          for (fni in 2:length(c_featnes)) {
            mv = cbind(mv,mvv+max(mv))
          }
        }
        # mv = matrix(c(1,1,2,2,3,4,5,6,7,8,8,9,9,10,11,12,13,14),ncol=2,byrow=F)
        png(paste0(root,"/pval/",data,"/gr_",classn,ptype,"-",adj,".png"), height=nrow(mv)*250, width=ncol(mv)*500)
        layout(mv) # scatterplot + histograms
        par(cex=1.2)
        for (featne in c_featnes) {
          grt = grclusthist = grhubhist = NULL
          for (feat in c_feats[grepl(featne,c_feats)]) {
            grname = paste0(feat,"_",uc)
            if (!grname%in%names(grs)) next
            
            # build graph of sig only
            gr0 = graph_from_data_frame(d=grs[[grname]]$e, vertices=grs[[grname]]$v, directed=T)
            all_sig = pvals[[data]][[feat]][[uc]]$p[[ptype]][[adj]]$all<pt
            all_sig_ = names(all_sig)[all_sig]
            if (featne%in%c("cell","group")) {
              gr = gr0 - setdiff(V(gr0)$name, all_sig_)
            } else {
              gr_e = as.data.frame(Reduce("rbind",str_split(all_sig_,"_")))
              colnames(gr_e) = c("from","to")
              gr = graph_from_data_frame(gr_e)
            }

           # get graph metrics
            con = components(gr) # membership (cluster id/feat), csize (cluster sizes), no (of clusers)
            grd = decompose.graph(gr)
            vertex_connectivity(grd[[1]])
            
            
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
        
        for (featne in c_featnes) {
          for (feat in c_feats[grepl(featne,c_feats)]) {
            
            grname = paste0(feat,"_",uc)
            if (!grname%in%names(grs)) next
            gr_e = grs[[grname]]$e
            gr_v = grs[[grname]]$v
            
            # colour palette
            palblue = colorRampPalette(brewer.pal(3, "Blues"))
            
            # edit layout manually
            gr_vxy_ = layout.reingold.tilford(gr = graph_from_data_frame(gr_e)) # layout.circle
            gr_vxy = as.data.frame(gr_vxy_)
            gys = sort(unique(gr_vxy[,2]))
            gxns = sapply(gys,function(y) length(gr_vxy[gr_vxy[,2]==y,1]))
            names(gxns) = gys
            gxnmax = max(gxns)
            gxnmaxl = which(gxns==gxnmax)
            minwidthmid = 2
            maxwidthmid = 4
            gxnmaxwidth = gxnmax*minwidthmid-1
            gxos = unlist(llply(gys, function(gy) {
              gxtf = which(gr_vxy[,2]==gy)
              gx = gr_vxy[gxtf,1]
              gxtf[order(gx)]
            }))
            gr_vxy[gxos,1] = unlist(llply(1:length(gys), function(gyi) {
              if (gyi%in%gxnmaxl) return( seq(0,gxnmaxwidth,by=minwidthmid) )
              if (gxns[gyi]==1) return( gxnmaxwidth/2 )
              by = min(maxwidthmid,(gxnmaxwidth+1)/(gxns[gyi]+1))
              a = seq(0,gxns[gyi]-1)*by
              a + (gxnmaxwidth-(a[length(a)]-1))/2
            }))
            
            # get edge
            gr_e$from.x <- gr_vxy$x[match(gr_e$from, gr_v$name)]
            gr_e$from.y <- gr_vxy$y[match(gr_e$from, gr_v$name)]
            gr_e$to.x <- gr_vxy$x[match(gr_e$to, gr_v$name)]
            gr_e$to.y <- gr_vxy$y[match(gr_e$to, gr_v$name)]
            
            if (featne%in%c("cell","group")) {
              # gr = delete.vertices(gr, V(gr)[degree(gr)<3]) # exclude low degree from graph
              # V(gr)$color = ifelse(V(gr)$name=='CA', 'blue', 'red')
              # V(gr)$size = degree(gr)/10
              cols = as.numeric(cut(append(V(gr)$mean_uc,V(gr)$mean_ctrl),breaks=breaks))
              
              ggplot() +
                geom_segment(data=gr_e,aes(x=from.x,xend = to.x, y=from.y,yend = to.y,size=weight),colour=ifelse(gr_e$p==0, "gray80", "grey")) +
                geom_point(data=gr_v,aes(x=x,y=y),size=8,colour="black") +  # adds a black border around the nodes
                geom_point(data=gr_v,aes(x=x,y=y),size=5,colour="lightgrey") +
                geom_text(data=gr_v,aes(x=x,y=y,label=name)) + # add the node labels
                scale_x_continuous(expand=c(0,1))+  # expand the x limits 
                scale_y_continuous(expand=c(0,1))+ # expand the y limits
                theme_bw()+  # use the ggplot black and white theme
                theme(
                  axis.text.x = element_blank(),  # remove x-axis text
                  axis.text.y = element_blank(), # remove y-axis text
                  axis.ticks = element_blank(),  # remove axis ticks
                  axis.title.x = element_blank(), # remove x-axis labels
                  axis.title.y = element_blank(), # remove y-axis labels
                  panel.background = element_blank(), 
                  panel.border =element_blank(), 
                  panel.grid.major = element_blank(),  #remove major-grid labels
                  panel.grid.minor = element_blank(),  #remove minor-grid labels
                  plot.background = element_blank())
              
              V(gr)$color = gr_vxy$color = palblue(breaks)[cols[1:(length(cols)/2)]]
              V(gr)$frame.color = gr_vxy$frame.color = palblue(breaks)[cols[(length(cols)/2+1):length(cols)]]
              V(gr)$size = gr_vxy$size = 5#ifelse(V(gr)$p<pt,3,1)
              # V(gr)$vertex.shape="none"
              V(gr)$label = gr_vxy$label = V(gr)$name
              V(gr)$label.dist = 1
              V(gr)$label.color = gr_vxy$label.color = ifelse(V(gr)$p<pt,"black","grey")
              V(gr)$label.font = 2
              V(gr)$label.cex = gr_vxy$label.cex = 1
              E(gr)$color = ifelse(E(gr)$p==0, "gray80", "grey")
              # E(net)$width <- E(net)$weight/6
              E(gr)$arrow.size <- .5
              # E(net)$width <- 1+E(net)$weight/12
              
              par(mai=c(0,0,1,0)) #this specifies the size of the margins. the default settings leave too much free space on all sides (if no axes are printed)
              plot(gr, layout=grlo,#, layout.circle
                   main=paste0(data," data; class=",uc, ", feature=",feat,", method=",ptype,".",adj)
              )
              par(mar=c(3,1,1,1))
              image.scale(cols, col=palblue(length(breaks)-1), breaks=breaks, horiz=TRUE)
              
            } else if (featne=="edge") {
              
            }

# Save and export the plot. The plot can be copied as a metafile to the clipboard, or it can be saved as a pdf or png (and other formats).
# For example, we can save it as a png:
png(filename="org_network.png", height=800, width=600) #call the png writer
#run the plot
dev.off() #dont forget to close the device
#And that's the end for now.





          }
        }
        
      }
    }
  }
}
time_output(start1,"graph plots")


