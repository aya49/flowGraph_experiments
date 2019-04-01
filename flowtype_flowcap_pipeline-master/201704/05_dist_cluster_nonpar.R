# Use distance metrics to cluster tubes
# aya43@sfu.ca 20170419

root = "~/projects/flowCAP-II"
result_dir = "result"
setwd(root)

options(stringsAsFactors=FALSE)
#options(device="cairo")
options(na.rm=T)

#Input
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
dist_dir = paste(result_dir, "/dist", sep="")
plot_dir = paste(result_dir, "/plots", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dir[i])) }
plot_dist_dir = paste(plot_dir, "/dist", sep=""); for(i in 1:length(plot_dist_dir)) { suppressWarnings(dir.create(plot_dist_dir[i])) }
plot_dist_source_dir = paste0(plot_dist_dir, "/source"); for(i in 1:length(plot_dist_source_dir)) { suppressWarnings(dir.create(plot_dist_source_dir[i])) }

doUnderflow = T #if numbers used to calculate scores --> 0, do work-around (slow)
KOonly = F #delete all WT samples
adjustD = T #divide by a number if the numbers are too big e.g. manhattan
pamtries = 10
overwrite = T #redo and overwrite all past scores

width = 700; height = 700 # of each plot in png
rowplot = 3

rwThres = c(0,.05,.1,.2,.5,.75)

scoretype = 1 #from adjustedRand, choose 1 only
scorename = c("AdjustedRandInd")


ignoredist = ".csv"


#Output
dist_score_dir = paste(dist_dir, "/dist_score", sep=""); for (i in 1:length(dist_score_dir)) { suppressWarnings(dir.create(dist_score_dir[i])) }
cl_score_result_dir = paste0(dist_score_dir,"/score_clusterlv_list")

libr(stringr)
libr(foreach)
libr(doMC)
libr(igraph)
libr(PerfMeas)
libr(clues)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")


interested = c("tube","aml") #sampleMeta columns to plot
splitby = c("none","tube") #match with above

start = Sys.time()

no_cores = 8#detectCores()-1
registerDoMC(no_cores)

# ftWTGT = "normal"

sampleMeta0 = get(load(sampleMeta_dir))
distmfile = list.files(dist_dir, recursive=F, full.names=T, pattern=".Rdata")
delpvalind = grepl(ignoredist,distmfile,ignore.case=T)
if (sum(delpvalind)>0) distmfile = distmfile[!delpvalind]
if (!overwrite) {
  f0 = get(load(paste0(cl_score_result_dir,".Rdata")))
  doneind = which(!is.na(match(names(f0),distmfile)))
  if (length(doneind)>0) distmfile = distmfile[-doneind]
}
distmfilenames = fileNames(distmfile)
distmfilelayer = str_extract(distmfile, "layer[0-9]+")

interestedCols = match(interested,colnames(sampleMeta0))
splitbyCols = match(splitby,colnames(sampleMeta0))

#dimension of plot
nplot = 0
for (i in 1:length(interested)) {
  if (splitby[i]=="none") {
    nplot = nplot + length(rwThres)
  } else {
    nplot = nplot + length(rwThres)*length(unique(sampleMeta0[,splitbyCols[i]]))
  }
}
colplot = ceiling(nplot/rowplot)

#calculate score_nac for each distance matrix & interested classes
#result = list()
#for(i in 1:length(distmfile)) {
f = f0 = foreach(i=order(distmfilelayer)) %dopar% {
  tryCatch({
    #loop.ind = 1:length(distmfile)
    #for(i in 1:length(distmfile)) {
    #result = foreach(i=loop.ind) %dopar% { 
    start1 = Sys.time()
    cat("\n", distmfile[i], ", loading dist matrix", sep="")
    #result[[i]] = list()
    d0 = as.matrix(get(load(distmfile[i])))
    dtsne = get(load(paste0(plot_dist_source_dir,"/tsne_",distmfilenames[i])))
    # if (length(d)==1) { d=d[[1]]; save(d,file=distmfile[i])}
    if (sum(is.na((d0)))>0 | sum((d0)==Inf)>0) return(NULL)
    drows = colnames(as.matrix(d0))
    
    
    #adjust sampleMeta
    dcol = rep(0,ncol(sampleMeta0))
    for (j in 1:ncol(sampleMeta0)) {
      l = sum(!is.na(match(drows,sampleMeta0[,j])))
      if (l>0) dcol[j] = l
    }
    dcol = which.max(dcol)
    if (length(unique(sampleMeta0[,dcol]))==nrow(sampleMeta0)) {
      sampleMeta <- sampleMeta0[match(drows,sampleMeta0[,dcol]),]
      dtsne = dtsne[match(sampleMeta$fileName,rownames(dtsne)),]
      # interestedCols = interestedCols0
      
      #distmfile[i] = paste0("noScore_",distmfile[i])
    } else {
      # interestedCols = 1
      # sampleMeta = as.data.frame(colnames(as.matrix(d)))
      return(NULL)
    }
    
    #adjust distance values
    if (adjustD) {
      d0 = (d0)
      if (min(d0[which(d0>0)])>1 & mean(d0)>10000) d0 = d0/mean(d0)
      d0 = as.dist(d0)
    }
    
    pngname = paste0(plot_dist_dir, "/tsne_", gsub(".Rdata","",fileNames(distmfile[i])),"_lvcluster.png")
    
    png(pngname, width=width*colplot, height=height*rowplot)
    par(mfrow=c(rowplot,colplot))
    
    cat(", clustering ", sep="")
    fm = list()
    for (col in 1:length(interestedCols)) {
      colnam = paste0("interested-",interested[col],"_splitby-",splitby[col])
      fm[[colnam]] = list()
      
      #assuming d is by filename, split by tube first
      d = list()
      sm = list()
      dt = list()
      if (splitby[col]!="none") {
        tubecol = which(colnames(sampleMeta)==splitby[col])
        for (tubei in 1:length(unique(sampleMeta[,tubecol]))) {
          index = which(sampleMeta[,tubecol]==unique(sampleMeta[,tubecol])[tubei])
          d[[tubei]] = as.matrix(d0)[index,index]
          sm[[tubei]] = sampleMeta[index,]
          dt[[tubei]] = dtsne[index,]
        }
        names(d) = names(sm) = names(dt) = unique(sampleMeta[,tubecol])
      } else {
        sm[[1]] = sampleMeta
        d[[1]] = as.matrix(d0)
        dt[[1]] = dtsne
        names(d) = names(sm) = names(dt) = "all"
      }
      
      for (dind in 1:length(d)) {
        la = sm[[dind]][,interestedCols[col]]
        d1 = get_graph(d[[dind]])

        fm[[colnam]][[dind]] = list()
        for (rwt in rwThres) {
          d2 = d1
          tops = quantile(as.vector(d1),rwt)
          d2[d2<tops] = 0
          gr = graph_from_adjacency_matrix(d2, weighted=T,mode='undirected', diag=F)
          
          cl0 = cluster_louvain(gr)$membership
          if (!length(unique(la))>1) {
            fm[[colnam]][[dind]] = rep(1,length(scorename))
            names(fm[[colnam]][[dind]]) = scorename
            cl = cl0
            clul = length(unique(cl0))
            break
          }
          
          #prepare for calculate f measure
          ar0 = adjustedRand(as.integer(factor(la)),as.integer(factor(cl0)))
          fm[[colnam]][[dind]][[as.character(rwt)]] = ar0[scoretype]
          names(fm[[colnam]][[dind]][[as.character(rwt)]]) = scorename
          cl = cl0
          clul = length(unique(cl0))
          
          # if (sum(ar0)>maxsum) {
          #   maxsum = sum(ar0)
          # }
          # }
          
          # gr = graph_from_adjacency_matrix(1/(max(as.matrix(d[[dind]]))-as.matrix(d[[dind]])), weighted=T)
          mst = spantree(d[[dind]])
          
          # #create edge list
          # node = 1:ncol(as.matrix(d[[dind]]))
          # node = colnames(as.matrix(d[[dind]]))
          # arcs0 = as.data.frame(foreach(from=1:length(mst$kid), .combine='rbind') %dopar% { return(c(mst$kid[from], from+1, mst$dist[from])) })
          # arcs0[,1] = node[arcs0[,1]]
          # arcs0[,2] = node[arcs0[,2]]
          # arcs0[,3] = as.numeric(arcs0[,3])
          # arcs = arcs0
          # am = arcs
          # am[,3] = 1/as.numeric(am[,3])
          # gr = graph_from_data_frame(am, directed=F)
          # wc = cluster_walktrap(gr)
          # 
          # #use min cut to cluster, don't use
          # mincut = NULL
          # for (clusteri in 1:(laul-1)) {
          #   mincut[[clusteri]] = findMinCut(node,arcs)$cut.set
          #   arcs = arcs[-which(arcs[,2]==mincut[[clusteri]]$cut.set[,2]),3]
          # }
          # cl0 = getComponents(node, arcs)
          # 
          # #make into hclust, too slow, don't use
          # depths = spandepth(mst)
          # hc <- as.hclust(mst)
          # plot(hc)
          # plot(hc, col = cutree(hc, 5), pch=16)
          
          
          
          
          
          
          #make into adjacency matrix, =0 if under threshold
          # wthres = 0 #[0,1] (1=cut all)
          # 
          # am0 = as.matrix(d[[dind]])
          # am0[am0==0] = min(am0>0)
          # am0 = 1/am0
          # 
          # am = am0
          # am[am<max(am)*wthres] = 0
          # gr = graph_from_adjacency_matrix(am, mode="undirected", weighted=T, diag=F)
          
          
          
          #rtsne for plot
          di = dt[[dind]]
          
          
          #plot
          
          main = paste0("(0=class, .=label) flowcap-II AML, ",splitby[col],"-",names(d)[dind],", cut edges under ",rwt,", clustering \n", 
                        paste0(paste("__",names(fm[[colnam]][[dind]][[as.character(rwt)]]),"-"), round(fm[[colnam]][[dind]][[as.character(rwt)]],2), collapse=""))
          
          colour = rainbow(max(length(unique(cl)),length(unique(la))))
          cp = expand.grid(c(1:length(unique(cl))),c(20,1:19,21:25))
          plot(di, t='n', main=main)
          lines(mst,di,col='grey')
          for (g in 1:max(length(unique(cl)),length(unique(la)))) {
            if (length(unique(cl))>=g) points(di[cl%in%unique(cl)[g],], col=colour[g], cex=2)
            if (length(unique(la))>=g) points(di[la%in%unique(la)[g],], col=colour[g], pch=16, cex=1)
          }
          legend("topleft", legend=unique(cl), fill=colour[rep(c(1:clul),g)], pch=rep(c(20,1:19,21:25),each=clul)[c(1:g)])
          
          TimeOutput(start1)
          
          
          cat("; Done ", sep="")
          TimeOutput(start1)
        }
      }
    }
    graphics.off()
    return(fm)
  }, error = function(e) {
    cat(paste("ERROR:  ",e)); graphics.off(); return(NULL)
  })
  
}

names(f) = distmfile
distmfile_working = which(vapply(f, Negate(is.null), NA))
if (!length(distmfile_working)>0) { f = NULL; break }
f = f[distmfile_working]
if (!overwrite) f = append(f0,f)
distmfile = names(f)
distmfilenames = names(f) = fileNames(distmfile)
save(f,file=paste0(cl_score_result_dir,".Rdata"))

# split scores from distance matrices on differenct sets of objects
score = NULL
rn = NULL
for (i in 1:length(f)) {
  for (j in 1:length(f[[i]])) {
    for (k in 1:length(f[[i]][[j]])) {
      for (l in 1:length(f[[i]][[j]][[k]]))
      rn = append(rn,paste0(names(f)[i],"__",names(f[[i]])[j],"__",names(f[[i]][[j]][[k]])[l]))
      score = rbind(score,f[[i]][[j]][[k]][[l]])
    }
  }
}
rownames(score) = rn
colnames(score) = scorename
write.csv(score,file=paste0(cl_score_result_dir,".csv"))

TimeOutput(start)








