# measure distance between distance matrices
# aya43@sfu.ca 20170419

root = "~/projects/flowCAP-II"
result_dir = "result"
setwd(root)

options(stringsAsFactors=FALSE)
#options(device="cairo")
options(na.rm=T)

#Input
# sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
# matrix_dir = paste(result_dir, "/matrix", sep="")
dist_dir = paste(result_dir, "/dist", sep="")
plot_dir = paste(result_dir, "/plots", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dir[i])) }
dist_score_dir = paste(dist_dir, "/dist_score", sep=""); for (i in 1:length(dist_score_dir)) { suppressWarnings(dir.create(dist_score_dir[i])) }
dist_clusterf_dir = paste(dist_score_dir, "/dist_cluster_f", sep=""); for (i in 1:length(dist_clusterf_dir)) { suppressWarnings(dir.create(dist_clusterf_dir[i])) }
cl_score_result_dir = paste0(dist_score_dir,"/score_cl_list")

matrix_type_features = c("CountAdj","Prop","LogFold","Pval","Child_entropy","Parent_entropy","Child_prop","Child_pnratio","Parent_contrib","Parent_effort","Freqp0.5")
plotcols = c("feature","layer","featend","dist","norm")
splitcols = c("none","rand","featend")

metrics = c("silmed","NCA","f","r","p","Rand","pearsongamma","dunn","dunn2","entropy","ch")
clustnumbefore = c(T,T,T,T,T,T,F,F,F,F,F) #before is accuracy, after is cluster quality
cltypes = c("distmatrix","knn","kmed","lv","spec1","spec","hc","dc1","dc")#,"spec","dc") #"rw1" (produces one cluster...)
clplotallpar = c("hc")#plot only best

interested = c("tube","aml") #sampleMeta column; choose one
# splitby = c("none","tube") #match with above
plotdist = c("feature","norm","dist","layer")
plotProp = "noProp" #"all"

matrix_count = c("CountAdj")

overwritedist = F
overwritedistd = F
overwritedistf = T

quantilnorm = T

trimonly = T
ignoredist = ".csv|simmatrix"

dis = c("euclidean", "maximum", "manhattan", "canberra", "binary", "bray", "kulczynski", "minkowski", "morisita", "horn", "binomial", "rbf","poly","vanilla","tanh","laplace","bessel","anova","spline")


maxcl = 110 # max number of clusters, else evaluation is slow
maxpl = 60 # max number of plots per png

#Output
plot_compare_dir = paste(plot_dir, "/dist_compare", sep=""); for(i in 1:length(plot_compare_dir)) { suppressWarnings(dir.create(plot_compare_dir[i])) }
dist_comp_dir = paste(result_dir, "/dist_compare", sep="")
dist_comp_dist_dir = paste(result_dir, "/dist_compare_dist", sep="")
dist_comp_score_dir = paste(dist_score_dir, "/dist_compare", sep="")



source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
libr("stringr")
libr("foreach")
libr("doMC")
libr("Rtsne")
libr("devtools")
libr("Biobase")
libr("preprocessCore")


start = Sys.time()

no_cores = detectCores()-1
registerDoMC(no_cores)

#list out dist files
# distmfilenames = list.files(dist_clusterf_dir, recursive=F, full.names=F, pattern=".Rdata")
# distmfile = paste0(dist_dir,"/",distmfilenames)
distmfile = list.files(dist_dir, recursive=F, full.names=T, pattern=".Rdata")
distmfile = c(distmfile[grep("Freqp",distmfile)],
              distmfile[grep("entropyTRIM",distmfile)],distmfile[grep("LogFoldTRIM",distmfile)],distmfile[grep("PvalTRIM",distmfile)],
              distmfile[grep("LogFold_",distmfile)],distmfile[grep("Pval_",distmfile)],distmfile[grep("entropy_",distmfile)],
              distmfile[grep("/CountAdj_",distmfile)],distmfile[grep("/Prop_",distmfile)],
              distmfile[grep("effortTRIM",distmfile)],distmfile[grep("contribTRIM",distmfile)],
              distmfile[grep("pnratioTRIM",distmfile)],distmfile[grep("propTRIM",distmfile)],
              distmfile[grep("effort_",distmfile)],distmfile[grep("contrib_",distmfile)],
              distmfile[grep("pnratio_",distmfile)],distmfile[grep("prop_",distmfile)])
distmfile0 = distmfile
distmfile = c(distmfile[!grepl("_rw_",distmfile)],distmfile[grepl("_rw_",distmfile)])
# distmfile = distmfile[grepl("Freqp0.5_orig",distmfile)]

distmfile = distmfile[!grepl(ignoredist,distmfile,ignore.case=T)]
distmfilenames = fileNames(distmfile)
distmfile = paste0(dist_dir,"/",intersect(distmfilenames,list.files(dist_clusterf_dir, recursive=F, full.names=F, pattern=".Rdata")))
distmfilenames = fileNames(distmfile)


#create distMeta
distMeta = distMetafun(distmfile=distmfile, dis=dis,features=matrix_type_features)
# distmfile = distmfile[distMeta$rand==0 & !distMeta$sim & !distMeta$rw]
# if (trimonly) distmfile = distmfile[grepl("TRIM",distMeta$type)]
# dm = distMeta[match(distmfile,distMeta$path),]





## Preparing distances ----------------------------------------------------------

if (overwritedist | !file.exists(paste0(dist_comp_dir,".Rdata"))) {
  cat("samples in common ")
  
  # samples0 = lapply(1:length(distmfile), function(i) {
  samples0 = foreach(i=1:length(distmfile)) %dopar% {
    #load/Prep files
    d0 = get(load(distmfile[i]))
    return(rownames(as.matrix(d0)))
    # })
  }
  samples1 = Reduce('intersect',samples0)
  TimeOutput(start)
  
  
  cat(", collating distances ")
  d = d0 = foreach(i=1:length(distmfile),.combine="rbind") %dopar% {
    # d0 = lapply(1:length(distmfile), function(i) {
    #load/Prep files
    d0 = as.matrix(get(load(distmfile[i])))
    ind = match(samples1,rownames(d0))
    d0 = d0[ind,ind]
    d = unlist(as.list(as.dist(d0)))
    return(d)
  }
  # })
  # d = d0 = Reduce('rbind',d0)
  save(d0,file=paste0(dist_comp_dir,".Rdata"))
} else {
  d0 = get(load(paste0(dist_comp_dir,".Rdata")))
}

d1 = t(apply(d0,1,function(x) exp(log(x,getlv(max(x),100)))/100 ))
TimeOutput(start)




cat(", normalizing ")
if (quantilnorm) {
  colramp = colorRampPalette(c(3,"white",2))(20)
  #x11();par(mfrow=c(2,1))
  #plot(density(log(d0+1)[,1]),col=colramp[1],lwd=3)
  d = normalize.quantiles(as.matrix(d0))
  #plot(density(log(d+1)[,1]),col=colramp[1],lwd=3)
}
TimeOutput(start)


cat(", distances ")
if (overwritedist | !file.exists(paste0(dist_comp_dist_dir,".Rdata"))) {
  did = dist(d)
  did1 = dist(d1)
} else {
  did0 = get(load(paste0(dist_comp_dist_dir,".Rdata")))
  did = did0$did
  did1 = did0$did1
}
TimeOutput(start)

cat(", Rtsne ")
tsne = Rtsne(as.matrix(did),is_distance=T, perplexity = 1)
tsne1 = Rtsne(as.matrix(did1),is_distance=T, perplexity = 1)
TimeOutput(start)








## Preparing scores -----------------------------------------------------------

if (!file.exists(paste0(dist_comp_score_dir,".Rdata")) | overwritedistf) {
  
  cat(", collating medsil scores ")
  # result = sapply(1:length(distmfile), function(i) {
  result = foreach(i=1:length(distmfile)) %dopar% {
    silm = list()
    fm = get(load(paste0(dist_clusterf_dir,"/",distmfilenames[i])))
    fm00 = unlist(fm)
    
    
    for (mind in 1:length(metrics)) {
      metric = metrics[mind]
      metricpos = which(grepl(paste0("[.]",metric,"$"), names(fm00),ignore.case=T) & !is.na(fm00))
      fm0 = fm00[metricpos]
      
      #find cluster number position and set to 0 if needed
      gaptoclno = grep("[.]cluster.number$",names(fm00))
      if (clustnumbefore[mind]) { gaptoclno0 = sapply(metricpos, function(x) min(gaptoclno[gaptoclno>x]))
      } else { gaptoclno0 = sapply(metricpos, function(x) max(gaptoclno[gaptoclno<x])) }
      fm0[fm00[gaptoclno0]==1] = 0
      
      fm1 = strsplit(names(fm0),"[.]") #col, dind, cltype, par, silmed
      fm1 = lapply(fm1, function(x) {
        if(length(x)>5) return(c(x[1:3],paste0(x[4:5],collapse="."),x[6]))
        return(x)
      })
      fm2 = sapply(fm1,function(x) paste(x,collapse="."))
      fm1 = Reduce('rbind',fm1); colnames(fm1) = c("col", "dind", "cltype", "par", "metric")
      # fm0 = fm0[!duplicated(fm2)]
      
      
      for (col in interested) {
        colind = grepl(paste0("^interested-",col),fm1[,"col"])
        if (sum(colind)==0) next
        cltypes = unique(fm1[colind,"cltype"])
        
        
        for (cltype in cltypes) {
          # if ((grepl("silmed",metric) & !cltype%in%c("distmatrix","spec")) | (!clustnumbefore[mind] & cltype%in%c("knn")) | grepl("NCA",metric) & !cltype%in%c("distmatrix","spec")) next
          clind = colind & fm1[,"cltype"]==cltype
          if (sum(clind)==0) next
          dindnames = unique(fm1[clind,"dind"])
          calcavg = F
          if (length(dindnames)>1) calcavg = T #calculate average score for everything excluding "all"
          
          
          silm[[paste0(col,"_",cltype,"_",metric)]] = NULL
          
          pars = unique(fm1[clind,"par"])
          if (!cltype%in%clplotallpar) pars = "none"
          
          for (par0 in pars) {
            if (par0 == "none") { parind = clind
            } else { parind = clind & fm1[,"par"]==par0 }
            
            for (dindname in dindnames) {
              ddind = parind & fm1[,"dind"]==dindname
              if (sum(ddind)==0) next
              
              fm01 = fm0[ddind]
              fm11 = fm1[ddind,]; if (is.null(dim(fm11))) { fm11 = matrix(fm11,nrow=1); colnames(fm11) = colnames(fm1) }
              fm21 = fm2[ddind]
              
              b = NULL
              if (par0=="none") { b = max(fm01)
              } else { b = max(fm01[fm11[,"par"]==par0]) }
              silm[[paste0(col,"_",cltype,"_",metric)]][paste0("dind-",dindname,"_",par0)] = b
            }
            if (length(dindnames)>1) { #average all dinds except "all"
              ddind0 = parind & fm1[,"dind"]!="all"
              if (sum(ddind0)==0) next
              
              aa = silm[[paste0(col,"_",cltype,"_",metric)]][!grepl("all",names(silm[[paste0(col,"_",cltype,"_",metric)]]))]
              if (length(aa)>0) { if (sum(!is.na(aa))>0) a = mean(aa[!is.na(aa)]) }
              silm[[paste0(col,"_",cltype,"_",metric)]][paste0("dind-avg","_",par0)] = a
              
            }
          }
        }
      }
    }
    return(silm)
    # })
  }
  TimeOutput(start)
  
  names(result) = distmfile
  silm1 = list()
  for (i in names(result[[which.max(sapply(result, function(x) length(x)))]])) {
    silm1[[i]] = foreach(x=result) %dopar% {
      if (!i%in%names(x)) return(NULL)
      a = x[[i]]
      names(a) = names(x[[i]])
      return(a)
    }
    goodind = (vapply(silm1[[i]], Negate(is.null), NA))
    if (sum(goodind)==0) next
    if (length(goodind)<length((distmfile))) goodind = append(goodind, rep(F,length(distmfile)-length(goodind)))
    colnamess = names(silm1[[i]][[which(goodind)[1]]])
    rownamess = distmfile[goodind]
    aa0 = silm1[[i]][goodind]
    aanames = sort(Reduce('union',lapply(aa0, function(x) names(x))))
    aa1 = sapply(aa0, function(x) x[match(aanames,names(x))]) ## some are missing because repeating clusterings are deleted -- set these to min!!
    aa1[is.na(aa1)] = min(aa1/2)
    if (is.null(dim(aa1))){ aa2 = matrix(aa1,nrow=1); rownames(aa2) = aanames; colnames(aa2) = names(aa1); aa1 = aa2 }
    silm1[[i]] = t(aa1)
    colnames(silm1[[i]]) = aanames
    rownames(silm1[[i]]) = rownamess
  }
  save(silm1,file=paste0(dist_comp_score_dir,".Rdata"))
  
} else {
  silm1 = get(load("result/dist/dist_score/dist_compare.Rdata"))
}


# silm0 = foreach(i=1:length(distmfile),.combine="c") %dopar% {
#   fm = get(load(paste0(dist_clusterf_dir,"/",distmfilenames[i])))
#   fm0 = unlist(fm)
#   fm0 = fm0[grepl("[.]silmed", names(fm0)) & grepl(interested, names(fm0))]
#   fm1 = lapply(strsplit(names(fm0),"[.]"),function(x) x[c(1,2,length(x))])
#   fm1 = sapply(fm1,function(x) paste(x,collapse="."))
#   fm0 = fm0[!duplicated(fm1)]
#   names(fm0) = fm1[!duplicated(fm1)]
#   silm = mean(fm0)
#   
#   return(silm)
# }
TimeOutput(start)





## Plotting ---------------------------------------------------------


cat(", plotting... ")

# colour = rainbow(max(sapply(plotdist,function(x)length(unique(distMeta[,x])))))
colour = rainbow(max(length(unique(distMeta[,"feature"])),length(unique(distMeta[,"layer"]))))
cp1 = c(15:18,7:14)[1:length(unique(distMeta[,"featend"]))]
cp2 = c(0,1,5,2,6)[length(unique(distMeta[,"dist"]))] #empty inside

#plotno = 4+length(plotcols)
plotno = length(plotcols)

a=foreach (colscore = names(silm1)) %dopar% {
  silm2 = silm1[[colscore]]
  dm2 = distMeta[distMeta$path%in%rownames(silm2),]
  dm2[dm2[,"weighted"],"dist"] = "weighted manhattan"
  if (plotProp=="noProp") {
    silm2 = silm2[!grepl("Prop",dm2[,"feature"]),]
    if (is.null(dim(silm2))) { silm2 = matrix(silm2,ncol=1); rownames(silm2) = rownames(silm1[[colscore]])[!grepl("Prop",dm2[,"feature"])]; colnames(silm2) = colnames(silm1[[colscore]])}
    dm2 = dm2[!grepl("Prop",dm2[,"feature"],ignore.case=F),]
  }
  
  for(splitcol in splitcols) {
    if (splitcol=="none") {
      splitcolc = list(); splitcolc[["all"]] = 1:nrow(silm2)
    } else { 
      splitcolc = lapply(unique(dm2[,splitcol]), function(x) dm2[,splitcol]==x)
      names(splitcolc) = unique(dm2[,splitcol])
      splitind = sapply(splitcolc, function(x) sum(x)>0)
      if (sum(sapply(splitcolc, function(x) sum(x)>0))==0) next
      splitcolc = splitcolc[sapply(splitcolc, function(x) sum(x)>0)]
    }
    for (splitcolind in 1:length(splitcolc)) {
      dm = dm2[splitcolc[[splitcolind]],]
      silm3 = silm2[splitcolc[[splitcolind]],]
      if (is.null(dim(silm3))) {
        silm3 = matrix(silm3,ncol=ncol(silm2))
        rownames(silm3) = rownames(silm2)[splitcolc[[splitcolind]]]
        colnames(silm3) = colnames(silm2)
      }
      
      if (all(silm3==silm3[1,1])) next
      
      png(filename=paste0(plot_compare_dir,"/", colscore,"__",splitcol,"-",names(splitcolc)[splitcolind],"_",plotProp,".png"), width=700*ncol(silm3), height=500*plotno)
      par(mfcol=c(plotno,ncol(silm3)), mar=c(13,5,5,5))
      main1 = paste0("Distance Parameter vs Quality Metric (split by ",splitcol,"-",names(splitcolc)[splitcolind],")\n")
      # for (i in 1:2) {
      #   if (i==1) { x = tsne; main1 = "distance measures\nsize = exp of avg med silhouette index of aml/healthy from each tube\nquantnormalized" }
      #   if (i==2) { x = tsne1; main1 = "distance measures\nsize = exp of avg med silhouette index of aml/healthy from each tube\n0-1 normalized" }
      #   plot(x$Y,t='n',main=main1, xlim=c(min(x$Y[,1])-.3*(max(x$Y[,1])-min(x$Y[,1])),max(x$Y[,1])))
      #   # for (col in plotdist) {
      #   # }
      #   points(x$Y[,1],x$Y[,2],col=colour[as.integer(factor(dm[,"feature"]))],pch=16,cex=exp(silm0/max(silm0)))
      #   points(x$Y[,1],x$Y[,2],col=colour[as.integer(factor(dm[,"layer"]))],cex=exp((silm0/max(silm0))+.5))
      #   if (sum(dm$weighted)>0) points(x$Y[dm$weighted,1],x$Y[dm$weighted,2],col="black",cex=exp((silm0/max(silm0))[dm$weighted]+1))
      #   if (sum(dm$weightedorig)>0) points(x$Y[dm$weightedorig,1],x$Y[dm$weightedorig,2],col="blue",cex=exp((silm0/max(silm0))[dm$weightedorig]+1))
      # 
      #   legend("topleft",
      #          legend=c(paste0("feature (1st ring): ",unique(dm[,"feature"])),
      #                   paste0("layer (2nd ring): ",unique(dm[,"layer"])), "weighted (3rd ring)"),
      #          col=c(colour[as.integer(factor(unique(dm[,"feature"])))],
      #                colour[as.integer(factor(unique(dm[,"layer"])))], "black"),
      #          pch =c(rep(16,length(unique(dm[,"feature"]))),
      #                 rep(1,length(unique(dm[,"layer"]))), 1) )
      # 
      #   plot(x$Y,t='n',main=main1, xlim=c(min(x$Y[,1])-.3*(max(x$Y[,1])-min(x$Y[,1])),max(x$Y[,1])))
      # 
      #   points(x$Y[,1],x$Y[,2],col=colour[as.integer(factor(dm[,"featend"]))],pch=16,cex=exp(silm0/max(silm0)))
      #   points(x$Y[,1],x$Y[,2],col=colour[as.integer(factor(dm[,"dist"]))],cex=exp((silm0/max(silm0))+.5))
      #   if (sum(dm$weighted)>0) points(x$Y[dm$weighted,1],x$Y[dm$weighted,2],col="black",cex=exp((silm0/max(silm0))[dm$weighted]+1))
      #   if (sum(dm$weightedorig)>0) points(x$Y[dm$weightedorig,1],x$Y[dm$weightedorig,2],col="blue",cex=exp((silm0/max(silm0))[dm$weightedorig]+1))
      #   legend("topleft",
      #          legend=c(paste0("feattype (1st ring): ",unique(dm[,"featend"])),
      #                   paste0("distmeasure (2nd ring): ",unique(dm[,"dist"])), "weighted (3rd ring)"),
      #          col=c(colour[as.integer(factor(unique(dm[,"featend"])))],
      #                colour[as.integer(factor(unique(dm[,"dist"])))], "black"),
      #          pch =c(rep(16,length(unique(dm[,"featend"]))),
      #                 rep(1,length(unique(dm[,"dist"]))), 1) )
      # }
      for (columnplot in 1:ncol(silm3)) {
        silm0 = silm3[,columnplot]
        
        for (testcol in plotcols) {
          dmc0 = as.integer(factor(dm[,testcol]))
          dmc = jitter(dmc0, factor=1)
          dmcl = lapply(sort(unique(dmc0)), function(x) silm0[dmc0==x])
          names(dmcl) = sort(unique(dm[,testcol]))
          #Jitter
          
          boxplot(dmcl, lwd = 1, 
                  main=paste0(main1,"_",testcol," vs ", colscore,"; Only distances of ",splitcol,"-",names(splitcolc)[splitcolind], "\n", colnames(silm3)[columnplot],"-par ;; ", testcol), 
                  xaxt ='n', ylab=colscore, cex.lab=1.5) #,xaxt ='n', xlab=testcol
          axis(1, at=as.integer(factor(unique(dm[,testcol]))), labels=unique(dm[,testcol]), cex.axis = 1.3, las=2)
          title(cex.main=1.5)
          points(dmc,silm0,col=colour[as.integer(factor(dm[,"feature"]))], pch=16,cex=exp((silm0/max(silm0)))/2)
          if (sum(dm$weighted)>0) points(dmc[dm$weighted],silm0[dm$weighted],col="black",cex=exp((silm0/max(silm0))[dm$weighted]+.75)/2)
          if (sum(dm$weightedorig)>0) points(dmc[dm$weightedorig],silm0[dm$weightedorig],col="blue",cex=((silm0/max(silm0))[dm$weightedorig]+.75)/2)
          
        } 
      }
      
      graphics.off()
    }
  }
}

TimeOutput(start)





cat(", PCA-ing... ")
pc = prcomp(d)
plotpc = 5 #number of pca pc to plot
silm0 = silm1[["silmed_aml_spec-none_avg"]]

png(filename=paste0(plot_compare_dir,"/distcomparePCA.png"), width=600, height=600*(1+plotpc))
par(mfrow=c(plotpc+1,1))
plot(pc$sdev^2/sum(pc$sdev^2), xlim = c(0, 15), type = "b", pch = 16, xlab = "principal components", ylab = "variance explained")
legend("topleft",
       legend=c(paste0("feature (1st ring): ",unique(distMeta[,"feature"])),
                paste0("feattype: ",unique(distMeta[,"featend"])),
                paste0("layer (2nd ring): ",unique(distMeta[,"layer"])),
                paste0("distmeasure: ",unique(distMeta[,"dist"])), "weighted (3rd ring)"),
       col=c(colour[as.integer(factor(unique(distMeta[,"feature"])))],
             rep("black",length(unique(distMeta[,"featend"]))), 
             colour[as.integer(factor(unique(distMeta[,"layer"])))],
             rep("black",length(unique(distMeta[,"dist"]))), "black"),
       pch =c(rep(16,length(unique(distMeta[,"feature"]))),
              cp1[as.integer(factor(unique(distMeta[,"featend"])))],
              rep(1,length(unique(distMeta[,"layer"]))),
              cp2[as.integer(factor(unique(distMeta[,"dist"])))], 1) )

for (i in 1:plotpc) {
  plot(pc$x[,i], pc$x[,i+1], t='n', main = "distance measures\nsize = exp of avg med silhouette index of aml/healthy from each tube PCA", xlab = paste0("PC_",i), ylab = paste0("PC_",i+1),xlim=c(min(pc$x[,i])-50,max(pc$x[,i])))
  points(pc$x[,1],pc$x[,2],col=colour[as.integer(factor(distMeta[,"feature"]))],pch=cp1[as.integer(factor(distMeta[,"featend"]))],cex=exp(1.5*silm0/max(silm0)))
  points(pc$x[,1],pc$x[,2],col=colour[as.integer(factor(distMeta[,"layer"]))],pch=cp2[as.integer(factor(distMeta[,"dist"]))],cex=exp((1.5*silm0/max(silm0))+.5))
  points(pc$x[distMeta$weighted,1],pc$x[distMeta$weighted,2],col="black",cex=exp((2*silm0/max(silm0))[dm$weighted]+1))
  points(pc$x[distMeta$weightedorig,1],pc$x[distMeta$weightedorig,2],col="bllue",cex=exp((2*silm0/max(silm0))[dm$weightedorig]+1))
  points(0, 0, pch = 3, cex = 4, lwd = 4)
}
graphics.off()




TimeOutput(start)








