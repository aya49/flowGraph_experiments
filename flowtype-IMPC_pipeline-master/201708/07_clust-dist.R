# Use distance metrics to cluster tubes FIX FOR IMPC, STILL CLUSTER FOR TUBES
# aya43@sfu.ca 20170419

#root = "~/projects/IMPC"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#,"Sanger_MLN","CIPHE","TCP","H")

options(stringsAsFactors=FALSE)
#options(device="cairo")
options(na.rm=T)

#Input
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
dist_dir = paste(result_dir, "/", panelL, "/", centreL, "/dist", sep="")
doUnderflow = T #if numbers used to calculate scores --> 0, do work-around (slow)

KOonly = F #delete all WT samples
adjustD = T #divide by a number if the numbers are too big e.g. manhattan
pamtries = 10
overwrite = T #redo and overwrite all past scores
tubecol = 2

width = 700; height = 700 # of each plot in png


#Output
dist_score_dir = paste(dist_dir, "/dist_score", sep=""); for (i in 1:length(dist_score_dir)) { suppressWarnings(dir.create(dist_score_dir[i])) }
cl_score_result_dir = paste0(dist_score_dir,"/score_cl_list")
plot_dir = paste(result_dir, "/", panelL, "/", centreL, "/plots", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dir[i])) }
plot_dist_dir = paste(plot_dir, "/dist", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dist_dir[i])) }

libr(stringr)
libr(foreach)
libr(doMC)
libr(lubridate) #if there are date variables
libr(fpc)
libr(cluster)
libr(Rtsne)
libr(PerfMeas)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

interested = c("tube") #sampleMeta columns to plot
dis = c("euclidean", "manhattan", "canberra")#, "mahalanobis") #assume after splitting dist filename by "_", distance is second element

start = Sys.time()

no_cores = 4#detectCores()-1
registerDoMC(no_cores)

# ftWTGT = "normal"

for (ci in 1:length(paste0(panelL,centreL))) {
  
  start1 = Sys.time()
  
  centre = paste0(panelL," ",centreL)[ci]
  cat("\n",centre,sep="")
  
sampleMeta0 = get(load(sampleMeta_dir[ci]))
distmfile = list.files(dist_dir[ci], recursive=F, full.names=T, pattern=".Rdata")
delpvalind = grep(".csv",distmfile,ignore.case=T)
delpvalind0 = grep("_original",distmfile,ignore.case=T)
delpvalind = append(delpvalind,delpvalind0)
if (length(delpvalind)>0) distmfile = distmfile[-delpvalind]
if (!overwrite) {
  f0 = get(load(paste0(cl_score_result_dir[ci],".Rdata")))
  doneind = which(!is.na(match(names(f0),distmfile)))
  if (length(doneind)>0) distmfile = distmfile[-doneind]
}
distmfilenames = fileNames(distmfile)
interestedCols0 = which(colnames(sampleMeta0)%in%interested)


#calculate score_nac for each distance matrix & interested classes
#result = list()
#for(i in 1:length(distmfile)) {
f = foreach(i=1:length(distmfile)) %dopar% {
  #loop.ind = 1:length(distmfile)
  #for(i in 1:length(distmfile)) {
  #result = foreach(i=loop.ind) %dopar% { 
  start1 = Sys.time()
  cat("\n", distmfile[i], ", loading dist matrix", sep="")
  #result[[i]] = list()
  d = get(load(distmfile[i]))
  # if (length(d)==1) { d=d[[1]]; save(d,file=distmfile[i])}
  if (sum(is.na(as.matrix(d)))>0 | sum(as.matrix(d)==Inf)>0) return(NULL)
  drows = colnames(as.matrix(d))
  
  
  #adjust sampleMeta
  dcol = rep(0,ncol(sampleMeta0))
  for (j in 1:ncol(sampleMeta0)) {
    l = sum(!is.na(match(drows,sampleMeta0[,j])))
    if (l>0) dcol[j] = l
  }
  dcol = which.max(dcol)
  if (length(unique(sampleMeta0[,dcol]))==nrow(sampleMeta0)) {
    sampleMeta <- sampleMeta0[match(drows,sampleMeta0[,dcol]),]
    interestedCols = interestedCols0
    #distmfile[i] = paste0("noScore_",distmfile[i])
  } else {
    interestedCols = 1
    sampleMeta = as.data.frame(colnames(as.matrix(d)))
    next
  }
  
  #adjust distance values
  if (adjustD) {
    d = as.matrix(d)
    if (min(d[which(d>0)])>1 & mean(d)>10000) d = d/mean(d)
    d = as.dist(d)
  }
  
  
  pp = NULL
  fm = NULL
  maxsum = 0
  fm0 = NULL
  cat(", clustering ", sep="")
  tryCatch({
    la = sampleMeta[,tubecol]
    for (j in 1:pamtries) {
      #result[[i]][[colnames(sampleMeta)[interestedCols[j]]]] = NCA_score(d, sampleMeta[,interestedCols[j]], doUnderflow=doUnderflow) # 15s for 2500*2500 d
      #pp[[j]] = pam(d,k=pamk(d,diss=T)$nc,diss=T)
      # plot(pp[[j]])
      # clusplot(pp[[j]])
      pp[[j]] = pam(d, k=length(unique(la)),diss=T) #2min, 21s
      cl0 = NULL
      for (tubei in unique(pp[[j]]$clustering)) {
        #clustering tubei is a part of real tubej
        tci = which(pp[[j]]$clustering==tubei)
        tubej = Mode(la[tci])
        cl0[tci] = tubej ## MAJORITY IN 2+ classes?
      }
      
      #prepare for calculate f measure
      clul0 = length(unique(cl0))
      cldf = foreach(clind=1:length(cl0),.combine='rbind') %dopar% {
        a = rep(0,clul0)
        a[cl0[clind]] = 1
        return(a)
      }
      laul = length(unique(la))
      ladf = foreach(laind=1:length(la),.combine='rbind') %dopar% {
        a = rep(0,laul)
        a[la[laind]] = 1
        return(a)
      }
      if (ncol(cldf)<ncol(ladf)) cldf = cbind(cldf,matrix(0,ncol=ncol(ladf)-ncol(cldf),nrow=nrow(cldf)))
      colnames(ladf) = colnames(cldf) = 1:laul
      rownames(ladf) = rownames(ladf) = names(pp[[j]]$clustering)
      
      #calculte f measure
      fm0[[j]] = F.measure.single.over.classes(ladf,cldf)
      if (sum(fm0[[j]]$average[-length(fm0[[j]]$average)])>maxsum) {
        fm = fm0[[j]]
        cl = cl0
        clul = clul0
        maxsum = sum(fm0[[j]]$average[-length(fm0[[j]]$average)])
      }
    }
    #rtsne for plot
    di = Rtsne(d,is_distance=T)$Y
    rownames(di) = cl
    
    #plot
    pngname = paste0(plot_dist_dir, "/tsne_", fileNames(distmfile[i]),"_cluster.png")
    main = paste0("flowcap-II AML tube cluster",
                  paste0(c("\nPrec-","_Recall/Sens-", "_Spec-", "_F-", "_0/1LossAcc-"), round(fm0$average[-length(fm0$average)],2), collapse=""))
    
    png(file=pngname , width=width, height=height)
    par(mar=(c(5,5,5,40) + 0.1))
    
    colour = rainbow(clul)
    cp = expand.grid(c(1:clul),c(20,1:19,21:25))
    plot(di, t='n', main=main)
    for (g in 1:length(unique(cl))) {
      points(di[which(cl%in%unique(cl)[g]),], col=colour[cp[g%%nrow(cp),1]], pch=cp[g%%nrow(cp),2])
    }
    legend("topleft", legend=unique(cl), fill=colour[rep(c(1:clul),g)], pch=rep(c(20,1:19,21:25),each=clul)[c(1:g)])
    graphics.off()
    
    TimeOutput(start1)
    
    
    cat("; Done ", sep="")
    TimeOutput(start1)
    return(fm)
  }, error = function(e) {
    cat(paste("ERROR:  ",e));
  })
}

names(f) = distmfile
distmfile_working = which(vapply(f, Negate(is.null), NA))
if (!length(distmfile_working)>0) { f = NULL; break }
f = f[distmfile_working]
if (!overwrite) f = append(f0,f)
distmfile = names(f)
distmfilenames = names(f) = fileNames(distmfile)
save(f,file=paste0(cl_score_result_dir[ci],".Rdata"))

# split scores from distance matrices on differenct sets of objects
score = foreach(i=1:length(f),.combine='rbind') %dopar% {
  return(f[[i]]$average[-length(f[[i]]$average)])
}
rownames(score) = names(f)
colnames(score) = c("Prec","Recall/Sens", "Spec", "F", "0/1LossAcc")
write.csv(score,file=paste0(cl_score_result_dir[ci],".csv"))

cat("\n centre ", centre, " ",TimeOutput(start1)," \n",sep="")

}
TimeOutput(start)








