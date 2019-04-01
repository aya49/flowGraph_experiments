# Use distance metrics to classify aml
# aya43@sfu.ca 20170419

root = "~/projects/flowCAP-II"
result_dir = "result"
setwd(root)

options(stringsAsFactors=FALSE)
#options(device="cairo")
options(na.rm=T)

#Input
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
sampleMetaTrain_dir = paste0("attachments/AMLTraining.csv")
matrix_dir = paste(result_dir, "/matrix", sep="")
dist_dir = paste(result_dir, "/dist", sep="")
plot_dir = paste(result_dir, "/plots", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dir[i])) }
plot_dist_dir = paste(plot_dir, "/dist", sep=""); for(i in 1:length(plot_dist_dir)) { suppressWarnings(dir.create(plot_dist_dir[i])) }
plot_dist_source_dir = paste0(plot_dist_dir, "/source"); for(i in 1:length(plot_dist_source_dir)) { suppressWarnings(dir.create(plot_dist_source_dir[i])) }
dist_dir = paste(result_dir, "/dist", sep="")
doUnderflow = T #if numbers used to calculate scores --> 0, do work-around (slow)

KOonly = F #delete all WT samples
adjustD = T #divide by a number if the numbers are too big e.g. manhattan
pamtries = 10
overwritef = T #redo and overwrite all past scores
overwritecl = T
overwriteplot = T
dof = T
docl = T
doplot = T

savekern = T

width = 700; height = 700 # of each plot in png
rowplot = 3

matrix_count = c("CountAdj")

cltypes = c("knn","kmed","lv","spec1","spec","hc","dc1","dc")#,"spec","dc") #"rw1" (produces one cluster...)
rwonly = T #for rw, ie cut edges, use rw distances only
cltypesclass = c("knn") # which is classification
cltypesorigm = c("spec","dc") #uses original matrix
cltypespar = list(knn=c(1:6), kmed=c(1:10), lv=c(0,.05,.1,.2,.5,.75), spec=list(methods=c("rbf"),kpar="automatic",tries=1), spec1=1, hc=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"), rw1=c(.01,.05,.1,.2,.3,.4,.5,.75))
# spec methods with no automatic kpar learning c("poly","vanilla","tanh","laplace","bessel","anova","spline")
clpar = c(knn="k",kmed="",lv="cutEunderQ",spec="kernel",spec1="",hc="link",dc="",dc1="")
clplotallpar = c("spec","hc") #else only plot best
ignoredist = ".csv|simmatrix"

dis = c("euclidean", "maximum", "manhattan", "canberra", "binary", "bray", "kulczynski", "minkowski", "morisita", "horn", "binomial", "rbf","poly","vanilla","tanh","laplace","bessel","anova","spline")

plot_sill_rp = c(T,T) #plot best sollhouette plot and/or recall precision plot

interested = c("tube","aml") #sampleMeta columns to plot
splitby = c("none","tube") #match with above

maxcl = 110 # max number of clusters, else evaluation is slow
maxpl = 60 # max number of plots per png


libr(stringr)
libr(foreach)
libr(doMC)
libr(lubridate) #if there are date variables
libr(FastKNN)
libr(fpc)
libr(cluster)
libr(mclust)
libr(igraph)
libr(kernlab)
libr(densitycut) #devtools::install_bitbucket("jerry00/densitycut_dev")
libr(vegan)
libr(Rtsne)
libr(fpc)
libr(PerfMeas)
libr(clues)
libr(clusteval)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")


#Output
plot_dist_eval_dir = paste(plot_dir, "/dist_eval", sep=""); for(i in 1:length(plot_dist_eval_dir)) { suppressWarnings(dir.create(plot_dist_eval_dir[i])) }
dist_score_dir = paste(dist_dir, "/dist_score", sep=""); for (i in 1:length(dist_score_dir)) { suppressWarnings(dir.create(dist_score_dir[i])) }
cl_score_result_dir = paste0(dist_score_dir,"/score_cl_list")
dist_clusterf_dir = paste(dist_score_dir, "/dist_cluster_f", sep=""); for (i in 1:length(dist_clusterf_dir)) { suppressWarnings(dir.create(dist_clusterf_dir[i])) }
dist_clustercl_dir = paste(dist_score_dir, "/dist_cluster_cl", sep=""); for (i in 1:length(dist_clustercl_dir)) { suppressWarnings(dir.create(dist_clustercl_dir[i])) }
dist_kern_dir = paste(dist_dir, "/dist_kern", sep=""); for (i in 1:length(dist_kern_dir)) { suppressWarnings(dir.create(dist_kern_dir[i])) }












start = Sys.time()

no_cores = 10
registerDoMC(no_cores)

#add 'label' col to sampleMeta0
m0 = get(load(paste0(matrix_dir, matrix_count,".Rdata")))
sampleMeta0 = get(load(sampleMeta_dir))
sampleMetaTrain = read.csv(sampleMetaTrain_dir)
ts = strsplit(sampleMeta0$fileName,"[A-Z]")
tu = sapply(ts, function(x) x[2])
sa = sapply(ts, function(x) x[3])
label = sapply(1:nrow(sampleMeta0), function(i) return(sampleMetaTrain[which(sampleMetaTrain[,"TubeNumber"]==tu[i] & sampleMetaTrain[,"SampleNumber"]==sa[i]),"Label"]))
sampleMeta0 = cbind(sampleMeta0,label)

#list out dist files
distmfile = list.files(dist_dir, recursive=F, full.names=T, pattern=".Rdata")
distmfile = c(distmfile[grep("entropyTRIM",distmfile)],distmfile[grep("LogFoldTRIM",distmfile)],distmfile[grep("PvalTRIM",distmfile)],
              distmfile[grep("LogFold_",distmfile)],distmfile[grep("Pval_",distmfile)],distmfile[grep("entropy_",distmfile)],
              distmfile[grep("/CountAdj_",distmfile)],distmfile[grep("/Prop_",distmfile)],
              distmfile[grep("effortTRIM",distmfile)],distmfile[grep("contribTRIM",distmfile)],
              distmfile[grep("pnratioTRIM",distmfile)],distmfile[grep("propTRIM",distmfile)],
              distmfile[grep("effort_",distmfile)],distmfile[grep("contrib_",distmfile)],
              distmfile[grep("pnratio_",distmfile)],distmfile[grep("prop_",distmfile)])
distmfile0 = distmfile
distmfile = c(distmfile[-grep("_rw_",distmfile)],distmfile[grep("_rw_",distmfile)])
delpvalind = grepl(ignoredist,distmfile,ignore.case=T)
if (sum(delpvalind)>0) distmfile = distmfile[!delpvalind]
distmfilenames = fileNames(distmfile)

distMeta = distMetafun(distmfile,dis)




a = foreach(i=1:length(distmfile),.combine="c") %dopar% {
  
  adf=F
  if (((grepl("effort|contrib|prop|pnratio",distmfile[i]) & distMeta$norm[i]=="none" & all(!distMeta[i,c("rw","sim","weighted")]) ) | 
                                      (!grepl("effort|contrib|prop|pnratio",distmfile[i]) & distMeta$norm[i]=="cellpop") & all(!distMeta[i,c("rw","sim","weighted")]))) {
    #load/Prep files
    start1 = Sys.time()
    cat("\n", distmfile[i], ": loading dist matrix", sep="")
    
    #upload dist, original matrix when needed, tsne corrdinates
    d0 = as.matrix(get(load(distmfile[i]))); if (sum(is.na(as.matrix(d0)))>0 | sum(as.matrix(d0)==Inf)>0) return(NULL)
    
    mmorig = NULL
    mresult = Loadintermatrices(paste0(matrix_dir,distMeta$type[i],".Rdata"),verbose=F)
    phenolevel = sapply(colnames(mresult$mml[[1]]), function(x){return(length(strsplit(x,"[+-]")[[1]])) } )
    mmlresult = trimMatrices(mresult$mml,m0,mresult$pt,mresult$gt,phenolevel,NULL,F,distMeta$count[i],distMeta$layer[i])
    mmorig = mmlresult$mml[[1]]
    pm = mmlresult$pm
    if (is.null(dim(mmorig))) mmorig = mmorig[[1]]
    
    
    #for each interested column
    for (col in 1:length(interested)) { 
      #split by another column first, assuming d is by filename
      split=T
      if (splitby[col]=="none") split=F
      d = splitmatrix(as.matrix(d0),sampleMeta[,colnames(sampleMeta)==splitby[col]],split=split)

      
      #for each split
      maxcltcol = 0
      for (dind in 1:length(d)) {
        if (!is.null(mmorig) & sum(is.na(match(rownames(d[[dind]]),rownames(mmorig))))>0) adf=T
      } 
    }
  }
  return(adf)
}
