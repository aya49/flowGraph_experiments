# Collates scores made in dist_score
# aya43@sfu.ca 20170419

root = "~/projects/flowCAP-II"
result_dir = "result"
setwd(root)

options(stringsAsFactors=FALSE)
#options(device="cairo")
options(include.rownames=F)
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="") #metadata for cell populations (phenotype)
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="") #metadata for FCM files
sampleMetaTrain_dir = paste0("attachments/AMLTraining.csv") #which FCM files are testing/training files
matrix_dir = paste(result_dir, "/matrix", sep="") #feature matrices
dist_dir = paste(result_dir, "/dist", sep="")
plot_dir = paste(result_dir, "/plots", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dir[i])) }
plot_dist_dir = paste(plot_dir, "/dist", sep=""); for(i in 1:length(plot_dist_dir)) { suppressWarnings(dir.create(plot_dist_dir[i])) }
plot_dist_source_dir = paste0(plot_dist_dir, "/source"); for(i in 1:length(plot_dist_source_dir)) { suppressWarnings(dir.create(plot_dist_source_dir[i])) }


adjustD = T #divide by a number if the numbers are too big e.g. manhattan

docl = T #do clustering
dof = T #do scoring
doplot = T #do plotting
overwritecl = T #redo and overwrite all past scores
overwritef = T
overwriteplot = T
overwritedist = T

deldup = F # delete duplicate clusterings (don't if want standardized columns)
doUnderflow = T #if numbers used to calculate scores --> 0, do work-around (slow)
maxcl = 110 # max number of clusters, else evaluation is slow
avgall = T
avgallcol = "specimen" #when averaging distances, match based on specimen
rwonly = T #for rw, ie cut edges, use rw distances only


width = 700; height = 700 # of each plot in png
maxpl = 60 # max number of plots per png
clplotallpar = c("spec","hc") #else only plot best parameter in cltype
plot_sill_rp = c(T,T) #plot best sollhouette plot and/or recall precision plot


cltypes = c("distmatrix","knn","kmed","lv","spec1","spec","hc","dc1","dc")#,"spec","dc") #"rw1" (produces one cluster...)
cltypesclass = c("knn") # which is classification
cltypesorigm = c("spec","dc") #uses original matrix
cltypespar = list(knn=c(1:6), kmed=c(1:6), lv=c(0,.05,.1,.2,.5,.75), spec=list(methods=c("rbf"),kpar="automatic",tries=1), spec1=1, hc=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"), rw1=c(.01,.05,.1,.2,.3,.4,.5,.75))
# spec methods with no automatic kpar learning c("poly","vanilla","tanh","laplace","bessel","anova","spline")
clpar = c(distmatrix="",knn="k",kmed="",lv="cutEunderQ",spec="kernel",spec1="",hc="link",dc="",dc1="")

sctypes_ivdist = c("NCA", "silmed") #metrics for distmatrix only


matrix_count = c("CountAdj")
matrix_type_features = c("CountAdj","Prop","LogFold","Pval","Child_entropy","Parent_entropy","Child_prop","Child_pnratio","Parent_contrib","Parent_effort","Freqp0.5")
ignoredist = ".csv|simmatrix"

dis = c("euclidean", "maximum", "manhattan", "canberra", "binary", "bray", "kulczynski", "minkowski", "morisita", "horn", "binomial", "rbf","poly","vanilla","tanh","laplace","bessel","anova","spline")

interested = c("tube","aml") #sampleMeta columns to plot
splitby = c("none","tube") #match with above


source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
libr("stringr")
libr("foreach")
libr("doMC")
libr("doSNOW")
libr("lubridate") #if there are date variables
libr("FastKNN")
libr("fpc")
libr("cluster")
libr("mclust")
libr("igraph")
libr("kernlab")
libr("densitycut") #devtools::install_bitbucket("jerry00/densitycut_dev")
libr("vegan")
libr("Rtsne")
libr("fpc")
libr("PerfMeas")
libr("clues")
libr("clusteval")
libr("xtable")




#Output
plot_dist_eval_dir = paste(plot_dir, "/dist_eval", sep=""); for(i in 1:length(plot_dist_eval_dir)) { suppressWarnings(dir.create(plot_dist_eval_dir[i])) }
dist_score_dir = paste(dist_dir, "/dist_score", sep=""); for (i in 1:length(dist_score_dir)) { suppressWarnings(dir.create(dist_score_dir[i])) }
cl_score_result_dir = paste0(dist_score_dir,"/score_cl_list")
dist_clusterf_dir = paste(dist_score_dir, "/dist_cluster_f", sep=""); for (i in 1:length(dist_clusterf_dir)) { suppressWarnings(dir.create(dist_clusterf_dir[i])) }
dist_clustercl_dir = paste(dist_score_dir, "/dist_cluster_cl", sep=""); for (i in 1:length(dist_clustercl_dir)) { suppressWarnings(dir.create(dist_clustercl_dir[i])) }
dist_kern_dir = paste(dist_dir, "/dist_kern", sep=""); for (i in 1:length(dist_kern_dir)) { suppressWarnings(dir.create(dist_kern_dir[i])) }
dist_score_all_dir = paste0(dist_score_dir,"/dist_score_all")


















## collating scores
start=Sys.time()

no_cores = detectCores()-1
registerDoMC(no_cores)


distmfile1 = list.files(dist_clusterf_dir, recursive=F, full.names=T)
fs = file.size(distmfile1)
distmfile10 = distmfile1 = distmfile1[order(fs)]; fs = sort(fs)
distmfile1names = fileNames(distmfile1)
#head(distmfile1names[fs<30000])
distMeta1 = distMetafun(distmfile1,c(dis,"other"), features = matrix_type_features)
# distMeta1$size = fs

f = foreach(di=distmfile1) %dopar% { get(load(di)) }
names(f) = fileNames(distmfile1)

result = foreach(i = names(f)) %dopar% { #dist
  b = NULL
  a = NULL
  aind = 1
  for (j in names(f[[i]])) { #col
    for (k in names(f[[i]][[j]])) { #dind
      for (l in names(f[[i]][[j]][[k]])) { #cltype
        for (m in names(f[[i]][[j]][[k]][[l]])) { #par
          b = rbind(b, c(i,j,k,l,m))
          aa = sapply(f[[i]][[j]][[k]][[l]][[m]], function(x) { #some are lists, so need to get them all to be vector
            if (!is.na(x) & !is.na(x)) return(x)
            return(NA)
          })
          names(aa) = names(f[[i]][[j]][[k]][[l]][[m]])
          a[[aind]] = aa[!is.na(aa)]
          aind = aind+1
        }
      }
    }
  }
  return(list(rn=b, score=a))
}

TimeOutput(start)


#collate scores and their attributes
scores = foreach(ii=1:length(result)) %dopar% { result[[ii]]$score }
scores = unlist(scores, recursive=F)
rn0 = foreach(ii=1:length(result),.combine="rbind") %dopar% { result[[ii]]$rn }
colnames(rn0) = c("fpath","colnam","dindname","cltype","par")

nrow(rn0)==length(scores) #should be true
nrow(rn0)==nrow(distMeta1) #should be true

#prepare attributes table
rn = cbind(distMeta1[match(rn0[,1],fileNames(distMeta1$path)),],rn0[,-1])

#collate scores into a table (NA for missing)
metricnames = Reduce('union',lapply(scores, function(x) names(x)))
score0 = foreach(ii=1:length(scores), .combine='rbind') %dopar% {
  a = rep(NA,length(metricnames))
  a[match(names(scores[[ii]]),metricnames)] = scores[[ii]]
  return(a)
}
rownames(score0) = rn[,1]
colnames(score0) = metricnames

score1 = cbind(rn[,-1],score0)
rownames(score1) = NULL#rn[,1]

save(score1,file=paste0(dist_score_all_dir,".Rdata"))
write.csv(score1,file=paste0(dist_score_all_dir,".csv"))



TimeOutput(start)









score1 = get(load(paste0(dist_score_all_dir,".Rdata")))
score2 = score1[!grepl("Prop|Freqp",score1$type,ignore.case=F),]
scorecol = 16

classind = score2$cltype%in%cltypesclass
distmind = score2$cltype%in%c("distmatrix")
clustind = !classind & !distmind
classcols = !sapply(scorecol:ncol(score2), function(x) any(is.na(score2[classind,x])))
distcols = !sapply(scorecol:ncol(score2), function(x) any(is.na(score2[distmind,x])))
clustcols = !sapply(scorecol:ncol(score2), function(x) any(is.na(score2[clustind,x])))

norand = score2$wrand==0
noweight = !score2$weighted

tubeind = score2$colnam=="interested-tube_splitby-none"
amlind = score2$colnam=="interested-aml_splitby-tube"
amlavgind = amlind & score2$dindname=="avg"








## create xtables






