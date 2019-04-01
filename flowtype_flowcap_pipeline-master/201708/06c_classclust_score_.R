## Input: distance metrics --> Output: classify/cluster, score, and plot
# aya43@sfu.ca 20170419

#Directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)
options(include.rownames=F)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="") #metadata for cell populations (phenotype)
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="") #metadata for FCM files
sampleMetaTrain_dir = paste0("attachments/AMLTraining.csv") #which FCM files are testing/training files
matrix_dir = paste(result_dir, "/matrix", sep="") #feature matrices
dist_dir = paste(result_dir, "/dist", sep="")
plot_dir = paste(result_dir, "/plots", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dir[i])) }
plot_dist_dir = paste(plot_dir, "/dist", sep=""); for(i in 1:length(plot_dist_dir)) { suppressWarnings(dir.create(plot_dist_dir[i])) }
plot_dist_source_dir = paste0(plot_dist_dir, "/source"); for(i in 1:length(plot_dist_source_dir)) { suppressWarnings(dir.create(plot_dist_source_dir[i])) }

#Output
plot_dist_eval_dir = paste(plot_dir, "/dist_eval", sep=""); for(i in 1:length(plot_dist_eval_dir)) { suppressWarnings(dir.create(plot_dist_eval_dir[i])) }
dist_score_dir = paste(dist_dir, "/dist_score", sep=""); for (i in 1:length(dist_score_dir)) { suppressWarnings(dir.create(dist_score_dir[i])) }
cl_score_result_dir = paste0(dist_score_dir,"/score_cl_list")
dist_clusterf_dir = paste(dist_score_dir, "/dist_cluster_f", sep=""); for (i in 1:length(dist_clusterf_dir)) { suppressWarnings(dir.create(dist_clusterf_dir[i])) }
dist_clustercl_dir = paste(dist_score_dir, "/dist_cluster_cl", sep=""); for (i in 1:length(dist_clustercl_dir)) { suppressWarnings(dir.create(dist_clustercl_dir[i])) }
dist_kern_dir = paste(dist_dir, "/dist_kern", sep=""); for (i in 1:length(dist_kern_dir)) { suppressWarnings(dir.create(dist_kern_dir[i])) }
dist_score_all_dir = paste0(dist_score_dir,"/dist_score_all")

#Libraries/Functions
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
libr(PerfMeas)
libr(clues)
libr(clusteval)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")

#Setup Cores
no_cores = 3#detectCores()-1
registerDoMC(no_cores)












#Options for script
adjustD = T #divide by a number if the numbers are too big e.g. manhattan

docl = T #do clustering
dof = T #do scoring
doplot = T #do plotting
overwritecl = T #redo and overwrite all clusterings
overwritef = T #redo and overwrite all scores
overwriteplot = T #redo and overwrite all plots

deldup = F # delete duplicate clusterings (don't if want to standardized columns)
doUnderflow = T #if numbers used to calculate scores are near 0, do work-around (slow) for NCA
maxcl = 110 # max number of clusters allowed, else evaluation is slow
avgall = T #avg all scores across variables
avgallcol = "specimen" #when averaging distances, match based on specimen
rwonly = T #for rw, ie cut edges, use rw distances only; not used for now

width = 700; height = 700 # of each plot in png
maxpl = 60 # max number of plots per png; not needed for flowCAP-II
clplotallpar = c("spec","hc") #else only plot best parameter in cltype
plot_sill_rp = c(T,T) #plot best sollhouette plot and/or recall precision plot

cltypes = c("distmatrix","knn","kmed","lv","spec1","spec","hc","dc1","dc")#,"spec","dc") #"rw1" (produces one cluster...)
cltypesclass = c("knn") # classification
cltypesorigm = c("spec","dc") #uses original matrix; see 06a_dist_score.R
sctypes_ivdist = c("NCA", "silmed") #metrics for distmatrix only
cltypespar = list(knn=c(1:6), kmed=c(1:6), lv=c(0,.05,.1,.2), spec=list(methods=c("rbf"),kpar="automatic",tries=1), spec1=1, hc=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"), rw1=c(.01,.05,.1,.2,.3,.4,.5,.75))
clpar = c(distmatrix="",knn="k",kmed="",lv="cutEunderQ",spec="kernel",spec1="",hc="link",dc="",dc1="")

matrix_count = c("CountAdj")
ignoredist = ".csv|simmatrix"

dis = c("euclidean", "maximum", "manhattan", "canberra", "binary", "bray", "kulczynski", "minkowski", "morisita", "horn", "binomial", "rbf","poly","vanilla","tanh","laplace","bessel","anova","spline")

interested = c("tube","aml") #sampleMeta columns to plot
splitby = c("none","tube") #match with above



#Prepare data
phenoMeta = get(load(phenoMeta_dir))
mcount = get(load(paste0(matrix_dir, matrix_count,".Rdata")))

#add 'label' col to sampleMeta0
sampleMeta0 = get(load(sampleMeta_dir))
sampleMetaTrain = read.csv(sampleMetaTrain_dir)
ts = strsplit(sampleMeta0$fileName,"[A-Z]")
tu = sapply(ts, function(x) x[2])
sa = sapply(ts, function(x) x[3])
label = sapply(1:nrow(sampleMeta0), function(i) return(sampleMetaTrain[which(sampleMetaTrain[,"TubeNumber"]==tu[i] & sampleMetaTrain[,"SampleNumber"]==sa[i]),"Label"]))
sampleMeta0 = cbind(sampleMeta0,label)

#list out dist files
distmfile = append(list.files(dist_dir, recursive=F, full.names=F, pattern=".Rdata"),list.files(dist_clustercl_dir, recursive=F, full.names=F, pattern=".Rdata")) 
# distmfile = c(distmfile[grep("Freqp",distmfile)],
#               distmfile[grep("entropyTRIM",distmfile)],distmfile[grep("LogFoldTRIM",distmfile)],distmfile[grep("PvalTRIM",distmfile)],
#               distmfile[grep("LogFold_",distmfile)],distmfile[grep("Pval_",distmfile)],distmfile[grep("entropy_",distmfile)],
#               distmfile[grep("/CountAdj_",distmfile)],distmfile[grep("/Prop_",distmfile)],
#               distmfile[grep("effortTRIM",distmfile)],distmfile[grep("contribTRIM",distmfile)],
#               distmfile[grep("pnratioTRIM",distmfile)],distmfile[grep("propTRIM",distmfile)],
#               distmfile[grep("effort_",distmfile)],distmfile[grep("contrib_",distmfile)],
#               distmfile[grep("pnratio_",distmfile)],distmfile[grep("prop_",distmfile)])
distmfile = c(distmfile[!grepl("_rw_",distmfile)],distmfile[grepl("_rw_",distmfile)])
distmfile = distmfile[!grepl(ignoredist,distmfile,ignore.case=T)]
distmfile = gsub("FULL-","",distmfile)
distmfile = distmfilenames = distmfile0 = unique(distmfile)
# distmfile = distmfile[grepl("Freqp0.5_orig",distmfile)]

# distmfile = distmfile[!distMeta[,c("rw")] & !distMeta[,c("sim")]& !distMeta[,c("weighted")] & distMeta[,"layer"]==max(distMeta[,"layer"]) &
#   ((grepl("effort|contrib|prop|pnratio",distmfile) & distMeta$norm=="none") | 
#      (!grepl("effort|contrib|prop|pnratio",distmfile) & distMeta$norm=="cellpop")) | 
#   !distmfilenames%in%list.files(dist_clusterf_dir, recursive=F, full.names=F, pattern=".Rdata")]
# distmfile = distmfile[!distmfile%in%list.files(dist_clusterf_dir, recursive=F, full.names=F, pattern=".Rdata")]

# distmfilenames = fileNames(distmfile)
distMeta = distMetafun(distmfile,append(dis,"other"))
# distMeta = distMeta[order(distMeta$layer,decreasing=T),]
# nonorm = ((grepl("effort|contrib|prop|pnratio",distMeta$type) & distMeta$norm=="none") | (!grepl("effort|contrib|prop|pnratio",distMeta$type) & distMeta$norm=="cellpop"))
# nopropfreq = !grepl("Prop",distMeta$type,ignore.case=F) & !grepl("Freqp",distMeta$type)
# norand = distMeta$rand==0
# distMeta = distMeta[nonorm & nopropfreq & norand,]#,distMeta[!nonorm,])
# distMeta = distMeta[!(nonorm & nopropfreq & norand),]
# weightednode = (!grepl("effort|contrib|prop|pnratio",distMeta$type) & distMeta$weighted)
# distMeta = distMeta[weightednode,]#,distMeta[!nonorm,])












start = Sys.time()

#cluster/classify, score, plot distance matrices
loop.ind = c(1:length(distMeta$path))#[!distMeta$path%in%list.files(dist_clusterf_dir,full.names=F)]
errors = foreach(i=loop.ind, .combine="c") %dopar% {
  #for (i in loop.ind) {
  tryCatch({

        
        
        
        
        
        
        
        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        ## score -------------------------------------------------------
        
        if (dof) {
          cltypes0 = cltypes
          if (i2==2) cltypes0 = cltypes[!grepl("distmatrix",cltypes)]
          
          cat("scoring (")
          for (cltype in cltypes0) {
            #to do or not to do
            if (cltype!="distmatrix") {
              if (!cltype%in%names(cm0[[colnam]][[dindname]])) next # if no clustering/classification is made available, skip
              if ((i2!=2 & cltype%in%cltypesorigm) | (i2==2 & !cltype%in%cltypesorigm)) next
              pars = colnames(clt1)
            } else { pars = "none" }
            
            #create score lists
            if (!cltype%in%names(fm[[colnam]][[dindname]]) | overwritef) fm[[colnam]][[dindname]][[cltype]] = list()
            
            start2 = Sys.time()
            cat(" ",cltype," ",sep="")
            
            #prepare clusterings/classifications
            if (!cltype%in%c("distmatrix")) {
              if (!cltype%in%cltypesclass) clt = cm0[[colnam]][[dindname]][[cltype]]$clt
              clt1 = cm0[[colnam]][[dindname]][[cltype]]$clt1
            }
            
            #prepare class labels & class label charts for direct F measure calculation
            if (cltype%in%cltypesclass) {
              la = la1
              ladf = la1df
            } else {
              la = la0
              ladf = la0df
            }
            
            #for each classification/clustering parameter
            for (pari in 1:length(pars)) {
              par = pars[pari]
              
              if (par%in%names(fm[[colnam]][[dindname]][[cltype]]) & !overwritef) next
              fm[[colnam]][[dindname]][[cltype]][[par]] = NULL
              
              
              if (!cltype%in%c("distmatrix")) {
                #prepare classification/clustering
                if (!cltype%in%c(cltypesclass)) cl = clt[,pari]
                cl1 = clt1[,pari]
                
                #prepare class label chart for direct F1 measure calculation
                cldf = t(sapply(cl1, function(clx) {
                  a = rep(0,length(unique(la)))
                  a[unique(la)==clx] = 1
                  return(a)
                }))
                colnames(cldf) = laname
                rownames(cldf) = names(cl1)
                
                
                
                ## external validation F1 (classification & clustering)
                F1 = F.measure.single.over.classes(ladf,cldf)$average[-6]
                f11 = f.measure.comembership(la,cl1); names(f11) = paste0(names(f11), "_co_1")
                #r1 = adjustedRand(la,cl1) 
                score = c(F1, f11)#, r1)
                if (length(unique(cl1))==1) score[1:length(score)] = NA #rand doesn't have a minimum of 0
                fm[[colnam]][[dindname]][[cltype]][[par]] = append(fm[[colnam]][[dindname]][[cltype]][[par]], score)
                
                
                
                ## external validation f1 (clustering)
                if (!cltype%in%cltypesclass) {
                  f1 = f.measure.comembership(la,cl); names(f1) = paste0(names(f1), "_co")
                  #r = adjustedRand(la,cl)
                  score = c(f1) #,r)
                  if (length(unique(cl))==1) score[1:length(score)] = NA
                  fm[[colnam]][[dindname]][[cltype]][[par]] = append(fm[[colnam]][[dindname]][[cltype]][[par]], score)
                }
              }
              
              
              
              if (!cltype%in%cltypesclass & !cltype%in%cltypesorigm) {
                
                ## internal validation NCA (distance)
                if (cltype%in%c("distmatrix")) {
                  cl = la0
                  if (!"NCA"%in%names(fm[[colnam]][[dindname]][[cltype]][[par]]) | overwritef) {
                    if (length(unique(cl))==1) { nca = NA
                    } else { nca = NCA_score(as.matrix(d[[dindname]]), cl)$p }
                    fm[[colnam]][[dindname]][[cltype]][[par]]["NCA"] = nca
                  }
                }
                
                
                
                ## internal validation silmed (distance & clustering)
                if (!"silmed"%in%names(fm[[colnam]][[dindname]][[cltype]][[par]]) | overwritef) {
                  if (length(unique(cl))==1) { sil = NA
                  } else { sil = median(silhouette(cl,d[[dindname]])[,3]) }
                  fm[[colnam]][[dindname]][[cltype]][[par]]["silmed"] = sil
                }
                
                if (length(unique(cl))==1) {
                  score = rep(NA,9)
                  names(score) = c("average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")
                } else { 
                  score = unlist(cluster.stats(as.dist(d[[dindname]]),cl))[c("average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")] 
                }
                fm[[colnam]][[dindname]][[cltype]][[par]] = append(fm[[colnam]][[dindname]][[cltype]][[par]], score)
                
                
                
                ## internal validation (adjusted clustering)
                if (!cltype%in%c("distmatrix")) {
                  if (length(unique(cl1))==1) { sil = NA
                  } else { sil = median(silhouette(cl1,d[[dindname]])[,3]) }
                  fm[[colnam]][[dindname]][[cltype]][[par]]["silmed_1"] = sil
                  
                  if (length(unique(cl1))==1) {
                    score = rep(NA,9)
                    names(score) = c("average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")
                  } else {
                    score = unlist(cluster.stats(as.dist(d[[dindname]]),cl1))[c("average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")] 
                  }
                  names(score) = paste0(names(score),"_1")
                  fm[[colnam]][[dindname]][[cltype]][[par]] = append(fm[[colnam]][[dindname]][[cltype]][[par]], score)
                }
              }
            }
            TimeOutput(start2)
          }
          cat(")")
          
          
          
          ## caculate avg score across variable split
          if (length(sm)>1) {
            for (cltype in names(fm[[colnam]][[dindname]])) {
              for (par in names(fm[[colnam]][[dindname]][[cltype]])) {
                for (sctype in names(fm[[colnam]][[dindname]][[cltype]][[par]])) {
                  if (!"avg"%in%names(fm[[colnam]])) fm[[colnam]][["avg"]] = list()
                  if (!cltype%in%names(fm[[colnam]][["avg"]])) fm[[colnam]][["avg"]][[cltype]] = list()
                  if (!par%in%names(fm[[colnam]][["avg"]][[cltype]])) fm[[colnam]][["avg"]][[cltype]][[par]] = NULL
                  meanvals = sapply(names(fm[[colnam]])[!names(fm[[colnam]])%in%c("all","avg")], function(x) fm[[colnam]][[x]][[cltype]][[par]][[sctype]])
                  fm[[colnam]][["avg"]][[cltype]][[par]][sctype] = mean(meanvals[!is.na(meanvals)])
                }
              }
            }
          }
          
        }
        
        
        
        
        
        
        
        

    #save scores and clusterings
    if (dof) save(fm,file=paste0(dist_clusterf_dir,"/",fileNames(distpath)))

    return(F)
  }, error = function(e) {
    cat(paste("ERROR:  ",e)); return(T)
  })
}

TimeOutput(start)












#which distance matrices don't have scores yet
distMeta$path[!distMeta$path%in%(list.files(dist_clusterf_dir,full.names=F))]
distMeta$path[errors]








