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
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
libr("stringr")
libr("foreach")
libr("doMC")
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
libr("PerfMeas")
libr("clues")
libr("clusteval")

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
    i2 = 1 #input: dist matrix
    if (distMeta$dist[i]=="other") i2 = 2 #input: clustering made directly from features
    if (distMeta$dist[i]=="rbf") i2 = 3 #input: rbf distance matrices, split by variables already
    
    distpath = distMeta$path[i]
    if (i2!=2) distpath = paste0(dist_dir,"/",distMeta$path[i])
    
    #upload distance matrix
    if (i2==3) { 
      d0 = get(load(distpath))
      
      #create 2D tsne display
      d0fullind = which.max(sapply(d0,function(x) nrow(as.matrix(x[[1]])))) #index with most FCM files
      dtsnefile = paste0(plot_dist_source_dir,"/tsne_",fileNames(distpath))
      if (!file.exists(dtsnefile) ) { 
        dtsne = Rtsne(as.dist(d0[[d0fullind]][[1]]),is_distance=T)$Y
        rownames(dtsne) = rownames(as.matrix(d0[[d0fullind]][[1]]))
        save(dtsne,file=dtsnefile)
      } else {
        dtsne = get(load(dtsnefile))
        if (sum(is.na(dtsne))>0 | is.null(rownames(dtsne))) {
          dtsne = Rtsne(as.dist(d0[[d0fullind]][[1]]),is_distance=T)$Y
          rownames(dtsne) = rownames(as.matrix(d0[[d0fullind]][[1]]))
          save(dtsne,file=dtsnefile)
        }
      }
      
    } else if (i2==1) {
      d0 = as.matrix(get(load(distpath))); if (sum(is.na(as.matrix(d0)))>0 | sum(as.matrix(d0)==Inf)>0) stop("distance is blank")
      
      #create 2D tsne display
      dtsnefile = paste0(plot_dist_source_dir,"/tsne_",fileNames(distpath))
      if (!file.exists(dtsnefile) ) { 
        dtsne = Rtsne(as.dist(d0),is_distance=T)$Y
        rownames(dtsne) = rownames(d0)
        save(dtsne,file=dtsnefile)
      } else {
        dtsne = get(load(dtsnefile))
        if (sum(is.na(dtsne))>0) {
          dtsne = Rtsne(as.dist(d0),is_distance=T)$Y
          rownames(dtsne) = rownames(d0)
          save(dtsne,file=dtsnefile)
        }
      }
      
      #adjust distance values
      if (adjustD) if (min(d0[which(d0>0)])>1 & mean(d0)>10000) d0 = d0/mean(d0)
      
      #adjust sampleMeta
      dcol = rep(0,ncol(sampleMeta0)) # won't use dcol here... assume distance matrices and clusterings are marked by fileName
      for (j in 1:ncol(sampleMeta0)) { l = sum(!is.na(match(colnames(as.matrix(d0)),sampleMeta0[,j]))); if (l>0) dcol[j] = l }
      dcol = which.max(dcol) #sampleMeta column corresponding to distance matrix rownames
      if (!length(unique(sampleMeta0[,dcol]))==nrow(sampleMeta0)) stop("distance row labels not found in sampleMeta")
      sampleMeta = sampleMeta0[match(colnames(d0),sampleMeta0[,dcol]),] #sampleMeta in the same order as d
    }
    
    #create score list
    fm = list()
    if (file.exists(paste0(dist_clusterf_dir,"/",fileNames(distpath)))) fm = get(load(paste0(dist_clusterf_dir,"/",fileNames(distpath))))
    
    #create classification/clustering list
    cm0 = list()
    if (file.exists(paste0(dist_clustercl_dir,"/",fileNames(distpath)))) { cm0 = get(load(paste0(dist_clustercl_dir,"/",fileNames(distpath))))
    } else if (i2==2) { next } #must have cm0 in mode 2
    
    #for each intersted column
    for (col in 1:length(interested)) { 
      #create score list
      colnam = paste0("interested-",interested[col],"_splitby-",splitby[col])
      if (!colnam%in%names(fm)) fm[[colnam]] = list()
      if (!colnam%in%names(cm0)) cm0[[colnam]] = list()
      
      #split distance and sampleMeta by a variable + do an average
      split=T
      if (splitby[col]=="none") split=F
      
      if (i2!=2) {
        if (i2==1) { #split distance matrix and sampleMeta
          d = splitmatrix(as.matrix(d0),sampleMeta[,colnames(sampleMeta)==splitby[col]],split=split)
          sm = splitmatrix(sampleMeta,sampleMeta[,colnames(sampleMeta)==splitby[col]],split=split,missym=F)
          dt = splitmatrix(dtsne,sampleMeta[,colnames(sampleMeta)==splitby[col]],split=split,missym=F)
        } else if (i2==3) { #use already split distance matrix and split sampleMeta
          d = d0[[colnam]]
          sm = lapply(d, function(x) sampleMeta0[match(rownames(as.matrix(x)),sampleMeta0$fileName),])
          dt = lapply(d, function(x) dtsne[match(rownames(as.matrix(x)),rownames(dtsne)),])
          names(sm) = names(dt) = names(d)
        }
        if (split & avgall) { #average all distance matrices
          specimenavg = Reduce('intersect',lapply(sm,function(x) x[,avgallcol]))
          dtemp0 = lapply(1:length(d),function(x) {
            smtemp = sm[[x]][match(rownames(as.matrix(d[[x]])),sm[[x]][,"fileName"]),] #order sm[[x]] --> d[[x]]
            dtemp = as.matrix(d[[x]])[match(specimenavg,smtemp[,avgallcol]), match(specimenavg,smtemp[,avgallcol])] #order sm/d[[x]] --> specimenavg
            return(dtemp)
          })
          d[["all"]] = Reduce("+", dtemp0) / length(dtemp0)
          sm[["all"]] = sm[[1]][match(specimenavg,sm[[1]][,avgallcol]),colnames(sm[[1]])!=splitby[col]]
        }
      } else if (i2==2) { #split sampleMeta only; clustering directly taken from cm0
        sm = lapply(cm0[[colnam]], function(x) sampleMeta0[match(rownames(x$spec$clt),sampleMeta0$fileName),]) #assume spec clustering done
        names(sm) = names(cm0[[colnam]])
      }
      
      #for each interested column variable
      for (dindname in names(sm)) {
        cat("\n",interested[col],": ",dindname,".", sep="")
        
        #create score list
        if (!dindname%in%names(fm[[colnam]])) fm[[colnam]][[dindname]] = list()
        if (!dindname%in%names(cm0[[colnam]])) cm0[[colnam]][[dindname]] = list()
        
        #list out class labels
        testind = which(is.na(sm[[dindname]]$label)) #testing files
        laname = sort(unique(sm[[dindname]][,interested[col]])) #true class labels
        la1 = sapply(sm[[dindname]][testind,interested[col]], function(x) which(laname==x)) #label for unknown only; sapply just to make sure indices match...
        la0 = sapply(sm[[dindname]][,interested[col]], function(x) which(laname==x)) #label for all samples
        la0 = sapply(sm[[dindname]][,interested[col]], function(x) which(laname==x)) #label for all samples
        
        #prepare class labels for direct F measure calculation
        la0df = t(sapply(la0, function(lax) {
          a = rep(0,length(unique(la0)))
          a[unique(la0)==lax] = 1
          return(a)
        }))
        rownames(la0df) = sm[[dindname]][,"fileName"]
        colnames(la0df) = laname
        la1df = t(sapply(la1, function(lax) {
          a = rep(0,length(unique(la1)))
          a[unique(la1)==lax] = 1
          return(a)
        }))
        rownames(la1df) = sm[[dindname]][testind,"fileName"]
        colnames(la1df) = unique(la1)
        
        
        
        
        
        
        
        
        
        
        
        
        
        ## cluster ---------------------------------------------------------
        if (docl) {
          cat(" (")
          for (cltype in cltypes) {
            #to do or not to do
            if (cltype%in%names(cm0[[colnam]][[dindname]]) & !overwritecl) next
            if (cltype%in%c("distmatrix")) next
            if (!cltype%in%cltypesorigm & i2==2) next
            if (cltype%in%cltypesorigm) { #for other distance, these are clusterings made directly from original features
              if (i2!=2) next
              if (!"clt"%in%names(cm0[[colnam]][[dindname]][[cltype]])) next
              clt = cm0[[colnam]][[dindname]][[cltype]]$clt
            }
            
            #create class/cluster list
            cm0[[colnam]][[dindname]][[cltype]] = list()
            
            cat(" ",cltype," ",sep="")
            start2 = Sys.time()
            
            ## knn for each k (classification)
            if (cltype=="knn") { 
              labelss = la0
              labelss[is.na(sm[[dindname]]$label)] = NA
              cm0[[colnam]][[dindname]][[cltype]]$clt1 = knntable(as.matrix(d[[dindname]]),cltypespar[["knn"]],labelss)
              next
            } #parameter k
            
            ## kmedoids
            if (cltype=="kmed") { clt = pamtable(d[[dindname]],length(unique(la0)),cltypespar[["kmed"]]) } #number of tries
            
            ## louvain
            if (cltype=="lv") { #input is a similarity matrix
              if (grepl("dist.Rdata",distpath)) {
                sim0 = get(load(paste0(dist_dir,"/",gsub("dist.Rdata","simmatrix.Rdata",distpath))))
                simind = match(rownames(as.matrix(d[[dindname]])),rownames(sim0))
                sim = sim0[simind,simind]
              } else { sim = get_graph(d[[dindname]]) }
              clt = lvtable(sim,cltypespar[["lv"]])
            }
            
            ## random walk refined distance matrices only
            if (cltype=="rw1") {
              if (grepl("dist.Rdata",distpath)) {
                sim0 = get(load(paste0(dist_dir,"/",gsub("dist.Rdata","simmatrix.Rdata",distpath))))
                simind = match(rownames(as.matrix(d[[dindname]])),rownames(sim0))
                sim = sim0[simind,simind]
              } else if (!rwonly) { sim = get_graph(d[[dindname]]) 
              } else { next }
              clt = rw1table(sim,cltypespar[["rw1"]])
            }
            
            ## spectral clustering via distance matrix
            if (cltype=="spec1") {
              if (grepl("dist.Rdata",distpath)) {
                sim0 = get(load(paste0(dist_dir,"/",gsub("dist.Rdata","simmatrix.Rdata",distpath))))
                simind = match(rownames(as.matrix(d[[dindname]])),rownames(sim0))
                sim = sim0[simind,simind]
              } else { sim = get_graph(d[[dindname]]) }
              clt = spec1table(sim,length(unique(la0)))
            }
            
            ## hierarchical clustering
            if (cltype=="hc") { clt = hctable(d[[dindname]],length(unique(la0)),cltypespar[["hc"]]) }
            
            ## densitycut via distance matrix (k=3)
            if (cltype=="dc1") { clt = dc1table(d[[dindname]]) }
            
            TimeOutput(start2)     
            
            
            
            #make an adjusted version of clustering, such that they match with actual labels
            clt1 = matrix(NA,ncol=ncol(clt),nrow=nrow(clt), dimnames=dimnames(clt))
            for (j in 1:ncol(clt)) {
              cl = clt[,j]
              la = la0 #for all metrics
              
              # cl1 s.t. each cluster is labeled by what majority of its real contents; if label is 
              cl1 = rep(NA,length(la))
              if (T | length(unique(la))>2) {
                tubes0 = unique(cl)[order(table(cl))]
                for (tubei in tubes0) {
                  tci = which(cl==tubei) #index of cluster tubei
                  tubej = Mode(la[tci]) #tubei is label taking up majority of cluster tubei
                  cl1[tci] = tubej ## MAJORITY IN 2+ classes?
                }
              } else if (F & length(unique(la))==2) {
                if (sum(la==unique(la)[1])>sum(la==unique(la)[2])) {
                  tci = which(la==unique(la)[1])
                  cl1[cl[tci]==Mode(cl[tci])] = unique(la)[1]
                  cl1[cl[tci]!=Mode(cl[tci])] = unique(la)[2]
                } else {
                  tci = which(la==unique(la)[2])
                  cl1[cl[tci]==Mode(cl[tci])] = unique(la)[2]
                  cl1[cl[tci]!=Mode(cl[tci])] = unique(la)[1]
                }
              }
              
              # cl such that class matches up with cl1; prioritize clusters to larger size labels
              tubes = unique(cl1)[order(table(cl1),decreasing=F)] #don't switch these again!
              for (k in 1:length(tubes)) {
                tubek = tubes[k]
                tck = which(cl1==tubek)
                tubei = Mode(cl[tck])
                stop = F
                if (k>1) {
                  tck = tck[cl[tck]!=tubei & cl[tck]>k]
                  while (length(tck)>0 & tubei%in%tubes[1:(k-1)] & !stop) {
                    tubei = Mode(cl[tck])
                    tck = tck[cl[tck]!=tubei]
                    if (!length(tck)>0) stop = T
                  }
                }
                if (tubei!=tubek & !stop) {
                  clk = which(cl==tubek)
                  cli = which(cl==tubei)
                  if (length(clk)>0) cl[clk] = tubei
                  cl[cli] = tubek
                }
              }
              clt1[,j] = cl1
              clt[,j] = cl
            }
            
            #done!
            cm0[[colnam]][[dindname]][[cltype]]$clt = clt
            cm0[[colnam]][[dindname]][[cltype]]$clt1 = clt1
            
          } #cltype
          cat("), ")
        }
        
        
        
        
        
        
        
        
        
        
        
        
        
        
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
        
        
        
        
        
        
        
        
        
        ## plot -------------------------------------------------------------
        try({ #never stop the loop because of plotting errors!
          if (doplot & i2!=2 & (length(cm0[[colnam]][[dindname]])>0 & length(fm[[colnam]][[dindname]])>0) & ((!dindname=="all" & length(sm)>1) | (dindname=="all" & length(sm)==1))) {
            
            cat(" plotting (",sep="")
            
            #make plot layout
            cltypes0 = names(cm0[[colnam]][[dindname]])
            cols = sapply(cm0[[colnam]][[dindname]],function(x)ncol(x$clt))
            cols[!names(cm0[[colnam]][[dindname]])%in%clplotallpar] = 1
            cols = cols+sum(plot_sill_rp)+1 #plus a silhouette plot, r/p plot, and original plot
            layoutm = matrix(0, ncol=max(cols), nrow=length(cm0[[colnam]][[dindname]]))
            for (cind in 1:length(cols)) { layoutm[cind,1:cols[cind]] = c(1:cols[cind])+suppressWarnings(max(layoutm)) }
            
            #split up layoutm if too much plots per png
            # layoutmm = layoutm
            # nopng = ceiling(max(layoutm)/maxpl) #number of png's to make
            # perpng = ceiling(nrow(layoutm)/nopng) #number of rows per png
            # layoutm0 = list()
            # for (p in 1:nopng) {
            #   layoutm0[[p]] = layoutm[1:min(nrow(layoutm),perpng),]
            #   layoutm0[[p]] = layoutm0[[p]][,!apply(layoutm0[[p]],2,function(x) all(x==0))]
            #   layoutm = layoutm[-c(1:min(nrow(layoutm),perpng)),]
            #   if (nrow(layoutm)==0) break
            #   layoutm[layoutm>0] = layoutm[layoutm>0]-min(layoutm[layoutm>0])+1
            # }
            
            #prepare png
            pngname = paste0(plot_dist_eval_dir, "/tsne_", gsub(".Rdata","",fileNames(distpath)),"_class_",interested[col],"-",splitby[col],"-",dindname,"__",paste(cltypes0,collapse="-"),".png")
            if (!file.exists(pngname) | overwriteplot) {
              graphics.off()
              png(pngname, width=width*ncol(layoutm), height=height*nrow(layoutm))
              par(mar=c(7,2,20,2))
              layout(layoutm)
              
              # plotcount = 0
              # pngcount = 1
              
              for (cltype in names(cm0[[colnam]][[dindname]])) {
                start2 = Sys.time()
                cat(" ",cltype,sep="")
                
                #list out class labels & clusterings/parameters
                if (cltype%in%cltypesclass) { la = la1
                } else { la = la0 }
                
                clt = cm0[[colnam]][[dindname]][[cltype]]$clt
                clt1 = cm0[[colnam]][[dindname]][[cltype]]$clt1
                par = clpar[cltype]
                
                #rtsne for plot
                # mst = spantree(d[[dindname]])
                di0 = dt[[dindname]]
                di = di0[testind,]
                if (nrow(clt)>length(testind)) di = di0
                
                #only plot one par (best) or all pars
                summetrics = 1
                if ("silmed"%in%names(fm[[colnam]][[dindname]][[cltype]])) summetrics = unlist(lapply(fm[[colnam]][[dindname]][[cltype]],function(x) return(x["silmed"])))
                if ("F"%in%names(fm[[colnam]][[dindname]][[cltype]])) summetrics = unlist(lapply(fm[[colnam]][[dindname]][[cltype]],function(x) return(x["F"])))
                
                mmax = which(summetrics==max(summetrics[!is.na(summetrics)]))[1]
                jbest = ceiling(mmax/2)
                jind = 1:ncol(clt)
                if (!cltype%in%clplotallpar) jind = jbest
                
                # pngname = paste0(plot_dist_eval_dir, "/tsne_", gsub(".Rdata","",fileNames(distmfile[i])),"_class_",interested[col],"-",splitby[col],names(d)[dind],"_",cltype,".png")
                # png(pngname, width=width*ceiling(sqrt(length(jind))), height=height*ceiling(sqrt(length(jind))))
                # par(mfrow=c(ceiling(sqrt(length(jind))),ceiling(sqrt(length(jind)))), mar=c(3,3,15,2), oma=c(2,2,5,2))
                
                main0 = paste0("flowcap-II AML, ",cltype,".",par,"=","parinsert","; ",interested[col],", ",splitby[col],"-",dindname," ;; (big0=cluster, 0=labeledcluster, .=label)")
                
                for (j in c(0,jind)) {
                  
                  #plot original plot in front of every row
                  if (j==0) {
                    colour = rainbow(length(unique(la0)))
                    plot(di0, t='n', main=main0)
                    # lines(mst,di0,col='grey')
                    for (g in 1:length(colour)) {
                      if (length(unique(la0))>=g) points(di0[la0%in%unique(la0)[g],1],di0[la0%in%unique(la0)[g],2], col=colour[g], pch=16, cex=1)
                    }
                    legend("topleft", legend=laname[unique(la0)], fill=colour[rep(c(1:g),g)], pch=rep(c(20,1:19,21:25),each=g)[c(1:g)])
                    next
                  }
                  
                  
                  if (par!="") { main1 = gsub("parinsert","",colnames(clt)[j],main0)
                  } else { main1 = gsub(paste0(".",par,"=","parinsert"),"",main0) }
                  main = paste0(main1,"\n\n",paste0(paste0(" ",names(fm[[colnam]][[dindname]][[cltype]][[as.character(colnames(clt)[j])]]),"-"),
                                                    signif(unlist(fm[[colnam]][[dindname]][[cltype]][[as.character(colnames(clt)[j])]]),2), collapse=" "))
                  main = paste(str_wrap(main,width=70),collapse="\n")
                  
                  colour = rainbow(max(length(unique(clt1[,j])),length(unique(clt[,j])),length(unique(la))))
                  
                  # plotcount = plotcount+1
                  # if (plotcount>max(layoutm0[[pngcount]])) {
                  #   plotcount = 1
                  #   pngcount = pngcount+1
                  #   graphics.off()
                  #   png(gsub(".png",paste0("_",pngcount,".png"),pngname), width=width*ncol(layoutm0[[pngcount]]), height=height*nrow(layoutm0[[pngcount]])) 
                  #   par(mar=c(5,2,15,2))
                  #   layout(layoutm0[[pngcount]])
                  # }
                  
                  plot(di0, t='n', main=main)
                  #lines(mst,di0,col='grey')
                  for (g in 1:length(colour)) {
                    if (length(unique(clt1[,j]))>=g & sum(is.na(clt1[,j]))==0) 
                      points(di[clt1[,j]%in%unique(clt1[,j])[g],1], di[clt1[,j]%in%unique(clt1[,j])[g],2], col=colour[g], cex=3)
                    if (length(unique(clt[,j]))>=g) points(di[clt[,j]%in%unique(clt[,j])[g],1],di[clt[,j]%in%unique(clt[,j])[g],2], col=colour[g], cex=2)
                    if (length(unique(la0))>=g) points(di0[la0%in%unique(la0)[g],1],di0[la0%in%unique(la0)[g],2], col=colour[g], pch=16, cex=1)
                  }
                  legend("topleft", legend=laname[unique(clt[,j])], fill=colour[rep(c(1:g),g)], pch=rep(c(20,1:19,21:25),each=g)[c(1:g)])
                  if (j==jbest & plot_sill_rp[1]) {
                    cls = clt1[,j]
                    if (mmax%%2==1) cls = clt[,j]
                    dls = as.dist(d[[dindname]])
                    if (length(cls)<nrow(as.matrix(d[[dindname]]))) dls = as.dist(as.matrix(d[[dindname]])[testind,testind])
                    
                    # plotcount = plotcount+1
                    if (length(unique(cls))>1) plot(silhouette(cls,dls))
                  }
                }
                
                #plot recall precision
                if (plot_sill_rp[2]) { 
                  clstats = list()
                  a = fm[[colnam]][[dindname]][[cltype]]
                  if (!cltype%in%cltypesclass) sil_orig = sapply(a, function(x) return(x["silmed"]))
                  
                  clstats$cl1P = sapply(a, function(x) return(x[grep("^P",names(x))]))
                  clstats$cl1R = sapply(a, function(x) return(x[grep("^R",names(x))]))
                  clstats$cl1p = sapply(a, function(x) return(x[grep("p_co1",names(x))[2]])); clstats$cl1p[is.na(clstats$cl1p)] = clstats$clp[is.na(clstats$cl1p)]
                  clstats$cl1r = sapply(a, function(x) return(x[grep("r_co1",names(x))[2]])); clstats$cl1r[is.na(clstats$cl1r)] = clstats$clr[is.na(clstats$cl1r)]
                  cl1nc = sapply(1:ncol(cm0[[colnam]][[dindname]][[cltype]]$clt1), function(x) length(unique(cm0[[colnam]][[dindname]][[cltype]]$clt1[,x])))
                  if (!cltype%in%cltypesclass) {
                    clstats$clp = sapply(a, function(x) return(x[grep("p_co",names(x))[1]])); clstats$clp[is.na(clstats$clp)] = clstats$cl1p[is.na(clstats$clp)]
                    clstats$clr = sapply(a, function(x) return(x[grep("r_co",names(x))[1]])); clstats$clr[is.na(clstats$clr)] = clstats$cl1r[is.na(clstats$clr)]
                    clnc = sapply(1:ncol(cm0[[colnam]][[dindname]][[cltype]]$clt), function(x) length(unique(cm0[[colnam]][[dindname]][[cltype]]$clt[,x])))
                  }
                  
                  # plotcount = plotcount+1
                  plot(1, t='n', main="recall & precision (label=cluster.number; recall=dashed)",ylim=c(0,1),xlim=c(1,ncol(clt)),xaxt="n", xlab=clpar[cltype])
                  colour = rainbow(1+length(clstats)/2)
                  for (clstatsi in 1:length(clstats)) {
                    if (grepl("r",names(clstats)[clstatsi],ignore.case=T)) {
                      lines(unlist(clstats[[clstatsi]]), col=colour[ceiling(clstatsi/2)],lty=2)
                      points(unlist(clstats[[clstatsi]]), col=colour[ceiling(clstatsi/2)])
                    } else { 
                      lines(unlist(clstats[[clstatsi]]), col=colour[ceiling(clstatsi/2)]) 
                      points(unlist(clstats[[clstatsi]]), col=colour[ceiling(clstatsi/2)]) 
                    }
                    texta = unlist(clstats[[clstatsi]])
                    if (!any(is.null(texta)) & !any(is.na(texta))) {
                      if (grepl("1",names(clstats)[clstatsi])) { text(unlist(clstats[[clstatsi]]), labels=cl1nc, cex=2, pos=3)
                      } else if (!cltype%in%cltypesclass) { text(unlist(clstats[[clstatsi]]), labels=clnc, cex= 2, pos=3) }
                    }
                  }
                  if (!cltype%in%cltypesclass) {
                    lines(unlist(sil_orig), col=colour[length(colour)])
                    points(unlist(sil_orig), col=colour[length(colour)])
                    legend("bottomleft", legend=c("cl1.standard","cl1","cl","silhouette med / truth"), fill=colour)
                  } else {
                    legend("bottomleft", legend=c("cl1.standard","cl1","cl"), fill=colour)
                  }
                  
                  if (clpar[cltype]!="") { axis(1, at=1:length(clstats[[1]]), labels=names(clstats$clr))
                  } else { axis(1, at=1:length(clstats[[1]]), labels=1:length(clstats[[1]])) }
                }
                #mtext(text=main0,side=3,line=0,outer=T)
              }
              cat(")",sep="")
              graphics.off()
            }
          }
        })
        
      } #dind
    } #col
    
    #save scores and clusterings
    if (dof) save(fm,file=paste0(dist_clusterf_dir,"/",fileNames(distpath)))
    if (docl) save(cm0,file=paste0(dist_clustercl_dir,"/",fileNames(distpath)))
    
    return(F)
  }, error = function(e) {
    cat(paste("ERROR:  ",e)); graphics.off(); return(T)
  })
}

TimeOutput(start)












#which distance matrices don't have scores yet
distMeta$path[!distMeta$path%in%(list.files(dist_clusterf_dir,full.names=F))]
distMeta$path[errors]








