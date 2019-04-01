# Use distance metrics to classify aml
# aya43@sfu.ca 20170419

root = "~/projects/flowCAP-II"
result_dir = "result"
setwd(root)

options(stringsAsFactors=FALSE)
#options(device="cairo")
options(include.rownames=F)
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
sampleMetaTrain_dir = paste0("attachments/AMLTraining.csv")
matrix_dir = paste(result_dir, "/matrix", sep="")
rw_dir = paste(result_dir, "/rw", sep="")

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

cltypes = c("distmatrix","knn","kmed","lv","spec1","spec","hc","dc1","dc")#,"spec","dc") #"rw1" (produces one cluster...)
rwonly = T #for rw, ie cut edges, use rw distances only
cltypesclass = c("knn") # which is classification
cltypesorigm = c("spec","dc") #uses original matrix
cltypespar = list(knn=c(1:6), kmed=c(1:10), lv=c(0,.05,.1,.2,.5,.75), spec=list(methods=c("rbf"),kpar="automatic",tries=1), spec1=1, hc=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"), rw1=c(.01,.05,.1,.2,.3,.4,.5,.75))
# spec methods with no automatic kpar learning c("poly","vanilla","tanh","laplace","bessel","anova","spline")
clpar = c(distmatrix="",knn="k",kmed="",lv="cutEunderQ",spec="kernel",spec1="",hc="link",dc="",dc1="")
clplotallpar = c("spec","hc") #else only plot best
ignoredist = ".csv|simmatrix"

dis = c("euclidean", "maximum", "manhattan", "canberra", "binary", "bray", "kulczynski", "minkowski", "morisita", "horn", "binomial", "rbf","poly","vanilla","tanh","laplace","bessel","anova","spline")

plot_sill_rp = c(T,T) #plot best sollhouette plot and/or recall precision plot

interested = c("tube","aml") #sampleMeta columns to plot
splitby = c("none","tube") #match with above

avgall = T
avgallcol = "specimen"

maxcl = 110 # max number of clusters, else evaluation is slow
deldup = F # delete duplicate clusterings (don't if want standardized columns)
maxpl = 60 # max number of plots per png


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


#Output
plot_dist_eval_dir = paste(plot_dir, "/dist_eval", sep=""); for(i in 1:length(plot_dist_eval_dir)) { suppressWarnings(dir.create(plot_dist_eval_dir[i])) }
dist_score_dir = paste(dist_dir, "/dist_score", sep=""); for (i in 1:length(dist_score_dir)) { suppressWarnings(dir.create(dist_score_dir[i])) }
cl_score_result_dir = paste0(dist_score_dir,"/score_cl_list")
dist_clusterf_dir = paste(dist_score_dir, "/dist_cluster_f", sep=""); for (i in 1:length(dist_clusterf_dir)) { suppressWarnings(dir.create(dist_clusterf_dir[i])) }
dist_clustercl_dir = paste(dist_score_dir, "/dist_cluster_cl", sep=""); for (i in 1:length(dist_clustercl_dir)) { suppressWarnings(dir.create(dist_clustercl_dir[i])) }
dist_kern_dir = paste(dist_dir, "/dist_kern", sep=""); for (i in 1:length(dist_kern_dir)) { suppressWarnings(dir.create(dist_kern_dir[i])) }












start = Sys.time()

no_cores = detectCores()-1
registerDoMC(no_cores)

# cl <- makeCluster(no_cores, outfile="")
# registerDoSNOW(cl)

#add 'label' col to sampleMeta0
phenoMeta = get(load(phenoMeta_dir))
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
distmfile = distmfile[grepl("/rw_",distmfile)]
# distmfile = c(distmfile[grep("Freqp",distmfile)],
#               distmfile[grep("entropyTRIM",distmfile)],distmfile[grep("LogFoldTRIM",distmfile)],distmfile[grep("PvalTRIM",distmfile)],
#               distmfile[grep("LogFold_",distmfile)],distmfile[grep("Pval_",distmfile)],distmfile[grep("entropy_",distmfile)],
#               distmfile[grep("/CountAdj_",distmfile)],distmfile[grep("/Prop_",distmfile)],
#               distmfile[grep("effortTRIM",distmfile)],distmfile[grep("contribTRIM",distmfile)],
#               distmfile[grep("pnratioTRIM",distmfile)],distmfile[grep("propTRIM",distmfile)],
#               distmfile[grep("effort_",distmfile)],distmfile[grep("contrib_",distmfile)],
#               distmfile[grep("pnratio_",distmfile)],distmfile[grep("prop_",distmfile)])
# distmfile0 = distmfile
# distmfile = c(distmfile[!grepl("_rw_",distmfile)],distmfile[grepl("_rw_",distmfile)])
# 
# distmfile = distmfile[!grepl(ignoredist,distmfile,ignore.case=T)]
distmfilenames = fileNames(distmfile)
distMeta = distMetafun(distmfile,dis)

# distmfile = distmfile[!distMeta[,c("rw")] & !distMeta[,c("sim")]& !distMeta[,c("weighted")] & distMeta[,"layer"]==max(distMeta[,"layer"]) &
#   ((grepl("effort|contrib|prop|pnratio",distmfile) & distMeta$norm=="none") | 
#      (!grepl("effort|contrib|prop|pnratio",distmfile) & distMeta$norm=="cellpop")) | 
#   !distmfilenames%in%list.files(dist_clusterf_dir, recursive=F, full.names=F, pattern=".Rdata")]
#distmfile = distmfile[!distmfilenames%in%list.files(dist_clusterf_dir, recursive=F, full.names=F, pattern=".Rdata")]

# distmfile = distmfile[grepl("_weighted",distmfile)]
# distmfilenames = fileNames(distmfile)
# distmfilefnames = list.files(dist_clusterf_dir, recursive=F, full.names=F, pattern=".Rdata")
# distmfile = distmfile[!distmfilenames%in%distmfilefnames]


distmfilenames = fileNames(distmfile)
distMeta = distMetafun(distmfile,dis)









# nplot = 0
# for (i in order(distmfilelayer,decreasing=T)) {
#   if (splitby[i]=="none") {
#     nplot = nplot+1
#   } else {
#     nplot = nplot + length(unique(sampleMeta0[,splitbyCols[i]]))
#   }
# }
# colplot = ceiling(nplot/rowplot)

#calculate score_nac for each distance matrix & interested classes
#result = list()
#for(i in 1:length(distmfile)) {
loop.ind=1:length(distmfile)
# pb <- txtProgressBar(max = max(loop.ind), style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress = progress) #, .options.snow = opts
redodist = foreach(i=loop.ind) %dopar% {
    
    # tryCatch({
    #   fm = get(load(paste0(dist_clusterf_dir,"/",distmfilenames[i])))
    #   return(F)
    # }, error = function(e) {
    #   return(T)
    # })
    
  tryCatch({
    redo = F
    
    #load/Prep files
    start1 = Sys.time()
    cat("\n", distmfile[i], ": loading dist matrix", sep="")
    
    #upload dist, original matrix when needed, tsne corrdinates
    d0 = as.matrix(get(load(distmfile[i]))); if (sum(is.na(as.matrix(d0)))>0 | sum(as.matrix(d0)==Inf)>0) return(NULL)
    
    #for matrices without rownames
    # if (is.null(rownames(d0))) {
    #   mresult = Loadintermatrices(paste0(matrix_dir,distMeta$type[i],".Rdata"),verbose=F)
    #   phenolevel = sapply(colnames(mresult$mml[[1]]), function(x){return(length(strsplit(x,"[+-]")[[1]])) } )
    #   mmlresult = trimMatrices(mresult$mml,m0,mresult$pt,mresult$gt,phenolevel,NULL,F,distMeta$count[i],distMeta$layer[i])
    #   m = mmlresult$mml[[1]]
    #   if (is.null(dim(m))) rownames(d0) = rownames(m[[1]])
    #   if (!is.null(dim(m))) rownames(d0) = rownames(m)
    #   save(as.dist(d0),file=distmfile[i])
    #   return(T)
    # }
    # return(F)
    
    mmorig = NULL
    mcp0 = distMeta$type[i]
    if (distMeta$rand[i]>0) mcp0 = paste0(distMeta$type[i],"_",distMeta$rand[i])
    if (grepl("rw_",mcp0)) {
      not_normalized = F
    } else if (any(cltypesorigm%in%cltypes) & docl & all(!distMeta[i,c("rw","sim","weighted")]) & distMeta[i,"layer"]==max(distMeta[,"layer"]) &
               ((grepl("effort|contrib|prop|pnratio",distmfile[i]) & distMeta$norm[i]=="none") | 
                (!grepl("effort|contrib|prop|pnratio",distmfile[i]) & distMeta$norm[i]=="cellpop"))) {
      not_normalized = T
    } else {
      not_normalized = F
    }
    if (not_normalized) {
      mresult = Loadintermatrices(paste0(matrix_dir,mcp0,".Rdata"),verbose=F)
      if (!is.null(dim(mresult$mml))) phenolevel = phenoMeta$phenolevel[match(rownames(mresult$mml),phenoMeta$phenotype)]
      if (is.null(dim(mresult$mml))) phenolevel = phenoMeta$phenolevel[match(rownames(mresult$mml[[1]]),phenoMeta$phenotype)]
      mmlresult = trimMatrices(mresult$mml,m0,mresult$pt,mresult$gt,phenolevel,NULL,F, distMeta$count[i],distMeta$layer[i])
      mmorig = mmlresult$mml[[1]]
      pm = mmlresult$pm
      if (is.null(dim(mmorig))) mmorig = Reduce("cbind",mmorig)
      if (sum(is.na(match(rownames(d0),rownames(mmorig))))>0) redo = T
    } 
    
    if (is.null(rownames(d0))) {
      if (is.null(mmorig)) {
        mresult = Loadintermatrices(paste0(matrix_dir,mcp0,".Rdata"),verbose=F)
        if (!is.null(dim(mresult$mml))) phenolevel = phenoMeta$phenolevel[match(rownames(mresult$mml),phenoMeta$phenotype)]
        if (is.null(dim(mresult$mml))) phenolevel = phenoMeta$phenolevel[match(rownames(mresult$mml[[1]]),phenoMeta$phenotype)]
        mmlresult = trimMatrices(mresult$mml,m0,mresult$pt,mresult$gt,phenolevel,NULL,F, distMeta$count[i],distMeta$layer[i])
        mmorig = mmlresult$mml[[1]]
        if (is.null(dim(mmorig))) mmorig = mmorig[[1]]
      }
      rownames(d0) = colnames(d0) = rownames(mmorig)
      mmorig = NULL
    }
    
    dtsnefile = paste0(plot_dist_source_dir,"/tsne_",distmfilenames[i])
    if (!file.exists(dtsnefile) ) { 
      dtsne = Rtsne(as.dist(d0),is_distance=T)$Y
      rownames(dtsne) = rownames(d0)
      save(dtsne,file=dtsnefile)
    } else {
      dtsne = get(load(dtsnefile))
      if (sum(is.na(dtsne))>0) {
        dtsne = Rtsne(as.dist(d0),is_distance=T)$Y;
        rownames(dtsne) = rownames(d0)
        save(dtsne,file=dtsnefile)
      }
    }
    
    #adjust distance values
    if (adjustD) if (min(d0[which(d0>0)])>1 & mean(d0)>10000) d0 = d0/mean(d0)
    
    #adjust sampleMeta
    dcol = rep(0,ncol(sampleMeta0))
    for (j in 1:ncol(sampleMeta0)) {
      l = sum(!is.na(match(colnames(as.matrix(d0)),sampleMeta0[,j])))
      if (l>0) dcol[j] = l
    }
    dcol = which.max(dcol)
    if (length(unique(sampleMeta0[,dcol]))==nrow(sampleMeta0)) {
      sampleMeta = sampleMeta0[match(colnames(d0),sampleMeta0[,dcol]),]
      dtsne = dtsne[match(sampleMeta$fileName,rownames(dtsne)),]
    } else { return(NULL) }
    
    
    
    # pngname = paste0(plot_dist_dir, "/tsne_", gsub(".Rdata","",fileNames(distmfile[i])),"_class.png")
    # png(pngname, width=width*colplot, height=height*rowplot)
    # par(mfrow=c(rowplot,colplot))
    
    
    fm = list()
    if (file.exists(paste0(dist_clusterf_dir,"/",distmfilenames[i])) & (!overwritef | !dof)) fm = get(load(paste0(dist_clusterf_dir,"/",distmfilenames[i])))
    
    cm0 = list()
    if (file.exists(paste0(dist_clustercl_dir,"/",distmfilenames[i])) & (!overwritecl | !docl)) cm0 = get(load(paste0(dist_clustercl_dir,"/",distmfilenames[i])))
    
    #for each interested column
    for (col in 1:length(interested)) { 
      colnam = paste0("interested-",interested[col],"_splitby-",splitby[col])
      if (!colnam%in%names(fm)) fm[[colnam]] = list()
      if (!colnam%in%names(cm0)) cm0[[colnam]] = list()
      
      #split by another column first, assuming d is by filename
      split=T
      if (splitby[col]=="none") split=F
      d = splitmatrix(as.matrix(d0),sampleMeta[,colnames(sampleMeta)==splitby[col]],split=split)
      sm = splitmatrix(sampleMeta,sampleMeta[,colnames(sampleMeta)==splitby[col]],split=split,missym=F)
      dt = splitmatrix(dtsne,sampleMeta[,colnames(sampleMeta)==splitby[col]],split=split,missym=F)
      if (split & avgall) {
        dtemp = d
        specimenavg = Reduce('intersect',lapply(sm,function(x) x[,avgallcol]))
        dtemp0 = lapply(1:length(d),function(x) {
          smtemp = sm[[x]][match(rownames(d[[x]]),sm[[x]][,"fileName"]),] #order sm[[x]] --> d[[x]]
          dtemp = d[[x]][match(specimenavg,smtemp[,avgallcol]), match(specimenavg,smtemp[,avgallcol])] #order sm/d[[x]] --> specimenavg
          return(dtemp)
        })
        d[["all"]] = Reduce("+", dtemp0) / length(dtemp0)
        sm[["all"]] = sm[[1]][match(specimenavg,sm[[1]][,avgallcol]),colnames(sm[[1]])!=splitby[col]]
        dt[["all"]] = Rtsne(as.dist(d[["all"]]),is_distance=T)$Y
        rownames(dt[["all"]]) = rownames(d[["all"]])
      }
      
      #for each split
      maxcltcol = 0
      for (dind in 1:length(d)) {
        start0 = Sys.time()
        dindname = names(d)[dind]
        cat("\n",interested[col],": ",names(d)[dind],".", sep="")
        testind = which(is.na(sm[[dind]]$label)) #testing files
        laname = sort(unique(sm[[dind]][,interested[col]])) #true class labels
        la1 = sapply(sm[[dind]][testind,interested[col]], function(x) which(laname==x)) #label for unknown only
        la0 = sapply(sm[[dind]][,interested[col]], function(x) which(laname==x)) #label for all samples
        
        #if no class/clustering needed
        if (!length(unique(la0))>1) next
        
        #for each class/clustering algorithm
        if (!dindname%in%names(cm0[[colnam]])) cm0[[colnam]][[dindname]] = list()
        if (!dindname%in%names(fm[[colnam]])) fm[[colnam]][[dindname]] = list()
        
        
        
        ## cluster ---------------------------------------------------------
        if (docl) {
          mm0ready = F
          cat(" (")
          for (cltype in cltypes) {
            start2 = Sys.time()
            if (is.null(mmorig) & cltype%in%cltypesorigm) next
            if (!(!cltype%in%names(cm0[[colnam]][[dindname]]) | overwritecl | cltype%in%c("distmatrix")) ) next
            
            clt = NULL
            par = clpar[cltype]
            if (cltype%in%cltypesclass) { la = la1
            } else { la = la0 }
            
              cm0[[colnam]][[dindname]][[cltype]] = list()
              
              cat(cltype," ",sep="")
              
              if (cltype=="distmatrix") {
                clt = matrix(la0,ncol=1)
                colnames(clt) = "1"
                rownames(clt) = rownames(d[[dind]]) 
                if ("sil"%in%names(fm[[colnam]][[dindname]])) {
                  cm0[[colnam]][[dindname]][[cltype]]$sil = fm[[colnam]][[dindname]]$sil
                  fm[[colnam]][[dindname]]$sil = NULL
                }
              }
              
              #knn for each k
              if (cltype=="knn") { clt = knntable(d[[dind]],cltypespar[["knn"]],as.integer(factor(sm[[dind]]$label))) } #parameter k
              
              #kmedoids
              if (cltype=="kmed") { clt = pamtable(d[[dind]],length(unique(la0)),cltypespar[["kmed"]]) } #number of tries
              
              #louvain
              if (cltype=="lv") {
                if (grepl("dist",distmfilenames[i])) {
                  sim0 = get(load(gsub("dist.Rdata","simmatrix.Rdata",distmfile[i])))
                  simind = match(rownames(d[[dind]]),rownames(sim0))
                  sim = sim0[simind,simind]
                } else { sim = get_graph(d[[dind]]) }
                clt = lvtable(sim,cltypespar[["lv"]])
              }
              
              #random walk only
              if (cltype=="rw1") {
                if (grepl("dist",distmfilenames[i])) {
                  sim0 = get(load(gsub("dist.Rdata","simmatrix.Rdata",distmfile[i])))
                  simind = match(rownames(d[[dind]]),rownames(sim0))
                  sim = sim0[simind,simind]
                } else if (!rwonly) { sim = get_graph(d[[dind]]) 
                } else { next }
                clt = rw1table(sim,cltypespar[["rw1"]])
              }
              
              #spectral; don't know why parent_contrib sometimes doesn't work, negative and 0 values? therefore skip if so
              if (cltype=="spec") { 
                
                if (!mm0ready) {
                  mm0 = as.matrix(mmorig)[match(rownames(d[[dind]]),rownames(mmorig)),]
                  mm0non0 = !sapply(1:ncol(mm0), function(x) all(mm0[,x]==0))
                  if (sum(mm0non0)==0) next
                  mm0 = mm0[,mm0non0]
                  mm0ready = T
                }
                if (savekern) { 
                  dis_dir_temp = distmfilenames[i]
                  dis_dir_temp = strsplit(distmfilenames[i],"_")[[1]]
                  dis_dir_temp[dis_dir_temp%in%dis] = "kerns"
                  dis_dir_temp = paste(dis_dir_temp,collapse="_")
                  kern_dir = paste0(dist_kern_dir,"/",gsub(".Rdata","",dis_dir_temp),"_",colnam,"-",dindname)
                  # clt = spectable(mm0,length(unique(la0)),la0,cltypespar[["spec"]]$methods,cltypespar[["spec"]]$tries, savedist=kern_dir, savesim=kern_dir,replace="kerns")
                  
                  nclass = length(unique(la0)); tries = cltypespar[["spec"]]$tries; replace = "kerns"; method = "rbf"
                  parlist = clt = parlist0 = sil = NCA = cl = sim = dist = NULL
                  ii = 1 # for (ii in 1:tries) {
                    sp = specc(as.matrix(mm0),kernel=paste0(method,"dot"),centers=nclass)
                    cl[[ii]] = sp@.Data
                    parlist0[[ii]] = unlist(sp@kernelf@kpar)
                    sim[[ii]] = kernelMatrix(kernel=match.fun(paste0(method,"dot"))(as.numeric(parlist0[[ii]])),x=as.matrix(mm0))
                    dist[[ii]] = get_graphd(sim[[ii]])
                    sil = append(sil,median(silhouette(la0,dist[[ii]])[,3]))
                    NCA = append(NCA,NCA_score(as.matrix(dist[[ii]]), la0, doUnderflow=doUnderflow)$p)
                  # }
                  maxsilind = which.max(sil)
                  clt = cbind(clt, cl[[maxsilind]])
                  parlist = c(parlist, paste0(paste0(names(sp@kernelf@kpar),"-",signif(parlist0[[maxsilind]],4),collapse="_"),"_silmed-",signif(max(sil),5),"_NCA-",signif(NCA[maxsilind],5),collapse=""))
                  sss = sim[[maxsilind]]
                  ddd = dist[[maxsilind]]
                  save(sss,file=paste0(gsub(replace,paste0("kern_",method,"_",paste0(names(sp@kernelf@kpar),"-",signif(parlist0[[maxsilind]],4),collapse="_")),kern_dir),"_simmatrix.Rdata"))
                  save(ddd,file=paste0(gsub(replace,paste0("kern_",method,"_",paste0(names(sp@kernelf@kpar),"-",signif(parlist0[[maxsilind]],4),collapse="_")),kern_dir),"_dist.Rdata"))
                  colnames(clt) = parlist
                  rownames(clt) = rownames(mm0)
                  
                } else { clt = spectable(mm0,length(unique(la0)),la0,cltypespar[["spec"]]$methods,cltypespar[["spec"]]$tries) }
                
              }
              
              if (cltype=="spec1") {
                if (grepl("dist",distmfilenames[i])) {
                  sim0 = get(load(gsub("dist.Rdata","simmatrix.Rdata",distmfile[i])))
                  simind = match(rownames(d[[dind]]),rownames(sim0))
                  sim = sim0[simind,simind]
                } else { sim = get_graph(d[[dind]]) }
                clt = spec1table(sim,length(unique(la0)))
              }
              
              
              #hierarchical clustering
              if (cltype=="hc") { clt = hctable(d[[dind]],length(unique(la0)),cltypespar[["hc"]]) }
              
              #density cut
              if (cltype=="dc") { 
                ## do only if distance ends with correct 'normalize'
                if (!mm0ready) {
                  mm0 = as.matrix(mmorig)[match(rownames(d[[dind]]),rownames(mmorig)),]
                  mm0non0 = !sapply(1:ncol(mm0), function(x) all(mm0[,x]==0))
                  if (sum(mm0non0)==0) next
                  mm0 = mm0[,mm0non0]
                  mm0ready = T
                }
                clt = dctable(mm0)
              }
              
              if (cltype=="dc1") { clt = dc1table(d[[dind]]) }
              
              
              
              
              
              
              
              
              
              TimeOutput(start2)     
              
              if (which(cltypes==cltype)==length(cltypes)) { cat("), ")
              } else { cat(" ") }
              
              #determine if valid clusterings exist
              # delcol_maxcl = apply(clt,2,function(x) return( length(unique(x))>maxcl))
              # delcol_dup = duplicated(t(clt))
              # leavcol = !delcol_maxcl & !delcol_dup #delete clustering if too many classes or duplicated
              # if (is.null(clt) | sum(leavcol)==0) {
              #   if (length(cm0[[colnam]][[dindname]])==1) { cm0[[colnam]][[dindname]] = list()
              #   } else { cm0[[colnam]][[dindname]][[cltype]] = NULL } 
              #   next 
              # }
              # #leave only valid clusterings
              # clt0 = clt
              # if (sum(leavcol)==1) {
              #   clt = matrix(clt[,leavcol],ncol=1)
              #   rownames(clt) = rownames(clt0)
              #   colnames(clt) = colnames(clt0)[leavcol]
              # } else { clt = clt[,leavcol] }
              
              if (cltype%in%c("distmatrix") | nrow(clt)==length(la1)) { clt1 = clt
              } else {
                #make an adjusted version of clustering, such that they match with actual labels
                clt1 = matrix(NA,ncol=ncol(clt),nrow=nrow(clt))
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
              }
              cm0[[colnam]][[dindname]][[cltype]]$clt = clt
              cm0[[colnam]][[dindname]][[cltype]]$clt1 = clt1
            
          } #cltype
        }
        
        
        
        
        
        
        
        ## score -------------------------------------------------------
        if (dof & length(cm0[[colnam]][[dindname]])>0) {
          cat(" scoring (",sep="")
          
          for (cltype in names(cm0[[colnam]][[dindname]])) { #[!grepl("sil",names(cm0[[colnam]][[dindname]]))]) {
            start2 = Sys.time()
            cat(" ",cltype,sep="")
            par = clpar[cltype]
            
            if (cltype%in%cltypesclass) { la = la1
            } else { la = la0 }
            
            if (length(cm0[[colnam]][[dindname]][[cltype]])==0) {                 
              cm0[[colnam]][[dindname]] = cm0[[colnam]][[dindname]][!names(cm0[[colnam]][[dindname]])%in%cltype] 
              if (cltype%in%names(fm[[colnam]][[dindname]])) fm[[colnam]][[dindname]] = fm[[colnam]][[dindname]][!names(cm0[[colnam]][[dindname]])%in%cltype]
              next 
            }
            
            clt = NULL
            clt = cm0[[colnam]][[dindname]][[cltype]]$clt
            clt1 = cm0[[colnam]][[dindname]][[cltype]]$clt1
            
            #record clusterings with too many clusters and duplicates
            delcol_maxcl = apply(clt,2,function(x) return( length(unique(x))>maxcl))
            delcol_dup = duplicated(t(clt))
            col_dup = duplicateindM(clt)

            if (!cltype%in%names(fm[[colnam]][[dindname]]) | overwritef)  fm[[colnam]][[dindname]][[cltype]] = list()
            
            for (j in 1:ncol(clt)) {
              if (as.character(colnames(clt)[j])%in%names(fm[[colnam]][[dindname]][[cltype]])) next
              
              if (delcol_maxcl[j]) { # if too many clusters
                fm[[colnam]][[dindname]][[cltype]][[as.character(colnames(clt)[j])]] = rep(NA,51)
              } else if (delcol_dup[j]) { # if duplicated
                fm[[colnam]][[dindname]][[cltype]][[as.character(colnames(clt)[j])]] = fm[[colnam]][[dindname]][[cltype]][[as.character(colnames(clt)[col_dup[j]])]]
              } else {
              
              cl = clt[,j]
              cl1 = clt1[,j]
              
              # prepare for f measure traditional; use cl1 here
              ladf = t(sapply(la, function(lax) {
                a = rep(0,length(unique(la)))
                a[unique(la)==lax] = 1
                return(a)
              }))
              cldf = t(sapply(cl1, function(clx) {
                a = rep(0,length(unique(la)))
                a[unique(la)==clx] = 1
                return(a)
              }))
              colnames(ladf) = colnames(cldf) = unique(la)
              rownames(ladf) = rownames(cldf) = rownames(clt)
              
              #calculte f measure
              dd0 = d[[dind]]
              if (length(cl1)==length(la1)) dd0 = dd0[testind,testind]
              if (length(unique(cl1))==1) {
                cs1 = c()
                cs1u = rep(NA,10)
              } else {
                cs1 = cluster.stats(as.dist(dd0),cl1)
                cs1u = unlist(cs1)[c("cluster.number","average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")]
              }
              if (length(unique(cl))==1) {
                cs = c()
                csu = rep(NA,10)
              } else {
                cs = cluster.stats(as.dist(dd0),cl)
                csu = unlist(cs)[c("cluster.number","average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")]
              }
              
              cl1eval = unlist(c(F.measure.single.over.classes(ladf,cldf)$average[-6],
                                 adjustedRand(la,cl1),
                                 f.measure.comembership(la,cl1),                
                                 cs1u))
              if (length(unique(cl1))==1) cl1eval["cluster.number"] = 1
                cleval = unlist(c(adjustedRand(la,cl),
                                  f.measure.comembership(la,cl),
                                  csu))
              if (length(unique(cl))==1) cleval["cluster.number"] = 1
              if (identical(cl,cl1)) cleval[1:length(cleval)] = NA
              
              cl1silmed = clsilmed = NCA = NA
              if (!cltype%in%cltypesclass) {
                if (length(unique(cl1))>1) cl1silmed = median(silhouette(cl1,d[[dind]])[,3])
                if (length(unique(cl))>1 & !identical(cl,cl1)) clsilmed = median(silhouette(cl,d[[dind]])[,3])
              }
              if (cltype=="distmatrix") NCA = NCA_score(d[[dind]], la0, doUnderflow=doUnderflow)$p
              
              cl.oneClass=0; if (length(unique(cl1))==1) cl.oneClass=1
              cl1.oneClass=0; if (length(unique(cl))==1) cl1.oneClass=1
              cl1Equalcl=0; if (identical(cl1,cl)) cl1Equalcl=1
              eval = c(NCA=NCA,silmed=cl1silmed,cl1eval,
                       cl1.oneClass=cl1.oneClass,cl1Equalcl=cl1Equalcl,cl.oneClass=cl.oneClass,
                       silmed=clsilmed,cleval)
              if (cltype=="spec") {
                a = strsplit(colnames(clt)[j],"_")[[1]]
                b = gsub("silmed-","",a[grep("^silmed-",a)])
                eval["silmed"] = as.numeric(b)
                
                c = gsub("NCA-","",a[grep("^NCA-",a)])
                eval["NCA"] = as.numeric(c)
              }
              fm[[colnam]][[dindname]][[cltype]][[as.character(colnames(clt)[j])]] = eval
            }
            }
          } #cltype score
          cat(") ",sep="")
        }
        
        
        
        
        
        
        
        
        
        
        
        ## plot -------------------------------------------------------------
        if (doplot | (length(cm0[[colnam]][[dindname]])>0 & length(fm[[colnam]][[dindname]])>0)) {
          cat(" plotting (",sep="")
          cols = sapply(cm0[[colnam]][[dindname]],function(x)ncol(x$clt))
          cols[!names(cm0[[colnam]][[dindname]])%in%clplotallpar] = 1
          cols = cols+sum(plot_sill_rp)
          layoutm = matrix(0, ncol=max(cols), nrow=length(cm0[[colnam]][[dindname]]))
          for (cind in 1:length(cols)) { layoutm[cind,1:cols[cind]] = c(1:cols[cind])+suppressWarnings(max(layoutm)) }
          
          #split up layoutm if plots inside of it is too much
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
          
          pngname = paste0(plot_dist_eval_dir, "/tsne_", gsub(".Rdata","",fileNames(distmfile[i])),"_class_",interested[col],"-",splitby[col],"-",dindname,"__",paste(cltypes,collapse="-"),".png")
          if (!file.exists(pngname) | overwriteplot) {
            png(pngname, width=width*ncol(layoutm), height=height*nrow(layoutm))
            par(mar=c(7,2,20,2))
            layout(layoutm)
            
            # plotcount = 0
            # pngcount = 1
            
            for (cltype in names(cm0[[colnam]][[dindname]])) {
              start2 = Sys.time()
              cat(" ",cltype,sep="")
              clt = NULL
              par = clpar[cltype]
              if (cltype%in%cltypesclass) { la = la1
              } else { la = la0 }
              
              clt = cm0[[colnam]][[dindname]][[cltype]]$clt
              clt1 = cm0[[colnam]][[dindname]][[cltype]]$clt1
              par = clpar[cltype]
              
              #rtsne for plot
              mst = spantree(d[[dindname]])
              di0 = dt[[dindname]]
              di = di0[testind,]
              if (nrow(clt)>length(testind)) di = di0
              
              #only plot one par or all pars
              summetrics = unlist(lapply(fm[[colnam]][[dindname]][[cltype]], function(x) return(c(sum(x[c(11,12)]),sum(x[c(34,35)])))))
              mmax = which(summetrics==max(summetrics[!is.na(summetrics)]))[1]
              jbest = ceiling(mmax/2)
              jind = 1:ncol(clt)
              if (!cltype%in%clplotallpar) jind = jbest
              
              # pngname = paste0(plot_dist_eval_dir, "/tsne_", gsub(".Rdata","",fileNames(distmfile[i])),"_class_",interested[col],"-",splitby[col],names(d)[dind],"_",cltype,".png")
              # png(pngname, width=width*ceiling(sqrt(length(jind))), height=height*ceiling(sqrt(length(jind))))
              # par(mfrow=c(ceiling(sqrt(length(jind))),ceiling(sqrt(length(jind)))), mar=c(3,3,15,2), oma=c(2,2,5,2))
              
              main0 = paste0("flowcap-II AML, ",cltype,".",par,"=","parinsert","; ",interested[col],", ",splitby[col],"-",dindname," ;; (big0=cluster, 0=labeledcluster, .=label)")
              
              for (j in jind) {
                if (par!="") { main1 = gsub("parinsert",colnames(clt)[j],main0)
                } else { main1 = gsub(paste0(".",par,"=","parinsert"),"",main0) }
                main = paste0(main1,"\n\n",paste0(paste0("_",names(fm[[colnam]][[dindname]][[cltype]][[as.character(colnames(clt)[j])]]),"-"),
                                                  signif(fm[[colnam]][[dindname]][[cltype]][[as.character(colnames(clt)[j])]],2), collapse=""))
                main = gsub("cl1Equalcl","\n\ncl1Equalcl",main)
                main = gsub("_p-","\np-",main)
                main = gsub("_Rand-","\nRand-",main)
                main = gsub("_cluster.number-","\ncluster.number-",main)
                main = gsub("_pearsongamma-","\npg-",main)
                main = gsub("_cl1silmed","\n\ncl1silmed",main)
                
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
                lines(mst,di0,col='grey')
                for (g in 1:length(colour)) {
                  if (length(unique(clt1[,j]))>=g & fm[[colnam]][[dindname]][[cltype]][[as.character(colnames(clt)[j])]]["cl1Equalcl"]==0 & sum(is.na(clt1[,j]))==0)
                    points(di[clt1[,j]%in%unique(clt1[,j])[g],1],di[clt1[,j]%in%unique(clt1[,j])[g],2], col=colour[g], cex=3)
                  if (length(unique(clt[,j]))>=g) points(di[clt[,j]%in%unique(clt[,j])[g],1],di[clt[,j]%in%unique(clt[,j])[g],2], col=colour[g], cex=2)
                  if (length(unique(la0))>=g) points(di0[la0%in%unique(la0)[g],1],di0[la0%in%unique(la0)[g],2], col=colour[g], pch=16, cex=1)
                }
                legend("topleft", legend=laname[unique(clt[,j])], fill=colour[rep(c(1:g),g)], pch=rep(c(20,1:19,21:25),each=g)[c(1:g)])
                if (j==jbest & plot_sill_rp[1]) {
                  cls = clt1[,j]
                  if (mmax%%2==1) cls = clt[,j]
                  dls = as.dist(d[[dindname]])
                  if (length(cls)<nrow(d[[dindname]])) dls = as.dist(d[[dindname]][testind,testind])
                  
                  # plotcount = plotcount+1
                  if (length(unique(cls))>1) plot(silhouette(cls,dls))
                }
                
              }
              if (plot_sill_rp[2]) { #plot recall precision
                clstats = list()
                a = fm[[colnam]][[dindname]][[cltype]]
                sil_orig = sapply(a, function(x) return(x[which(names(x)=="silmed")]))
                cl1nc = sapply(a, function(x) return(x[which(names(x)=="cluster.number")[1]]))
                clstats$cl1P = sapply(a, function(x) return(x[which(names(x)=="P")]))
                clstats$cl1R = sapply(a, function(x) return(x[which(names(x)=="R")]))
                clstats$cl1p = sapply(a, function(x) return(x[which(names(x)=="p")[1]])); clstats$cl1p[is.na(clstats$cl1p)] = clstats$clp[is.na(clstats$cl1p)]
                clstats$cl1r = sapply(a, function(x) return(x[which(names(x)=="r")[1]])); clstats$cl1r[is.na(clstats$cl1r)] = clstats$clr[is.na(clstats$cl1r)]
                clnc = sapply(a, function(x) return(x[which(names(x)=="cluster.number")[2]]))
                clstats$clp = sapply(a, function(x) return(x[which(names(x)=="p")[2]])); clstats$clp[is.na(clstats$clp)] = clstats$cl1p[is.na(clstats$clp)]
                clstats$clr = sapply(a, function(x) return(x[which(names(x)=="r")[2]])); clstats$clr[is.na(clstats$clr)] = clstats$cl1r[is.na(clstats$clr)]
                
                # plotcount = plotcount+1
                plot(1, t='n', main="recall & precision (label=cluster.number; recall=dashed)",ylim=c(0,1),xlim=c(1,ncol(clt)),xaxt="n", xlab=clpar[cltype])
                colour = rainbow(1+length(clstats)/2)
                for (clstatsi in 1:length(clstats)) {
                  if (grepl("r",names(clstats)[clstatsi],ignore.case=T)) {
                    lines(clstats[[clstatsi]], col=colour[ceiling(clstatsi/2)],lty=2)
                    points(clstats[[clstatsi]], col=colour[ceiling(clstatsi/2)])
                  } else { 
                    lines(clstats[[clstatsi]], col=colour[ceiling(clstatsi/2)]) 
                    points(clstats[[clstatsi]], col=colour[ceiling(clstatsi/2)]) 
                  }
                  if (grepl("1",names(clstats)[clstatsi])) { text(clstats[[clstatsi]], labels=cl1nc, cex=2, pos=3)
                  } else { text(clstats[[clstatsi]], labels=clnc, cex= 2, pos=3) }
                }
                lines(sil_orig, col=colour[length(colour)])
                points(sil_orig, col=colour[length(colour)])
                
                legend("bottomleft", legend=c("cl1.standard","cl1","cl","silhouette med / truth"), fill=colour)
                if (clpar[cltype]!="") { axis(1, at=1:length(clstats[[1]]), labels=names(clstats$clr))
                } else { axis(1, at=1:length(clstats[[1]]), labels=1:length(clstats[[1]])) }
              }
              #mtext(text=main0,side=3,line=0,outer=T)
            }
            cat(")",sep="")
            graphics.off()
          }
        }
        
        
      } #dind
    } #col
    
    save(fm,file=paste0(dist_clusterf_dir,"/",distmfilenames[i]))
    save(cm0,file=paste0(dist_clustercl_dir,"/",distmfilenames[i]))
    TimeOutput(start1)
    
    return(list(redo=redo,error=F))
  }, error = function(e) {
    cat(paste("ERROR:  ",e)); graphics.off(); return(list(redo=redo,error=T))
  })
}
TimeOutput(start)

redo0 = sapply(redodist,function(x) x$redo)
errors = sapply(redodist,function(x) x$error)












## collate scores ---------------------------------------------------------------
distmfile1 = list.files(dist_clusterf_dir, recursive=F, full.names=T)
fs = file.size(distmfile1)
distmfile10 = distmfile1 = distmfile1[order(fs)]; fs = sort(fs)
distmfile1names = fileNames(distmfile1)
#head(distmfile1names[fs<30000])
dm = distMetafun(distmfile1,dis)
#distmfile1 = paste0(dist_dir,"/",distmfile1names)
dm$size = fs
# f = lapply(distmfile1, function(x) {
#   tryCatch({
#     get(load(x));
#   }, error = function(e) {
#     cat(paste(x," ERROR: ",e)); return(NULL)
#   })
# })

redodist = foreach(di=distmfile1,.combine='c') %dopar% {
  tryCatch({
    fm = get(load(di))
    return(F)
  }, error = function(e) {
    return(T)
  })
}

f = foreach(di=distmfile1) %dopar% { get(load(di)) }

#make a list of failed files
# failind = sapply(distmfile1, function(x) {
#   tryCatch({ get(load(x)); return(F)
#   }, error = function(e) { cat(paste(x," ERROR: ",e)); return(T) })
# })
# distmfile1names = fileNames(distmfile1[failind])
# distmfile1 = paste0(dist_dir,"/",distmfile1names)

distmfile1names = names(f) = fileNames(distmfile1)

#create distMeta
distMeta1 = distMetafun(distmfile1names,dis)

result = foreach(i = names(f)) %dopar% { #dist
  b = NULL
  a = NULL
  aind = 1
  for (j in names(f[[i]])) { #col
    for (k in names(f[[i]][[j]])) { #dind
      for (l in names(f[[i]][[j]][[k]])) { #cltype / remove sil
        if (l=="sil") next
        for (m in names(f[[i]][[j]][[k]][[l]])) {
          b = rbind(b, c(i,j,k,l,m))
          a[[aind]] = (f[[i]][[j]][[k]][[l]][[m]])
          aind = aind+1
        }
      }
    }
  }
  al = sapply(a,function(x) length(x))
  an = sapply(a,function(x) sum(names(x)=="")>0)
  alm = Mode(al[!an])
  al!=alm
  a = lapply(a,function(x) {
    if (length(x)==alm) return(x)
    if (length(x)<alm) return(c(x[1:which(names(x)=="" & is.na(x))[1]], rep(NA,alm-length(x)),  x[(which(names(x)=="")[1]+1):length(x)]))
    if (length(x)>alm) return(x[-which(names(x)=="" & is.na(x))[1:(length(x)-alm)]])
  })
  return(list(rn=b, score=Reduce('rbind',a), colnamea=names(f[[i]][[j]][[k]][[l]][[m]])))
}

rn = foreach(ii=1:length(result),.combine="rbind") %dopar% { result[[ii]]$rn }
score = foreach(ii=1:length(result),.combine="rbind") %dopar% { result[[ii]]$score }
aa = foreach(ii=1:length(result)) %dopar% { result[[ii]]$colnamea }
a = aa[[which(sapply(aa, function(x) length(x))==ncol(score))[1]]]

rownames(rn) = NULL
colnames(score) = gsub("average.","avg ",a)
colnames(score) = gsub("avg.","avg ",colnames(score))
colnames(score) = gsub("pearsongamma","pg",colnames(score))
score0 = cbind(distMeta1[match(rn[,1],distMeta1$path),-1],rn[,-1],score)
rownames(score0) = NULL
write.csv(score0,file=paste0(cl_score_result_dir,"_",paste(cltypes,collapse="-"),".csv"))
save(score0,file=paste0(cl_score_result_dir,"_",paste(cltypes,collapse="-"),".Rdata"))

TimeOutput(start)



## Convert to latex


require(xtable)

#pure distance only
pureind = !score0[,"rw"] & !score0[,"sim"] & !score0[,"weighted"] & 
  score0[,"layer"]==7 & score0[,"dist"]=="manhattan" &  score0[,"rand"]==0 & score0[,"3"]=="distmatrix" &
  ((grepl("effort|contrib|prop|pnratio",score0[,"type"]) & score0[,"norm"]=="none") | 
     (!grepl("effort|contrib|prop|pnratio",score0[,"type"]) & score0[,"norm"]=="cellpop"))

tubeind = grepl("^interested-tube",score0[,"1"])
amlind = grepl("^interested-aml",score0[,"1"])
noProp = grepl("Prop",score0[,"type"])
phenodev = grepl("effort|contrib|contrib|effort|TRIM",score0[,"type"])

distcolind = c("layer", "norm", "count", "weighted", "weightedorig", "dist","type", "rand", "rw", "sim")

#tube: Non-cluster scores
xt = xtable( score0[pureind & tubeind ,c("type","NCA","silmed","avg between","avg within","avg silwidth","pg","dunn","dunn2","entropy","wb.ratio","ch")] )
xt = xtable( score0[pureind & amlind & grepl("all",score0[,"2"]) ,c("type","NCA","silmed","avg between","avg within","avg silwidth","pg","dunn","dunn2","entropy","wb.ratio","ch")] )





# 
# #delete what's not in dist and is in cl and f; temporary
# distmfile = list.files(dist_dir, recursive=F, full.names=T, pattern=".Rdata$")
# distmfilenames = fileNames(distmfile)
# distmfilef = list.files(dist_clusterf_dir, recursive=F, full.names=T, pattern=".Rdata$")
# distmfilefnames = fileNames(distmfilef)
# file.remove(distmfilef[!distmfilefnames%in%distmfilenames])
# distmfilecl = list.files(dist_clustercl_dir, recursive=F, full.names=T, pattern=".Rdata$")
# distmfileclnames = fileNames(distmfilecl)
# file.remove(distmfilecl[!distmfileclnames%in%distmfilenames])
# 



# # delete all but most recent kern files
# kerndistfiles = list.files("result/dist/dist_kern",full.names=T)
# kerndistinfo = file.info(kerndistfiles,full.names=T)
# 
# kerndistfiles = kerndistfiles[order(kerndistinfo[,"mtime"])]
# kerndistinfo = kerndistinfo[order(kerndistinfo[,"mtime"]),]
# 
# a = str_split(kerndistfiles,"_")
# a = sapply(a, function(x) paste(x[!grepl("sigma",x)],collapse="_"))
# file.remove(kerndistfiles[!duplicated(a)])
# 
# # adupind = duplicateind(a)
