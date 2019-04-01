# measure distance between distance matrices
# aya43@sfu.ca 20170419

root = "~/projects/flowCAP-II"
result_dir = "result"
setwd(root)

options(stringsAsFactors=FALSE)
#options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/matrix", sep="")
dist_dir = paste(result_dir, "/dist", sep="")
plot_dir = paste(result_dir, "/plots", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dir[i])) }
dist_score_dir = paste(dist_dir, "/dist_score", sep=""); for (i in 1:length(dist_score_dir)) { suppressWarnings(dir.create(dist_score_dir[i])) }
dist_clusterf_dir = paste(dist_score_dir, "/dist_cluster_f", sep=""); for (i in 1:length(dist_clusterf_dir)) { suppressWarnings(dir.create(dist_clusterf_dir[i])) }
cl_score_result_dir = paste0(dist_score_dir,"/score_cl_list")

matrix_type_features = c("CountAdj","Prop","LogFold","Pval","Child_entropy","Parent_entropy","Child_prop","Child_pnratio","Parent_contrib","Parent_effort","Freqp0.5")
plotcols = c("feature","layer","featend","dist","norm")
splitcols = c("none","rand")

metrics = c("silmed","NCA","f","r","p","Rand","pearsongamma","dunn","dunn2","entropy","ch")
clustnumbefore = c(T,T,T,T,T,T,F,F,F,F,F) #before is accuracy, after is cluster quality
cltypes = c("distmatrix","knn","kmed","lv","spec1","spec","hc","dc1","dc")#,"spec","dc") #"rw1" (produces one cluster...)
clplotallpar = c("hc")#plot only best

interested = c("tube","aml") #sampleMeta column; choose one
# splitby = c("none","tube") #match with above
plotdist = c("feature","norm","dist","layer")
plotProp = "_noProp" #"all"
plotNorm = "_noNorm"
plotFreqp = "_noFreqp"
plotOrig = "_orig"

matrix_count = c("CountAdj")

overwritedist = F
overwritedistd = F
overwritedistf = F
recalcNCA = F

quantilnorm = T

trimonly = T
ignoredist = ".csv|simmatrix"

dis = c("euclidean", "maximum", "manhattan", "canberra", "binary", "bray", "kulczynski", "minkowski", "morisita", "horn", "binomial", "rbf","poly","vanilla","tanh","laplace","bessel","anova","spline")


maxcl = 110 # max number of clusters, else evaluation is slow
maxpl = 60 # max number of plots per png
maxcol = 9

#Output
plot_compare_dir = paste(plot_dir, "/dist_compare", sep=""); for(i in 1:length(plot_compare_dir)) { suppressWarnings(dir.create(plot_compare_dir[i])) }
dist_comp_dir = paste(result_dir, "/dist_compare", sep="")
dist_comp_dist_dir = paste(result_dir, "/dist_compare_dist", sep="")
dist_comp_score_dir = paste(dist_score_dir, "/dist_compare", sep="")
dist_kern_dir = paste(dist_dir, "/dist_kern", sep=""); for (i in 1:length(dist_kern_dir)) { suppressWarnings(dir.create(dist_kern_dir[i])) }



source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
libr("stringr")
libr("foreach")
libr("doMC")
libr("Rtsne")
libr("devtools")
libr("Biobase")
libr("preprocessCore")
libr("kernlab")


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


sampleMeta = get(load(sampleMeta_dir))
phenoMeta = get(load(phenoMeta_dir))



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
  
  metrics = c("silmed","NCA","f","r","p","Rand","pearsongamma","dunn","dunn2","entropy","ch")
  start = Sys.time()
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
              if (length(aa)>0) { if (sum(!is.na(aa))>0) {
                a = mean(aa[!is.na(aa)]) 
                silm[[paste0(col,"_",cltype,"_",metric)]][paste0("dind-avg","_",par0)] = a
              }}
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
  silm1 = lapply(names(result[[which.max(sapply(result, function(x) length(x)))]]), function(i) {
    aa3 = foreach(x=result) %dopar% {
      if (!i%in%names(x)) return(NULL)
      a = x[[i]]
      names(a) = names(x[[i]])
      return(a)
    }
    goodind = (vapply(aa3, Negate(is.null), NA))
    if (sum(goodind)==0) next
    if (length(goodind)<length((distmfile))) goodind = append(goodind, rep(F,length(distmfile)-length(goodind)))
    colnamess = names(aa3[[which(goodind)[1]]])
    rownamess = distmfile[goodind]
    aa0 = aa3[goodind]
    aanames = sort(Reduce('union',lapply(aa0, function(x) names(x))))
    aa1 = sapply(aa0, function(x) x[match(aanames,names(x))]) ## some are missing because repeating clusterings are deleted -- set these to min!!
    aa1[is.na(aa1)] = min(aa1[!is.na(aa1)])/2
    if (is.null(dim(aa1))){ aa2 = matrix(aa1,nrow=1); rownames(aa2) = aanames; colnames(aa2) = names(aa1); aa1 = aa2 }
    aa3 = t(aa1)
    colnames(aa3) = aanames
    rownames(aa3) = rownamess
    return(aa3)
  })
  save(silm1,file=paste0(dist_comp_score_dir,".Rdata"))
  
} else {
  silm1 = get(load(paste0(dist_comp_score_dir,".Rdata")))
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

#recalculate NCA

if (recalcNCA) {
  start=Sys.time()
  
  # kernfilenames0 = list.files(dist_dir,full.names=F)
  # kernfilenames0 = kernfilenames0[!grepl("simmatrix",kernfilenames0)]
  # kernfilenames = str_split(kernfilenames0,"_")
  # kernfilenames = lapply(kernfilenames, function(x) c(paste(x[!grepl("sigma|kern|rbf|interested|splitby|Rdata",x)],collapse="_"), gsub("interested-","",x[grepl("interested",x)]), gsub("splitby-tube-|splitby-none-","",x[grepl("splitby",x)])  ))
  # kernfilenames = Reduce('rbind',kernfilenames)
  
  loop.ind=unique(distMeta$path)
  a = foreach (i = loop.ind, .combine="rbind") %dopar% {
    dnca = get(load(i))
    tubescore = NCA_score(dnca, sampleMeta$tube[match(colnames(as.matrix(dnca)),sampleMeta$fileName)])$p
    amlscore = NULL
    distt2 = list()
    for (x in 1:7) {
      aa = sampleMeta[match(rownames(as.matrix(dnca)),sampleMeta$fileName),c("tube","aml")]
      distt2[[x]] = as.dist(as.matrix(dnca)[aa[,"tube"]==x,aa[,"tube"]==x])
      amlscore[x] = (NCA_score(distt2[[x]], aa[aa[,"tube"]==x,"aml"])$p)
    }
    dtemp = distt2
    specimenavg = Reduce(intersect,lapply(distt2,function(x) gsub("T[0-9]S|FT","",colnames(as.matrix(x)))))
    dtemp0 = lapply(1:length(distt2),function(x) {
      smtemp = sampleMeta[match(rownames(as.matrix(distt2[[x]])),sampleMeta[sampleMeta$tube==x,"fileName"]),] #order sm[[x]] --> d[[x]]
      inds = match(as.numeric(specimenavg),smtemp[,"specimen"])
      dtemp = as.matrix(distt2[[x]])[inds,inds] #order sm/d[[x]] --> specimenavg
      return(dtemp)
    })
    distt2[["all"]] = Reduce("+", dtemp0) / length(dtemp0)
    
    amlscore = append(amlscore,NCA_score(distt2[["all"]], sampleMeta[match(rownames(as.matrix(distt2[["all"]])),sampleMeta[,"fileName"]),"aml"])$p)
    amlscore = append(amlscore,mean(amlscore[1:7]))
    
    return(c(tubescore,amlscore))
  }
  amlscores = a[,2:ncol(a)]
  tubescores = matrix(a[,1],ncol=1)
  rownames(tubescores) = rownames(amlscores) = unique(distMeta$path)
  colnames(tubescores) = "dind-all_none"
  colnames(amlscores) = paste0("dind-",c(1:7,"all","avg"),"_none")
  silm1[["aml_distMatrix_NCA"]] = amlscores
  silm1[["tube_distMatrix_NCA"]] = tubescores
  
  save(silm1,file=paste0(dist_comp_score_dir,".Rdata"))
  TimeOutput(start)
  
  
  
  
  
  
  start=Sys.time()
  
  m0 = get(load(paste0(matrix_dir, matrix_count,".Rdata")))
  
  matrix_type = list.files(path=result_dir,pattern=glob2rx("matrix*.Rdata"))
  matrix_type = matrix_type[!matrix_type=="matrixCount.Rdata"]
  matrix_type = gsub("matrix","",matrix_type)
  matrix_type = gsub(".Rdata","",matrix_type)
  matrix_type = matrix_type[!grepl("KO|Max",matrix_type)]
  matrix_type = matrix_type[!grepl("Freqp_orig",matrix_type)]
  loop.ind = matrix_type
  b = foreach (mcp = loop.ind, .combine="rbind") %dopar% {
    bb = NULL
    if (!file.exists(paste0(matrix_dir, mcp,".Rdata"))) {cat("doesn't exist"); return(NULL)}
    
    mresult = Loadintermatrices(paste0(matrix_dir, mcp,".Rdata"))
    mml0 = mresult$mml
    mmlname = names(mml0)
    pt = mpi = mresult$pt
    gt = mresult$gt
    
    #get to-delete low count phenotype indices; CountAdj should be first one
    countThres = 1200
    if (grepl("TRIM",mcp)) { k0=c(4,7)
    } else { k0=c(1,4,7) }
    #get to-delete high no of marker phenotypes
    for (k in k0) { cat(" level",k," ",sep="")
      dnamee = paste(dist_kern_dir, "/", mcp, "_rbf_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_normalize-none",sep="")
      dname = paste(dnamee,"_interested-", sep = "" )
      
      mmlresult = trimMatrices(mml0,m0,pt,gt,phenoMeta,NULL,F, countThres,k)
      if (is.null(mmlresult)) return(NULL)
      m = mmlresult$mml[[1]]
      if (is.null(dim(m))) m = Reduce('cbind',m)
      m = as.matrix(m)
      m = m[!apply(m[,-1], 1, function(x) all(x==0)),!apply(m[-1,], 2, function(x) all(x==0))]
      pm = mmlresult$pm
      
      labelss = sampleMeta$tube[match(rownames(m),sampleMeta$fileName)]
      labelssaml = sampleMeta$aml[match(rownames(m),sampleMeta$fileName)]
      labelssaml = lapply(unique(labelss),function(x) (labelssaml[labelss==x]))
      
      
      if (length(unique(labelss))<7) { next
      } else if (min(sapply(labelssaml,function(x)length(unique(x))))<2) { next  #if not both aml and healthy in each tube
      } else if (length(Reduce('intersect',lapply(labelssaml,function(x) x[duplicated(x)])))<2) { next } #each tube has at least two heathy and two aml
      
      sp = specc(x=m,kernel="rbfdot",centers=7)
      parlist0 = unlist(sp@kernelf@kpar)
      sim = kernelMatrix(kernel=rbfdot(as.numeric(parlist0)), x=m)
      distt1 = get_graphd(sim)
      dname0 = paste0(dname,"tube_splitby-none-all")
      save(sim,file=paste0(dname0,"_simmatrix.Rdata"))
      save(distt1,file=paste0(dname0,"_dist.Rdata"))
      tubescorespec = NCA_score(distt1, labelss)$p
      amlscorespec = NULL
      distt2 = list()
      
      for (t in 1:7) {
        mind = rownames(m)%in%sampleMeta$fileName[sampleMeta$tube==t]
        labelss = sampleMeta$aml[match(rownames(as.matrix(m))[mind],sampleMeta$fileName)]
        sp = specc(x=m[mind,],kernel="rbfdot",centers=7)
        parlist0 = unlist(sp@kernelf@kpar)
        sim = kernelMatrix(kernel=rbfdot(as.numeric(parlist0)), x=m[mind,])
        distt2t = distt2[[t]] = get_graphd(sim)
        dname0 = paste0(dname,"aml_splitby-tube-",t)
        save(sim,file=paste0(dname0,"_simmatrix.Rdata"))
        save(distt2t,file=paste0(dname0,"_dist.Rdata"))
        amlscorespec[t] = NCA_score(distt2t, labelss)$p
      }
      
      dtemp = distt2
      specimenavg = Reduce(intersect,lapply(distt2,function(x) gsub("T[0-9]S|FT","",colnames(as.matrix(x)))))
      dtemp0 = lapply(1:length(distt2),function(x) {
        smtemp = sampleMeta[match(rownames(as.matrix(distt2[[x]])),sampleMeta[sampleMeta$tube==x,"fileName"]),] #order sm[[x]] --> d[[x]]
        inds = match(as.numeric(specimenavg),smtemp[,"specimen"])
        dtemp = as.matrix(distt2[[x]])[inds,inds] #order sm/d[[x]] --> specimenavg
        return(dtemp)
      })
      distt2[["all"]] = Reduce("+", dtemp0) / length(dtemp0)
      
      #all, avg
      amlscorespec[8] = NCA_score(distt2[["all"]], sampleMeta$tube[match(colnames(as.matrix(distt2[["all"]])),sampleMeta$fileName)])$p
      amlscorespec[9] = mean(amlscorespec[1:7])
      
      bb = rbind(bb,c(dnamee,tubescorespec,amlscorespec))
    }
    return(bb)
  }
  
  tubescorespecs = as.matrix(b[,2],ncol=1)
  amlscorespecs = b[,3:ncol(b)]
  rownames(tubescorespecs) = rownames(amlscorespecs) = b[,1]
  colnames(tubescorespecs) = "dind-all_none"
  colnames(amlscorespecs) = paste0("dind-",c(1:7,"all","avg"),"_none")
  silm1[["aml_spec_NCA"]] = amlscorespecs
  silm1[["tube_spec_NCA"]] = tubescorespecs
  
  save(silm1,file=paste0(dist_comp_score_dir,".Rdata"))
  TimeOutput(start)
  
}









## Plotting ---------------------------------------------------------


cat(", plotting... ")
start=Sys.time()

silm1 = silm11 = get(load(paste0(dist_comp_score_dir,".Rdata")))

leaverows = rep(T,nrow(distMeta))
if (plotOrig=="_orig") leaverows = leaverows & !distMeta[,"rw"] & !distMeta[,"sim"] & !distMeta[,"weighted"] & distMeta[,"layer"]==7 & distMeta[,"dist"]=="manhattan" # & distMeta[,"featend"]==""
if (plotProp=="_noProp") leaverows = leaverows & !grepl("Prop",distMeta[,"type"],ignore.case=F)
if (plotFreqp=="_noFreqp") leaverows = leaverows & !grepl("Freqp",distMeta[,"feature"],ignore.case=F)
if (plotNorm=="_noNorm") leaverows = leaverows & ((grepl("effort|contrib|prop|pnratio",distMeta[,"feature"]) & distMeta[,"norm"]=="none") | 
                                                    (!grepl("effort|contrib|prop|pnratio",distMeta[,"feature"]) & distMeta[,"norm"]=="cellpop"))
silm1 = lapply(1:length(silm11),function(x) {
  indss = leaverows[match(rownames(silm11[[x]]),distMeta$path)]
  if (length(indss)==0) return(NULL)
  y = silm11[[x]][indss[!is.na(indss)],]
  if (is.null(dim(y))) { y = matrix(y,ncol=1); rownames(y) = rownames(silm11[[x]])[indss[!is.na(indss)]]; colnames(y) = colnames(silm11[[x]])}
  return(y)
})
names(silm1) = names(silm11)[1:length(silm1)]
silm1 = silm1[!sapply(silm1, is.null)]


# colour = rainbow(max(sapply(plotdist,function(x)length(unique(distMeta[,x])))))
colour = rainbow(max(length(unique(distMeta[leaverows,"feature"])),length(unique(distMeta[leaverows,"layer"]))))
cp1 = c(15:18,7:14)[1:length(unique(distMeta[leaverows,"featend"]))]
cp2 = c(0,1,5,2,6)[length(unique(distMeta[leaverows,"dist"]))] #empty inside


metrics = sapply(str_split(names(silm1),"_"), function(x) {
  if (x[3]!="NCA") return(x[3])
  return(paste0(x[1],"_",x[3]))
})
maxvals = sapply(unique(metrics), function(x) max(sapply(silm1[metrics==x], function(y) max(y[!is.na(y) & y!=Inf]))))
minvals = sapply(unique(metrics), function(x) min(sapply(silm1[metrics==x], function(y) min(y[!is.na(y) & y!=Inf]))))
#plotno = 4+length(plotcols)
plotno = length(plotcols)

a=foreach (colscorei = 1:length(silm1)) %dopar% {
  colscore = names(silm1)[colscorei]
  metric = metrics[colscorei]
  silm2 = silm1[[colscore]]; silm2[is.na(silm2)] = min(silm2[!is.na(silm2)])/2
  dm2 = distMeta[match(rownames(silm2),distMeta$path),]
  dm2[dm2[,"weighted"],"dist"] = "weighted manhattan"
  
  
  
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
      
      if (length(unique(as.vector(silm3)))<2) next
      
      #if too many columns, split into seperate png's
      silm30 = silm3
      silm4 = list()
      pngno = ""
      if (ncol(silm3)>maxcol) {
        for (silm3col in 1:ceiling(ncol(silm3)/maxcol)) {
          b = silm3[,1:min(maxcol,ncol(silm3))]; if (is.null(dim(b))) { b = matrix(b,ncol=1); rownames(b) = rownames(silm3); colnames(b) = colnames(silm3)[1:min(maxcol,ncol(silm3))]}
          silm4[[silm3col]] = b
          
          if (!ncol(silm3)>maxcol) break
          a = silm3[,-c(1:min(maxcol,ncol(silm3)))]; if (is.null(dim(a))) { a = matrix(a,ncol=1); rownames(a) = rownames(silm3); colnames(a) = colnames(silm3)[-c(1:min(maxcol,ncol(silm3)))] }
          silm3 = a
        }
        pngno = paste0("_",1:length(silm4))
      } else { silm4[[1]] = silm3 }
      
      
      
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
      
      
      for (pngind in 1:length(silm4)) {
        silm5 = silm4[[pngind]]
        png(filename=paste0(plot_compare_dir,"/", colscore,"__",splitcol,"-",names(splitcolc)[splitcolind],pngno[pngind],plotProp,plotNorm,plotFreqp,plotOrig,".png"), width=700*(ncol(silm5)), height=700*plotno)
        par(mfcol=c(plotno,ncol(silm5)), mar=c(10,5,5,5))
        main1 = paste0("Distance Parameter vs Quality Metric (split by ",splitcol,"-",names(splitcolc)[splitcolind],")\n")
        for (columnplot in 1:ncol(silm5)) {
          plotall = 1
          #if end of a set of tubes, plot on all 7 tubes together
          if (columnplot%in%grep("dind-avg",colnames(silm5))) {
            silm3col = which(colnames(silm3)==colnames(silm5)[columnplot])
            silm0 = unlist(sapply((silm3col-9):(silm3col-2), function(x) silm3[,x]))
            dm00 = do.call("rbind", replicate(7, dm, simplify=F))
          }
          # if (sum(silm0==Inf)>0) print(paste0(plot_compare_dir,"/", colscore,"__",splitcol,"-",names(splitcolc)[splitcolind],pngno[pngind],plotProp,plotNorm,".png"))
          silm0plot = silm0-min(silm0[!is.na(silm0)])
          silm0plot = silm0plot/max(silm0plot[!is.na(silm0plot)]) +.1
          
          for (testcol in plotcols) {
            dmc0 = as.integer(factor(dm00[,testcol]))
            dmc = jitter(dmc0, factor=1)
            dmcl = lapply(sort(unique(dmc0)), function(x) silm0[dmc0==x])
            names(dmcl) = sort(unique(dm00[,testcol]))
            #Jitter
            
            boxplot(dmcl, lwd = 1, outline=T, ylim=c(minvals[metric],maxvals[metric]),
                    main=paste0(main1,"_",testcol," vs ", colscore,"; Only distances of ",splitcol,"-",names(splitcolc)[splitcolind], "\n", colnames(silm3)[columnplot],"-par ;; ", testcol), 
                    xaxt ='n', yaxt='n', ylab=colscore, cex.lab=1.5) #,xaxt ='n', xlab=testcol
            axis(1, at=as.integer(factor(unique(dm00[,testcol]))), labels=unique(dm00[,testcol]), cex.axis = 1.5, las=2)
            axis(2, cex.axis = 1.5)
            title(cex.main=1.5)
            points(dmc,silm0,col=colour[as.integer(factor(dm00[,"feature"]))], pch=16,cex=exp(silm0plot))
            if (sum(dm00$weighted)>0) points(dmc[dm00$weighted],silm0[dm00$weighted],col="black",cex=exp(silm0plot[dm00$weighted]+.75))
            if (sum(dm00$weightedorig)>0) points(dmc[dm00$weightedorig],silm0[dm00$weightedorig],col="blue",cex=exp(silm0plot[dm00$weightedorig]+.75))
            
          } 
          
          
        }
        
        graphics.off()
      }
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
  points(pc$x[distMeta$weightedorig,1],pc$x[distMeta$weightedorig,2],col="blue",cex=exp((2*silm0/max(silm0))[dm$weightedorig]+1))
  points(0, 0, pch = 3, cex = 4, lwd = 4)
}
graphics.off()




TimeOutput(start)








