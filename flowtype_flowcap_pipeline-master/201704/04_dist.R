# Uses different distance measures to calculate distance & plot samples
# aya43@sfu.ca 20161220

root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/matrix", sep="")
# mcp_types = c("/original", "/child_pn","/child_prop", "/trim") 
# matrix_type1 = c("CountAdj") #, "Prop", "Count",  #countAdj comes first, to set which columns are deleted based on cellcountthres
# matrix_type2 = c("Child_prop")
# matrix_type3 = c("Child_pnratio", "Child_entropy")
# matrix_type4 = c("Child_pnratio", "Child_entropy")
matrix_count = c("CountAdj")
# matrix_type = c("Child_entropyTRIM_CountAdj", "Parent_entropyTRIM_CountAdj", "LogFoldTRIM_CountAdj", "PvalTRIM_CountAdj",
#                 "Child_entropyTRIM_Prop", "Parent_entropyTRIM_Prop", "LogFoldTRIM_Prop", "PvalTRIM_Prop",
#                 "Parent_effortTRIM_CountAdj", "Parent_contribTRIM_CountAdj", "Child_pnratioTRIM_CountAdj", "Child_propTRIM_CountAdj",
#                 "Parent_effortTRIM_Prop", "Parent_contribTRIM_Prop", "Child_pnratioTRIM_Prop", "Child_propTRIM_Prop",
#                 "LogFold_CountAdj", "CountAdj", "Pval_CountAdj", "Parent_entropy", "Child_entropy",
#                 "Parent_effort_CountAdj", "Parent_effort_Prop", "Parent_contrib_CountAdj", "Parent_contrib_Prop", "Child_pnratio", "Child_prop")
matrix_type = c("Child_entropyTRIM_CountAdj", "Child_entropyTRIM_Prop", "Parent_entropyTRIM_CountAdj", "Parent_entropyTRIM_Prop", 
                "LogFoldTRIM_CountAdj", "LogFoldTRIM_Prop", "PvalTRIM_CountAdj", "PvalTRIM_Prop",
                "CountAdj", "Parent_entropy", "Child_entropy",
                "LogFold_CountAdj", "LogFold_Prop", "Pval_CountAdj", "Pval_Prop",
                "Parent_effortTRIM_CountAdj", "Parent_effortTRIM_Prop", "Parent_contribTRIM_CountAdj", "Parent_contribTRIM_Prop", 
                "Child_pnratioTRIM_CountAdj", "Child_pnratioTRIM_Prop", "Child_propTRIM_CountAdj", "Child_propTRIM_Prop",
                "Parent_effort_CountAdj", "Parent_effort_Prop", "Parent_contrib_CountAdj", "Parent_contrib_Prop")
#"Child_pnratio", "Child_prop")
# matrix_type = c("Parent_effortTRIM", "Parent_contribTRIM", "Child_pnratioTRIM", "Child_propTRIM")
phenoChild_dir = paste(result_dir, "/phenoChild.Rdata",sep="")
phenoChild_ind_dir = paste(result_dir, "/phenoChild_ind.Rdata",sep="")
phenoChildpn_dir = paste(result_dir, "/phenoChildpn.Rdata",sep="")
phenoChildpn_ind_dir = paste(result_dir, "/phenoChildpn_ind.Rdata",sep="")
cellCountThres = c(1200) #insignificant if count under

overwrite=T

#Output
dist_dir = paste(result_dir, "/dist", sep=""); for (i in 1:length(dist_dir)) { suppressWarnings(dir.create (dist_dir[i])) }
dist_type_dir = NULL
# for(i in 1:length(mcp_types)) { for(j in 1:length(dist_dir)) { dist_type_dir = append(dist_type_dir, paste0(dist_dir[j], mcp_types[i]))} }
# for (i in 1:length(dist_type_dir)) { suppressWarnings(dir.create (dist_type_dir[i])) }

libr(stringr)
libr(colorspace)
libr(vegan) # libr(proxy)
libr(foreach)
libr(doMC)
libr(lubridate) #if there are date variables
source("~/projects/IMPC/code/_funcAlice.R")


dis = c("manhattan")#,"euclidean", "canberra")#, "mahalanobis") #
#dis = c("euclidean", "canberra", "binomial", "cao") # "jaccard","mahalanobis", "gower","altGower","morisita","horn", "manhattan", "maximum", "binary", "minkowski", "raup", "mountford" (all at same level),)
disnoavg = c("euclidean", "manhattan", "canberra", "bray", "kulczynski", "morisita", "horn", "binomial")  #dis measures that don't average over number of features
disindist = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski") #use dist function when possible, instead of vegdist
disnoneg = c("canberra")

no_cores = 8#detectCores()-1
registerDoMC(no_cores)








start = Sys.time()


#for (ci in 4:1) {

phenoMeta = get(load(phenoMeta_dir))
#k0 = c(1,3,max(phenoMeta$phenolevel)) # how many markers to consider i.e. k=max(phenolevel) only
k0 = c(1,4,max(phenoMeta$phenolevel)) # how many markers to consider i.e. k=max(phenolevel) only
m0 = get(load(paste0(matrix_dir, matrix_count,".Rdata")))

#load different matrices
# for (mcp in matrix_type_) {
a = foreach(mcp=matrix_type) %dopar% {
  cat("\n", mcp, " ",sep="")
  start2 = Sys.time()
  leavePhenotype = list()
  doneAll = F
  if (!file.exists(paste0(matrix_dir, mcp,".Rdata"))) {cat("doesn't exist"); return(NULL)}
  # mm = get(load(paste0(matrix_dir, mcp,".Rdata")))
  # 
  # if (!is.null(dim(mm))) {
  #   mpi = colnames(mm)
  #   gt = rownames(mm)
  #   # print(mm[50:60,20:30])
  # } else {
  #   mpi = union(names(mm), Reduce('union',lapply(mm,function(y) return(colnames(y)))) )
  #   gt = rownames(mm[[1]])
  #   # print(head(mm[[25]]))
  # }
  
  mresult = Loadintermatrices(paste0(matrix_dir, mcp,".Rdata"))
  mml0 = mresult$mml
  mmlname = names(mml0)
  pt = mpi = mresult$pt
  gt = mresult$gt
  
  #get to-delete low count phenotype indices; CountAdj should be first one
  for (countThres in cellCountThres) {cat("\ncountThres: ",countThres," >",sep="")
    
    #get to-delete high no of marker phenotypes
    for (k in k0) { cat(" level",k," ",sep="")
      
      mmlresult = trimMatrices(mml0,m0,pt,gt,phenoMeta,leavePhenotype,doneAll, countThres,k)
      if (is.null(mmlresult)) return(NULL)
      m = mmlresult$mml[[1]]
      pm = mmlresult$pm
      leavePhenotype = mmlresult$leavePhenotype
      doneAll = mmlresult$doneAll
      
      if (is.null(m)) return(NULL)
      
      
      # For every distance type --------------------------------------------------------
      
      loop.ind = 1:length(dis)
      a = match(disnoneg,dis)
      if (sum(!is.na(a))>0) loop.ind = loop.ind[-a]
      
      tryCatch({
        
        #for (i in 1:length(dis)) {
        # foreach(i=loop.ind) %dopar% { #for each phenotype
        for (i in loop.ind) {
          cat(", ", length(dis)-i+1, ":", dis[i], " ", sep="")
          
          normalize = c("none","cellpop", "layer") # by none (for child matrices only), cell pop, layer
          
          #assume after splitting dist filename by "_", distance is second element
          if (doneAll) dname = paste(dist_dir, "/", mcp, "_", dis[i], "_FULL_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_normalize-", sep="")
          dname = paste(dist_dir, "/", mcp, "_", dis[i], "_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_normalize-", sep = "" )
          
          
          
          #try({ if (overwrite | !file.exists(dname)) {
          if (is.null(dim(m))) { #give each layer same weight, m is a list
            if ("none"%in%normalize) { #none
              dname0 = paste0(dname, "none.Rdata", sep="")
              if (overwrite | !file.exists(dname0)) {
                a = Reduce('cbind',m)
                if (sum(disindist%in%dis[i])>0) {
                  d = dist(a, method=dis[i])
                } else {
                  d = vegdist(a, method=dis[i])
                }
                save(d, file=dname0)
                write.csv(as.matrix(d), file=paste0(checkm(d,dname0),".csv"))
              }
            } 
            if ("cellpop"%in%normalize | "layer"%in%normalize) { #cell pop
              if ("cellpop"%in%normalize) {
                dc = matrix(0,nrow=nrow(m[[1]]),ncol=nrow(m[[1]]))
                dname1 = paste0(dname, "cellpop.Rdata", sep="")
              }
              if ("layer"%in%normalize) {
                dl = lapply(unique(pm$phenolevel), function(y) return(matrix(0, nrow=nrow(m[[1]]), ncol=nrow(m[[1]]))) )
                names(dl) = as.character(unique(pm$phenolevel))
                dname2 = paste0(dname, "layer.Rdata", sep="")
              }
              if (overwrite | !file.exists(dname1) | !file.exists(dname2)) {
                for (x in 1:length(m)) {
                  if (sum(disindist%in%dis[i])>0) {
                    dx = as.matrix(dist(m[[x]],method=dis[i]))
                  } else {
                    dx = as.matrix(vegdist(m[[x]],method=dis[i]))
                  }
                  if (sum(disnoavg%in%dis[i])>0) dx = dx/ncol(m[[x]])
                  if ("cellpop"%in%normalize) dc = dc+dx
                  if ("layer"%in%normalize) dl[[as.character(pm$phenolevel[x])]] = dl[[as.character(pm$phenolevel[x])]]+dx #1 indexing (layer starts at 0)
                }
                if ("cellpop"%in%normalize) {
                  save(dc, file=dname1)
                  write.csv(as.matrix(dc), file=paste0(checkm(d,dname1),".csv"))
                }
                if ("layer"%in%normalize) {
                  dl = lapply(unique(pm$phenolevel), function(x) {
                    a = dl[[as.character(x)]]
                    a[a>0] = a[a>0]/sum(pm$phenolevel==x)
                    return(a)
                  })
                  a = Reduce('+',dl)
                  a[a>0] = a[a>0]/length(dl)
                  dl = as.dist(a)
                  save(dl, file=dname2)
                  write.csv(as.matrix(dl), file=paste0(checkm(d,dname2),".csv"))
                }
              }
            }
            
          } else {
            if ("cellpop"%in%normalize) {
              dname1 = paste0(dname, "cellpop.Rdata")
              if (overwrite | !file.exists(dname1)) {
                if (sum(disindist%in%dis[i])>0) {
                  d = dist(m, method=dis[i])
                } else {
                  d = vegdist(m, method=dis[i])
                }
                save(d, file=dname1)
                write.csv(as.matrix(d), file=paste0(checkm(d,dname1),".csv"))
              }
            }
            if ("layer"%in%normalize) { #cell pop
              dname2 = paste0(dname, "layer.Rdata", sep="")
              if (overwrite | !file.exists(dname2)) {
                layers = unique(pm$phenolevel)
                #if (sum(layers==0)>0) layers = layers[-which(layers==0)]
                if (sum(disindist%in%dis[i])>0) {
                  dd = lapply(layers, function(x) return(as.matrix(dist(m[,which(pm$phenolevel==x)], method=dis[i]))) )
                } else {
                  dd = lapply(layers, function(x) return(as.matrix(vegdist(m[,which(pm$phenolevel==x)], method=dis[i]))) )
                }
                names(dd) = as.character(layers)
                #dis measures that don't average over number of features
                if (sum(disnoavg %in% dis[i]) > 0) dd = lapply(layers, function(x) {
                  a = dd[[as.character(x)]]
                  a[a>0] = a[a>0]/sum(pm$phenolevel==x)
                  return(a)
                })
                a = Reduce('+',dd)
                a[a>0] = a[a>0]/length(dd)
                d = as.dist(a)
                save(d,file=dname2)
                write.csv(as.matrix(d), file=paste0(checkm(d,dname2),".csv"))
              }
            }
          }
        }
      }, error = function(err) {
        cat(paste("ERROR:  ",err))
      })
      
    }
  }
  TimeOutput(start2)
}



TimeOutput(start)