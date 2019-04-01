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

sampletimes = 5

matrix_count = c("CountAdj")
matrix_type = c(paste0("Child_entropyTRIM_CountAdj_",c(1:sampletimes)), paste0("Parent_entropyTRIM_CountAdj_",c(1:sampletimes)), 
                paste0("LogFoldTRIM_CountAdj_",c(1:sampletimes)), paste0("PvalTRIM_CountAdj_",c(1:sampletimes)), 
                paste0("LogFold_CountAdj_",c(1:sampletimes)), paste0("Pval_CountAdj_",c(1:sampletimes)), 
                paste0("Parent_effortTRIM_CountAdj_",c(1:sampletimes)), paste0("Parent_contribTRIM_CountAdj_",c(1:sampletimes)), 
                paste0("Child_pnratioTRIM_CountAdj_",c(1:sampletimes)), paste0("Child_propTRIM_CountAdj_",c(1:sampletimes)), 
                paste0("Parent_effort_CountAdj_",c(1:sampletimes)), paste0("Parent_contrib_CountAdj_",c(1:sampletimes)))
# matrix_type = c("Child_entropyTRIM_CountAdj", "Child_entropyTRIM_Prop", "Parent_entropyTRIM_CountAdj", "Parent_entropyTRIM_Prop", 
#                 "LogFoldTRIM_CountAdj", "LogFoldTRIM_Prop", "PvalTRIM_CountAdj", "PvalTRIM_Prop",
#                 "CountAdj", "Parent_entropy", "Child_entropy",
#                 "LogFold_CountAdj", "LogFold_Prop", "Pval_CountAdj", "Pval_Prop",
#                 "Parent_effortTRIM_CountAdj", "Parent_effortTRIM_Prop", "Parent_contribTRIM_CountAdj", "Parent_contribTRIM_Prop", 
#                 "Child_pnratioTRIM_CountAdj", "Child_pnratioTRIM_Prop", "Child_propTRIM_CountAdj", "Child_propTRIM_Prop",
#                 "Parent_effort_CountAdj", "Parent_effort_Prop", "Parent_contrib_CountAdj", "Parent_contrib_Prop")
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

source("~/projects/IMPC/code/_funcAlice.R")
libr("stringr")
libr("colorspace")
libr("vegan") # libr(proxy)
libr("foreach")
libr("doMC")
libr("lubridate") #if there are date variables


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
                m = Reduce('cbind',m)
                
                save(d, file=dname0)
                write.csv(as.matrix(d), file=paste0(checkm(d,dname0),".csv"))
              }
            } 

          } else {
            if (dis[i]=="rbf") d = as.matrix(kernelMatrix(rbfdot(0.0838810518463767 ),m))
            simmatrix = kernelMatrix()
            sp = specc(m,centers=2,kernel="rbfdot")
            sp0 = specc(as.kernelMatrix(d),centers=2)
            rts = Rtsne(d, is_distance=T)$Y
            sm = sampleMeta[match(rownames(m),sampleMeta$fileName),]
            plot(rts,col=as.integer(factor(sm$aml)), pch=16, cex=1)
            plot(rts,col=sp@.Data, cex=2)
            plot(rts,col=sp0@.Data, cex=2)
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