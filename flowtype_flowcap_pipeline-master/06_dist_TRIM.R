# Dist TRIM: for trimmed only
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
matrix_count = c("CountAdj")
matrix_type = c("Child_entropyTRIM_CountAdj", "Child_entropyTRIM_Prop", "Parent_entropyTRIM_CountAdj", "Parent_entropyTRIM_Prop", 
                "LogFoldTRIM_CountAdj", "LogFoldTRIM_Prop", "PvalTRIM_CountAdj", "PvalTRIM_Prop",
                #"CountAdj", "Parent_entropy", "Child_entropy",
                #"LogFold_CountAdj", "LogFold_Prop", "Pval_CountAdj", "Pval_Prop",
                "Parent_effortTRIM_CountAdj", "Parent_effortTRIM_Prop", "Parent_contribTRIM_CountAdj", "Parent_contribTRIM_Prop", 
                "Child_pnratioTRIM_CountAdj", "Child_pnratioTRIM_Prop", "Child_propTRIM_CountAdj", "Child_propTRIM_Prop")
#matrix_type = c("PvalTRIM_CountAdj", "LogFoldTRIM_CountAdj", "Parent_effort", "Parent_contrib")
weight_matrix = c("MaxCountAdj_CountAdj")

phenoChild_dir = paste(result_dir, "/phenoChild.Rdata",sep="")
phenoChild_ind_dir = paste(result_dir, "/phenoChild_ind.Rdata",sep="")
phenoChildpn_dir = paste(result_dir, "/phenoChildpn.Rdata",sep="")
phenoChildpn_ind_dir = paste(result_dir, "/phenoChildpn_ind.Rdata",sep="")
cellCountThres = c(1200) #insignificant if count under



#Output
dist_dir = paste(result_dir, "/dist", sep=""); for (i in 1:length(dist_dir)) { suppressWarnings(dir.create (dist_dir[i])) }
dist_type_dir = NULL
# for(i in 1:length(mcp_types)) { for(j in 1:length(dist_dir)) { dist_type_dir = append(dist_type_dir, paste0(dist_dir[j], mcp_types[i]))} }
# for (i in 1:length(dist_type_dir)) { suppressWarnings(dir.create (dist_type_dir[i])) }

source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
libr("stringr")
libr("colorspace")
libr("vegan") # libr(proxy)
libr("foreach")
libr("doMC")
libr("lubridate") #if there are date variables
libr("pracma") #for dot product


dis = c("euclidean", "manhattan", "canberra", "mahalanobis") #
#dis = c("euclidean", "canberra", "binomial", "cao") # "jaccard","mahalanobis", "gower","altGower","morisita","horn", "manhattan", "maximum", "binary", "minkowski", "raup", "mountford" (all at same level),)
disnoavg = c("euclidean", "manhattan", "canberra", "bray", "kulczynski", "morisita", "horn", "binomial")  #dis measures that don't average over number of features
disindist = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski") #use dist function when possible, instead of vegdist


no_cores = 5#detectCores()-1
registerDoMC(no_cores)



start = Sys.time()


#for (ci in 4:1) {
phenoMeta = get(load(phenoMeta_dir))
sampleMeta = get(load(sampleMeta_dir))
#k0 = c(1,3,max(phenoMeta$phenolevel)) # how many markers to consider i.e. k=max(phenolevel) only
k0 = c(1,max(phenoMeta$phenolevel)) # how many markers to consider i.e. k=max(phenolevel) only

# delPhenoInd = list()
# doneAll = F #will mark as TRUE once we are done calculating d of the full mm, don't repeat!
# delPhenoInd = list()
# delPhenoInd_child = list()
# delPhenoInd_childpn = list()
# doneAll = F #will mark as TRUE once we are done calculating d of the full mm, don't repeat!
# doneAll_child = F
# doneAll_child_pn = F

start2 = Sys.time()
leavePhenotype = list()
doneAll = F

mw0 = get(load(paste0(matrix_dir, weight_matrix,".Rdata")))
m0 = get(load(paste0(matrix_dir, matrix_count,".Rdata")))

for (mcp in matrix_type) {
  cat("\n", mcp, " ",sep="")
  start2 = Sys.time()
  leavePhenotype = list()
  doneAll = F
  if (!file.exists(paste0(matrix_dir, mcp,".Rdata"))) {cat("doesn't exist"); next}
  mm = get(load(paste0(matrix_dir, mcp,".Rdata")))
  
  if (!is.null(dim(mm))) {
    mpi = colnames(mm)
    # print(mm[50:60,20:30])
  } else {
    mpi = names(mm)
    # print(head(mm[[25]]))
  }
  
  #get to-delete low count phenotype indices; CountAdj should be first one
  for (countThres in cellCountThres) {
    cat("\ncountThres: ",countThres," > ",sep="")
    lowCountpt = colnames(m0)[apply(m0,2,function(x) all(x<=countThres))]
    if (!length(lowCountpt)>0) lowCountpt = c()
    
    #get to-delete high no of marker phenotypes
    for (k in k0) {
      cat("level",k," ",sep="")
      highLevelpt = phenoMeta[phenoMeta$phenolevel>k,"phenotype"]
      if (!length(highLevelpt)>0) highLevelpt = c()
      
      
      # cat("length of countthres: ", length(lowCountInd))
      # cat("\nlength of highlevelind: ", length(highLevelInd))
      ## Load & fix cell count/countAdj/proportion matrix -----------------------------------------------------
      dpi0 = union(lowCountpt,highLevelpt)
      # cat("\nlength of dpi: ", length(dpi0), "\n")
      
      m = mm
      
      lpi = setdiff(mpi,dpi0) #phenotypes to leave in analysis
      lpisplit = lapply(lpi, function(x) return(unlist(strsplit(x, "(?<=[+-])", perl=TRUE))))
      
      #check if indices already calculated for on this matrix, if yes, skip, else trim and do
      if (Position(function(x) identical(x, lpi), leavePhenotype, nomatch = 0) > 0) { cat("; child skipped", sep=""); next }
      leavePhenotype[[paste0(k="k.",k,"_countThres.",countThres)]] = lpi
      
      if (!doneAll & length(lpi)==length(mpi)) doneAll = T
      
      if (!is.null(dim(mm))) { m = m[,match(lpi,mpi)]
      } else { m = m[match(lpi,mpi)] }
      
      pm = phenoMeta[match(lpi,phenoMeta$phenotype),]
      
      
      #seperate out the weight matrix
      #seperate out the weight matrix
      mw = mw0[match(rownames(m),rownames(mw0)),match(lpi,colnames(mw0))]
      mw[mw<1]=1; mwl = log(mw,getlv(max(mw),100))/100
      mwl = mwl[,match(colnames(m),colnames(mwl))] #no normalization, arrange like m
      
      
      
      # Calculate distance (normalize by cell population) --------------------------------------------------------
      
      if (doneAll) dname = paste(dist_dir, "/overlap_", paste0(matrix_type,"-",weight_matrix, collapse="_"), "_FULL_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, sep="")
      dname = paste(dist_dir, "/overlap_", paste0(matrix_type,"-",weight_matrix, collapse="_"), "_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, sep = "" )
      if (file.exists(dname)) { cat(" exists, skipped."); next }
      d = matrix(0,nrow=nrow(m),ncol=nrow(m))
      dsim = matrix(Inf,nrow=nrow(m),ncol=nrow(m))
      ddif = matrix(0,nrow=nrow(m),ncol=nrow(m)) #if no overlap, contains weights to nonoverlapping
      
      for (i in 2:nrow(m)) {
        dsim[1:(i-1),i] = dsim[i,1:(i-1)] = foreach(j = 1:(i-1), .combine="c") %dopar% {
          col0 = apply(m[c(i,j),], 2, function(x) all(x>0))==T
          if (sum(col0)>0) {
            x = sum( pmax(mwl[i,col0],mwl[j,col0]) * abs(m[i,col0]-m[j,col0]) ) #smallDifference * largeImportance = small is good
          } else { x = Inf }
          return(x)
        }
        ddif[1:(i-1),i] = ddif[i,1:(i-1)] = foreach(j = 1:(i-1), .combine="c") %dopar% {
          colij = apply(m[c(i,j),], 2, function(x) sum(x>0)==1)==T ## use frequent pattern mining... compare frequent itemsets
          if (sum(colij)>0) {
            y = sum( pmax(mwl[i,colij],mwl[j,colij]) * pmax(m[i,colij], m[j,colij]) ) /2 #value * largeImportance = small is good
          } else { y = 0 }
          return(y)
        }
      }
      d = dsim/ddif
      for (i in 1:nrow(d)) { d[i,i] = 0 }
      d[ddif==0] = 0 #super close
      d[(dsim==0 & ddif>0) | dsim==Inf] = max(d[d!=Inf])*2
      d[dsim==Inf] = 
      dsimdiff = list(dsim=dsim,ddif=ddif)
      
      
      # if (absweight) save(d,file=paste0(dname,"_absweight_original.Rdata"))
      # if (logweight) save(dl,file=paste0(dname,"_logweight_original.Rdata"))
      save(d,file=paste0(dname,".Rdata"))
      save(dsimdiff,file=paste0(dname,"/simdiff.Rdata"))
      
      #weigh each feature matrix; scalendev = 0 (uhh...), T (scaling is done by dividing the centred columns by their SD), F (none)
      for (s in 1:3) {
        # for (l in 1:2) {
        # if (l==1 & logweight) dfinal = dl
        # if (l==2 & absweight) dfinal = d
        
        # if (s==1) { scalendev=0; dn = dfinal
        if (s==1) { scalendev=0; dn = d
        } else {
          if (s==2) scalendev=T
          if (s==3) scalendev=F
          dn = foreach(i=1:length(d)) %dopar% { 
            a = scale(as.vector(d[[i]]), scale=scalendev)
            a = a-min(a)
            return(matrix(a,nrow=nrow(d[[i]])))
          }
        }
        d0 = foreach(i=1:length(dn)) %dopar% { return(matrix_weights[i]*dn[[i]]) }
        d1 = Reduce("+",d0)
        colnames(d1) = rownames(d1) = gt
        
        # if (l==1 & logweight) {
        #   save(d1, file=paste0(dname,"_logweight_scaleSD-",scalendev,".Rdata"))
        #   write.csv(d1, file=paste0(dname,"_logweight_scaleSD-",scalendev,".csv"))
        # }
        # if (l==2 & absweight) {
        #   save(d1, file=paste0(dname,"_absweight_scaleSD-",scalendev,".Rdata"))
        #   write.csv(d1, file=paste0(dname,"_absweight_scaleSD-",scalendev,".csv"))
        # }
        save(d1, file=paste0(dname,"_scaleSD-",scalendev,".Rdata"))
        write.csv(d1, file=paste0(dname,"_scaleSD-",scalendev,".csv"))
      }
    }
    
    TimeOutput(start2)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  }
  TimeOutput(start2)
}

}




TimeOutput(start)



#require(RDRToolbox)
#libr(rgl)
## Isomap ##



# libr(dbscan)
# oc = optics(d, eps=10, minPts=3) # no results (gravitate towards each other)
