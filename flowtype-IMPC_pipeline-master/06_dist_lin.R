# aya43@sfu.ca 20161220
# Uses different distance measures to calculate distance & plot samples

#root = "~/projects/IMPC"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#,"Sanger_MLN","CIPHE","TCP","H")

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrix", sep="")
# mcp_types = c("/original", "/child_pn","/child_prop", "/trim") 
# matrix_type1 = c("CountAdj") #, "Prop", "Count",  #countAdj comes first, to set which columns are deleted based on cellcountthres
# matrix_type2 = c("Child_prop")
# matrix_type3 = c("Child_pnratio", "Child_entropy")
# matrix_type4 = c("Child_pnratio", "Child_entropy")

# node only
# matrix_type = c("Child_entropy", "Child_pnratio", "Child_prop")
# matrix_weights = c(.3,.3,.)
# weight_matrix = c("CountAdj")

# with pvalued features
# matrix_count = c("CountAdj")
# matrix_type = c("PvalTRIM_CountAdj", "LogFoldTRIM_CountAdj", "Child_entropy", "Child_pnratio", "Child_prop")#, "Parent_effortTRIM", "Parent_contribTRIM") #--blank :( fix!
# matrix_weights = c(.3,.3,.15,.15,.1)
# weight_matrix = c("MaxCountAdj_CountAdj")
# logweight = T #log weight first?

matrix_count = c("CountAdj")
matrix_type = c("PvalTRIM_CountAdj", "LogFoldTRIM_CountAdj", "Parent_effort", "Parent_contrib")
matrix_weights = c(.25,.25,.25,.25)
weight_matrix = c("MaxCountAdj_CountAdj")
logweight = T #log weight first?

phenoChild_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChild.Rdata",sep="")
phenoChild_ind_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChild_ind.Rdata",sep="")
phenoChildpn_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChildpn.Rdata",sep="")
phenoChildpn_ind_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChildpn_ind.Rdata",sep="")
cellCountThres = c(200,200,500,500,500) #insignificant if count under



#Output
dist_dir = paste(result_dir, "/", panelL, "/", centreL, "/dist", sep=""); for (i in 1:length(dist_dir)) { suppressWarnings(dir.create (dist_dir[i])) }
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


dis = c("euclidean", "manhattan", "canberra", "mahalanobis") #
#dis = c("euclidean", "canberra", "binomial", "cao") # "jaccard","mahalanobis", "gower","altGower","morisita","horn", "manhattan", "maximum", "binary", "minkowski", "raup", "mountford" (all at same level),)
disnoavg = c("euclidean", "manhattan", "canberra", "bray", "kulczynski", "morisita", "horn", "binomial")  #dis measures that don't average over number of features
disindist = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski") #use dist function when possible, instead of vegdist


no_cores = 7#detectCores()-1
registerDoMC(no_cores)



start = Sys.time()


#for (ci in 4:1) {
for (ci in 1:length(paste0(panelL,centreL))) {
  
  start1 = Sys.time()
  
  centre = paste0(panelL," ",centreL)[ci]
  cat("\n",centre,sep="")
  
  phenoMeta = get(load(phenoMeta_dir[ci]))
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
  
  m0 = get(load(paste0(matrix_dir[ci], matrix_count,".Rdata")))
  
  #load different matrices
  mind = 1
  mml = list()
  mmlname = c()
  cat(" loading feature matrices ",sep="")
  for (mcp in c(weight_matrix,matrix_type)) {
    
    #pt = overlapping phenotypes & samples on all matrices
    if (!file.exists(paste0(matrix_dir[ci], mcp,".Rdata"))) {cat("doesn't exist"); next}
    mml[[mind]] = get(load(paste0(matrix_dir[ci], mcp,".Rdata")))
    mmlname[mind] = mcp
    if (mind==1) {
      if (is.null(dim(mml[[1]]))) {
        pt = names(mml[[1]])
        gt = rownames(mml[[1]][[1]])
      } else {
        pt = colnames(mml[[1]])
        gt = rownames(mml[[1]])
      }
    } else {
      if (is.null(dim(mml[[mind]]))) {
        pt0 = names(mml[[mind]])
        gt0 = rownames(mml[[mind]][[1]])
      } else {
        pt0 = colnames(mml[[mind]])
        gt0 = rownames(mml[[mind]])
      }
      pt = intersect(pt, pt0)
      gt = intersect(gt, gt0)
    }
    mind = mind+1
  }
  
  #trim matrices (sample)
  mml = foreach (i = 1:length(mml)) %dopar% {
    if (is.null(dim(mml[[i]]))) {
      gt0 = rownames(mml[[i]][[i]])
      ind = match(gt,gt0)
      mmlnames_temp = names(mml[[i]])
      a = lapply(1:length(mml[[i]]), function(j) {
        if (is.null(dim(mml[[i]][[j]]))) { return(as.matrix(mml[[i]][[j]],ncol=1)[ind,])
        } else { return( mml[[i]][[j]][ind,]) }
      })
      names(a) = mmlnames_temp
    } else {
      gt0 = rownames(mml[[i]])
      a = mml[[i]][match(gt,gt0),]
    }
    return(a)
  }
  names(mml) = mmlname
  mml0 = mml
  
  #get to-delete low count phenotype indices; CountAdj should be first one
  for (countThres in cellCountThres[ci]) {
    cat("\ncountThres: ",countThres," > ",sep="")
    lowCountpt = colnames(m0)[apply(m0,2,function(x) all(x<=countThres))]
    if (!length(lowCountpt)>0) lowCountpt = c()
    
    #get to-delete high no of marker phenotypes
    for (k in k0) {
      cat("level",k," ",sep="")
      highLevelpt = colnames(m0)[phenoMeta$phenolevel>k]
      if (!length(highLevelpt)>0) highLevelpt = c()
      
      
      # cat("length of countthres: ", length(lowCountInd))
      # cat("\nlength of highlevelind: ", length(highLevelInd))
      ## Load & fix cell count/countAdj/proportion matrix -----------------------------------------------------
      dpi0 = union(lowCountpt,highLevelpt)
      # cat("\nlength of dpi: ", length(dpi0), "\n")
      
      
      lpi = setdiff(pt,dpi0)
      
      #check if indices already calculated for on this matrix
      if (Position(function(x) identical(x, lpi), leavePhenotype, nomatch = 0) > 0) {cat("-skipped ", sep=""); next}
      leavePhenotype[[paste0(k="k.",k,"_countThres.",countThres)]] = lpi
      
      if (!doneAll & length(lpi)==length(pt)) doneAll = T
      
      
      #trim matrices (phenotype)
      mml = foreach (i = 1:length(mml0)) %dopar% {
        if (is.null(dim(mml0[[i]]))) {
          a = mml0[[i]][match(lpi,names(mml0[[i]]))]
        } else {
          a = mml0[[i]][,match(lpi,colnames(mml0[[i]]))]
        }
        return(a)
      }
      names(mml) = mmlname
      
      pm = phenoMeta[match(lpi,phenoMeta$phenotype),]
      
      
      #seperate out the weight matrix
      mw = matrix(1,ncol=length(lpi),nrow=length(gt),dimnames=(list(gt,lpi)))
      if (!is.null(weight_matrix)) { mw = mml[[1]]; mml = mml[seq(1,length(mml))[-1]] }
      if(logweight) mw = log(mw)
      
      
      
      
      # Calculate distance (normalize by cell population) --------------------------------------------------------
      
      if (doneAll) dname = paste(dist_dir[ci], "/linear_", paste0(matrix_type,"-",matrix_weights, collapse="_"), "_FULL_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_normalize-cellpop.Rdata", sep="")
      dname = paste(dist_dir[ci], "/linear_", paste0(matrix_type,"-",matrix_weights, collapse="_"), "_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_normalize-cellpop.Rdata", sep = "" )
      d = list()
      if (file.exists(dname)) { cat(" exists, skipped."); next }
      
      
      #calculate distance for every matrix/list
      
      for (matrix_ind in length(mml):1) {
        start3 = Sys.time()
        cat(" ", matrix_ind," ",sep="")
        dm = matrix(0,nrow=length(gt),ncol=length(gt))
        m = mml[[matrix_ind]]
        
        if (is.null(dim(m))) {
          for (i in 2:length(gt)) {
            dm[1:(i-1),i] = dm[i,1:(i-1)] = foreach(j = 1:(i-1), .combine="c") %dopar% {
              
              #distance edge metric (mean manhattan)
              x = sapply(1:length(m), function(h) {
                if (is.null(dim(m[[h]]))) {
                  a = matrix(m[[h]],ncol=1)
                } else {
                  a = m[[h]]
                }
                jx = a[j,]
                ix = a[i,]
                dist_mx = mean(abs(jx-ix))
                return(dist_mx)
              })
              
              #distance metric (sum)
              a = pmax(mw[i,],mw[j,]) * x
              a = sqrt(a[a!=0])
              dist_m = sum(a)
              return(dist_m)
            }
          }
        } else {
          for (i in 2:length(gt)) {
            dm[1:(i-1),i] = dm[i,1:(i-1)] = foreach(j = 1:(i-1), .combine="c") %dopar% {
              
              #distance metric (maxweight * manhattan)
              a = pmax(mw[i,],mw[j,]) * (abs(m[j,]-m[i,]))
              a = sqrt(a[a!=0])
              dist_m = sum(a)
              return(dist_m)
              
            }
          }
          
        }
        d[[matrix_ind]] = dm
        TimeOutput(start3)
      }
      #weigh each feature matrix
      d0 = foreach(i=1:length(d)) %dopar% { return(matrix_weights[i]*d[[i]]) }
      d1 = Reduce("+",d0)
      colnames(d1) = rownames(d1) = gt
      
      save(d1, file=dname)
      
      TimeOutput(start2)
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
    }
  }
  cat("\n centre ", centre, " ",TimeOutput(start1)," \n",sep="")
  
}



TimeOutput(start)



#require(RDRToolbox)
#libr(rgl)
## Isomap ##



# libr(dbscan)
# oc = optics(d, eps=10, minPts=3) # no results (gravitate towards each other)
