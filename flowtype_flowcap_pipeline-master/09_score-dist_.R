## Input: original features --> Output: clusterings and distance matrices
# aya43@sfu.ca 20170419

#Directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="") #metadata for cell populations (phenotype)
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="") #metadata for FCM files
matrix_dir = paste(result_dir, "/matrix", sep="") #feature matrices
dist_dir = paste(result_dir, "/dist", sep="")

#Output
dist_score_dir = paste(dist_dir, "/dist_score", sep=""); for (i in 1:length(dist_score_dir)) { suppressWarnings(dir.create(dist_score_dir[i])) }
dist_clustercl_dir = paste(dist_score_dir, "/dist_cluster_cl", sep=""); for (i in 1:length(dist_clustercl_dir)) { suppressWarnings(dir.create(dist_clustercl_dir[i])) }

#Libraries/Functions
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
libr("stringr")
libr("foreach")
libr("doMC")
libr("FastKNN")
libr("cluster")
libr("mclust")
libr("kernlab")
libr("densitycut") #devtools::install_bitbucket("jerry00/densitycut_dev")

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)













#Options for script
overwritecl = T #redo and overwrite all clusters
overwritedist = T #redo and overwrite all distances

avgallcol = "specimen" #when averaging distances, match based on specimen

matrix_count = c("CountAdj")

dis = c("euclidean", "maximum", "manhattan", "canberra", "binary", "bray", "kulczynski", "minkowski", "morisita", "horn", "binomial", "rbf","poly","vanilla","tanh","laplace","bessel","anova","spline")

cltypesorigm = c("spec","dc") #uses original matrix

interested = c("tube","aml") #sampleMeta columns to plot
splitby = c("none","tube") #match with above



#Prepare data
phenoMeta = get(load(phenoMeta_dir))
sampleMeta0 = get(load(sampleMeta_dir))
mcount = get(load(paste0(matrix_dir, matrix_count,".Rdata")))

cellCountThres = c(1200)
k0=c(1,4,max(phenoMeta$phenolevel))

#list out feature matrices
matrix_type = list.files(path=result_dir,pattern=glob2rx("matrix*.Rdata"))
matrix_type = matrix_type[!matrix_type=="matrixCount.Rdata"]
matrix_type = gsub("matrix","",matrix_type)
matrix_type = gsub(".Rdata","",matrix_type)
matrix_type = matrix_type[!grepl("KO|Max",matrix_type)]
# matrix_type = matrix_type[!grepl("Freqp_orig",matrix_type)]
# nonrand = is.na(as.integer(sapply(matrix_type, function(x) substrRight(x,1))))
# nopropfreq = !grepl("Prop|Freqp",matrix_type,ignore.case=F)
# matrix_type = matrix_type[nonrand & nopropfreq]
# matrix_type = matrix_type[!(nonrand & nopropfreq)]











start = Sys.time()

#calculate rbf dist matrices & clustering for ones whose input is the original feature
start=Sys.time()

loop.ind = matrix_type
errorsb = foreach (mcp = loop.ind) %dopar% {
  tryCatch({
    if (!file.exists(paste0(matrix_dir, mcp,".Rdata"))) {cat("doesn't exist"); return(NULL)}
    
    #load feature matrix
    mresult = Loadintermatrices(paste0(matrix_dir, mcp,".Rdata"))
    mml0 = mresult$mml
    norm = "cellpop"
    if (is.null(dim(mml0[[1]])))norm = "none"
    mmlname = names(mml0)
    pt = mpi = mresult$pt
    gt = mresult$gt
    
    #get to-delete low count phenotype indices; CountAdj should be first one
    for (countThres in cellCountThres) { cat("\ncountThres: ",countThres," >",sep="")

      #get to-delete high no of marker phenotypes
      for (k in ifelse(grepl("TRIM",mcp), k0[!k0==1], k0)) { cat(" level",k," ",sep="")
        
        #trim feature matrix
        mmlresult = trimMatrices(mml0,mcount,pt,gt,phenoMeta,NULL,F, countThres,k)
        if (is.null(mmlresult)) next
        m0 = mmlresult$mml[[1]]
        if (norm=="none") m0 = Reduce('cbind',m0)
        m0 = as.matrix(m0)
        m0 = m0[!apply(m0[,-1], 1, function(x) all(x==0)),!apply(m0[-1,], 2, function(x) all(x==0))]
        pm = mmlresult$pm
        
        dnamee = paste(dist_dir, "/", mcp, "_rbf_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_normalize-",norm,".Rdata", sep="")
        dnameo = paste(dist_clustercl_dir, "/", mcp, "_other_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_normalize-",norm,".Rdata", sep="")
        
        sampleMeta = sampleMeta0[match(rownames(m0),sampleMeta0$fileName),]
        
        #create score list
        cm0 = list()
        if (file.exists(dnameo) & !overwritecl) cm0 = get(load(dnameo))
        
        #create distance/similarity matrix list
        drbf = list()
        simrbf = list()
        if (file.exists(dnamee) & !overwritedist) next
        
        #for each intersted column
        for (col in 1:length(interested)) { 
          colnam = paste0("interested-",interested[col],"_splitby-",splitby[col])
          
          #create lists
          if (!colnam%in%names(cm0)) cm0[[colnam]] = list()
          drbf[[colnam]] = list()
          simrbf[[colnam]] = list()
          
          #split distance and sampleMeta by a variable + do an average
          split=T
          if (splitby[col]=="none") split=F
          
          sm = splitmatrix(sampleMeta,sampleMeta[,colnames(sampleMeta)==splitby[col]],split=split,missym=F)
          m = splitmatrix(m0,sampleMeta[,colnames(sampleMeta)==splitby[col]],split=split,missym=F)
          
          #ensure every matrices have values and more than one class
          if (length(Reduce('intersect',lapply(sm, function(x) unique(x[,interested[col]]))))<length(unique(sampleMeta0[,interested[col]]))) next
          if (length(Reduce('intersect',lapply(sm, function(x) x[duplicated(x[,interested[col]]),interested[col]])))<length(unique(sampleMeta0[,interested[col]]))) next  #if not both aml and healthy in each tube
          
          #for each interested column variable
          for (dindname in names(m)) {
            #create score list
            if (!dindname%in%names(cm0[[colnam]])) cm0[[colnam]][[dindname]] = list()
            
            #list out class labels
            la0 = sm[[dindname]][,interested[col]]
            
            ## spectral clustering and rbf kernel parameter learning
            sp = NULL
            tr <- try({ sp = specc(x=m[[dindname]],kernel="rbfdot",centers=length(unique(la0))) })
            if ("try-error" %in% class(tr)) {
              mmm = (m[[dindname]]-min(m[[dindname]])) / (max(m[[dindname]])-min(m[[dindname]]))
              sp = specc(x=mmm,kernel="rbfdot",centers=length(unique(la0)))
            }
            
            ## calculate rbf distance
            parlist0 = unlist(sp@kernelf@kpar)
            simrbf[[colnam]][[dindname]] = kernelMatrix(kernel=rbfdot(as.numeric(parlist0)), x=m[[dindname]])
            drbf[[colnam]][[dindname]] = get_graphd(simrbf[[colnam]][[dindname]])
            
            ## get cluster results for spectral and densitycut
            for (cltype in cltypesorigm) {
              if (cltype%in%names(cm0[[colnam]][[dindname]]) & !overwritecl) next
              if (is.null(sp) & cltype=="spec") next
              cm0[[colnam]][[dindname]][[cltype]] = list()
              
              if (cltype=="spec") {
                clt = matrix(sp@.Data,ncol=1)
              } 
              if (cltype=="dc" & (!cltype%in%names(cm0[[colnam]][[dindname]]) | overwritecl)) {
                alpha=.85; nu=seq(0.0, 1.0, by=0.05)
                clt = DensityCut(X=m[[dindname]], alpha=alpha, nu=nu, show.plot=F)$cluster
                clt = matrix(clt,ncol=1)
              }
              rownames(clt) = rownames(m[[dindname]])
              colnames(clt) = "none"
              cm0[[colnam]][[dindname]][[cltype]]$clt = clt
            }
          }
        }
        save(cm0,file=dnameo)
        
        save(simrbf,file=gsub(".Rdata","_simmatrix.Rdata",dnamee))
        save(drbf,file=dnamee)
      }
    }
    return(T)
  }, error = function(e) {
    cat(paste("ERROR:  ",e)); graphics.off(); return(T)
  })
}

TimeOutput(start)









#check which features weren't finished
aa = distMetafun(list.files(dist_clustercl_dir,pattern="other"),dis=c(dis,"other"))
mattype = sapply(1:nrow(aa), function(x) {
  if (aa$rand[x]==0) { add = ""
  } else { add = paste0("_",aa$rand[x]) }
  return(paste0(aa$type[x],add))
})
mattype2 = sapply(unique(mattype), function(x) {
  al = aa$layer[which(mattype==x)]
  if (grepl("TRIM",x) & length(al)>1) return(T)
  if (!grepl("TRIM",x) & length(al)>2) return(T)
  return(F)
})
# unique(mattype)[mattype2]
matrix_type[!matrix_type%in%unique(mattype)[mattype2]] #not done



