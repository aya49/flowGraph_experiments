## Input: original/trimmed features --> Output: trimmed features + rchyoptimyx additional nodes
# aya43@sfu.ca 20161220

#Directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/matrix", sep="")

#Output
dist_dir = paste(result_dir, "/dist", sep=""); for (i in 1:length(dist_dir)) { suppressWarnings(dir.create (dist_dir[i])) }
rchy_dir = paste(result_dir, "/rchy", sep="")


#Libraries/Functions
source("~/projects/IMPC/code/_funcAlice.R")
libr("stringr")
libr("colorspace")
libr("vegan") # libr(proxy)
libr("foreach")
libr("doMC")
libr("lubridate") #if there are date variables
libr("kernlab")

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)








#Options for script
sampletimes = 5
overwrite = T #overwrite distances?
writecsv = T

cellCountThres = c(1200) #insignificant if count under

kpaths = 5

matrix_count = c("CountAdj")



#Prepare data
weight_matrix = c("MaxCountAdj_CountAdj") #assume by cell pop

matrix_type = c("PvalTRIM_CountAdj")
matrix_all_type = c("Pval_CountAdj")
matrix_edge_type = c("Parent_contrib_CountAdj") #plot only

# matrix_type = list.files(path=result_dir,pattern=glob2rx("matrix*.Rdata"))
# matrix_type = matrix_type[!matrix_type=="matrixCount.Rdata"]
# matrix_type = gsub("matrix|.Rdata","",matrix_type)
# matrix_type = matrix_type[!grepl("KO|Max",matrix_type)]
# matrix_type = matrix_type[!grepl("Freqp_orig",matrix_type)]

# matrix_type = c("CountAdj", "Parent_entropy", "Child_entropy", 
#                 paste0("Child_entropyTRIM_CountAdj_",c(1:sampletimes)), "Child_entropyTRIM_CountAdj", "Child_entropyTRIM_Prop", 
#                 paste0("Parent_entropyTRIM_CountAdj_",c(1:sampletimes)), "Parent_entropyTRIM_CountAdj", "Parent_entropyTRIM_Prop", 
#                 paste0("LogFoldTRIM_CountAdj_",c(1:sampletimes)), "LogFoldTRIM_CountAdj", "LogFoldTRIM_Prop", 
#                 paste0("PvalTRIM_CountAdj_",c(1:sampletimes)), "PvalTRIM_CountAdj", "PvalTRIM_Prop",
#                 paste0("LogFold_CountAdj_",c(1:sampletimes)), "LogFold_CountAdj", "LogFold_Prop", 
#                 paste0("Pval_CountAdj_",c(1:sampletimes)), "Pval_CountAdj", "Pval_Prop",
#                 "Child_pnratio", "Child_prop",
#                 paste0("Parent_effortTRIM_CountAdj_",c(1:sampletimes)), "Parent_effortTRIM_CountAdj", "Parent_effortTRIM_Prop", 
#                 paste0("Parent_contribTRIM_CountAdj_",c(1:sampletimes)), "Parent_contribTRIM_CountAdj", "Parent_contribTRIM_Prop", 
#                 paste0("Child_pnratioTRIM_CountAdj_",c(1:sampletimes)), "Child_pnratioTRIM_CountAdj", "Child_pnratioTRIM_Prop", 
#                 paste0("Child_propTRIM_CountAdj_",c(1:sampletimes)), "Child_propTRIM_CountAdj", "Child_propTRIM_Prop",
#                 paste0("Parent_effort_CountAdj_",c(1:sampletimes)), "Parent_effort_CountAdj", "Parent_effort_Prop", 
#                 paste0("Parent_contrib_CountAdj_",c(1:sampletimes)), "Parent_contrib_CountAdj", "Parent_contrib_Prop")

m0 = get(load(paste0(matrix_dir, matrix_count,".Rdata")))
phenoMeta = get(load(phenoMeta_dir))

k0 = c(1,4,max(phenoMeta$phenolevel)) # how many markers to consider i.e. k=max(phenolevel) only











start = Sys.time()

#load different features matrix and calculate distance matrix
#a = foreach(mt=1:length(matrix_type)) %dopar% {
for (mt in 1:length(matrix_type)) {
  tryCatch({
    mcp = matrix_type[mt]
    map = matrix_all_type[mt]
    mep = matrix_edge_type[mt]
    cat("\n", mcp, " ",sep="")
    start2 = Sys.time()
    
    #start a list of phenotypes included in each distance matrix calculation such that none are repeated
    leavePhenotype = list()
    doneAll = F
    
    #load feature matrix
    mresult = Loadintermatrices(c(paste0(matrix_dir, mcp,".Rdata"),paste0(matrix_dir, map,".Rdata"),paste0(matrix_dir, mcp,".Rdata")))
    mml0 = mresult$mml
    mmlname = names(mml0)
    pt = mpi = mresult$pt
    gt = mresult$gt
    
    #get to-delete low count phenotype indices; CountAdj should be first one
    for (countThres in cellCountThres) { cat("\ncountThres: ",countThres," >",sep="")
      
      #get to-delete high no of marker phenotypes
      for (k in k0) { cat(" level",k," ",sep="")
        
        #trim feature matrix
        mmlresult = trimMatrices(mml0,m0,pt,gt,phenoMeta,leavePhenotype,doneAll, countThres,k)
        if (is.null(mmlresult)) return(NULL)
        m = mmlresult$mml[[1]]
        ma = mmlresult$mml[[2]]
        e = mmlresult$mml[[3]]
        pm = mmlresult$pm
        leavePhenotype = mmlresult$leavePhenotype
        doneAll = mmlresult$doneAll
        
        if (is.null(m)) return(NULL)
        if (is.null(dim(m))) {
          a = Reduce('cbind',m); if (all(a==0)) next
        } else { if (all(m==0)) next }
        
        fname = paste0(rchy_dir,"_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres); suppressWarnings(dir.create(fname))
        
        foreach(i=rownames(m)) %dopar% {
          rchy = NULL
          phenocodes = phenoMeta$phenocode[match(colnames(ma),phenoMeta$phenotype)]
          
          # positive change pvalues
          if (sum(m[i,]>0)>0) {
            phenoscore = ma[i,]+min(ma[i,])
            startpheno = phenoMeta$phenocode[match(colnames(m)[m[i,]>0],phenoMeta$phenotype)]
            rchy0 = lapply(startpheno, function(sp) RchyOptimyx(pheno.codes=phenocodes, phenotypeScores=phenoscore, startPhenotype=sp, pathCount=kpaths, trimPaths=F))
            rchy = rchy0[[1]]; if (length(rchy0)>1) for (ri in 2:length(rchy0)) { rchy = merge(rchy,rchy0[[ri]]) }
          }
          
          # negative change pvalues
          if (sum(m[i,]<0)>0) {
            phenoscore = (-ma[i,])+min(-ma[i,])
            startpheno = phenoMeta$phenocode[match(colnames(m)[m[i,]<0],phenoMeta$phenotype)]
            rchy1 = lapply(startpheno, function(sp) RchyOptimyx(pheno.codes=phenocodes, phenotypeScores=phenoscore, startPhenotype=sp, pathCount=kpaths, trimPaths=F))
            if (!sum(m[i,]>0)>0) rchy = rchy1[[1]]
            rchy = merge(rchy,rchy1[[1]]); if (length(rchy1)>1) for (ri in 2:length(rchy1)) { rchy = merge(rchy,rchy1[[ri]]) }
          }
          
          if (!is.null(rchy)) save(rchy,file=paste0(fname,"/",i))
        }
        
      } #layer
    } #countThres
    
    TimeOutput(start2)
  }, error = function(err) { cat(paste("ERROR:  ",err)) })
}

TimeOutput(start)












#delete distance matrices with all 0
distmfile = list.files(dist_dir, recursive=F, full.names=T, pattern=".Rdata$")
a = foreach (i = 1:length(distmfile),.combine="c") %dopar% {
  a = F
  d = get(load(distmfile[i]))
  if (any(is.na(as.matrix(d)))) { rm(d); return(T) }
  if (all(as.matrix(d)==as.matrix(d)[1,1])) { rm(d); return(T) }
  rm(d)
  return(a)
}
file.remove(distmfile[a])




