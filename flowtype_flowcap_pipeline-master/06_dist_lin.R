## Input: distance matrices --> Output: distance matrix that is linear combo of input distance matrices
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
phenoChild_dir = paste(result_dir, "/phenoChild.Rdata",sep="")
phenoChild_ind_dir = paste(result_dir, "/phenoChild_ind.Rdata",sep="")
phenoChildpn_dir = paste(result_dir, "/phenoChildpn.Rdata",sep="")
phenoChildpn_ind_dir = paste(result_dir, "/phenoChildpn_ind.Rdata",sep="")
cellCountThres = c(200) #insignificant if count under

#Output
dist_dir = paste(result_dir, "/dist", sep=""); for (i in 1:length(dist_dir)) { suppressWarnings(dir.create (dist_dir[i])) }
dist_type_dir = NULL

#Libraries/Functions
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
libr("stringr")
libr("colorspace")
libr("vegan") # libr(proxy)
libr("foreach")
libr("doMC")
libr("lubridate") #if there are date variables
libr("pracma") #for dot product


#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)












#Options for script
overwrite = T

ignoredist = ".csv"
matrix_count = c("CountAdj")
matrix_type0 = list(c(273,471,460,505)) #list according to distmfilenames
matrix_weights0 = list(c(.25,.25,.25,.25))

countThres = c(1200) #just one is enough



#Prepare data
distmfile = list.files(dist_dir, recursive=F, full.names=T, pattern=".Rdata")
delpvalind = grepl(ignoredist,distmfile,ignore.case=T)
if (sum(delpvalind)>0) distmfile = distmfile[!delpvalind]
distmfilenames = fileNames(distmfile)

sampleMeta0 = get(load(sampleMeta_dir))
phenoMeta = get(load(phenoMeta_dir))
k0 = c(1,4,max(phenoMeta$phenolevel)) # how many markers to consider i.e. k=max(phenolevel) only

m0 = get(load(paste0(matrix_dir, matrix_count,".Rdata")))










start = Sys.time()

for (mi in 1:length(matrix_type0)) {
  leavePhenotype = list()
  doneAll = F
  
  start2 = Sys.time()
  matrix_type. = distmfile[matrix_type0[[mi]]]
  matrix_weights = matrix_weights0[[mi]]
  cat("\n"); cat(paste(matrix_weights,matrix_type., sep=" x "),sep=",  ")
  
  #load different matrices
  mresult = Loadintermatrices(matrix_type.)
  mml = mresult$mml
  mmlname = names(mml)
  pt = mresult$pt
  gt = mresult$gt
  
  #get to-delete low count phenotype indices; CountAdj should be first one
  for (countThres in cellCountThres) {cat("\ncountThres: ",countThres," >",sep="")
    
    #get to-delete high no of marker phenotypes
    for (k in k0) { cat(" level",k," ",sep="")
      
      #trim feature matrix
      mmlresult = trimMatrices(mml,m0,pt,gt,phenoMeta,leavePhenotype,doneAll, countThres,k)
      if (is.null(mmlresult)) next
      mml = mmlresult$mml
      pm = mmlresult$pm
      leavePhenotype = mmlresult$leavePhenotype
      doneAll = mmlresult$doneAll
      
      # Calculate distance (normalize by cell population) --------------------------------------------------------
      
      dname = paste(dist_dir, "/linear_", paste0(matrix_type.,"-",matrix_weights, collapse="_"), "_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_normalize-cellpop", sep = "" )
      if (grepl(dname,distmfilenames) & overwrite) { cat(" exists, skipped."); next }
      
      #weigh each feature matrix; scalendev = 0 (uhh...), T (scaling is done by dividing the centred columns by their SD), F (none)
      for (s in 1:3) {
        if (s==1) { scalendev=0; dn = mml
        } else {
          if (s==2) scalendev=T
          if (s==3) scalendev=F
          dn = foreach(i=1:length(mml)) %dopar% { 
            a = scale(as.vector(mml[[i]]), scale=scalendev)
            a = a-min(a)
            return(matrix(a,nrow=nrow(mml[[i]])))
          }
        }
        d0 = foreach(i=1:length(dn)) %dopar% { return(matrix_weights[i]*dn[[i]]) }
        d1 = Reduce("+",d0)
        colnames(d1) = rownames(d1) = gt
        
        save(d1, file=paste0(dname,"_scaleSD-",scalendev,".Rdata"))
        write.csv(d1, file=paste0(dname,"_scaleSD-",scalendev,".csv"))
      }
    }
  }
  TimeOutput(start2)
}
TimeOutput(start)



