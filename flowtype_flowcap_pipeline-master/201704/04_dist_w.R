# Dist weighted: Manhattan distances
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
matrix_type = c("Child_entropyTRIM_CountAdj", "Child_entropyTRIM_Prop", "Parent_entropyTRIM_CountAdj", "Parent_entropyTRIM_Prop", 
                "LogFoldTRIM_CountAdj", "LogFoldTRIM_Prop", "PvalTRIM_CountAdj", "PvalTRIM_Prop",
                "LogFold_CountAdj", "LogFold_Prop", "Pval_CountAdj", "Pval_Prop",
                "Child_entropy", "Parent_entropy",
                "Parent_effortTRIM_CountAdj", "Parent_effortTRIM_Prop", "Parent_contribTRIM_CountAdj", "Parent_contribTRIM_Prop", 
                "Child_pnratioTRIM_CountAdj", "Child_pnratioTRIM_Prop", "Child_propTRIM_CountAdj", "Child_propTRIM_Prop",
                "Parent_effort_CountAdj", "Parent_effort_Prop", "Parent_contrib_CountAdj", "Parent_contrib_Prop", "Child_pnratio", "Child_prop")

weight_matrix = c("MaxCountAdj_CountAdj") #assume by cell pop
# logweight = T #log weight first?
# absweight = T #pick at least one of logweight or absweight, else just use dist ==> manhatta, so much faster

phenoChild_dir = paste(result_dir, "/phenoChild.Rdata",sep="")
phenoChild_ind_dir = paste(result_dir, "/phenoChild_ind.Rdata",sep="")
phenoChildpn_dir = paste(result_dir, "/phenoChildpn.Rdata",sep="")
phenoChildpn_ind_dir = paste(result_dir, "/phenoChildpn_ind.Rdata",sep="")
cellCountThres = c(200) #insignificant if count under

overwrite = T

ignoredist = ".csv"


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
libr(pracma) #for dot product
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")




no_cores = 7#detectCores()-1
registerDoMC(no_cores)



start = Sys.time()


phenoMeta = get(load(phenoMeta_dir))
#k0 = c(1,3,max(phenoMeta$phenolevel)) # how many markers to consider i.e. k=max(phenolevel) only
k0 = c(1,max(phenoMeta$phenolevel)) # how many markers to consider i.e. k=max(phenolevel) only
m0 = get(load(paste0(matrix_dir, matrix_count,".Rdata")))

#load different matrices
# for (mcp in matrix_type_) {
for (mcp in matrix_type) {
  cat("\n", mcp, " ",sep="")
  start2 = Sys.time()
  leavePhenotype = list()
  doneAll = F
  if (!file.exists(paste0(matrix_dir, mcp,".Rdata"))) {cat("doesn't exist"); next}
  mm = get(load(paste0(matrix_dir, mcp,".Rdata")))
  
  start2 = Sys.time()
  matrix_type. = append(paste0(matrix_dir,weight_matrix,".Rdata"), paste0(matrix_dir,mcp,".Rdata"))
  
  #load different matrices
  mresult = Loadintermatrices(matrix_type.)
  mml0 = mresult$mml
  mmlname = names(mml0)
  pt = mresult$pt
  gt = mresult$gt
  
  #get to-delete low count phenotype indices; CountAdj should be first one
  for (countThres in cellCountThres) {cat("\ncountThres: ",countThres," >",sep="")
    
    #get to-delete high no of marker phenotypes
    for (k in k0) { cat(" level",k," ",sep="")
      
      mmlresult = trimMatrices(mml0,m0,pt,gt,phenoMeta,leavePhenotype,doneAll, countThres,k)
      if (is.null(mmlresult)) next
      mml = mmlresult$mml
      pm = mmlresult$pm
      leavePhenotype = mmlresult$leavePhenotype
      doneAll = mmlresult$doneAll
      
      #seperate out the weight matrix
      m = NULL
      mw = matrix(1,ncol=nrow(pm),nrow=length(gt),dimnames=(list(gt,pm$phenotype)))
      if (!is.null(weight_matrix)) {
        mw = mml[[1]]; if (length(mml)>1) m = mml[[2]] 
      } else { if (length(mml)>0) m = mml[[1]] }
      #if(logweight) { mw[mw<1]=1; mw = log(mw) }
      # mwl = mw; mwl[mwl<1]=1; mwl = log(mwl)
      if (is.null(m) | !length(m)>0) next
      mw[mw<1]=1; mwl = exp(log(mw,getlv(max(mw),100)))/100
      
      if (doneAll) dname = paste(dist_dir, "/", mcp, "_manhattan_FULL_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_normalize-cellpop", sep="")
      dname = paste(dist_dir, "/", mcp, "_manhattan_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_normalize-cellpop", sep = "" )
      if (file.exists(paste0(gsub("FULL_","",dname),"_weighted.Rdata")) & !overwrite) { cat(" exists, skipped."); next }
      
      
      #calculate distance manhattan

      if (is.null(dim(m))) {
        dm = matrix(0,nrow=nrow(m[[1]]),ncol=nrow(m[[1]]))
        for (i in 2:nrow(dm)) {
          
          dm[1:(i-1),i] = dm[i,1:(i-1)] = foreach(j=1:(i-1),.combine='c') %dopar% {
            x = sapply(m, function(a) {
              if (is.null(dim(a))) a = matrix(a,ncol=1)
              return(mean(abs(a[j,]-a[i,])))
            })
            #distance metric (sum)
            weight = pmax(mwl[i,],mwl[j,])
            dist_m = sum( weight[x!=0] * x[x!=0] )
            return(dist_m)
          }
          
        }
        
      } else {
        dm = matrix(0,nrow=nrow(m),ncol=nrow(m))
        for (i in 2:nrow(dm)) {
          dm[1:(i-1),i] = dm[i,1:(i-1)] = foreach(j = 1:(i-1), .combine="c") %dopar% {
            x = abs(m[j,]-m[i,])
            weight = pmax(mwl[i,],mwl[j,])
            dist_m = sum(weight[x!=0]*x[x!=0])
            return(dist_m)
          }
        }
      }
      rownames(dm) = colnames(dm) = gt
      save(dm,file=paste0(dname,"_weighted.Rdata"))
      #weigh each feature matrix; scalendev = 0 (uhh...), T (scaling is done by dividing the centred columns by their SD), F (none)
      TimeOutput(start2)
    }
    
    
  }
  TimeOutput(start2)

}

TimeOutput(start)

