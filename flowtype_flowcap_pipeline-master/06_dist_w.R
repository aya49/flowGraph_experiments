## Input: original features --> Output: weighted manhattan distance matrices
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
libr("progress")
libr("pracma") #for dot product

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)










#Options for script
overwrite = F
cellCountThres = c(200) #insignificant if count under

weightings = "" #"orig-"

matrix_count = matrix_weight = c("CountAdj") #matrix_weight = NULL if no weight matrix



#Prepare data
matrix_type = list.files(path=result_dir,pattern=glob2rx("matrix*.Rdata"))
matrix_type = matrix_type[!matrix_type=="matrixCount.Rdata"]
matrix_type = gsub("matrix","",matrix_type)
matrix_type = gsub(".Rdata","",matrix_type)
matrix_type = matrix_type[!grepl("KO|Max",matrix_type)]
# matrix_type = c(matrix_type[!grepl("pnratio|prop|effort|contrib",matrix_type)],matrix_type[grepl("pnratio|prop|effort|contrib",matrix_type)])
# matrix_type = c(matrix_type[grepl("TRIM|Freqp",matrix_type) & is.na(as.integer(sapply(matrix_type,function(x)substr(x,nchar(x),nchar(x)))))])#, matrix_type[!grepl("TRIM",matrix_type)])

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

phenoMeta = get(load(phenoMeta_dir))
k0 = c(1,4,max(phenoMeta$phenolevel)) # how many markers to consider i.e. k=max(phenolevel) only

m0 = get(load(paste0(matrix_dir, matrix_count,".Rdata")))










start = Sys.time()

#load different matrices
for (mcp in matrix_type) {
  cat("\n", mcp, " ",sep="")
  start2 = Sys.time()
  
  #start a list of phenotypes included in each distance matrix calculation such that none are repeated
  leavePhenotype = list()
  doneAll = F
  
  #load feature matrix
  if (!file.exists(paste0(matrix_dir, mcp,".Rdata"))) {cat("doesn't exist"); next}
  mm = get(load(paste0(matrix_dir, mcp,".Rdata")))
  
  matrix_type. = append(paste0(matrix_dir,matrix_weight,".Rdata"), paste0(matrix_dir,mcp,".Rdata"))
  
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
      
      #trim feature matrix
      mmlresult = trimMatrices(mml0,m0,pt,gt,phenoMeta,leavePhenotype,doneAll, countThres,k)
      if (is.null(mmlresult)) next
      mml = mmlresult$mml
      pm = mmlresult$pm
      leavePhenotype = mmlresult$leavePhenotype
      doneAll = mmlresult$doneAll
      
      #seperate out the weight matrix
      m = NULL
      mwl = matrix(1,ncol=nrow(pm),nrow=length(gt),dimnames=(list(gt,pm$phenotype)))
      if (!is.null(matrix_weight)) {
        mw = mml[[1]]; if (length(mml)>1) m = mml[[2]] 
      } else { if (length(mml)>0) m = mml[[1]] }
      if (is.null(m) | !length(m)>0) next
      
      dname = paste(dist_dir, "/", mcp, "_manhattan_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_normalize-cellpop", sep = "" )

      for (weighting in weightings) {
        
        dname1 = paste0(dname,"_",weighting,"weighted.Rdata")
        mw[mw<1]=1
        if (weighting=="") { mwl = exp(log(mw,getlv(max(mw),100)))/100 
        } else { mwl = mw/max(mw) }
        
        if (file.exists(dname1) & !overwrite) { cat(" exists, skipped."); next }

        ## calculate weighted distance manhattan
        if (is.null(dim(m))) { #edge feature
          dm = matrix(0,nrow=nrow(m[[1]]),ncol=nrow(m[[1]]))
          pb <- progress_bar$new( format = "  calculating [:bar] :percent eta: :eta", total = nrow(dm), clear = FALSE, width= 60)
          for (i in 2:nrow(dm)) {
            pb$tick();Sys.sleep(1 / nrow(dm))
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
          
        } else { #node feature
          dm = matrix(0,nrow=nrow(m),ncol=nrow(m))
          pb <- progress_bar$new( format = "  calculating [:bar] :percent eta: :eta", total = nrow(dm), clear = FALSE, width= 60)
          for (i in 2:nrow(dm)) {
            pb$tick();Sys.sleep(1 / nrow(dm))
            dm[1:(i-1),i] = dm[i,1:(i-1)] = foreach(j = 1:(i-1), .combine="c") %dopar% {
              x = abs(m[j,]-m[i,])
              weight = pmax(mwl[i,],mwl[j,])
              dist_m = sum(weight[x!=0]*x[x!=0])
              return(dist_m)
            }
          }
        }
        rownames(dm) = colnames(dm) = rownames(mwl)
        save(dm,file=dname1)
        #save(dm,file=paste0(dname,"_weighted.Rdata"))
        #weigh each feature matrix; scalendev = 0 (uhh...), T (scaling is done by dividing the centred columns by their SD), F (none)
        TimeOutput(start2)
      }
    }
    
  }
  TimeOutput(start2)
}
TimeOutput(start)

