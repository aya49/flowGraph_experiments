## Input: original normalized count matrix + p values --> Output: a combined RchyOptimyx plot for each file
# aya43@sfu.ca 20170131

#Directory
root = "~/projects/IMPC"
result_dir = "result"
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#, "Sanger_MLN","CIPHE","TCP","H")

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrix", sep="")

#Output
rchy_dir = paste(result_dir, "/", panelL, "/", centreL, "/rchy", sep=""); for (ri in 1:length(rchy_dir)) { suppressWarnings(dir.create(rchy_dir[ri])) }

#Libraries/Functions
source("~/projects/IMPC/code/_funcAlice.R")
libr("Matrix")
libr("stringr")
libr("foreach")
libr("doMC")
libr("RchyOptimyx")

#Setup Cores
no_cores = detectCores()-2
registerDoMC(no_cores)




#Options for script
cellCountThres = c(200,750) #insignificant if count under

kpaths = 1#4
doall = F #if true, will get an rchyoptimyx for all significant nodes, if not do all, will do in order of most to least significant until significant nodes are covered by paths made thus far

matrix_count = c("CountAdj")
# matrix_TRIM = c("matrixPvalTRIM_CountAdj.Rdata") #matrix with 0s for insignificant nodes
# matrix_rchy = c("matrixPval_CountAdj.Rdata") #matrix as input into rchy

weight_matrix = c("MaxCountAdj_CountAdj") #assume by cell pop

matrix_type = c("PvalTRIM_CountAdj", "PvalTRIM_CountAdjBH", "PvalTRIM_CountAdjBY")
matrix_all_type = c(NA,NA,NA) #NA if just the above without TRIM
matrix_edge_type = c("Parent_contrib_CountAdj","Parent_contrib_CountAdj","Parent_contrib_CountAdj") #plot only










start = Sys.time()

for (ci in length(paste0(panelL,centreL)):1) {
  centre = paste0(panelL," ",centreL)[ci]
  cat("\n",centre,sep="")
  
  start1 = Sys.time()
  
  
  
  
  
  m0 = get(load(paste0(matrix_dir[ci], matrix_count,".Rdata")))
  phenoMeta = get(load(phenoMeta_dir))
  
  k0 = c(max(phenoMeta$phenolevel))#,1,4) # how many markers to consider i.e. k=max(phenolevel) only
  
  
  #get list of children for each non-leaf node & save
  
  #Prepare data
  for (mti in 1:length(matrix_type)) {
    tryCatch({
      mcp = matrix_type[mti]
      map = gsub("TRIM","",mcp)
      if (!is.na(matrix_all_type[mti])) map = matrix_all_type[mti]
      #mep = matrix_edge_type[mti]
      cat("\n", mcp, " ",sep="")
      start2 = Sys.time()
      
      #start a list of phenotypes included in each distance matrix calculation such that none are repeated
      leavePhenotype = list()
      doneAll = F
      
      #load feature matrix
      mresult = Loadintermatrices(c(paste0(matrix_dir[ci], mcp,".Rdata"),paste0(matrix_dir[ci], map,".Rdata"),paste0(matrix_dir[ci], mep,".Rdata")))
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
          
          fname = paste0(rchy_dir[ci],"/doAll-",doall,"_dir_k=",kpaths,"_",matrix_all_type[mti],"_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres); suppressWarnings(dir.create(fname))
          fname0 = gsub("dir_k","all_k",fname); suppressWarnings(dir.create(fname0))
          
          a = foreach(ii=1:nrow(m)) %dopar% {
            i = rownames(m)[ii]
            phenocodes = phenoMeta$phenocode[match(colnames(ma),phenoMeta$phenotype)]
            
            
            # either change pvalues
            rchy0 = NULL
            if (any(m[ii,]!=0)) {
              phenoscore = abs(ma[i,])
              startphenoscore = c(m[ii,-1])[m[ii,-1]!=0]
              startpheno = phenoMeta$phenocode[match(colnames(m)[m[ii,]!=0],phenoMeta$phenotype)]
              startpheno = setdiff(startpheno,paste0(rep(0,nchar(phenoMeta$phenocode[1]))))
              rchy0 = list()
              if (doall) {
                rchy0 = lapply(startpheno, function(sp) {
                  rr = RchyOptimyx(pheno.codes=phenocodes, phenotypeScores=phenoscore, startPhenotype=sp, pathCount=kpaths, trimPaths=F)
                  rr@nodes[2,] = ma[ii,match(rr@nodes[1,],phenoMeta$phenocode[match(colnames(ma),phenoMeta$phenotype)])]
                  rownames(rr@nodes) = c("phenocode",matrix_type[mti],"colour")
                  return(rr)
                })
                names(rchy0) = startpheno
                
              } else {
                startpheno1 = startpheno[order(startphenoscore,decreasing=T)]
                while (length(startpheno1)>0) {
                  sp = startpheno1[1]
                  rr = RchyOptimyx(pheno.codes=phenocodes, phenotypeScores=phenoscore, startPhenotype=sp, pathCount=kpaths, trimPaths=F)
                  rr@nodes[2,] = ma[ii,match(rr@nodes[1,],phenoMeta$phenocode[match(colnames(ma),phenoMeta$phenotype)])]
                  rownames(rr@nodes) = c("phenocode",matrix_type[mti],"colour")
                  rchy0[[sp]] = rr
                  startpheno1 = setdiff(startpheno1[-1],rr@nodes["phenocode",])
                }
              }
              # rchy = rchy1[[1]]; if (length(rchy1)>1) for (ri in 2:length(rchy1)) { rchy = merge(rchy,rchy1[[ri]]) }
            }
            
            rchy0_m = NULL
            if (!is.null(rchy0)) {
              rchy0_m = rchy0[[1]]
              if (length(rchy0)>1) {
                for (ri in 2:length(rchy0)) {
                  rchy0_m = merge(rchy0_m,rchy0[[ri]])
                }
              }
            }
            
            
            
            rchy = NULL
            # positive change pvalues
            rchy1 = NULL
            if (sum(m[ii,]>0)>0) {
              phenoscore = ma[ii,]
              startphenoscore = c(m[ii,-1])[m[ii,-1]!=0]
              if (min(phenoscore)<0) {
                phenoscore = phenoscore+min(-ma[ii,])
                startphenoscore = startphenoscore+min(-ma[ii,])
              }
              startpheno = phenoMeta$phenocode[match(colnames(m)[m[i,]>0],phenoMeta$phenotype)]
              startpheno = setdiff(startpheno,paste0(rep(0,nchar(phenoMeta$phenocode[1]))))
              rchy1 = list()
              if (doall) {
                rchy1 = lapply(startpheno, function(sp) {
                  rr = RchyOptimyx(pheno.codes=phenocodes, phenotypeScores=phenoscore, startPhenotype=sp, pathCount=kpaths, trimPaths=F)
                  rr@nodes[2,] = ma[ii,match(rr@nodes[1,],phenoMeta$phenocode[match(colnames(ma),phenoMeta$phenotype)])]
                  rownames(rr@nodes) = c("phenocode",matrix_type[mti],"colour")
                  return(rr)
                })
                names(rchy1) = startpheno
                # rchy = rchy1[[1]]; if (length(rchy1)>1) for (ri in 2:length(rchy1)) { rchy = merge(rchy,rchy1[[ri]]) }
              } else {
                startpheno1 = startpheno[order(startphenoscore,decreasing=T)]
                while (length(startpheno1)>0) {
                  sp = startpheno1[1]
                  rr = RchyOptimyx(pheno.codes=phenocodes, phenotypeScores=phenoscore, startPhenotype=sp, pathCount=kpaths, trimPaths=F)
                  rr@nodes[2,] = ma[ii,match(rr@nodes[1,],phenoMeta$phenocode[match(colnames(ma),phenoMeta$phenotype)])]
                  rownames(rr@nodes) = c("phenocode",matrix_type[mti],"colour")
                  rchy1[[sp]] = rr
                  startpheno1 = setdiff(startpheno1[-1],rr@nodes["phenocode",])
                }
              }
            }
            # negative change pvalues
            rchy2 = NULL
            if (sum(m[ii,]<0)>0) {
              phenoscore = -ma[ii,]
              startphenoscore = -c(m[ii,-1])[m[ii,-1]!=0]
              if (min(phenoscore)<0) {
                phenoscore = phenoscore+min(-ma[ii,])
                startphenoscore = startphenoscore+min(-ma[ii,])
              } 
              startpheno = phenoMeta$phenocode[match(colnames(m)[m[i,]<0],phenoMeta$phenotype)]
              startpheno = setdiff(startpheno,paste0(rep(0,nchar(phenoMeta$phenocode[1])),collapse=""))
              rchy2 = list()
              if (doall) {
                rchy2 = lapply(startpheno, function(sp) {
                  rr = RchyOptimyx(pheno.codes=phenocodes, phenotypeScores=phenoscore, startPhenotype=sp, pathCount=kpaths, trimPaths=F)
                  rr@nodes[2,] = ma[ii,match(rr@nodes[1,],phenoMeta$phenocode[match(colnames(ma),phenoMeta$phenotype)])]
                  rownames(rr@nodes) = c("phenocode",matrix_type[mti],"colour")
                  return(rr)
                })
                names(rchy2) = startpheno
                # if (!sum(m[i,]>0)>0) rchy = rchy2[[1]]
                # rchy = merge(rchy,rchy2[[1]]); if (length(rchy2)>1) for (ri in 2:length(rchy2)) { rchy = merge(rchy,rchy2[[ri]]) }
              } else {
                startpheno1 = startpheno[order(startphenoscore,decreasing=T)]
                while (length(startpheno1)>0) {
                  sp = startpheno1[1]
                  rr = RchyOptimyx(pheno.codes=phenocodes, phenotypeScores=phenoscore, startPhenotype=sp, pathCount=kpaths, trimPaths=F)
                  rr@nodes[2,] = ma[ii,match(rr@nodes[1,],phenoMeta$phenocode[match(colnames(ma),phenoMeta$phenotype)])]
                  rownames(rr@nodes) = c("phenocode",matrix_type[mti],"colour")
                  rchy2[[sp]] = rr
                  startpheno1 = setdiff(startpheno1[-1],rr@nodes["phenocode",])
                }
              }
            }
            rchy = append(rchy1,rchy2)
            
            rchy_m = NULL
            if (!is.null(rchy)) {
              rchy_m = rchy[[1]]
              if (length(rchy)>1) {
                for (ri in 2:length(rchy)) {
                  rchy_m = merge(rchy_m,rchy[[ri]])
                }
              }
            }
            
            if (!is.null(rchy0)) save(rchy0,file=paste0(fname0,"/",i,".Rdata"))
            if (!is.null(rchy)) save(rchy,file=paste0(fname,"/",i,".Rdata"))
            if (!is.null(rchy0_m)) save(rchy0,file=paste0(fname0,"/",i,"_merged.Rdata"))
            if (!is.null(rchy_m)) save(rchy,file=paste0(fname,"/",i,"_merged.Rdata"))
            
          }
          
        } #layer
      } #countThres
      
      TimeOutput(start2)
    }, error = function(err) { cat(paste("ERROR:  ",err)) })
    
  } #matrix_type
  
}

TimeOutput(start)