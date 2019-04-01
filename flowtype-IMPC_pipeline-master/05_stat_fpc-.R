# aya43@sfu.ca 20161220
# combine different features and discriminant plot

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
matrix_count = c("CountAdj")
matrix_type = c("CountAdj", "Child_entropy", 
                "LogFold_CountAdj", "Pval_CountAdj", "LogFoldTRIM_CountAdj", "PvalTRIM_CountAdj", "Child_entropy_CountAdj",
                "Parent_effort", "Parent_contrib", "Child_pnratio", "Child_prop",
                "Parent_effortTRIM", "Parent_contribTRIM", "Child_pnratioTRIM", "Child_propTRIM")
phenoChild_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChild.Rdata",sep="")
phenoChild_ind_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChild_ind.Rdata",sep="")
phenoChildpn_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChildpn.Rdata",sep="")
phenoChildpn_ind_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChildpn_ind.Rdata",sep="")

cellCountThres = c(200,200,500,500,500) #insignificant if count under
pclustermethods = c("dc","mvdc","nc","wnc")#,"bc","vbc","adc","awc","arc","anc") #assume two clusters

interested = c("gene", "gender", "date", "colony", "strain", "fur", "birth_date", "specimen", "sample") #sampleMeta columns to plot


#Output
plot_dir = paste(result_dir, "/", panelL, "/", centreL, "/plots", sep="")
pc_dir = paste(plot_dir, "/plotcluster", sep=""); for(i in 1:length(pc_dir)) { suppressWarnings(dir.create(pc_dir[i])) }

source("~/projects/IMPC/code/_funcAlice.R")
libr("stringr")
libr("colorspace")
libr("fpc") # libr(proxy)
libr("foreach")
libr("doMC")
libr("lubridate") #if there are date variables


no_cores = detectCores()-1
registerDoMC(no_cores)








start = Sys.time()


#for (ci in 4:1) {
for (ci in 1:length(paste0(panelL,centreL))) {
  start1 = Sys.time()
  centre = paste0(panelL," ",centreL)[ci]
  cat("\n",centre,sep="")
  
  phenoMeta = get(load(phenoMeta_dir[ci]))
  sampleMeta = get(load(sampleMeta_dir[ci]))
  #k0 = c(1,3,max(phenoMeta$phenolevel)) # how many markers to consider i.e. k=max(phenolevel) only
  k0 = c(max(phenoMeta$phenolevel)) # how many markers to consider i.e. k=max(phenolevel) only
  m0 = get(load(paste0(matrix_dir[ci], matrix_count,".Rdata")))
  interestedCols = which(colnames(sampleMeta)%in%interested)
  
  for (mcp in matrix_type) {
    cat("\n", mcp, " ",sep="")
    start2 = Sys.time()
    leavePhenotype = list()
    doneAll = F
    if (!file.exists(paste0(matrix_dir[ci], mcp,".Rdata"))) {cat("doesn't exist"); next}
    mm = get(load(paste0(matrix_dir[ci], mcp,".Rdata")))
    
    #get to-delete low count phenotype indices; CountAdj should be first one
    for (countThres in cellCountThres[ci]) {
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
        if (!is.null(dim(m))) { mpi = colnames(m)
        } else { mpi = names(m) }
        
        lpi = setdiff(mpi,dpi0) #phenotypes to leave in analysis
        
        #check if indices already calculated for on this matrix, if yes, skip, else trim and do
        if (Position(function(x) identical(x, lpi), leavePhenotype, nomatch = 0) > 0) {cat("; child skipped", sep=""); next}
        leavePhenotype[[paste0(k="k.",k,"_countThres.",countThres)]] = lpi
        
        if (!doneAll & length(lpi)==length(mpi)) doneAll = T
        
        if (!is.null(dim(mm))) { m = m[,match(lpi,mpi)]
        } else {
          ml = m[match(lpi,mpi)]
          m = foreach(k=1:length(m),.combine="cbind") %dopar% { return(m[[k]]) }
        }
        
        pm = phenoMeta[match(lpi,phenoMeta$phenotype),]
        sm = sampleMeta[match(rownames(m),sampleMeta$fileName),]

        foreach(pcm=pclustermethods) {
          for (coli in interestedCols) {
            
          }
          pcname = paste0(pc_dir[ci],"/plotcluster", mcp, "_", pcm, "_col-",colnames(sm)[coli],"_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, ".png")
          plotcluster(m,as.numeric(factor(sm[,coli])))
          
          plotcluster(face,grface, pclustermethods[1])
          plotcluster(face,grface, pclustermethods[2])
          plotcluster(face,grface, pclustermethods[3])
          plotcluster(face,grface, pclustermethods[4])
        }
      }
    }
  }
}

TimeOutput(start)


