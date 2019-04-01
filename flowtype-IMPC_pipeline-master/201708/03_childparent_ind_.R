# aya43@sfu.ca 20170131
# creates list of children for each non-leaf node and creates children/parent list[[matrices]] (-/+ are only for phenotypes where both -,+ data exists)

#root = "~/projects/IMPC"
result_dir = "result"
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#, "Sanger_MLN","CIPHE","TCP","H")

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoMeta", sep="")
# matrixCount_dir = paste(result_dir, "/matrixCount.Rdata", sep="")
pheno_type = c("")#, "Prop")


# createChildEntropy = T
# createParentEntropy = T

#Output
phenoChild_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChild",sep="")
phenoChild_ind_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChild_ind",sep="")
phenoChildpn_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChildpn",sep="")
phenoChildpn_ind_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChildpn_ind",sep="")
phenoParent_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoParent",sep="")
phenoParent_ind_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoParent_ind",sep="")
phenoParentpn_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoParentpn",sep="")
phenoParentpn_ind_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoParentpn_ind",sep="")

libr(foreach)
libr(doMC)
source("~/projects/IMPC/code/_funcAlice.R")



start = Sys.time()

no_cores = detectCores() - 1
registerDoMC(no_cores)


for (ci in length(paste0(panelL,centreL)):1) {
  centre = paste0(panelL," ",centreL)[ci]
  
  start1 = Sys.time()
  
  for (pcp in pheno_type) {
    cat("\n",pheno_type,sep="")
    
    phenoMeta = get(load(paste0(phenoMeta_dir[ci], pcp,".Rdata")))
    # colnames(phenoMeta) = c("phenotype","phenocode","phenolevel")
    # save(phenoMeta,file=paste0(phenoMeta_dir[ci], pcp,".Rdata"))
    
    #get list of children for each non-leaf node & save
    cat("\ncreating children indices ")
    pc = getphenoChild(phenoMeta, no_cores=no_cores)
    
    phenoChild = pc$phenoChild
    phenoChild_ind = pc$phenoChild_ind
    phenoChildpn = pc$phenoChildpn
    phenoChildpn_ind = pc$phenoChildpn_ind
    
    save(phenoChild, file=paste0(phenoChild_dir[ci], pcp,".Rdata"))
    save(phenoChild_ind, file=paste0(phenoChild_ind_dir[ci], pcp,".Rdata"))
    save(phenoChildpn, file=paste0(phenoChildpn_dir[ci], pcp,".Rdata"))
    save(phenoChildpn_ind, file=paste0(phenoChildpn_ind_dir[ci], pcp,".Rdata"))
    
    TimeOutput(start1)
    
    cat("\ncreating parent indices ")
    pp = getphenoParent(phenoMeta, phenoChildpn, phenoChildpn_ind, no_cores)
    
    phenoParent = pp$phenoParent
    phenoParent_ind = pp$phenoParent_ind
    phenoParentpn = pp$phenoParentpn #only children with a positive/negative counterpart information available
    phenoParentpn_ind = pp$phenoParentpn_ind
    
    save(phenoParent, file=paste0(phenoParent_dir[ci], pcp,".Rdata"))
    save(phenoParent_ind, file=paste0(phenoParent_ind_dir[ci], pcp,".Rdata"))
    save(phenoParentpn, file=paste0(phenoParentpn_dir[ci], pcp,".Rdata")) 
    save(phenoParentpn_ind, file=paste0(phenoParentpn_ind_dir[ci], pcp,".Rdata"))
    
    TimeOutput(start1)
  }
}

TimeOutput(start)




