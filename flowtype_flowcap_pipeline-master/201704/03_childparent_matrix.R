# creates matrices
# aya43@sfu.ca 20170131

root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
matrixCountAdj_dir = paste(result_dir, "/matrixCountAdj.Rdata", sep="") #the entorpy and child proportions are based on original count matrix
phenoChild_dir = paste(result_dir, "/phenoChild.Rdata",sep="")
phenoChild_ind_dir = paste(result_dir, "/phenoChild_ind.Rdata",sep="")
phenoChildpn_dir = paste(result_dir, "/phenoChildpn.Rdata",sep="")
phenoChildpn_ind_dir = paste(result_dir, "/phenoChildpn_ind.Rdata",sep="")
phenoParent_dir = paste(result_dir, "/phenoParent.Rdata",sep="")
phenoParent_ind_dir = paste(result_dir, "/phenoParent_ind.Rdata",sep="")
phenoParentpn_dir = paste(result_dir, "/phenoParentpn.Rdata",sep="")
phenoParentpn_ind_dir = paste(result_dir, "/phenoParentpn_ind.Rdata",sep="")
matrixLogFold_dir = paste(result_dir, "/matrixLogFold",sep="")
matrixKOCountAdj_dir = paste(result_dir, "/matrixKOCountAdj",sep="")
matrixPvalTRIM_dir = paste(result_dir, "/matrixPvalTRIM",sep="")

pheno_type = c("CountAdj", "Prop")

trimOnly = F

#Output
matrixChild_pnratio_dir = paste(result_dir, "/matrixChild_pnratio",sep="")
matrixChild_prop_dir = paste(result_dir, "/matrixChild_prop",sep="")
matrixChild_entropy_dir = paste(result_dir, "/matrixChild_entropy",sep="")
matrixParent_entropy_dir = paste(result_dir, "/matrixParent_entropy",sep="")
matrixParent_contrib_dir = paste(result_dir, "/matrixParent_contrib",sep="")
matrixParent_effort_dir = paste(result_dir, "/matrixParent_effort",sep="")

libr(Matrix)
libr(stringr)
libr(entropy)
libr(foreach)
libr(doMC)
source("~/projects/IMPC/code/_funcAlice.R")



start = Sys.time()

no_cores = detectCores() - 1
registerDoMC(no_cores)



#get list of children for each non-leaf node & save

m = get(load(matrixCountAdj_dir))
phenoMeta = get(load(paste0(phenoMeta_dir,".Rdata")))
sampleMeta = get(load(sampleMeta_dir))

phenoChild = get(load(phenoChild_dir))
phenoChild_ind = get(load(phenoChild_ind_dir))
phenoChildpn = get(load(phenoChildpn_dir))
phenoChildpn_ind = get(load(phenoChildpn_ind_dir))

mpt = list()
mlf = list()
mkc = list()
for (pcpi in 1:length(pheno_type)) {
  pcp = pheno_type[pcpi]
  mpt[[pcpi]] = get(load(paste0(matrixPvalTRIM_dir, "_",pcp,".Rdata")))
  mlf[[pcpi]] = get(load(paste0(matrixLogFold_dir, "_",pcp,".Rdata")))
  mkc[[pcpi]] = get(load(paste0(matrixKOCountAdj_dir, "_",pcp,".Rdata")))
}
phenoParent = get(load(phenoParent_dir))
phenoParent_ind = get(load(phenoParent_ind_dir))
phenoParentpn = get(load(phenoParentpn_dir))
phenoParentpn_ind = get(load(phenoParentpn_ind_dir))







## Childentropy ---------------------------------------------------------------
cat(", childentropy")

if (!trimOnly) {
  matrixChild_entropy = foreach(i = 1:length(phenoChildpn), .combine='cbind') %dopar% { #for each phenotype
    en = rep(0,nrow(m))
    
    parent = m[,phenoChildpn_ind[i]]
    children = m[,unlist(phenoChildpn[[i]])]
    if (length(phenoChildpn[[i]])==1) children = matrix(children,ncol=1)
    
    no_child = length(phenoChildpn[[i]][[1]])
    non0parents = which(parent>0)
    en[non0parents] = sapply(non0parents, function(x) return(entropy(children[x,]/parent[x])/no_child)) #average entropy over # of markers added
    names(en) = rownames(m)
    
    #mlist[[i]] = list(pnratio=pnratio, childprop=childprop)
    return(entropy=en) #ratio = +child_prop / -child_prop
  }
  colnames(matrixChild_entropy) = phenoMeta$phenotype[phenoChildpn_ind]
  save(matrixChild_entropy, file=paste0(matrixChild_entropy_dir,".Rdata"))
  write.csv(matrixChild_entropy, file=paste0(matrixChild_entropy_dir,".csv"))
} else {
  load(paste0(matrixChild_entropy_dir,".Rdata"))
  colnames(matrixChild_entropy) = phenoMeta$phenotype[phenoChildpn_ind]
  save(matrixChild_entropy, file=paste0(matrixChild_entropy_dir,".Rdata"))
  write.csv(matrixChild_entropy, file=paste0(matrixChild_entropy_dir,".csv"))
}

#trim based on pvalue trim
for (pcpi in 1:length(pheno_type)) { pcp = pheno_type[pcpi]
  cat(" ",pheno_type[pcpi],sep="")
  lpi = intersect(colnames(matrixChild_entropy), colnames(mpt[[pcpi]]))
  matrixChild_entropyTRIM = matrixChild_entropy[!is.na(match(rownames(matrixChild_entropy),rownames(mpt[[pcpi]]))),match(lpi,colnames(matrixChild_entropy))]
  matrixChild_entropyTRIM[as.matrix(mpt[[pcpi]][,!is.na(match(colnames(mpt[[pcpi]]),colnames(matrixChild_entropyTRIM)))])==0] = 0
  matrixChild_entropyTRIM = Matrix(matrixChild_entropyTRIM, sparse=T)
  save(matrixChild_entropyTRIM, file=paste0(matrixChild_entropy_dir,"TRIM_",pcp,".Rdata"))
  write.csv(as.matrix(matrixChild_entropyTRIM), file=paste0(matrixChild_entropy_dir,"TRIM_",pcp,".csv"))
}
# for (pcpi in 1:length(pheno_type)) { pcp = pheno_type[pcpi]
#   cat("\n",pheno_type,sep="")
#   pm = get(load(paste0(phenoMeta_dir,pcp,".Rdata")))
#   phenind = match(names(pnratio), pm$phenotype)
#   
#   pnratio_trim <- pnratio[-is.na(phenind)] #ind of pm
#   matrixChild_entropy_trim <- matrixChild_entropy[,-is.na(phenind)]
#   
#   save(pnratio_trim, file=paste0(matrixChild_pnratio_dir,pcp,".Rdata"))
#   save(matrixChild_entropy, file=paste0(matrixChild_entropy_dir,pcp,".Rdata"))
#   write.csv(Reduce('cbind',pnratio_trim), file=paste0(matrixChild_pnratio_dir,pcp,".csv"))
#   write.csv(matrixChild_entropy, file=paste0(matrixChild_entropy_dir,pcp,".csv"))
# }


TimeOutput(start)







## Parententropy ---------------------------------------------------------------


cat(", parententropy")
if (!trimOnly) {
  matrixParent_entropy = foreach(i = 1:length(phenoParentpn), .combine='cbind') %dopar% { #for each phenotype
    en = rep(0,nrow(m))
    
    childr = m[,phenoParentpn_ind[i]]
    parent = m[,phenoParentpn[[i]]]
    if (length(phenoParentpn[[i]])==1) parent = matrix(parent,ncol=1)
    
    no_parent = length(phenoParentpn[[i]])
    non0childrs = which(childr>0)
    en[non0childrs] = sapply(non0childrs, function(x) return(entropy(childr[x]/parent[x,])/no_parent)) #average entropy over # of markers added
    names(en) = rownames(m)
    return(en) #ratio = +child_prop / -child_prop
  }
  
  colnames(matrixParent_entropy) = phenoMeta$phenotype[phenoParentpn_ind]
  save(matrixParent_entropy, file=paste0(matrixParent_entropy_dir,".Rdata"))
  write.csv(matrixParent_entropy, file=paste0(matrixParent_entropy_dir,".csv"))
} else {
  load(paste0(matrixParent_entropy_dir,".Rdata"))
  colnames(matrixParent_entropy) = phenoMeta$phenotype[phenoParentpn_ind]
  save(matrixParent_entropy, file=paste0(matrixParent_entropy_dir,".Rdata"))
  write.csv(matrixParent_entropy, file=paste0(matrixParent_entropy_dir,".csv"))
}

#trim based on pvalue trim
for (pcpi in 1:length(pheno_type)) { pcp = pheno_type[pcpi]
  lpi = intersect(colnames(matrixParent_entropy), colnames(mpt[[pcpi]]))
  matrixParent_entropyTRIM = matrixParent_entropy[!is.na(match(rownames(matrixParent_entropy),rownames(mpt[[pcpi]]))), match(lpi,colnames(matrixParent_entropy))]
  matrixParent_entropyTRIM[as.matrix(mpt[[pcpi]][,!is.na(match(colnames(mpt[[pcpi]]),colnames(matrixParent_entropyTRIM)))])==0] = 0
  save(matrixParent_entropyTRIM, file=paste0(matrixParent_entropy_dir,"TRIM_",pcp,".Rdata"))
  write.csv(as.matrix(matrixParent_entropyTRIM), file=paste0(matrixParent_entropy_dir,"TRIM_",pcp,".csv"))
}
# for (pcpi in 1:length(pheno_type)) { pcp = pheno_type[pcpi]
#   cat("\n",pheno_type,sep="")
#   pm = get(load(paste0(phenoMeta_dir,pcp,".Rdata")))
#   phenind = match(colnames(matrixParent_entropy), pm$phenotype)
#   
#   matrixParent_entropy_trim <- matrixParent_entropy[,-is.na(phenind)]
#   
#   save(matrixParent_entropy, file=paste0(matrixParent_entropy_dir,pcp,".Rdata"))
#   write.csv(matrixParent_entropy, file=paste0(matrixParent_entropy_dir,pcp,".csv"))
# }
TimeOutput(start)














## Childprop ---------------------------------------------------------------



#mlist=list()
cat("; childprop")
if (!trimOnly) {
  childprop = foreach(i = 1:length(phenoChild)) %dopar% { #for each phenotype
    #for (i in 1:length(phenoChild)) {
    parent = m[,phenoChild_ind[i]]
    
    # Child proportion matrix
    children = m[,unlist(phenoChild[[i]])]
    children[which(children<1)] = 1
    parent0 = parent
    parent0[which(parent<1)] = 1
    childprop = exp(children/parent0)
    rownames(childprop) = rownames(m)
    return(childprop)
  }
  names(childprop) = phenoMeta$phenotype[phenoChild_ind]
  save(childprop, file=paste0(matrixChild_prop_dir,".Rdata"))
  childpropREDUCE = foreach(k=1:length(childprop),.combine="cbind") %dopar% { return(childprop[[k]]) }
  write.csv(childpropREDUCE, file=paste0(matrixChild_prop_dir,".csv"))
} else {
  load(paste0(matrixChild_prop_dir,".Rdata"))
  names(childprop) = phenoMeta$phenotype[phenoChild_ind]
  save(childprop, file=paste0(matrixChild_prop_dir,".Rdata"))
  childpropREDUCE = foreach(k=1:length(childprop),.combine="cbind") %dopar% { return(childprop[[k]]) }
  write.csv(childpropREDUCE, file=paste0(matrixChild_prop_dir,".csv"))
}

#trim based on pvalue (phenotype)
for (pcpi in 1:length(pheno_type)) { pcp = pheno_type[pcpi]
  cat(" ",pheno_type[pcpi],sep="")
  childpropTRIM = trimlistmatrix(mpt[[pcpi]],refdel=0, childprop,targetdel=0, ncore=no_cores)
  save(childpropTRIM, file=paste0(matrixChild_prop_dir,"TRIM_",pcp,".Rdata"))
  childpropTRIMREDUCE = foreach(k=1:length(childpropTRIM),.combine="cbind") %dopar% { return(childpropTRIM[[k]]) }
  write.csv(as.matrix(childpropTRIMREDUCE), file=paste0(matrixChild_prop_dir,"TRIM_",pcp,".csv"))
}
TimeOutput(start)











## Childpnratio ---------------------------------------------------------------
cat(", childpnratio")

if (!trimOnly) {
  pnratio = foreach(i = 1:length(phenoChildpn)) %dopar% { #for each phenotype
    parent = m[,phenoChildpn_ind[i]]
    
    # P/N ratio matrix
    pnratio = NULL
    
    pos = m[,phenoChildpn_ind[i]]-m[,phenoChildpn[[i]][[1]]]
    neg = m[,phenoChildpn[[i]][[1]]]
    pnratio = pos/neg ##get rid of 0, Inf
    pnratio[pos==0 & neg==0] = 1
    pnratio[pos==0 & neg!=0] = 1/neg[pos==0 & neg!=0]
    pnratio[pos!=0 & neg==0] = pos[pos!=0 & neg==0]
    pnratio = matrix(log(pnratio),nrow=nrow(m))
    
    if (sum(parent==0)>0) {
      pnratio[which(parent==0),] = rep(0,length(phenoChildpn[[i]][[1]]))
    } 
    rownames(pnratio) = rownames(m)
    colnames(pnratio) = colnames(pos)
    
    return(pnratio) #ratio = +child_prop / -child_prop
  }
  names(pnratio) = phenoMeta$phenotype[phenoChildpn_ind]
  save(pnratio, file=paste0(matrixChild_pnratio_dir,".Rdata"))
  pnratioREDUCE = foreach(k=1:length(pnratio),.combine="cbind") %dopar% { return(pnratio[[k]]) }
  write.csv(pnratioREDUCE, file=paste0(matrixChild_pnratio_dir,".csv"))
} else {
  load(paste0(matrixChild_pnratio_dir,".Rdata"))
  names(pnratio) = phenoMeta$phenotype[phenoChildpn_ind]
  save(pnratio, file=paste0(matrixChild_pnratio_dir,".Rdata"))
  pnratioREDUCE = foreach(k=1:length(pnratio),.combine="cbind") %dopar% { return(pnratio[[k]]) }
  write.csv(pnratioREDUCE, file=paste0(matrixChild_pnratio_dir,".csv"))
}

#trim based on pvalue trim
for (pcpi in 1:length(pheno_type)) { pcp = pheno_type[pcpi]
  cat(" ",pheno_type[pcpi],sep="")
  pnratioTRIM = trimlistmatrix(mpt[[pcpi]],refdel=0, pnratio,targetdel=0, ncore=no_cores)
  save(pnratioTRIM, file=paste0(matrixChild_pnratio_dir,"TRIM_",pcp,".Rdata"))
  save(matrixChild_entropyTRIM, file=paste0(matrixChild_entropy_dir,"TRIM_",pcp,".Rdata"))
  pnratioTRIMREDUCE = foreach(k=1:length(pnratioTRIM),.combine="cbind") %dopar% { return(pnratioTRIM[[k]]) }
}
# for (pcpi in 1:length(pheno_type)) { pcp = pheno_type[pcpi]
#   cat("\n",pheno_type,sep="")
#   pm = get(load(paste0(phenoMeta_dir,pcp,".Rdata")))
#   phenind = match(names(pnratio), pm$phenotype)
#   
#   pnratio_trim <- pnratio[-is.na(phenind)] #ind of pm
#   matrixChild_entropy_trim <- matrixChild_entropy[,-is.na(phenind)]
#   
#   save(pnratio_trim, file=paste0(matrixChild_pnratio_dir,pcp,".Rdata"))
#   save(matrixChild_entropy, file=paste0(matrixChild_entropy_dir,pcp,".Rdata"))
#   write.csv(Reduce('cbind',pnratio_trim), file=paste0(matrixChild_pnratio_dir,pcp,".csv"))
#   write.csv(matrixChild_entropy, file=paste0(matrixChild_entropy_dir,pcp,".csv"))
# }


TimeOutput(start)











## Parent Contribution ---------------------------------------------------------------
cat(", parentcontribution")
for (pcpi in 1:length(pheno_type)) { pcp = pheno_type[pcpi]
  cat(" ",pcp)
  if (!trimOnly) {
    #build original values of wt and ko
    matrixWTCountAdjFULL = foreach(i=1:ncol(mlf[[pcpi]]),.combine="cbind") %dopar% { return(mkc[[pcpi]][,i]/exp(mlf[[pcpi]][,i])) }
    colnames(matrixWTCountAdjFULL) = colnames(mkc[[pcpi]])
    contrib = foreach(i = 1:length(phenoParent)) %dopar% { #for each phenotype
      change = mkc[[pcpi]][,phenoParent_ind[i]] - matrixWTCountAdjFULL[,phenoParent_ind[i]]
      koparents = mkc[[pcpi]][,phenoParent[[i]]]
      wtparents = matrixWTCountAdjFULL[,phenoParent[[i]]]
      if (length(phenoParent[[i]])==1) {
        koparents = matrix(koparents,ncol=1)
        wtparents = matrix(wtparents,ncol=1)
      }
      changeparents = koparents - wtparents
      if (length(which(changeparents==0))>0) changeparents[changeparents==0] = 1
      
      contrib = change / changeparents
      
      rownames(contrib) = rownames(mkc[[pcpi]])
      colnames(contrib) = colnames(mkc[[pcpi]])[phenoParent[[i]]]
      return(contrib)
    }
    
    names(contrib) <- phenoMeta$phenotype[phenoParent_ind]
    save(contrib, file=paste0(matrixParent_contrib_dir,"_",pcp,".Rdata"))
    contribREDUCE = foreach(k=1:length(contrib),.combine="cbind") %dopar% { return(contrib[[k]]) }
    write.csv(contribREDUCE, file=paste0(matrixParent_contrib_dir,"_",pcp,".csv"))
    
  } else {
    load(paste0(matrixParent_contrib_dir,"_",pcp,".Rdata"))
    names(contrib) <- phenoMeta$phenotype[phenoParent_ind]
    save(contrib, file=paste0(matrixParent_contrib_dir,"_",pcp,".Rdata"))
    contribREDUCE = foreach(k=1:length(contrib),.combine="cbind") %dopar% { return(contrib[[k]]) }
    write.csv(contribREDUCE, file=paste0(matrixParent_contrib_dir,"_",pcp,".csv"))
  } 
  
  
  #trim based on pvalue trim
  cat(" ",pheno_type[pcpi],sep="")
  contribTRIM = trimlistmatrix(mpt[[pcpi]],refdel=0, contrib,targetdel=0, ncore=no_cores)
  save(contribTRIM, file=paste0(matrixParent_contrib_dir,"TRIM_",pcp,".Rdata"))
  a = foreach(k=1:length(contribTRIM),.combine="cbind") %dopar% { return(contribTRIM[[k]]) }; rownames(a) = sampleMeta$gene[match(rownames(a),sampleMeta$fileName)]
  write.table(as.matrix(a), file=paste0(matrixParent_contrib_dir,"TRIM_",pcp,".csv"), sep=",", col.names=F)
}
TimeOutput(start)










## Parent effort -------------------------------------------------------
cat(", parenteffort")
for (pcpi in 1:length(pheno_type)) { pcp = pheno_type[pcpi]
  cat(" ",pcp)
  if (!trimOnly) {
    if (!exists("matrixWTCountAdjFULL")) {
      matrixWTCountAdjFULL = foreach(i=1:ncol(mlf[[pcpi]]),.combine="cbind") %dopar% { return(mkc[[pcpi]][,i]/exp(mlf[[pcpi]][,i])) }
      colnames(matrixWTCountAdjFULL) = colnames(mkc[[pcpi]])
    }
    effort = foreach(i = 1:length(phenoParent)) %dopar% {
      effort = mlf[[pcpi]][,phenoParent_ind[i]] - mlf[[pcpi]][,phenoParent[[i]]]
      if (length(phenoParent[[i]])==1) effort = matrix(effort,ncol=1)
      rownames(effort) = rownames(mkc[[pcpi]])
      colnames(effort) = colnames(mkc[[pcpi]])[phenoParent[[i]]]
      return(effort)
    }
    names(effort) <- phenoMeta$phenotype[phenoParent_ind]
    save(effort, file=paste0(matrixParent_effort_dir,"_",pcp,".Rdata"))
    effortREDUCE = foreach(k=1:length(effort),.combine="cbind") %dopar% { return(effort[[k]]) }
    write.csv(effortREDUCE, file=paste0(matrixParent_effort_dir,"_",pcp,".csv"))
  } else {
    load(paste0(matrixParent_effort_dir,"_",pcp,".Rdata"))
    names(effort) <- phenoMeta$phenotype[phenoParent_ind]
    save(effort, file=paste0(matrixParent_effort_dir,"_",pcp,".Rdata"))
    effortREDUCE = foreach(k=1:length(effort),.combine="cbind") %dopar% { return(effort[[k]]) }
    write.csv(effortREDUCE, file=paste0(matrixParent_effort_dir,"_",pcp,".csv"))
  }
  
  #trim based on pvalue trim
  cat(" ",pheno_type[pcpi],sep="")
  effortTRIM = trimlistmatrix(mpt[[pcpi]],refdel=0, effort,targetdel=0, ncore=no_cores)
  save(effortTRIM, file=paste0(matrixParent_effort_dir,"TRIM_",pcp,".Rdata"))
  b = foreach(k=1:length(effortTRIM),.combine="cbind") %dopar% { return(effortTRIM[[k]]) }; rownames(b) = sampleMeta$gene[match(rownames(b),sampleMeta$fileName)]
  write.table(as.matrix(b), file=paste0(matrixParent_effort_dir,"TRIM_",pcp,".csv"), sep=",", col.names=F)
}

TimeOutput(start)



