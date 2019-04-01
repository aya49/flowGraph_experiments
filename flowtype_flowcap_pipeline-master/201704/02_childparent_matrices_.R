# aya43@sfu.ca 20170131
# creates list of children for each non-leaf node and creates children/parent list[[matrices]] (-/+ are only for phenotypes where both -,+ data exists)

root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
matrixCountAdj_dir = paste(result_dir, "/matrixCountAdj", sep="")
phenoChild_dir = paste(result_dir, "/phenoChild.Rdata",sep="")
phenoChild_ind_dir = paste(result_dir, "/phenoChild_ind.Rdata",sep="")
phenoChildpn_dir = paste(result_dir, "/phenoChildpn.Rdata",sep="")
phenoChildpn_ind_dir = paste(result_dir, "/phenoChildpn_ind.Rdata",sep="")
phenoParent_dir = paste(result_dir, "/phenoParent.Rdata",sep="")
phenoParent_ind_dir = paste(result_dir, "/phenoParent_ind.Rdata",sep="")
phenoParentpn_dir = paste(result_dir, "/phenoParentpn.Rdata",sep="")
phenoParentpn_ind_dir = paste(result_dir, "/phenoParentpn_ind.Rdata",sep="")

createChildEntropy = T
createParentEntropy = T

#Output
matrixChild_pnratio_dir = paste(result_dir, "/matrixChild_pnratio.Rdata",sep="")
matrixChild_prop_dir = paste(result_dir, "/matrixChild_prop.Rdata",sep="")
matrixChild_entropy_dir = paste(result_dir, "/matrixChild_entropy.Rdata",sep="")
matrixParent_entropy_dir = paste(result_dir, "/matrixParent_entropy.Rdata",sep="")
matrixChild_pnratio_dircsv = paste(result_dir, "/matrixChild_pnratio.csv",sep="")
matrixChild_prop_dircsv = paste(result_dir, "/matrixChild_prop.csv",sep="")
matrixChild_entropy_dircsv = paste(result_dir, "/matrixChild_entropy.csv",sep="")
matrixParent_entropy_dircsv = paste(result_dir, "/matrixParent_entropy.csv",sep="")

libr(stringr)
libr(entropy)
libr(foreach)
libr(doMC)
source("code/_funcAlice.R")



start = Sys.time()

no_cores = detectCores() - 1
registerDoMC(no_cores)

start1 = Sys.time()

#get list of children for each non-leaf node & save
cat("\ncreating child matrix")

m = get(load(matrixCount_dir))
phenoMeta = get(load(phenoMeta_dir))
phenoChild = get(load(phenoChild_dir))
phenoChild_ind = get(load(phenoChild_ind_dir))


#mlist=list()
cat("; childprop")
childprop = foreach(i = 1:length(phenoChild)) %dopar% { #for each phenotype
  #for (i in 1:length(phenoChild)) {
  parent = m[,phenoChild_ind[i]]

  # Child proportion matrix
  children = m[,unlist(phenoChild[[i]])]
  children[which(children<1)] = 1
  parent0 = parent
  parent0[which(parent<1)] = 1
  childprop = exp(children/parent0)
  
  return(childprop)
}
names(childprop) = phenoMeta$phenotype[phenoChild_ind]
save(childprop, file=matrixChild_prop_dir)
write.csv(Reduce('cbind',childprop), file=matrixChild_prop_dircsv)


cat(", childratio + childentropy")
phenoChildpn = get(load(phenoChildpn_dir))
phenoChildpn_ind = get(load(phenoChildpn_ind_dir))

mlist = foreach(i = 1:length(phenoChildpn)) %dopar% { #for each phenotype
  parent = m[,phenoChildpn_ind[i]]
  
  # P/N ratio matrix
  pnratio = NULL
  
  pos = m[,phenoChildpn_ind[i]]-m[,phenoChildpn[[i]][[1]]]
  neg = m[,phenoChildpn[[i]][[1]]]
  pnratio = pos/neg ##get rid of 0, Inf
  pnratio[pos==0 & neg==0] = 1
  pnratio[pos==0 & neg!=0] = 1/neg[pos==0 & neg!=0]
  pnratio[pos!=0 & neg==0] = pos[pos!=0 & neg==0]
  pnratio = log(pnratio)
  
  if (sum(parent==0)>0) {
    if (is.null(dim(pnratio))) {
      pnratio = matrix(pnratio,ncol=1)
    }
    pnratio[which(parent==0),] = rep(0,length(phenoChildpn[[i]][[1]]))
  } 
  
  # Entropy matrix
  en = rep(0,nrow(m))
  
  parent = m[,phenoChildpn_ind[i]]
  children = m[,unlist(phenoChildpn[[i]])]
  if (length(phenoChildpn[[i]])==1) children = matrix(children,ncol=1)
  
  
  no_child = length(phenoChildpn[[i]][[1]])
  non0parents = which(parent>0)
  en[non0parents] = sapply(non0parents, function(x) return(entropy(children[x,]/parent[x])/no_child)) #average entropy over # of markers added
  
  #mlist[[i]] = list(pnratio=pnratio, childprop=childprop)
  return(list(entropy=en, pnratio=pnratio)) #ratio = +child_prop / -child_prop
}

# childprop = list()
# for(i in 1:length(mlist)) { childprop[[i]] = mlist[[i]] }

  
pnratio = list()
for(i in 1:length(mlist)) { pnratio[[i]] = mlist[[i]]$pnratio }

en = list()
for(i in 1:length(mlist)) { en[[i]] = mlist[[i]]$en }
matrixChild_entropy = Reduce('cbind',en)

rm(mlist)

names(pnratio) = colnames(matrixChild_entropy) = phenoMeta$phenotype[phenoChildpn_ind]
save(pnratio, file=matrixChild_pnratio_dir)
save(matrixChild_entropy, file=matrixChild_entropy_dir)
write.csv(Reduce('cbind',pnratio), file=matrixChild_pnratio_dircsv)
write.csv(matrixChild_entropy, file=matrixChild_entropy_dircsv)





cat(", parententropy")
matrixParent_entropy = foreach(i = 1:length(phenoParentpn), .combine='cbind') %dopar% { #for each phenotype
  # Entropy matrix
  en = rep(0,nrow(m))
  
  childr = m[,phenoParentpn_ind[i]]
  parent = m[,phenoParentpn[[i]]]
  if (length(phenoParentpn[[i]])==1) parent = matrix(parent,ncol=1)
  
  no_parent = length(phenoParentpn[[i]])
  non0childrs = which(childr>0)
  en[non0childrs] = sapply(non0childrs, function(x) return(entropy(childr[x]/parent[x,])/no_parent)) #average entropy over # of markers added
  
  #mlist[[i]] = list(pnratio=pnratio, childprop=childprop)
  return(en) #ratio = +child_prop / -child_prop
}

colnames(matrixParent_entropy) = phenoMeta$phenotype[phenoParentpn_ind]
save(matrixParent_entropy, file=matrixParent_entropy_dir)
write.csv(matrixParent_entropy, file=matrixParent_entropy_dircsv)







TimeOutput(start1)



TimeOutput(start)




