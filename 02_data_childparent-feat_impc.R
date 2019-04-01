## Input: original normalized count matrix --> Output: abasolute feature + trimmed versions
# aya43@sfu.ca 20170131

## root directory
root = "~/projects/IMPC"
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#,"Sanger_MLN","CIPHE","TCP","H")
controlL = c("+_+|+_Y","+_+|+_Y","WildType","WildType","WildType") #control value in target_col column
ci = 1; panel = panelL[ci]; centre = centreL[ci]

result_dir = paste0("result/", panelL, "/", centreL); suppressWarnings(dir.create (result_dir))

## input directories
meta_dir = paste0(result_dir,"/meta")
meta_cell_dir = paste(meta_dir, "/cell", sep="")
meta_cell_child_dir = paste(meta_dir, "/cell_child",sep="") #specifies a phenotypes children
meta_cell_child_names_dir = paste(meta_dir, "/cell_child_names",sep="") #specifies a phenotypes children
meta_cell_child_ind_dir = paste(meta_dir, "/cell_child_ind",sep="")
meta_cell_childpn_dir = paste(meta_dir, "/cell_childpn",sep="") #specifies a phenotypes children and splits them into +/- (only for when both -/+ exists)
meta_cell_childpn_names_dir = paste(meta_dir, "/cell_childpn_names",sep="") #specifies a phenotypes children and splits them into +/- (only for when both -/+ exists)
meta_cell_childpn_ind_dir = paste(meta_dir, "/cell_childpn_ind",sep="")
meta_cell_parent_dir = paste(meta_dir, "/cell_parent",sep="") #specifies a phenotypes parents
meta_cell_parent_names_dir = paste(meta_dir, "/cell_parent_names",sep="") #specifies a phenotypes parents
meta_cell_parent_ind_dir = paste(meta_dir, "/cell_parent_ind",sep="")
meta_cell_parentpn_dir = paste(meta_dir, "/cell_parentpn",sep="") #specifies a phenotypes parents and splits them into +/- (only for when both -/+ exists)
meta_cell_parentpn_names_dir = paste(meta_dir, "/cell_parentpn_names",sep="")
meta_cell_parentpn_ind_dir = paste(meta_dir, "/cell_parentpn_ind",sep="")

feat_count = "file-cell-countAdj"

feat_dir = paste(result_dir, "/feat", sep=""); dir.create(feat_dir, showWarnings=F)
feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")
feat_file_cell_countAdjKO_dir = paste(feat_dir, "/file-cell-countAdjKO.",feat_count,".FULL",sep="")
feat_file_cell_logfold_dir = paste(feat_dir, "/file-cell-logfold.",feat_count,".FULL", sep="")

## output directories
feat_file_edge_pnratio_dir = paste(feat_dir, "/file-edge-pnratio",sep="")
feat_file_edge_prop_dir = paste(feat_dir, "/file-edge-prop",sep="")
feat_file_cell_entropychild_dir = paste(feat_dir, "/file-cell-entropychild",sep="")
feat_file_cell_entropyparent_dir = paste(feat_dir, "/file-cell-entropyparent",sep="")

feat_file_edge_contrib_dir = paste(feat_dir, "/file-cell-contrib.",feat_count,sep="")
feat_file_edge_effort_dir = paste(feat_dir, "/file-cell-effort.",feat_count,sep="")

## libraries
source("code/_funcAlice.R")
libr("stringr")
libr("entropy")
libr("foreach")
libr("doMC")





## cores
no_cores = detectCores() - 1
registerDoMC(no_cores)




## options
options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

create_child_entropy = T
create_parent_entropy = T



writecsv = F




start = Sys.time()








start1 = Sys.time()

#get list of children for each non-leaf node & save
cat("\ncreating child matrix")

m = get(load(paste0(feat_file_cell_count_dir,".Rdata")))
meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))

feat_file_cell_logfold = get(load(paste0(feat_file_cell_logfold_dir,".Rdata")))
feat_file_cell_countAdjKO = get(load(paste0(feat_file_cell_countAdjKO_dir,".Rdata")))




meta_cell_child = get(load(paste0(meta_cell_child_dir,".Rdata")))
meta_cell_child_names = get(load(paste0(meta_cell_child_names_dir,".Rdata")))
meta_cell_child_ind = get(load(paste0(meta_cell_child_ind_dir,".Rdata")))

meta_cell_childpn = get(load(paste0(meta_cell_childpn_dir, ".Rdata")))
meta_cell_childpn_names = get(load(paste0(meta_cell_childpn_names_dir, ".Rdata")))
meta_cell_childpn_ind = get(load(paste0(meta_cell_childpn_ind_dir, ".Rdata")))

meta_cell_parent = get(load(paste0(meta_cell_parent_dir, ".Rdata")))
meta_cell_parent_names = get(load(paste0(meta_cell_parent_names_dir, ".Rdata")))
meta_cell_parent_ind = get(load(paste0(meta_cell_parent_ind_dir, ".Rdata")))




## child proportion --------------------------------------------

start1 = Sys.time()

#mlist=list()
cat("; childprop")
childprop0 = foreach(i = 1:length(meta_cell_child)) %dopar% { #for each phenotype
  #for (i in 1:length(meta_cell_child)) {
  mpind = meta_cell_child_ind[i]
  parent = m[,mpind]
  
  # Child proportion matrix
  cnames = unlist(meta_cell_child_names[[i]])
  children = m[,cnames]
  children[which(children<1)] = 1
  parent0 = parent
  parent0[which(parent<1)] = 1
  childprop = exp(children/parent0)
  if (is.null(dim(childprop))) childprop = matrix(childprop,ncol=1)
  
  pnames = names(meta_cell_child)[i]
  colnames(childprop) = paste0(pnames, "_", cnames)
  
  return(childprop)
}
# for (i in 1:length(childprop)) {
#   colnames(childprop[[i]]) = paste0(names(childprop)[i], "_", colnames(childprop[[i]]))
# }
childprop = Reduce("cbind",childprop0)
save(childprop, file=paste0(feat_file_edge_prop_dir,".Rdata"))
if (writecsv) write.csv(childprop, file=paste0(feat_file_edge_prop_dir,".csv"))


TimeOutput(start1)




## child pn ratio + child entropy --------------------------------------------

start1 = Sys.time()


cat(", childratio + childentropy")

yesind = which(!sapply(meta_cell_childpn_names, is.null))
mlist = foreach(i = yesind) %dopar% { #for each phenotype
  mpind = meta_cell_childpn_ind[i]
  pnames = names(meta_cell_childpn)[i]
  parent = m[,mpind]
  
  # P/N ratio matrix
  pnratio = NULL
  
  cnegnames = meta_cell_childpn_names[[i]][[1]]
  pos = m[,cnegnames]-m[,mpind]
  neg = m[,cnegnames]
  pnratio = pos/neg ##get rid of 0, Inf
  pnratio[pos==0 & neg==0] = 1
  pnratio[pos==0 & neg!=0] = 1/neg[pos==0 & neg!=0]
  pnratio[pos!=0 & neg==0] = pos[pos!=0 & neg==0]
  pnratio = log(pnratio)
  
  if (is.null(dim(pnratio))) pnratio = matrix(pnratio,ncol=1)
  if (sum(parent==0)>0) pnratio[which(parent==0),] = rep(0,length(cnegnames))
  colnames(pnratio) = paste0(pnames, "_", cnegnames)
  
  
  # Entropy matrix
  en = rep(0,nrow(m))
  
  # parent = m[,meta_cell_childpn_ind[i]]
  cnames = unlist(meta_cell_child_names[[i]])
  children = m[,cnames]
  if (length(cnames)==1) children = matrix(children,ncol=1)
  
  no_child = length(cnegnames)
  non0parents = which(parent>0)
  en[non0parents] = sapply(non0parents, function(x) return(entropy(children[x,]/parent[x])/no_child)) #average entropy over # of markers added
  
  #mlist[[i]] = list(pnratio=pnratio, childprop=childprop)
  return(list(entropy=en, pnratio=pnratio)) #ratio = +child_prop / -child_prop
}


pnratio = lapply(mlist, function(x) x$pnratio)
en = lapply(mlist, function(x) x$en)

pnratio = Reduce("cbind",pnratio)
feat_file_cell_entropychild = Reduce('cbind',en)
colnames(feat_file_cell_entropychild) = meta_cell_childpn[yesind]

save(pnratio, file=paste0(feat_file_edge_pnratio_dir, ".Rdata"))
save(feat_file_cell_entropychild, file=paste0(feat_file_cell_entropychild_dir, ".Rdata"))
if (writecsv) write.csv(pnratio, file=paste0(feat_file_edge_pnratio_dir, ".csv"))
if (writecsv) write.csv(feat_file_cell_entropychild, file=paste0(feat_file_cell_entropychild_dir, ".csv"))




TimeOutput(start1)





## parent entropy --------------------------------------------

start1 = Sys.time()


cat(", parententropy")


feat_file_cell_entropyparent = foreach(i = 1:length(meta_cell_parent_ind), .combine='cbind') %dopar% { #for each phenotype
  # Entropy matrix
  en = rep(0,nrow(m))
  
  mpind = meta_cell_parent_ind[i]
  childr = m[,mpind]
  
  pnames = meta_cell_parent_names[[i]]
  parent = m[,colnames(m)%in%pnames]
  if (length(pnames)==1) parent = matrix(parent,ncol=1)
  
  no_parent = length(pnames)
  non0childrs = which(childr>0)
  en[non0childrs] = sapply(non0childrs, function(x) return(entropy(childr[x]/parent[x,])/no_parent)) #average entropy over # of markers added
  
  #mlist[[i]] = list(pnratio=pnratio, childprop=childprop)
  return(en) #ratio = +child_prop / -child_prop
}

colnames(feat_file_cell_entropyparent) = meta_cell$phenotype[meta_cell_parent_ind]
feat_file_cell_entropyparent = feat_file_cell_entropyparent[,!apply(feat_file_cell_entropyparent,2,function(x) all(x==0))]
save(feat_file_cell_entropyparent, file=paste0(feat_file_cell_entropyparent_dir, ".Rdata"))
if (writecsv) write.csv(feat_file_cell_entropyparent, file=paste0(feat_file_cell_entropyparent_dir, ".csv"))













## Parent contribution ---------------------------------------------------------------
cat(", parentcontribution")

feat_file_cell_countAdjWT = foreach(i=1:ncol(feat_file_cell_logfold),.combine="cbind") %dopar% { 
  return(feat_file_cell_countAdjKO[,i]/exp(feat_file_cell_logfold[,i])) 
}
colnames(feat_file_cell_countAdjWT) = colnames(feat_file_cell_countAdjKO)
contrib0 = foreach(i = 1:length(meta_cell_parent)) %dopar% { #for each phenotype
  cnames = names(meta_cell_parent)[i]
  cind = which(colnames(feat_file_cell_countAdjKO)==cnames)
  
  change = feat_file_cell_countAdjKO[,cind] - feat_file_cell_countAdjWT[,cind]
  
  pnames = meta_cell_parent_names[[i]]
  pind = which(colnames(feat_file_cell_countAdjKO)%in%pnames)
  
  koparents = feat_file_cell_countAdjKO[,pind]
  wtparents = feat_file_cell_countAdjWT[,pind]
  if (length(pnames)==1) {
    koparents = matrix(koparents,ncol=1)
    wtparents = matrix(wtparents,ncol=1)
  }
  
  changeparents = koparents - wtparents
  if (length(which(changeparents==0))>0) changeparents[changeparents==0] = 1
  
  contrib = change / changeparents
  
  rownames(contrib) = rownames(feat_file_cell_countAdjKO)
  colnames(contrib) = paste0(pnames,"_",cnames)
  return(contrib)
}
contrib = Reduce("cbind",contrib0)
save(contrib, file=paste0(feat_file_edge_contrib_dir,".Rdata"))
if (writecsv) write.csv(contrib, file=paste0(feat_file_edge_contrib_dir,".csv"))

TimeOutput(start)










## Parent effort -------------------------------------------------------
cat(", parenteffort")
start = Sys.time()


effort0 = foreach(i = 1:length(meta_cell_parent)) %dopar% {
  cnames = names(meta_cell_parent)[i]
  cind = which(colnames(feat_file_cell_logfold)==cnames)
  pnames = meta_cell_parent_names[[i]]
  effort1 = feat_file_cell_logfold[,cind] - feat_file_cell_logfold[,colnames(feat_file_cell_logfold)%in%pnames]
  if (length(meta_cell_parent[[i]])==1) effort1 = matrix(effort1,ncol=1)
  rownames(effort1) = rownames(feat_file_cell_logfold)
  colnames(effort1) = paste0(pnames,"_",cnames)
  return(effort1)
}

effort = Reduce("cbind",effort0)
save(effort, file=paste0(feat_file_edge_effort_dir,".Rdata"))
if (writecsv) write.csv(effort, file=paste0(feat_file_edge_effort_dir,".csv"))


TimeOutput(start)







TimeOutput(start1)