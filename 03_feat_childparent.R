## input: meta_cell & base features e.g. countadj, props (each cell population's counts as a percentage of the total count), p values (for if we want to delete rows/columns with no significant values)
## output: features
## - child_prop (proportion of child count over parent count)
## - child_pnratio (ratio of positive over negative child e.g. a+b+/a+b-)
## - child/parent_entropy (entropy of values of a cell population's children/parents)
## - parent_contrib (ratio of the difference between the avg control fcm file's and experiment fcm file's child population over that of the parent)
## - parent_effort (same as above except instead of difference, use log fold, and instead of ratio use different as values are in log scale)
## - lnpropexpect (product of parent props over product of grandparent props -- both normalized via an exponent 1/(number of parents/grandparents))

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr(c("stringr", "Matrix", "entropy", "plyr",
       "foreach", "doMC"))


## cores
no_cores = detectCores() - 3
registerDoMC(no_cores)


## options
options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

create_child_entropy = T
create_parent_entropy = T

writecsv = F

feat_count = "file-cell-countAdj"

start = Sys.time()
for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)[-16]) {
  # result_dir = paste0(root, "/result/flowcap_panel6") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta")
  meta_cell_dir = paste(meta_dir, "/cell", sep="")
  # meta_cell_child_dir = paste(meta_dir, "/cell_child",sep="") #specifies a phenotypes children
  # meta_cell_child_names_dir = paste(meta_dir, "/cell_child_names",sep="") #specifies a phenotypes children
  # meta_cell_child_ind_dir = paste(meta_dir, "/cell_child_ind",sep="")
  # meta_cell_childpn_dir = paste(meta_dir, "/cell_childpn",sep="") #specifies a phenotypes children and splits them into +/- (only for when both -/+ exists)
  meta_cell_childpn_names_dir = paste(meta_dir, "/cell_childpn_names",sep="") #specifies a phenotypes children and splits them into +/- (only for when both -/+ exists)
  # meta_cell_childpn_ind_dir = paste(meta_dir, "/cell_childpn_ind",sep="")
  # meta_cell_parent_dir = paste(meta_dir, "/cell_parent",sep="") #specifies a phenotypes parents
  meta_cell_parent_names_dir = paste(meta_dir, "/cell_parent_names",sep="") #specifies a phenotypes parents
  # meta_cell_parent_ind_dir = paste(meta_dir, "/cell_parent_ind",sep="")
  # meta_cell_parentpn_dir = paste(meta_dir, "/cell_parentpn",sep="") #specifies a phenotypes parents and splits them into +/- (only for when both -/+ exists)
  # meta_cell_parentpn_names_dir = paste(meta_dir, "/cell_parentpn_names",sep="")
  # meta_cell_parentpn_ind_dir = paste(meta_dir, "/cell_parentpn_ind",sep="")
  
  
  feat_dir = paste(result_dir, "/feat", sep=""); dir.create(feat_dir, showWarnings=F)
  feat_file_cell_countAdj_dir = paste(feat_dir, "/file-cell-countAdj", sep="")
  feat_file_cell_prop_dir = paste(feat_dir, "/file-cell-prop", sep="")
  feat_file_cell_countAdjKO_dir = paste(feat_dir, "/file-cell-countAdjKOFULL.",feat_count,"",sep="")
  feat_file_cell_logfold_dir = paste(feat_dir, "/file-cell-logfoldFULL.",feat_count, sep="")
  
  ## output directories
  feat_file_edge_pnratio_dir = paste(feat_dir, "/file-edge-pnratio",sep="")
  feat_file_edge_prop_dir = paste(feat_dir, "/file-edge-prop",sep="")
  feat_file_cell_entropychild_dir = paste(feat_dir, "/file-cell-entropychild",sep="")
  feat_file_cell_entropyparent_dir = paste(feat_dir, "/file-cell-entropyparent",sep="")
  feat_file_cell_lnpropexpect_dir = paste(feat_dir, "/file-cell-lnpropexpect",sep="")
  
  feat_file_edge_contrib_dir = paste(feat_dir, "/file-edge-contrib-",feat_count,sep="")
  feat_file_edge_effort_dir = paste(feat_dir, "/file-edge-effort-",feat_count,sep="")
  
  
  
  start1 = Sys.time()
  
  #get list of children for each non-leaf node & save
  cat("\ncreating child matrix")
  
  m = Matrix(get(load(paste0(feat_file_cell_countAdj_dir,".Rdata"))))
  
  meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
  mp = Matrix(get(load(paste0(feat_file_cell_prop_dir,".Rdata"))))
  # meta_cell_child = get(load(paste0(meta_cell_child_dir,".Rdata")))
  # meta_cell_child_names = get(load(paste0(meta_cell_child_names_dir,".Rdata")))
  # meta_cell_child_ind = get(load(paste0(meta_cell_child_ind_dir,".Rdata")))
  
  # meta_cell_childpn = get(load(paste0(meta_cell_childpn_dir, ".Rdata")))
  meta_cell_childpn_names = get(load(paste0(meta_cell_childpn_names_dir, ".Rdata")))
  # meta_cell_childpn_ind = get(load(paste0(meta_cell_childpn_ind_dir, ".Rdata")))
  
  # meta_cell_parent = get(load(paste0(meta_cell_parent_dir, ".Rdata")))
  meta_cell_parent_names = get(load(paste0(meta_cell_parent_names_dir, ".Rdata")))
  # meta_cell_parent_ind = get(load(paste0(meta_cell_parent_ind_dir, ".Rdata")))
  
  
  
  ## child proportion --------------------------------------------
  
  start1 = Sys.time()
  
  #mlist=list()
  cat("; childprop")
  childprop0 = foreach(i = 1:length(meta_cell_childpn_names)) %dopar% { #for each phenotype
    #for (i in 1:length(meta_cell_childpn_names)) {
    pnames = names(meta_cell_childpn_names)[i]
    parent = m[,colnames(m)%in%pnames]
    cnames = unlist(meta_cell_childpn_names[[i]])
    childr = m[,cnames,drosp=F]
    childr[childr<1] = parent[parent<1] = 1
    
    # tryCatch({
    childprop = exp(childr/parent)
    # }, error = function(err) { 
    #   parent = unlist(parent)
    #   childprop = childr
    #   for (j in 1:ncol(childr)) childprop[,j] = exp(unlist(childr[,j])/parent)
    # })
    colnames(childprop) = paste0(pnames, "_", cnames)
    return(childprop)
  }
  # for (i in 1:length(childprop)) {
  #   colnames(childprop[[i]]) = paste0(names(childprop)[i], "_", colnames(childprop[[i]]))
  # }
  childprop = Reduce("cbind",childprop0)
  save(childprop, file=paste0(feat_file_edge_prop_dir,".Rdata"))
  if (writecsv) write.csv(childprop, file=paste0(feat_file_edge_prop_dir,".csv"))
  
  time_output(start1)
  
  
  
  ## child pn ratio --------------------------------------------
  
  start1 = Sys.time()
  cat(", childratio")
  
  pnratio0 = foreach(i = 1:length(meta_cell_childpn_names)) %dopar% { #for each phenotype
    pnames = names(meta_cell_childpn_names)[i]
    parent = m[,colnames(m)%in%pnames]
    cnamep = unlist(meta_cell_childpn_names[[i]]$pos)
    cnamen = unlist(meta_cell_childpn_names[[i]]$neg)
    cnamep = cnamep[1:min(length(cnamep),length(cnamen))]
    cnamen = cnamen[1:length(cnamep)]
    childp = m[,cnamep,drop=F]
    childn = m[,cnamen,drop=F]
    
    # P/N ratio matrix
    pnratio = childp/childn #get rid of 0, Inf
    pnratio[childp==0 & childn==0] = 1
    pnratio[childp==0 & childn!=0] = 1 / childn[childp==0 & childn!=0]
    pnratio[childp!=0 & childn==0] = childp[childp!=0 & childn==0]
    pnratio = log(pnratio)
    
    if (is.null(dim(pnratio))) pnratio = matrix(pnratio,ncol=1) #safeguard... probs don't need it but, just want a matrix output not vector
    if (sum(parent==0)>0) pnratio[which(parent==0),] = rep(0,length(cnamen))
    colnames(pnratio) = paste0(pnames, "_", cnamen)
    
    #mlist[[i]] = list(pnratio=pnratio, childprop=childprop)
    return(pnratio) #ratio = +child_prop / -child_prop
  }
  
  pnratio = Reduce("cbind",pnratio0)
  
  save(pnratio, file=paste0(feat_file_edge_pnratio_dir, ".Rdata"))
  if (writecsv) write.csv(pnratio, file=paste0(feat_file_edge_pnratio_dir, ".csv"))
  
  time_output(start1)
  
  
  
  ## child entropy --------------------------------------------
  
  start1 = Sys.time()
  cat(", childentropy")
  
  meta_cell_childpn_names_ = meta_cell_childpn_names
  for (i in 1:length(meta_cell_childpn_names)) {
    if (length(meta_cell_childpn_names[[i]])<2) meta_cell_childpn_names_[[i]] = NULL
  }
  feat_file_cell_entropychild = foreach(i = 1:length(meta_cell_childpn_names_), .combine="cbind") %dopar% { #for each phenotype
    if (length(meta_cell_childpn_names_[[i]])<2) return(NULL)
    pnames = names(meta_cell_childpn_names_)[i]
    parent = m[,colnames(m)%in%pnames]
    cnames = unlist(meta_cell_childpn_names_[[i]])
    childr = m[,cnames,drop=F]
    # childr[childr<1] = parent[parent<1] = 1
    
    # Entropy matrix
    en = rep(0,nrow(m))
    no_child = length(cnames)
    non0parents = parent>0
    if (sum(non0parents>0)>0) en[non0parents] = sapply(which(non0parents), function(x) entropy(childr[x,]/parent[x])/no_child) #average entropy over # of markers added
    
    #mlist[[i]] = list(pnratio=pnratio, childprop=childprop)
    return(en) #ratio = +child_prop / -child_prop
  }
  rownames(feat_file_cell_entropychild) = rownames(m)
  colnames(feat_file_cell_entropychild) = names(meta_cell_childpn_names_)
  
  save(feat_file_cell_entropychild, file=paste0(feat_file_cell_entropychild_dir, ".Rdata"))
  if (writecsv) write.csv(feat_file_cell_entropychild, file=paste0(feat_file_cell_entropychild_dir, ".csv"))
  
  time_output(start1)
  
  
  
  ## parent entropy --------------------------------------------
  
  start1 = Sys.time()
  cat(", parententropy")
  
  meta_cell_parent_names_ = meta_cell_parent_names
  for (i in 1:length(meta_cell_parent_names))
    if (length(meta_cell_parent_names[[i]])<2) meta_cell_parent_names_[[i]] = NULL
  feat_file_cell_entropyparent = foreach(i = 1:length(meta_cell_parent_names_), .combine='cbind') %dopar% { #for each phenotype
    # Entropy matrix
    en = rep(0,nrow(m))
    
    pnames = meta_cell_parent_names_[[i]]
    parent = m[,colnames(m)%in%pnames,drop=F]
    cnames = names(meta_cell_parent_names_)[i]
    childr = m[,cnames]
    
    non0childrs = childr>0
    if (sum(non0childrs)>0) en[non0childrs] = sapply(which(non0childrs), function(x) return(entropy(childr[x]/parent[x,])/length(pnames))) #average entropy over # of markers added
    
    #mlist[[i]] = list(pnratio=pnratio, childprop=childprop)
    return(en) #ratio = +child_prop / -child_prop
  }
  colnames(feat_file_cell_entropyparent) = names(meta_cell_parent_names_)
  rownames(feat_file_cell_entropyparent) = rownames(m)
  feat_file_cell_entropyparent = feat_file_cell_entropyparent[,!apply(feat_file_cell_entropyparent,2,function(x) all(x==0))]
  save(feat_file_cell_entropyparent, file=paste0(feat_file_cell_entropyparent_dir, ".Rdata"))
  if (writecsv) write.csv(feat_file_cell_entropyparent, file=paste0(feat_file_cell_entropyparent_dir, ".csv"))
  
  time_output(start1)
  
  ## prop/expected -----------------------------------------------------
  
  start1 = Sys.time()
  cat(", ln prop/expected")
  
  cellis = meta_cell$phenotype[meta_cell$phenolevel>2]
  cellis = cellis[cellis%in%names(meta_cell_parent_names)]
  mpe = mp
  mpe[mpe==0] = min(mp[mp!=0])
  lnpropexpect = foreach(i=cellis, .combine="cbind") %dopar% {
    pnames = meta_cell_parent_names[[i]]
    parent = mpe[,colnames(mpe)%in%pnames,drop=F]
    numrtr = apply(parent, 1, prod)^ncol(parent)
    
    gnames = unique(unlist(meta_cell_parent_names[pnames]))
    grprnt = mpe[,colnames(mpe)%in%gnames,drop=F]
    denmtr = apply(grprnt, 1, prod)
    
    expect = (numrtr/denmtr)^(1/ncol(grprnt))
    childr = mpe[,i]
    lnpropexpecti = log(childr/expect)
    
    return(lnpropexpecti)
  }
  rownames(lnpropexpect) = rownames(mpe)
  colnames(lnpropexpect) = cellis
  save(lnpropexpect, file=paste0(feat_file_cell_lnpropexpect_dir,".Rdata"))
  if (writecsv) write.csv(lnpropexpect, file=paste0(feat_file_cell_lnpropexpect_dir,".csv"))
  
  time_output(start1)
  
  
  
  
  ## features that require p value script -----------------------
  if (file.exists(paste0(feat_file_cell_logfold_dir,".Rdata"))) {
    feat_file_cell_logfold = Matrix(get(load(paste0(feat_file_cell_logfold_dir,".Rdata"))))
    feat_file_cell_countAdjKO = Matrix(get(load(paste0(feat_file_cell_countAdjKO_dir,".Rdata"))))
    
    
    ## Parent contribution ---------------------------------------------------------------
    
    start1 = Sys.time()
    cat(", parentcontribution")
    
    feat_file_cell_countAdjWT = foreach(i=1:ncol(feat_file_cell_logfold),.combine="cbind") %dopar% { 
      return(feat_file_cell_countAdjKO[,i]/exp(feat_file_cell_logfold[,i])) 
    }
    colnames(feat_file_cell_countAdjWT) = colnames(feat_file_cell_countAdjKO)
    contrib0 = llply(1:length(meta_cell_parent_names), function(i) { #for each phenotype
      cnames = names(meta_cell_parent_names)[i]
      change = feat_file_cell_countAdjKO[,colnames(feat_file_cell_countAdjKO)%in%cnames] - feat_file_cell_countAdjWT[,colnames(feat_file_cell_countAdjWT)%in%cnames]
      
      pnames = meta_cell_parent_names[[i]]
      koparents = feat_file_cell_countAdjKO[,colnames(feat_file_cell_countAdjKO)%in%pnames, drop=F]
      wtparents = feat_file_cell_countAdjWT[,colnames(feat_file_cell_countAdjWT)%in%pnames, drop=F]
      
      changeparents = koparents - wtparents
      changeparents[changeparents<1] = 1
      
      contrib = change / changeparents
      
      rownames(contrib) = rownames(feat_file_cell_countAdjKO)
      colnames(contrib) = paste0(pnames,"_",cnames)
      return(contrib)
    }, .parallel=T)
    contrib = Reduce("cbind",contrib0)
    save(contrib, file=paste0(feat_file_edge_contrib_dir,".Rdata"))
    if (writecsv) write.csv(contrib, file=paste0(feat_file_edge_contrib_dir,".csv"))
    
    time_output(start1)
    
    
    
    ## Parent effort -------------------------------------------------------
    
    start1 = Sys.time()
    cat(", parenteffort")
    
    effort0 = llply(1:length(meta_cell_parent_names), function(i) {
      pnames = meta_cell_parent_names[[i]]
      parent = feat_file_cell_logfold[,colnames(m)%in%pnames,drop=F] 
      cnames = names(meta_cell_parent_names)[i]
      childr = feat_file_cell_logfold[,cnames]
      
      effort1 = childr - parent
      if (length(pnames)==1) effort1 = matrix(effort1,ncol=1) #prob don't need, was a failsafe
      rownames(effort1) = rownames(feat_file_cell_logfold)
      colnames(effort1) = paste0(pnames,"_",cnames)
      return(effort1)
    }, .parallel=T)
    
    effort = Reduce("cbind",effort0)
    save(effort, file=paste0(feat_file_edge_effort_dir,".Rdata"))
    if (writecsv) write.csv(effort, file=paste0(feat_file_edge_effort_dir,".csv"))
    
    time_output(start1)
    
  }
  
}
time_output(start)

