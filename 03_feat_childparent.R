## input: meta_cell & base features e.g. countadj, props (each cell population's counts as a percentage of the total count), p values (for if we want to delete rows/columns with no significant values)
## output: features
## - child_prop (proportion of child count over parent count)
## - child_pnratio (ratio of positive over negative child e.g. a+b+/a+b-)
## - child/parent_entropy (entropy of values of a cell population's children/parents)
## - group_entropy (entropy of duplicate marker groups after removing -+ e.g. a+b+ a+b- a-b+ a-b-)
## - lnpropexpect (product of parent props over product of grandparent props -- both normalized via an exponent 1/(number of parents/grandparents))
## - lnpropexpectshort (lnpropexpect but without any - markers, so all + marker combinations only)
## - 

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr(c("stringr", "Matrix", "entropy", "plyr",
       "foreach", "doMC"))


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

feat_count = "file-cell-countAdj"

start = Sys.time()
for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  print(result_dir)
  # result_dir = paste0(root, "/result/flowcap_panel6") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta")
  meta_cell_dir = paste(meta_dir, "/cell", sep="")
  meta_cell_childpn_names_dir = paste(meta_dir, "/cell_childpn_names",sep="") #specifies a phenotypes children and splits them into +/- (only for when both -/+ exists)
  meta_cell_parent_names_dir = paste(meta_dir, "/cell_parent_names",sep="") #specifies a phenotypes parents
  
  
  feat_dir = paste(result_dir, "/feat", sep=""); dir.create(feat_dir, showWarnings=F)
  feat_file_cell_countAdj_dir = paste(feat_dir, "/file-cell-countAdj", sep="")
  feat_file_cell_prop_dir = paste(feat_dir, "/file-cell-prop", sep="")
  
  ## output directories
  feat_file_edge_pnratio_dir = paste(feat_dir, "/file-edge-pnratio",sep="")
  feat_file_edge_prop_dir = paste(feat_dir, "/file-edge-prop",sep="")
  feat_file_cell_entropychild_dir = paste(feat_dir, "/file-cell-entropychild",sep="")
  feat_file_cell_entropyparent_dir = paste(feat_dir, "/file-cell-entropyparent",sep="")
  feat_file_group_entropy_dir = paste(feat_dir, "/file-group-entropy",sep="")
  feat_file_group_var_dir = paste(feat_dir, "/file-group-var",sep="")
  
  feat_file_cell_lnpropexpect_dir = paste(feat_dir, "/file-cell-lnpropexpect",sep="")
  feat_file_cell_lnpropexpectshort_dir = paste(feat_dir, "/file-cell-lnpropexpectshort",sep="")
  
  
  m = as.matrix(Matrix(get(load(paste0(feat_file_cell_countAdj_dir,".Rdata")))))
  mp = foreach(xi=1:ncol(m), .combine='cbind') %dopar% { return(m[,xi]/m[,1]) }
  dimnames(mp) = dimnames(m)
  save(mp, file=paste0(feat_file_cell_prop_dir,".Rdata"))
  
  meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
  
  meta_cell_childpn_names = get(load(paste0(meta_cell_childpn_names_dir, ".Rdata")))
  meta_cell_parent_names = get(load(paste0(meta_cell_parent_names_dir, ".Rdata")))
  
  
  start = Sys.time()
  
  
  ## child proportion --------------------------------------------
  
  start1 = Sys.time()
  
  #mlist=list()
  cat("childprop")
  childprop = foreach(ii = loopInd(1:length(meta_cell_childpn_names),no_cores), .combine="cbind") %dopar% { #for each phenotype
    childprop1 = foreach(i=ii, .combine="cbind") %do% {
      #for (i in 1:length(meta_cell_childpn_names)) {
      pnames = names(meta_cell_childpn_names)[i]
      pnames = pnames[pnames%in%colnames(m)]
      parent = m[,match(pnames,colnames(m)),drop=F]
      cnames = unlist(meta_cell_childpn_names[[i]])
      cnames = cnames[cnames%in%colnames(m)]
      childr = m[,match(cnames,colnames(m)),drop=F]
      if (ncol(parent)==0 | ncol(childr)==0) return(NULL)
      childr[childr<1] = parent[parent<1] = 1
      
      # tryCatch({
      childprop2 = exp(childr/as.vector(parent))
      # }, error = function(err) { 
      #   parent = unlist(parent)
      #   childprop = childr
      #   for (j in 1:ncol(childr)) childprop[,j] = exp(unlist(childr[,j])/parent)
      # })
      colnames(childprop2) = paste0(pnames, "_", cnames)
      return(childprop2)
    }
    return(childprop1)
  }
  # for (i in 1:length(childprop)) {
  #   colnames(childprop[[i]]) = paste0(names(childprop)[i], "_", colnames(childprop[[i]]))
  # }
  save(childprop, file=paste0(feat_file_edge_prop_dir,".Rdata"))
  if (writecsv) write.csv(childprop, file=paste0(feat_file_edge_prop_dir,".csv"))
  
  time_output(start1)
  
  
  
  ## child pn ratio --------------------------------------------
  
  start1 = Sys.time()
  cat("childratio")
  
  pnratio = foreach(ii = loopInd(1:length(meta_cell_childpn_names),no_cores), .combine="cbind") %dopar% { #for each phenotype
    pnratio1 = foreach(i=ii, .combine="cbind") %do% {
      pnames = names(meta_cell_childpn_names)[i]
      pnames = pnames[pnames%in%colnames(m)]
      cnamen = unlist(meta_cell_childpn_names[[i]]$neg)
      cnamen = cnamen[cnamen%in%colnames(m)]
      cnamep = unlist(meta_cell_childpn_names[[i]]$pos)
      cnamep = cnamep[cnamep%in%colnames(m)]
      if (length(cnamen)==0 | length(cnamep)==0) return(NULL)
      
      cnamen = cnamen[match(gsub("[+-]","",cnamep),gsub("[+-]","",cnamen))]
      cnamen = cnamen[!is.na(cnamen)]
      cnamep = cnamep[match(gsub("[+-]","",cnamen),gsub("[+-]","",cnamep))]
      cnamep = cnamep[!is.na(cnamep)]
      
      
      parent = m[,match(pnames,colnames(m)),drop=F]
      childp = m[,match(cnamep,colnames(m)),drop=F]
      childn = m[,match(cnamen,colnames(m)),drop=F]
      if (ncol(parent)==0 | ncol(childp)==0 | ncol(childn)==0) return(NULL)
      
      # P/N ratio matrix
      pnratio2 = childp/childn #get rid of 0, Inf
      pnratio2[childp==0 & childn==0] = 1
      pnratio2[childp==0 & childn!=0] = 1 / childn[childp==0 & childn!=0]
      pnratio2[childp!=0 & childn==0] = childp[childp!=0 & childn==0]
      pnratio2 = log(pnratio2)
      
      if (is.null(dim(pnratio2))) pnratio2 = matrix(pnratio2,ncol=1) #safeguard... probs don't need it but, just want a matrix output not vector
      if (sum(parent==0)>0) pnratio2[which(parent==0),] = rep(0,length(cnamen))
      colnames(pnratio2) = paste0(pnames, "_", cnamen)
      
      #mlist[[i]] = list(pnratio=pnratio, childprop=childprop)
      return(pnratio2) #ratio = +child_prop / -child_prop
    }
    return(pnratio1)
  }
  
  save(pnratio, file=paste0(feat_file_edge_pnratio_dir, ".Rdata"))
  if (writecsv) write.csv(pnratio, file=paste0(feat_file_edge_pnratio_dir, ".csv"))
  
  time_output(start1)
  
  
  
  ## child entropy --------------------------------------------
  
  start1 = Sys.time()
  cat("childentropy")
  
  meta_cell_childpn_names_ = meta_cell_childpn_names
  for (i in 1:length(meta_cell_childpn_names)) {
    if (length(meta_cell_childpn_names[[i]])<2) meta_cell_childpn_names_[[i]] = NULL
  }
  feat_file_cell_entropychild = a = unlist(llply(loopInd(1:length(meta_cell_childpn_names_),no_cores), function(ii) { #for each phenotype
    llply(ii, function(i) {
      if (length(meta_cell_childpn_names_[[i]])<2) return(NULL)
      # i = which(names(meta_cell_childpn_names)==ii)
      pnames = names(meta_cell_childpn_names)[i]
      pnames = pnames[pnames%in%colnames(m)]
      parent = m[,match(pnames,colnames(m)),drop=F]
      cnames = unlist(meta_cell_childpn_names[[i]])
      cnames = cnames[cnames%in%colnames(m)]
      childr = m[,match(cnames,colnames(m)),drop=F]
      # childr[childr<1] = parent[parent<1] = 1
      
      # Entropy matrix
      en = rep(0,nrow(m))
      if (ncol(parent)==0 | ncol(childr)==0) return(en)
      
      no_child = length(cnames)
      non0parents = parent>0
      if (sum(non0parents>0)>0) en[non0parents] = sapply(which(non0parents), function(x) entropy(childr[x,]/parent[x])/no_child) #average entropy over # of markers added
      
      #mlist[[i]] = list(pnratio=pnratio, childprop=childprop)
      return(en) #ratio = +child_prop / -child_prop
    }, .parallel=F)
  }, .parallel=T), recursive=F)
  feat_file_cell_entropychild = Reduce("cbind",a)
  colnames(feat_file_cell_entropychild) = names(meta_cell_childpn_names_)[sapply(a, function(x) !is.null(x))]
  rownames(feat_file_cell_entropychild) = rownames(m)
  feat_file_cell_entropychild = feat_file_cell_entropychild[,apply(feat_file_cell_entropychild,2,function(x) any(x!=0))]
  
  save(feat_file_cell_entropychild, file=paste0(feat_file_cell_entropychild_dir, ".Rdata"))
  if (writecsv) write.csv(feat_file_cell_entropychild, file=paste0(feat_file_cell_entropychild_dir, ".csv"))
  
  time_output(start1)
  
  
  
  ## parent entropy --------------------------------------------
  
  start1 = Sys.time()
  cat("parententropy")
  
  meta_cell_parent_names_ = meta_cell_parent_names
  for (i in names(meta_cell_parent_names))
    if (length(meta_cell_parent_names[[i]])<2) meta_cell_parent_names_[[i]] = NULL
  feat_file_cell_entropyparent = foreach(ii = loopInd(1:length(meta_cell_parent_names_), no_cores), .combine='cbind') %dopar% { #for each phenotype
    # Entropy matrix
    a1 = foreach(i=ii, .combine='cbind') %do% { #for each phenotype
      
      pnames = meta_cell_parent_names_[[i]]
      pnames = pnames[pnames%in%colnames(m)]
      parent = m[,pnames,drop=F]
      cnames = names(meta_cell_parent_names_)[i]
      cnames = cnames[cnames%in%colnames(m)]
      childr = m[,match(cnames,colnames(m)),drop=F]
      
      en = rep(0,nrow(m))
      if (ncol(parent)==0 | ncol(childr)==0) return(en)
      childr[childr<1] = parent[parent<1] = 1
      
      non0childrs = childr>0
      if (sum(non0childrs)>0) en[non0childrs] = sapply(which(non0childrs), function(x) return(entropy(childr[x]/parent[x,])/length(pnames))) #average entropy over # of markers added
      
      #mlist[[i]] = list(pnratio=pnratio, childprop=childprop)
      return(en) #ratio = +child_prop / -child_prop
    }
    return(a1)
  }
  colnames(feat_file_cell_entropyparent) = names(meta_cell_parent_names_)
  rownames(feat_file_cell_entropyparent) = rownames(m)
  feat_file_cell_entropyparent = feat_file_cell_entropyparent[,!apply(feat_file_cell_entropyparent,2,function(x) all(x==0))]
  save(feat_file_cell_entropyparent, file=paste0(feat_file_cell_entropyparent_dir, ".Rdata"))
  if (writecsv) write.csv(feat_file_cell_entropyparent, file=paste0(feat_file_cell_entropyparent_dir, ".csv"))
  
  time_output(start1)
  
  ## group entropy & variance --------------------------------------------
  
  start1 = Sys.time()
  cat("groupentropy")
  
  allplus = which(!grepl("[-]",colnames(m)) & grepl("[+]",colnames(m)))
  nonp = gsub("[-|+]","",colnames(m))
  nonpu = unique(nonp)
  groupi = match(nonp, nonp)
  groups = llply(allplus, function(i) {
    ii = which(groupi==groupi[i])
    if (length(ii)>1) return(ii)
    return(NULL)
  })
  names(groups) = colnames(m)[allplus]
  groups = plyr::compact(groups)
  
  a = llply(loopInd(1:length(groups),no_cores), function(gi) 
    sapply(groups[gi],function(i) apply(m[,i],1,entropy)), .parallel=T)
  feat_file_group_entropy = Reduce("cbind",a)
  b = llply(loopInd(1:length(groups),no_cores), function(gi) 
    sapply(groups[gi],function(i) apply(m[,i],1,var)), .parallel=T)
  feat_file_group_var = Reduce("cbind",b)
  e = llply(loopInd(1:length(groups),no_cores), function(gi) 
    sapply(groups[gi],function(i) apply(mp[,i],1,entropy)),.parallel=T)
  feat_file_group_entropyp = Reduce("cbind",e)
  f = llply(loopInd(1:length(groups),no_cores), function(gi) 
    sapply(groups[gi],function(i) apply(mp[,i],1,var)), .parallel=T)
  feat_file_group_varp = Reduce("cbind",f)
  
  colnames(feat_file_group_entropy) = colnames(feat_file_group_var) = colnames(feat_file_group_entropyp) = colnames(feat_file_group_varp) =  names(groups)
  rownames(feat_file_group_entropy) = rownames(feat_file_group_var) = rownames(feat_file_group_entropyp) = rownames(feat_file_group_varp) = rownames(m)
  
  feat_file_group_entropy = feat_file_group_entropy[,apply(feat_file_group_entropy,2,function(x) !all(x==x[1]))]
  feat_file_group_var = feat_file_group_var[,apply(feat_file_group_var,2,function(x) !all(x==x[1]))]
  feat_file_group_entropyp = feat_file_group_entropyp[,apply(feat_file_group_entropyp,2,function(x) !all(x==x[1]))]
  feat_file_group_varp = feat_file_group_varp[,apply(feat_file_group_varp,2,function(x) !all(x==x[1]))]
  
  save(feat_file_group_entropy, file=paste0(feat_file_group_entropy_dir, ".Rdata"))
  if (writecsv) write.csv(feat_file_group_entropy, file=paste0(feat_file_group_entropy_dir, ".csv"))
  save(feat_file_group_var, file=paste0(feat_file_group_var_dir, ".Rdata"))
  if (writecsv) write.csv(feat_file_group_var, file=paste0(feat_file_group_var_dir, ".csv"))
  save(feat_file_group_entropyp, file=paste0(feat_file_group_entropy_dir, "p.Rdata"))
  if (writecsv) write.csv(feat_file_group_entropyp, file=paste0(feat_file_group_entropy_dir, "p.csv"))
  save(feat_file_group_varp, file=paste0(feat_file_group_var_dir, "p.Rdata"))
  if (writecsv) write.csv(feat_file_group_varp, file=paste0(feat_file_group_var_dir, "p.csv"))
  
  time_output(start1)
  
  
  
  ## prop/expected -----------------------------------------------------
  
  start1 = Sys.time()
  cat("ln prop/expected")
  
  markern = str_length(meta_cell[1,2])
  
  cellis1 = meta_cell$phenotype[meta_cell$phenolevel==1]
  cellis1 = append("",cellis1[cellis1%in%names(meta_cell_parent_names)])
  expe1 = matrix(.5,nrow=nrow(mp), ncol=length(cellis1), dimnames=list(rownames(mp),cellis1))
  expe1[,1] = 1
  
  cellis = meta_cell$phenotype[meta_cell$phenolevel>1]
  cellis = cellis[cellis%in%names(meta_cell_parent_names)]
  cellisn = str_count(cellis,"[+-]")
  
  mpe = mp
  mpe[mpe==0] = min(mp[mp!=0])
  mpe = mpe[,match(append(cellis1,cellis),colnames(mpe))]
  
  meta_cell_parent_names_ = meta_cell_parent_names[names(meta_cell_parent_names)%in%colnames(mpe)]
  meta_cell_parent_names_[cellis] = llply(meta_cell_parent_names_[cellis], function(x) {
    a = x[x%in%colnames(mpe)]
    if (length(a)==0) return(NULL)
    return(a)
  })
  meta_cell_parent_names_ = Filter(Negate(is.null), meta_cell_parent_names_)
  
  expec = llply(loopInd(1:length(cellis),no_cores), function(ii) {
    # expect1 = foreach(ic=ii, .combine="cbind") %do% {
    #   i = cellis[ic]
    #   il = cellisn[ic]
    #   pnames = meta_cell_parent_names_[[i]]
    #   parent = mpe[,pnames,drop=F]
    #   numrtr = apply(parent, 1, prod)^(il-1)
    #   
    #   gnames = unique(unlist(meta_cell_parent_names_[pnames]))
    #   if (gnames=="") gnames = colnames(mpe)==""
    #   grprnt = mpe[,gnames,drop=F]
    #   denmtr = apply(grprnt, 1, prod)
    #   
    #   expect = (numrtr/denmtr)^(1/choose(il,2))
    #   return(expect)
    # }
    expect2 = foreach(ic=ii, .combine="cbind") %do% {
      i = cellis[ic]
      il = cellisn[ic]
      pnames = meta_cell_parent_names_[[i]]
      parent = mpe[,pnames,drop=F]
      numrtr = apply(parent, 1, function(x) prod(x^(il-1)))
      
      gnames = unique(unlist(meta_cell_parent_names_[pnames]))
      if (gnames=="") gnames = colnames(mpe)==""
      grprnt = mpe[,gnames,drop=F]
      denmtr = apply(grprnt, 1, prod)
      
      expect = (numrtr/denmtr)^(1/choose(il,2))
      return(expect)
    }
    return(expect2)
  }, .parallel=T)
  expec = Reduce("cbind", expec)
  colnames(expec) = cellis
  exp1 = cbind(expe1,expec)
  lnpropexpect1 = log(mpe/exp1)
  lnpropexpect1[exp1==0] = log(mpe[exp1==0])

  save(exp1, file=paste0(feat_file_cell_lnpropexpect_dir,"_.Rdata"))
  save(lnpropexpect1, file=paste0(feat_file_cell_lnpropexpect_dir,".Rdata"))

  
  
  time_output(start1)
  
  
  ## prop/expected short -----------------------------------------------------
  
  start1 = Sys.time()
  cat("ln prop/expected short")
  
  # markern = str_length(meta_cell[1,2])
  
  pmnegi = !grepl("[-]",meta_cell$phenotype)
  cellis1 = meta_cell$phenotype[meta_cell$phenolevel==1 & pmnegi]
  cellis1 = append("",cellis1[cellis1%in%names(meta_cell_parent_names)])
  expe1 = matrix(.5,nrow=nrow(mp), ncol=length(cellis1), dimnames=list(rownames(mp),cellis1))
  expe1[,1] = 1
  
  cellis = meta_cell$phenotype[meta_cell$phenolevel>1 & pmnegi]
  cellis = cellis[cellis%in%names(meta_cell_parent_names)]
  cellisn = str_count(cellis,"[+-]")
  
  mpe = mp
  mpe[mpe==0] = min(mp[mp!=0])
  mpe = mpe[,match(append(cellis1,cellis),colnames(mpe))]
  
  meta_cell_parent_names_ = meta_cell_parent_names[names(meta_cell_parent_names)%in%colnames(mpe)]
  meta_cell_parent_names_[cellis] = llply(meta_cell_parent_names_[cellis], function(x) {
    a = x[x%in%colnames(mpe)]
    if (length(a)==0) return(NULL)
    return(a)
  })
  meta_cell_parent_names_ = Filter(Negate(is.null), meta_cell_parent_names_)
  
  expec = llply(loopInd(1:length(cellis),no_cores), function(ii) {
    expect1 = foreach(ic=ii, .combine="cbind") %do% {
      i = cellis[ic]
      il = cellisn[ic]
      pnames = meta_cell_parent_names_[[i]]
      parent = mpe[,pnames,drop=F]
      numrtr = apply(parent, 1, function(x) prod(x^(il-1)))
      
      gnames = unique(unlist(meta_cell_parent_names_[pnames]))
      if (gnames=="") gnames = colnames(mpe)==""
      grprnt = mpe[,gnames,drop=F]
      denmtr = apply(grprnt, 1, prod)
      
      expect = (numrtr/denmtr)^(1/choose(il,2))
      return(expect)
    }
    return(expect1)
  }, .parallel=T)
  expec1 = Reduce("cbind", expec)
  colnames(expec1) = colnames(expec) = cellis
  exp1 = cbind(expe1,expec)
  lnpropexpect1 = log(mpe/exp1)
  lnpropexpect1[exp1==0] = log(mpe[exp1==0])

  save(exp1, file=paste0(feat_file_cell_lnpropexpect_dir,"_-short.Rdata"))
  save(lnpropexpect1, file=paste0(feat_file_cell_lnpropexpect_dir,"-short.Rdata"))

  
  
  time_output(start1)
  
}
time_output(start)

