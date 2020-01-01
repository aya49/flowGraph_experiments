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
libr(c("stringr", "Matrix", "entropy", "plyr", "gsubfn",
       "foreach", "doMC"))


## cores
no_cores = detectCores() - 1
registerDoMC(no_cores)


## options
options(stringsAsFactors=FALSE)
options(na.rm=T)

writecsv = F

feat_count = "file-cell-countAdj"



start = Sys.time()
result_dirs = list.dirs(paste0(root, "/result"), full.names=T, recursive=F)
for (result_dir in result_dirs) {
  if (grepl("paired",result_dir)) next
  # if (!grepl("flowcap",result_dir)) next
  print(result_dir)
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta")
  meta_cell_parent_names_dir = paste(meta_dir, "/cell_parent_names",sep="") #specifies a phenotypes parents
  feat_dir = paste(result_dir, "/feat", sep=""); dir.create(feat_dir, showWarnings=F)
  feat_file_cell_countAdj_dir = paste(feat_dir, "/file-cell-countAdj", sep="")
  
  
  ## output directories
  meta_cell_dir = paste(meta_dir, "/cell", sep="")
  meta_cell_childpn_names_dir = paste(meta_dir, "/cell_childpn_names",sep="") #specifies a phenotypes children and splits them into +/- (only for when both -/+ exists)
  meta_cell_parent_names_dir = paste(meta_dir, "/cell_parent_names",sep="") #specifies a phenotypes parents
  meta_cell_graph_dir = paste(meta_dir, "/cell_graph",sep="") #specifies a phenotypes parents
  meta_cell_graphpos_dir = paste(meta_dir, "/cell_graphpos",sep="") #specifies a phenotypes parents
  feat_file_cell_prop_dir = paste(feat_dir, "/file-cell-prop", sep="")
  feat_file_edge_pnratio_dir = paste(feat_dir, "/file-edge-pnratio",sep="")
  feat_file_edge_prop_dir = paste(feat_dir, "/file-edge-prop",sep="")
  feat_file_cell_entropychild_dir = paste(feat_dir, "/file-cell-entropychild",sep="")
  feat_file_cell_entropyparent_dir = paste(feat_dir, "/file-cell-entropyparent",sep="")
  feat_file_group_entropy_dir = paste(feat_dir, "/file-group-entropy",sep="")
  feat_file_group_var_dir = paste(feat_dir, "/file-group-var",sep="")
  feat_file_cell_lnpropexpect_dir = paste(feat_dir, "/file-cell-lnpropexpect",sep="")
  feat_file_cell_lncountexpect_dir = paste(feat_dir, "/file-cell-lncountexpect",sep="")
  
  
  
  
  start = Sys.time()
  
  
  ## load data
  m = get(load(paste0(feat_file_cell_countAdj_dir,".Rdata")))
  
  ## make prop
  mp = m/m[,colnames(m)==""]
  dimnames(mp) = dimnames(m)
  save(mp, file=paste0(feat_file_cell_prop_dir,".Rdata"))
  if (writecsv) write.csv(mp, file=paste0(feat_file_cell_prop_dir,".csv"), row.names=T)
  
  ## make cell meta
  meta_cell = getPhen(colnames(m))
  save(meta_cell, file=paste0(meta_cell_dir,".Rdata"))
  
  pccell = getPhenCP(meta_cell=meta_cell,no_cores=no_cores)
  meta_cell_childpn_names = pchild = pccell$pchild
  meta_cell_parent_names = pparen = pccell$pparen
  gr = pccell$gr
  grp = pccell$grp
  save(pchild, file=paste0(meta_cell_childpn_names_dir, ".Rdata"))
  save(pparen, file=paste0(meta_cell_parent_names_dir, ".Rdata"))
  save(gr, file=paste0(meta_cell_graph_dir, ".Rdata"))
  save(grp, file=paste0(meta_cell_graphpos_dir, ".Rdata"))
  # graphs are saved; one with all cell populations, one with only positive ones a+b+... these can be converted to igraph with graph_from_data_frame
  
  
  
  
  # prepare parent list and list of cell populations
  markern = str_length(meta_cell[1,2])
  
  # for expected proportions only
  cellis1 = meta_cell$phenotype[meta_cell$phenolevel==1]
  cellis1 = append("",cellis1[cellis1%in%names(meta_cell_parent_names)])
  cellis = meta_cell$phenotype[meta_cell$phenolevel>1]
  cellis = cellis[cellis%in%names(meta_cell_parent_names)]
  
  cellis_ = append(cellis1,cellis)
  
  expe1 = matrix(.5,nrow=nrow(mp), ncol=length(cellis1), dimnames=list(rownames(mp),cellis1))
  expe1[,1] = 1
  
  mpe = mp
  # mpe[mpe==0] = min(mp[mp!=0])
  mpe = mpe[,match(append(cellis1,cellis),colnames(mpe))]
  
  meta_cell_parent_names_ = meta_cell_parent_names[names(meta_cell_parent_names)%in%colnames(mpe)]
  meta_cell_parent_names_[cellis] = llply(meta_cell_parent_names_[cellis], function(x) {
    a = x[x%in%colnames(mpe)]
    if (length(a)==0) return(NULL)
    return(a)
  })
  meta_cell_parent_names_ = plyr::compact(meta_cell_parent_names_)
  
  
  # goodind = cellis%in%names(meta_cell_parent_names_)
  # cellis = cellis[goodind]
  # mpe = mpe[,goodind]
  cellisn = str_count(cellis,"[+-]")
  
  
  
  
  ## child proportion --------------------------------------------
  
  start1 = Sys.time()
  
  #mlist=list()
  cat("childprop")
  childprop_ = foreach(ii = loopInd(2:length(cellis_),no_cores), .combine="cbind") %dopar% { #for each phenotype
    childprop1 = foreach(ic=ii, .combine="cbind") %do% {
      #for (i in 1:length(meta_cell_childpn_names)) {
      i = cellis_[ic]
      pnames = meta_cell_parent_names_[[i]]
      if (pnames[1]=="") {
        parent = mpe[,colnames(mpe)=="",drop=F]
      } else {
        parent = mpe[,pnames,drop=F]
      }
      # tryCatch({
      childprop2 = mpe[,i]/parent
      # }, error = function(err) { 
      #   parent = unlist(parent)
      #   childprop = childr
      #   for (j in 1:ncol(childr)) childprop[,j] = exp(unlist(childr[,j])/parent)
      # })
      colnames(childprop2) = paste0(pnames, "_", i)
      return(childprop2)
    }
    return(childprop1)
  }
  childprop_[is.nan(childprop_)] = 0
  childprop = exp(as.matrix(childprop_))
  # for (i in 1:length(childprop)) {
  #   colnames(childprop[[i]]) = paste0(names(childprop)[i], "_", colnames(childprop[[i]]))
  # }
  save(childprop, file=paste0(feat_file_edge_prop_dir,".Rdata"))
  if (writecsv) write.csv(childprop, file=paste0(feat_file_edge_prop_dir,".csv"))
  
  time_output(start1)
  
  
  
  
  if (F) {
    
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
    pnratio = as.matrix(pnratio)
    save(pnratio, file=paste0(feat_file_edge_pnratio_dir, ".Rdata"))
    if (writecsv) write.csv(pnratio, file=paste0(feat_file_edge_pnratio_dir, ".csv"))
    
    time_output(start1)
    
    
    
    ## child entropy --------------------------------------------
    
    start1 = Sys.time()
    cat("childentropy")
    
    meta_cell_childpn_names_ = meta_cell_childpn_names
    for (i in names(meta_cell_childpn_names)) {
      if (length(meta_cell_childpn_names[[i]])<2) meta_cell_childpn_names_[[i]] = NULL
    }
    enc = a = as.matrix(foreach(ii = loopInd(1:length(meta_cell_childpn_names_),no_cores), .combine="cbind") %dopar% { #for each phenotype
      return(foreach(i=ii, .combine="cbind") %do% {
        if (length(meta_cell_childpn_names_[[i]])<2) return(NULL)
        # i = which(names(meta_cell_childpn_names_)==ii)
        pnames = names(meta_cell_childpn_names_)[i]
        pnames = pnames[pnames%in%colnames(m)]
        parent = m[,match(pnames,colnames(m)),drop=F]
        cnames = unlist(meta_cell_childpn_names_[[i]])
        cnames = cnames[cnames%in%colnames(m)]
        childr = m[,match(cnames,colnames(m)),drop=F]
        # childr[childr<1] = parent[parent<1] = 1
        
        # Entropy matrix
        en = rep(0,nrow(m))
        if (ncol(parent)==0 | ncol(childr)==0) return(NULL)
        
        no_child = length(cnames)
        non0parents = parent>0
        non0childs = Reduce("&",llply(1:ncol(childr),function(x) childr[,x]>0))
        if (sum(non0parents & non0childs)>0) en[non0parents & non0childs] = sapply(which(non0parents & non0childs), function(x) entropy(childr[x,]/parent[x])/no_child) #average entropy over # of markers added
        
        en = matrix(en, ncol=1)
        colnames(en) = names(meta_cell_childpn_names_)[i]
        
        #mlist[[i]] = list(pnratio=pnratio, childprop=childprop)
        return(en) #ratio = +child_prop / -child_prop
      })
    })
    rownames(enc) = rownames(m)
    enc = enc[,apply(enc,2,function(x) any(x!=0))]
    
    enc = as.matrix(enc)
    save(enc, file=paste0(feat_file_cell_entropychild_dir, ".Rdata"))
    if (writecsv) write.csv(enc, file=paste0(feat_file_cell_entropychild_dir, ".csv"))
    
    time_output(start1)
    
    
    
    ## parent entropy --------------------------------------------
    
    start1 = Sys.time()
    cat("parententropy")
    
    meta_cell_parent_names_ = meta_cell_parent_names
    for (i in names(meta_cell_parent_names))
      if (length(meta_cell_parent_names[[i]])<2) meta_cell_parent_names_[[i]] = NULL
    enp = as.matrix(foreach(ii = loopInd(1:length(meta_cell_parent_names_), no_cores), .combine='cbind') %dopar% { #for each phenotype
      # Entropy matrix
      enpi = foreach(i=ii, .combine='cbind') %do% { #for each phenotype
        
        pnames = meta_cell_parent_names_[[i]]
        pnames = pnames[pnames%in%colnames(m)]
        parent = m[,pnames,drop=F]
        
        cnames = names(meta_cell_parent_names_)[i]
        cnames = cnames[cnames%in%colnames(m)]
        childr = m[,match(cnames,colnames(m)),drop=F]
        
        en = rep(0,nrow(m))
        if (ncol(parent)==0 | ncol(childr)==0) return(NULL)
        childr[childr<1] = parent[parent<1] = 1
        
        no_paren = length(cnames)
        non0childrs = childr>0
        non0paren = Reduce("&",llply(1:ncol(parent),function(x) parent[,x]>0))
        if (sum(non0paren & non0childrs)>0) en[non0paren & non0childrs] = sapply(which(non0paren & non0childrs), function(x) entropy(rep(childr[x],ncol(parent))/parent[x,])/no_paren) #average entropy over # of markers added
        
        
        en = matrix(en, ncol=1)
        colnames(en) = names(meta_cell_parent_names_)[i]
        
        return(en) #ratio = +child_prop / -child_prop
      }
      return(enpi)
    })
    rownames(enp) = rownames(m)
    enp = enp[,!apply(enp,2,function(x) all(x==0))]
    
    enp = as.matrix(enp)
    save(enp, file=paste0(feat_file_cell_entropyparent_dir, ".Rdata"))
    if (writecsv) write.csv(enp, file=paste0(feat_file_cell_entropyparent_dir, ".csv"))
    
    time_output(start1)
    
    
    
    ## group entropy & variance ---------------------------------------
    
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
      sapply(groups[gi],function(i) apply(m[,i],1, function(x) {
        if (all(x==0)) return(0)
        return(entropy(x))
      } )), .parallel=T)
    feat_file_group_entropy = Reduce("cbind",a)
    b = llply(loopInd(1:length(groups),no_cores), function(gi) 
      sapply(groups[gi],function(i) apply(m[,i],1,var)), .parallel=T)
    feat_file_group_var = Reduce("cbind",b)
    e = llply(loopInd(1:length(groups),no_cores), function(gi) 
      sapply(groups[gi],function(i) apply(mp[,i],1, function(x) {
        if (all(x==0)) return(0)
        return(entropy(x))
      } )),.parallel=T)
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
    
    feat_file_group_entropy = as.matrix(feat_file_group_entropy)
    save(feat_file_group_entropy, file=paste0(feat_file_group_entropy_dir, ".Rdata"))
    if (writecsv) write.csv(feat_file_group_entropy, file=paste0(feat_file_group_entropy_dir, ".csv"))
    
    feat_file_group_var = as.matrix(feat_file_group_var)
    save(feat_file_group_var, file=paste0(feat_file_group_var_dir, ".Rdata"))
    if (writecsv) write.csv(feat_file_group_var, file=paste0(feat_file_group_var_dir, ".csv"))
    
    feat_file_group_entropyp = as.matrix(feat_file_group_entropyp)
    save(feat_file_group_entropyp, file=paste0(feat_file_group_entropy_dir, "p.Rdata"))
    if (writecsv) write.csv(feat_file_group_entropyp, file=paste0(feat_file_group_entropy_dir, "p.csv"))
    
    feat_file_group_varp = as.matrix(feat_file_group_varp)
    save(feat_file_group_varp, file=paste0(feat_file_group_var_dir, "p.Rdata"))
    if (writecsv) write.csv(feat_file_group_varp, file=paste0(feat_file_group_var_dir, "p.csv"))
    
    time_output(start1)
    
    
    
  }
  
  
  
  
  ## prop/expected ---------------------------------------------------
  
  start1 = Sys.time()
  cat("ln prop/expected")
  
  
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
      il = cellisn[ic] # layer
      pnames = meta_cell_parent_names_[[i]]
      parent = mpe[,pnames,drop=F]
      numrtr = apply(parent, 1, function(x) prod(x^(il-1)))
      
      gnames = unique(unlist(meta_cell_parent_names_[pnames]))
      if (gnames[1]=="") gnames = colnames(mpe)==""
      grprnt = mpe[,gnames,drop=F]
      denmtr = apply(grprnt, 1, prod)
      
      expect = (numrtr/denmtr)^(1/ncol(grprnt))
      expect[numrtr==0 & denmtr==0] = 0
      return(expect)
    }
    return(expect2)
  }, .parallel=T)
  expec = Reduce("cbind", expec)
  # which(is.na(expec),arr.ind=T)
  colnames(expec) = cellis
  exp1 = as.matrix(cbind(expe1,expec))
  expc1 = exp1*m[,colnames(m)==""]
  lnpropexpect1 = log(mpe/exp1)
  mei = m[,match(colnames(expc1),colnames(m))]
  lncountexpect1 = log(mei/expc1)
  lnpropexpect1[exp1==0] = log(mpe[exp1==0])
  lncountexpect1[expc1==0] = log(mei[expc1==0])
  lnpropexpect1[mpe==0] = lncountexpect1[mei==0] = 0
  
  save(exp1, file=paste0(feat_file_cell_lnpropexpect_dir,"-raw.Rdata"))
  save(expc1, file=paste0(feat_file_cell_lncountexpect_dir,"-raw.Rdata"))
  save(lnpropexpect1, file=paste0(feat_file_cell_lnpropexpect_dir,".Rdata"))
  save(lncountexpect1, file=paste0(feat_file_cell_lncountexpect_dir,".Rdata"))
  
  time_output(start1)
  
  
  
  
  start1 = Sys.time()
  cat("ln prop/expected2")
  
  
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
      # mpe0 = mpe
      # mpe = matrix(colMeans(mpe[501:1000,]),nrow=1)
      # colnames(mpe) = colnames(mpe0)
      
      i = cellis[ic]
      il = cellisn[ic]
      
      # pnames = meta_cell_parent_names_[[i]]
      # if (il>3)
      #   for (z in 1:(il-3)) {
      #     pnames = unique(unlist(meta_cell_parent_names_[pnames]))
      #     parent = mpe[,pnames,drop=F]
      #     numrtr = apply(parent, 1, function(x) prod(x)^(il-1))
      #     
      #     gnames = unique(unlist(meta_cell_parent_names_[pnames]))
      #     if (gnames[1]=="") gnames = colnames(mpe)==""
      #     grprnt = mpe[,gnames,drop=F]
      #     denmtr = apply(grprnt, 1, function(x) prod(x)^2)
      #   }
      
      pnames = meta_cell_parent_names_[[i]]
      parent = mpe[,pnames,drop=F]
      numrtr = apply(parent, 1, function(x) prod(x)^(il-1))
      
      gnames = unique(unlist(meta_cell_parent_names_[pnames]))
      if (gnames[1]=="") gnames = colnames(mpe)==""
      grprnt = mpe[,gnames,drop=F]
      denmtr = apply(grprnt, 1, function(x) prod(x)^2)
      
      expect = (numrtr/denmtr)^(1/(il-1))
      expect[numrtr==0 & denmtr==0] = 0
      
      return(expect)
    }
    return(expect2)
  }, .parallel=T)
  expec = Reduce("cbind", expec)
  # which(is.na(expec),arr.ind=T)
  colnames(expec) = cellis
  exp1 = as.matrix(cbind(expe1,expec))
  expc1 = exp1*m[,colnames(m)==""]
  lnpropexpect1 = log(mpe/exp1)
  mei = m[,match(colnames(expc1),colnames(m))]
  lncountexpect1 = log(mei/expc1)
  lnpropexpect1[exp1==0] = log(mpe[exp1==0])
  lncountexpect1[expc1==0] = log(mei[expc1==0])
  lnpropexpect1[mpe==0] = lncountexpect1[mei==0] = 0
  
  save(exp1, file=paste0(feat_file_cell_lnpropexpect_dir,"2-raw.Rdata"))
  save(expc1, file=paste0(feat_file_cell_lncountexpect_dir,"2-raw.Rdata"))
  save(lnpropexpect1, file=paste0(feat_file_cell_lnpropexpect_dir,"2.Rdata"))
  save(lncountexpect1, file=paste0(feat_file_cell_lncountexpect_dir,"2.Rdata"))
  
  time_output(start1)
  
  
  
  start1 = Sys.time()
  cat("ln prop/expected3")
  
  childprop_names = Reduce("rbind",str_split(colnames(childprop_),"_"))
  colnames(childprop_names) = c("from","to")
  childprop_names = as.data.frame(childprop_names)
  rownames(childprop_names) = 1:nrow(childprop_names)
  
  cpind = which(!grepl("[-]",cellis))
  expecp = llply(loopInd(cpind,no_cores), function(ii) {
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
      # mpe0 = mpe
      # mpe = matrix(colMeans(mpe[501:1000,]),nrow=1)
      # colnames(mpe) = colnames(mpe0)
      
      cpop = i = cellis[ic]
      il = cellisn[ic]
      
      # ## gets all ancestors
      # ge_ = childprop_names[childprop_names$to==cpop,,drop=F]
      # cpopl = ge_$from
      # while (cpopl[1]!="") {
      #   etemp = childprop_names[childprop_names$to%in%cpopl,,drop=F]
      #   ge_ = rbind(ge_, etemp)
      #   cpopl = unique(etemp$from)
      # }
      # ge = childprop_[,as.numeric(rownames(ge_)),drop=F]
      # # gem = colMeans(ge[501:1000,]) # check positive controls with this
      # ge_ = cbind(ge_,ge[501,])
      
      pnames = meta_cell_parent_names_[[i]]
      parent = mpe[,pnames,drop=F]
      # parentmin_ = apply(parent,1,min)
      # parentmax_ = apply(parent,1,max)
      # parentmean_ = rowMeans(cbind(parentmax_,parentmin_))
      # parentsd_ = apply(cbind(parentmax_,parentmin_),1,sd)
      parento = laply(1:nrow(parent), function(xi) order(parent[xi,]))#, decreasing=T))
      
      
      pedges = Reduce(cbind,llply(1:length(pnames), function(gi) {
        gname = meta_cell_parent_names_[[pnames[gi]]]
        cns = paste0(pnames[gi],"_",gname)
        if (gname[1]=="") gname = colnames(mpe)==""
        grprnt = mpe[,gname,drop=F]
        edges = parent[,gi]/grprnt # matrix
        colnames(edges) = cns
        edges
      }))
      
      ## take the mean of min and max
      # parentr = apply(parent,1,function(x) max(x)/min(x))
      expect1 = expect2 = rep(0, nrow(parent))
      for (colj in unique(parento[,1])) {
        jind = parento[,1]==colj
        
        parentjmin = parent[jind,colj]
        # parentjmax = parentmax_[jind]
        # parentjmean = parentmean_[jind]
        # parentsd = parentsd_[jind]
        
        pnamej = gsub("[+]","[+]",pnames[colj])
        pnamej = gsub("[-]","[-]",pnamej)
        pedgesnj = pedges[jind,!grepl(pnamej,colnames(pedges)),drop=F]
        
        # shrink ratio, value too extreme otherwise
        # pedgesnjr = apply(pedgesnj,1,function(x) max(x)/min(x))
        # 
        # pedgesnjmax = apply(pedgesnj, 1, max)
        # pedgesnjmin = apply(pedgesnj, 1, min)
        # pedgesncbind = cbind(pedgesnjmin,pedgesnjmax)
        # pedgesnjmean = rowMeans(pedgesncbind)
        # pedgesnjsd = apply(pedgesncbind,1,sd)
        # 
        # parentjr = parentr[jind]
        # ebig = pedgesnjsd>parentsd
        # ebig[is.na(ebig)] = F
        # if (sum(ebig)>0) {
        #   pedgesnjmax[ebig] = pedgesnjmean[ebig] + (pedgesnjmax[ebig]-pedgesnjmean[ebig])*(parentsd[ebig]/pedgesnjsd[ebig])
        # }
        # pbig = pedgesnjsd<=parentsd
        # pbig[is.na(pbig)] = F
        # if (sum(pbig)>0) {
        #   parentjmin[pbig] = parentjmean[pbig] + (parentjmin[pbig]-parentjmean[pbig])*(pedgesnjsd[pbig]/parentsd[pbig])
        # }
        # expect[which(jind)[(ebig | pbig)]] = parentjmin[(ebig | pbig)] * pedgesnjmax[(ebig | pbig)]
        
        expect1[jind] = apply(pedgesnj,1,max) * parentjmin
      }
      
      # for (colj in unique(parento[,ncol(parento)])) {
      #   jind = parento[,ncol(parento)]==colj
      #   
      #   parentjmax = parent[jind,colj]
      # 
      #   pnamej = gsub("[+]","[+]",pnames[colj])
      #   pnamej = gsub("[-]","[-]",pnamej)
      #   pedgesnj = pedges[jind,!grepl(pnamej,colnames(pedges)),drop=F]
      #   
      #   expect2[jind] = apply(pedgesnj,1,min) * parentjmax
      # }
      
      expect = expect1 #rowMeans(cbind(expect1,expect2))
      return(expect)
      
    }
    return(expect2)
  }, .parallel=T)
  expecp = Reduce("cbind", expecp)
  expecp[is.nan(expecp)] = 0
  # which(is.na(expec),arr.ind=T)
  colnames(expecp) = cellis[cpind]
  
  # infer rest of cell populations expected proportion
  expec = as.matrix(cbind(expe1,expecp))
  cpopneg = setdiff(cellis,colnames(expec)) # cell pops to do
  cpopnegl = str_count(cpopneg,"[-+]")
  p = proto(i = 1, j = 1, fun = function(this, x) if (count >= i && count <= j) "+" else x) # replaces whatever pattern only once with "+";   gsubfn("_", p, "A_B_C_")

  for (lev in sort(unique(cpopnegl))) {
    cpopl = cpopneg[cpopnegl==lev]
    cpopnegno = str_count(cpopl,"[-]") # number of negatives
    for (negno in sort(unique(cpopnegno))) {
      cpopi = cpopl[cpopnegno==negno]
      sibs = gsubfn("-", p, cpopi)
      for (sib in unique(sibs)) {
        sibpars = meta_cell_parent_names_[[sib]]
        sibis = which(sibs==sib)
        for (sibi in sibis) {
          parsin = intersect(sibpars, meta_cell_parent_names_[[cpopi[sibi]]])
          # parsin = parsin[which(parsin%in%colnames(mpe))[1]]
          ex = mpe[,parsin] - mpe[,sib]
          expec = cbind(expec, ex)
          colnames(expec)[ncol(expec)] = cpopi[sibi]
        }
      }
    }
    expec_temp = Reduce(cbind,llply(cpopi, function(cpi) {
      cpig = gsub("[+-]","",cpi)
      parnames = intersect(meta_cell_parent_names_[[cpi]],colnames(mpe))
      a = NULL
      for (parname in parnames) {
        chilnames = unlist(pchild[[parname]])
        chilnames = chilnames[chilnames!=cpi]
        chilnamesg = gsub("[+-]","",chilnames)
        chilgreat = chilnamesg==cpig & chilnames%in%colnames(expec)
        if (any(chilgreat)) {
          a = matrix(mpe[,parname]-expec[,chilnames[chilgreat]],ncol=1)
          break
        }
      }
      if (is.null(a)) {
        a = b = mpe[,parnames,drop=F]
        if (ncol(b)>1) a = matrix(apply(b,1,min), ncol=1)
      }
      colnames(a) = cpi
      return(a)
    }))
    if (is.null(dim(expec_temp))) expec_temp = matrix(expec_temp,ncol=1)
    # colnames(expec_temp) = cpopi
    expec = cbind(expec, expec_temp)
  }
  exp1 = cbind(expe1,expec[,match(cellis,colnames(expec))])
    
  expc1 = exp1*m[,1]
  a = mpe/exp1
  a[is.infinite(a)] = max(a[is.finite(a)])
  a[is.nan(a)] = 0
  lnpropexpect1 = log(a)
  lnpropexpect1[is.nan(lnpropexpect1)] = 0
  mei = m[,match(colnames(expc1),colnames(m))]
  a = mei/expc1
  a = mpe/exp1
  a[is.infinite(a)] = max(a[is.finite(a)])
  a[is.nan(a)] = 0
  lncountexpect1 = log(a)
  lncountexpect1[is.nan(lncountexpect1)] = 0
  lnpropexpect1[exp1==0] = log(mpe[exp1==0])
  lncountexpect1[expc1==0] = log(mei[expc1==0])
  lnpropexpect1[mpe==0] = lncountexpect1[mei==0] = 0
  
  if (grepl("/pos",result_dir)) {
    # plot to compare prop and expected prop
    gr0 = get(load(paste0(result_dir,"/meta/cell_graph.Rdata")))
    gr = layout_gr(gr0$e,gr0$v)
    gr = gpdf(gr)
    grorder = match(gr$v$name,colnames(exp1))
    mfmc = colMeans(mpe[1:500,grorder])
    mfmu = colMeans(mpe[501:1000,grorder])
    memu = colMeans(exp1[501:1000,grorder])
    gr$v$color = ifelse(mfmc<mfmu,"_increase","decrease") # node colour
    gr$v$label = paste0(gr$v$name,":",round(mfmu,3),"/",round(memu,3))
    main0 = paste0(result_dir,"\n specenrV3 = calc all positive cell pops first, others can be inferred \nlabel=prop / expected prop")
    gr$v$size = abs(mfmu-memu)/mfmu
    gr$v$label_ind = gr$v$v_ind = gr$v$size>.05
    gr$e$e_ind = gr$e[,1]%in%gr$v$name[gr$v$v_ind] &gr$e[,2]%in%gr$v$name[gr$v$v_ind]
    main = paste0(main0,"\nsize=(actual-expected)/actual")
    
    gp = gggraph(gr, main=main)
    
    ggsave(paste0(result_dir,"/meta/actualVSexpected3.png"), plot=gp, scale=1, width=9, height=9, units="in", dpi=500, limitsize=T)
    
    
    # ## comment out, this is only for one population
    # cpop = "A+B+C+D+"
    # gr_ = list()
    # gr_$e = gr$e[gr$e$to==cpop,,drop=F]
    # cpopl = gr_$e$from
    # while (cpopl[1]!="") {
    #   etemp = gr$e[gr$e$to%in%cpopl,,drop=F]
    #   gr_$e = rbind(gr_$e, etemp)
    #   cpopl = unique(etemp$from)
    # }
    # gr_$v = gr$v[gr$v$name%in%unique(unlist(gr_$e[,c(1,2)])),,drop=F]
    # gr_$v$label_ind = T
    # gp = gggraph(gr_, main=cpop)
    # plot(gp)
    
  }
  
  save(exp1, file=paste0(feat_file_cell_lnpropexpect_dir,"3-raw.Rdata"))
  save(expc1, file=paste0(feat_file_cell_lncountexpect_dir,"3-raw.Rdata"))
  save(lnpropexpect1, file=paste0(feat_file_cell_lnpropexpect_dir,"3.Rdata"))
  save(lncountexpect1, file=paste0(feat_file_cell_lncountexpect_dir,"3.Rdata"))
  
  time_output(start1)
  
  
  
  
  ## expected count
  
  
  ## prop/expected short ---------------------------------------------
  
  # start1 = Sys.time()
  # cat("ln prop/expected short")
  # 
  # # markern = str_length(meta_cell[1,2])
  # 
  # pmnegi = !grepl("[-]",meta_cell$phenotype)
  # cellis1 = meta_cell$phenotype[meta_cell$phenolevel==1 & pmnegi]
  # cellis1 = append("",cellis1[cellis1%in%names(meta_cell_parent_names)])
  # cellis = meta_cell$phenotype[meta_cell$phenolevel>1 & pmnegi]
  # cellis = cellis[cellis%in%names(meta_cell_parent_names)]
  # 
  # expe1 = matrix(.5,nrow=nrow(mp), ncol=length(cellis1), dimnames=list(rownames(mp),cellis1))
  # expe1[,1] = 1
  # 
  # mpe = mp
  # # mpe[mpe==0] = min(mp[mp!=0])
  # mpe = mpe[,match(append(cellis1,cellis),colnames(mpe))]
  # 
  # meta_cell_parent_names_ = meta_cell_parent_names[names(meta_cell_parent_names)%in%colnames(mpe)]
  # meta_cell_parent_names_[cellis] = llply(meta_cell_parent_names_[cellis], function(x) {
  #   a = x[x%in%colnames(mpe)]
  #   if (length(a)==0) return(NULL)
  #   return(a)
  # })
  # meta_cell_parent_names_ = plyr::compact(meta_cell_parent_names_)
  # 
  # # goodind = cellis%in%names(meta_cell_parent_names_)
  # # cellis = cellis[goodind]
  # # mpe = mpe[,goodind]
  # cellisn = str_count(cellis,"[+-]")
  # 
  # expec = llply(loopInd(1:length(cellis),no_cores), function(ii) {
  #   # expect1 = foreach(ic=ii, .combine="cbind") %do% {
  #   #   i = cellis[ic]
  #   #   il = cellisn[ic]
  #   #   pnames = meta_cell_parent_names_[[i]]
  #   #   parent = mpe[,pnames,drop=F]
  #   #   numrtr = apply(parent, 1, prod)^(il-1)
  #   #   
  #   #   gnames = unique(unlist(meta_cell_parent_names_[pnames]))
  #   #   if (gnames=="") gnames = colnames(mpe)==""
  #   #   grprnt = mpe[,gnames,drop=F]
  #   #   denmtr = apply(grprnt, 1, prod)
  #   #   
  #   #   expect = (numrtr/denmtr)^(1/choose(il,2))
  #   #   return(expect)
  #   # }
  #   expect2 = foreach(ic=ii, .combine="cbind") %do% {
  #     i = cellis[ic]
  #     il = cellisn[ic]
  #     pnames = meta_cell_parent_names_[[i]]
  #     parent = mpe[,pnames,drop=F]
  #     numrtr = apply(parent, 1, function(x) prod(x^(il-1)))
  #     
  #     gnames = unique(unlist(meta_cell_parent_names_[pnames]))
  #     if (gnames[1]=="") gnames = colnames(mpe)==""
  #     grprnt = mpe[,gnames,drop=F]
  #     denmtr = apply(grprnt, 1, prod)
  #     
  #     expect = (numrtr/denmtr)^(1/choose(il,2))
  #     expect[numrtr==0 & denmtr==0] = 0
  #     return(expect)
  #   }
  #   return(expect2)
  # }, .parallel=T)
  # expec = Reduce("cbind", expec)
  # which(is.na(expec),arr.ind=T)
  # colnames(expec) = cellis
  # exp1 = as.matrix(cbind(expe1,expec))
  # lnpropexpect1 = log(mpe/exp1)
  # lnpropexpect1[exp1==0] = log(mpe[exp1==0])
  # lnpropexpect1[mpe==0] = 0
  
  
  
  
  
  # ## USE THIS ONE
  # exp11 = exp1[,!grepl("[-]",colnames(exp1))]
  # save(exp11, file=paste0(feat_file_cell_lnpropexpect_dir,"-short-raw.Rdata"))
  # lnpropexpect11 = lnpropexpect1[,!grepl("[-]",colnames(lnpropexpect1))]
  # save(lnpropexpect11, file=paste0(feat_file_cell_lnpropexpect_dir,"-short.Rdata"))
  # 
  # expc11 = expc1[,!grepl("[-]",colnames(expc1))]
  # save(expc11, file=paste0(feat_file_cell_lncountexpect_dir,"-short-raw.Rdata"))
  # lncountexpect11 = lncountexpect1[,!grepl("[-]",colnames(lncountexpect1))]
  # save(lncountexpect11, file=paste0(feat_file_cell_lncountexpect_dir,"-short.Rdata"))
  
  
  
  # time_output(start1)
}

time_output(start)

