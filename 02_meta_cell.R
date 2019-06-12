## input: meta_cell e.g. a+b+, a+b- ...
## output: cell_childpn_names & cell_parent_names; lists labelled with each cell population, the content of which are their cildren/parent (for children, they're split into positive and negative e.g. a+ --> a+b+ a+c+, a+b- a+c-)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)


## libraries
source("source/_func.R")
libr(c("foreach", "doMC","stringr","plyr"))

## cores
no_cores = 6# detectCores()-4
registerDoMC(no_cores)



## options
options(stringsAsFactors=F)
options(device="cairo")
options(na.rm=T)


for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  # result_dir = paste0(root, "/result/flowcap_panel6") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  
  start = Sys.time()
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta")
  meta_cell_dir = paste(meta_dir, "/cell", sep="")
  feat_dir = paste0(result_dir,"/feat")
  file_cell_count_dir = paste(feat_dir, "/file-cell-countAdj", sep="")
  
  ## output directories
  meta_cell_childpn_names_dir = paste(meta_dir, "/cell_childpn_names",sep="")
  meta_cell_parent_names_dir = paste(meta_dir, "/cell_parent_names",sep="") #specifies a phenotypes parents

  ## prepare data
  mc = get(load(paste0(file_cell_count_dir,".Rdata")))
  meta_cell = getPhen(colnames(mc))
  save(meta_cell, file=paste0(meta_cell_dir,".Rdata"))
  meta_cell_grid = ldply(1:nrow(meta_cell), function(i) {
    pc = meta_cell$phenocode[i]
    pc = as.numeric(laply(seq(1, nchar(pc), 1), function(i) substr(pc, i, i)))
  }, .parallel=T)
  colnames(meta_cell_grid) = unique(gsub("[+-]","",meta_cell$phenotype[meta_cell$phenolevel==1]))
  rownames(meta_cell_grid) = meta_cell$phenotype
  mcsum = rowSums(meta_cell_grid)
  mcpls = llply(unique(meta_cell$phenolevel), function(x) meta_cell$phenolevel==x)
  names(mcpls) = unique(meta_cell$phenolevel)
  maxl = max(meta_cell$phenolevel)
  minl = min(meta_cell$phenolevel)
  
  ilevel = meta_cell$phenolevel==minl
  iparen = NULL; ichild = meta_cell$phenolevel==minl+1
  res = NULL
  for (pl in minl:maxl) {
    start2 = Sys.time()
    cat(sum(ilevel)," pops ")
    
    if(!is.null(ichild)) ccand = meta_cell_grid[ichild,,drop=F]
    if(!is.null(iparen)) pcand = meta_cell_grid[iparen,,drop=F]
                  
    loop_ind = loopInd(which(ilevel),no_cores)
    result = llply(loop_ind, function(ii) {
      resulti = NULL
      if (pl==0) {
        pos_ind = apply(ccand,1,function(x) 2%in%x)
        resulti[[""]]$pos = meta_cell$phenotype[which(ichild)[pos_ind]]
        resulti[[""]]$neg = meta_cell$phenotype[which(ichild)[!pos_ind]]
      } else if(!is.null(ichild)) {
        for (i in ii) {
          ci = cn = cp = NULL
          mci = meta_cell_grid[i,,drop=T]
          ci_ind = llply(which(mci>0), function(j) ccand[,j]==mci[j])
          ci = which(ichild)[Reduce("&",ci_ind)]
          names(ci) = meta_cell$phenotype[ci]
          if (length(ci)>0) {
            cn = ci[mcsum[ci]-1==mcsum[i]]; if(length(cn)==0) cn = NULL
            cp = ci[mcsum[ci]-2==mcsum[i]]; if(length(cp)==0) cp = NULL
          }
          resulti[[meta_cell$phenotype[i]]]$pos = names(cp)
          resulti[[meta_cell$phenotype[i]]]$neg = names(cn)
          
        }
      } else {
        for (i in ii) resulti[[meta_cell$phenotype[i]]]$pos = NULL
        for (i in ii) resulti[[meta_cell$phenotype[i]]]$neg = NULL
      }
      
      if (pl==1) {
        for (i in ii) resulti[[meta_cell$phenotype[i]]]$parent = ""
      } else if(!is.null(iparen)) {
        for (i in ii) {
          pi = NULL
          pi_ind = sapply(which(mci>0), function(j) pcand[,j]==mci[j])
          pi = which(iparen)[rowSums(pi_ind)==(sum(mci>0)-1)]
          names(pi) = meta_cell$phenotype[pi]
          if (length(pi)>0)
            resulti[[meta_cell$phenotype[i]]]$parent = names(pi)
        }
      } else {
        for (i in ii) resulti[[meta_cell$phenotype[i]]]$parent = NULL
      }
    
      return(resulti)
    }, .parallel=T)
    result = unlist(result, recursive=F)
    names(result) = meta_cell$phenotype[ilevel]
    res = append(res, result)
    
    iparen = ilevel
    ilevel = ichild
    ichild = NULL; if ((pl+1)!=maxl) ichild = meta_cell$phenolevel==pl+2
    time_output(start2,paste0("layer",pl))
  }
  
  pchild = llply(res, function(x) {
    a = NULL
    if (!is.null(x$neg)) a$neg = x$neg
    if (!is.null(x$pos)) a$pos = x$pos
    return(a)
  })
  pchild = Filter(Negate(is.null), pchild)
  
  pparen = llply(res, function(x) return(x$parent))
  pparen = Filter(Negate(is.null), pparen)
  
  
  # start1 = Sys.time()
  
  # cat("\ncreating children indices ")
  # 
  # maxl = max(meta_cell$phenolevel)
  # minl = min(meta_cell$phenolevel)
  # 
  # phenolevel_ind = lapply(unique(meta_cell$phenolevel), function(x) meta_cell$phenolevel==x)
  # names(phenolevel_ind) = unique(meta_cell$phenolevel)
  # 
  # phenotype_split = str_extract_all(meta_cell$phenotype,"[a-zA-Z0-9]+[-|+]")
  # phenotype_split = lapply(phenotype_split, function(x) gsub("[+]","[+]",x))
  # phenotype_split = lapply(phenotype_split, function(x) gsub("[-]","[-]",x))
  # if (length(phenotype_split[[1]])==0) phenotype_split[[1]] = ""
  # 
  # loop.ind = loopInd(1:nrow(meta_cell), no_cores)
  # pcpp0 = foreach (ptii = loop.ind) %dopar% {
  #   pclist = list()
  #   pplist = list()
  #   for (pti in ptii) {
  #     pl = meta_cell$phenolevel[pti]
  #     pt = phenotype_split[[pti]]
  #     pc = pp = "NA"
  #     
  #     if (pl!=maxl) {
  #       pc_cand = meta_cell$phenotype[phenolevel_ind[[as.character(pl+1)]]]
  #       pc_cand_s = phenotype_split[ phenolevel_ind[[as.character(pl+1)]] ]
  #       pc_ind = Reduce("&",lapply(pt, function(x) grepl(x,pc_cand)))
  #       pc_ = pc_cand[pc_ind]
  #       pc_s = pc_cand_s[pc_ind]
  #       pc_pos_ind = sapply(pc_s, function(x) grepl("[+]",x[!x%in%pt]) )
  #       pc_neg = pc_[!pc_pos_ind]
  #       pc_pos = pc_[pc_pos_ind]
  #       pc = list(neg=pc_neg, pos=pc_pos)
  #       #order
  #       pc_inter = intersect(gsub("[-]","+",pc_neg), pc_pos)
  #       pc_neg = append(gsub("[+]","-",pc_inter), pc_neg[!gsub("[-]","+",pc_neg)%in%pc_inter])
  #       pc_pos = append(pc_inter, pc_pos[!pc_pos%in%pc_inter])
  #     }
  #     
  #     if (pl==1) {
  #       pp = ""
  #     } else if (pl!=minl) {
  #       pp_cand = meta_cell$phenotype[ phenolevel_ind[[as.character(pl-1)]] ]
  #       pp_cand_s = phenotype_split[ phenolevel_ind[[as.character(pl-1)]] ]
  #       pc_ind = sapply(pp_cand_s, function(x) all(x%in%pt))
  #       pp = pp_cand[pc_ind]
  #     }
  #     
  #     pclist[[length(pclist)+1]] = pc
  #     pplist[[length(pplist)+1]] = pp
  #   }
  #   return(list(pclist=pclist,pplist=pplist))
  # }
  # 
  # pchild0 = lapply(pcpp0, function(x) x$pclist)
  # pparen0 = lapply(pcpp0, function(x) x$pplist)
  # pchild1 = unlist(pchild0, recursive=F)
  # pparen1 = unlist(pparen0, recursive=F)
  # names(pchild1) = names(pparen1) = meta_cell$phenotype
  # pchild = pchild1[sapply(pchild1, function(x) unlist(x)[1]!="NA")]
  # pparen = pparen1[sapply(pparen1, function(x) x[1]!="NA")]
  
  
  #save
  save(pchild, file=paste0(meta_cell_childpn_names_dir, ".Rdata"))
  save(pparen, file=paste0(meta_cell_parent_names_dir, ".Rdata"))
  
  time_output(start)
}



