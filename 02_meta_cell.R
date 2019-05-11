## input: meta_cell e.g. a+b+, a+b- ...
## output: cell_childpn_names & cell_parent_names; lists labelled with each cell population, the content of which are their cildren/parent (for children, they're split into positive and negative e.g. a+ --> a+b+ a+c+, a+b- a+c-)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metrics"
setwd(root)
for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  # result_dir = paste0(root, "/result/flowcap_panel6") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen


## input directories
meta_dir = paste0(result_dir,"/meta")
meta_cell_dir = paste(meta_dir, "/cell", sep="")

## output directories
# meta_cell_child_dir = paste(meta_dir, "/cell_child",sep="") #specifies a phenotypes children
# meta_cell_child_ind_dir = paste(meta_dir, "/cell_child_ind",sep="")
# meta_cell_child_names_dir = paste(meta_dir, "/cell_child_names",sep="")
# meta_cell_childpn_dir = paste(meta_dir, "/cell_childpn",sep="") #specifies a phenotypes children and splits them into +/- (only for when both -/+ exists)
# meta_cell_childpn_ind_dir = paste(meta_dir, "/cell_childpn_ind",sep="")
meta_cell_childpn_names_dir = paste(meta_dir, "/cell_childpn_names",sep="")
# meta_cell_parent_dir = paste(meta_dir, "/cell_parent",sep="") #specifies a phenotypes parents
meta_cell_parent_names_dir = paste(meta_dir, "/cell_parent_names",sep="") #specifies a phenotypes parents
# meta_cell_parent_ind_dir = paste(meta_dir, "/cell_parent_ind",sep="")
# meta_cell_parentpn_dir = paste(meta_dir, "/cell_parentpn",sep="") #specifies a phenotypes parents and splits them into +/- (only for when both -/+ exists)
# meta_cell_parentpn_ind_dir = paste(meta_dir, "/cell_parentpn_ind",sep="")

## libraries
source("source/_funcAlice.R")
source("source/_funcdist.R")
libr(c("foreach", "doMC"))

## cores
no_cores = detectCores()-1
registerDoMC(no_cores)



## options
options(stringsAsFactors=F)
options(device="cairo")
options(na.rm=T)

## prepare data
meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))



start = Sys.time()

# start1 = Sys.time()

cat("\ncreating children indices ")

maxl = max(meta_cell$phenolevel)
minl = min(meta_cell$phenolevel)

phenolevel_ind = lapply(unique(meta_cell$phenolevel), function(x) meta_cell$phenolevel==x)
names(phenolevel_ind) = unique(meta_cell$phenolevel)

phenotype_split = str_extract_all(meta_cell$phenotype,"[a-zA-Z0-9]+[-|+]")
phenotype_split = lapply(phenotype_split, function(x) gsub("[+]","[+]",x))
phenotype_split = lapply(phenotype_split, function(x) gsub("[-]","[-]",x))
if (length(phenotype_split[[1]])==0) phenotype_split[[1]] = ""

loop.ind = loopInd(1:nrow(meta_cell), no_cores)
pcpp0 = foreach (ptii = loop.ind) %dopar% {
  pclist = list()
  pplist = list()
  for (pti in ptii) {
    pl = meta_cell$phenolevel[pti]
    pt = phenotype_split[[pti]]
    pc = pp = "NA"
    
    if (pl!=maxl) {
      pc_cand = meta_cell$phenotype[phenolevel_ind[[as.character(pl+1)]]]
      pc_cand_s = phenotype_split[ phenolevel_ind[[as.character(pl+1)]] ]
      pc_ind = Reduce("&",lapply(pt, function(x) grepl(x,pc_cand)))
      pc_ = pc_cand[pc_ind]
      pc_s = pc_cand_s[pc_ind]
      pc_pos_ind = sapply(pc_s, function(x) grepl("[+]",x[!x%in%pt]) )
      pc_neg = pc_[!pc_pos_ind]
      pc_pos = pc_[pc_pos_ind]
      pc = list(neg=pc_neg, pos=pc_pos)
      #order
      pc_inter = intersect(gsub("[-]","+",pc_neg), pc_pos)
      pc_neg = append(gsub("[+]","-",pc_inter), pc_neg[!gsub("[-]","+",pc_neg)%in%pc_inter])
      pc_pos = append(pc_inter, pc_pos[!pc_pos%in%pc_inter])
    }
    
    if (pl==1) {
      pp = ""
    } else if (pl!=minl) {
      pp_cand = meta_cell$phenotype[ phenolevel_ind[[as.character(pl-1)]] ]
      pp_cand_s = phenotype_split[ phenolevel_ind[[as.character(pl-1)]] ]
      pc_ind = sapply(pp_cand_s, function(x) all(x%in%pt))
      pp = pp_cand[pc_ind]
    }
    
    pclist[[length(pclist)+1]] = pc
    pplist[[length(pplist)+1]] = pp
  }
  return(list(pclist=pclist,pplist=pplist))
}

pchild0 = lapply(pcpp0, function(x) x$pclist)
pparen0 = lapply(pcpp0, function(x) x$pplist)
pchild1 = unlist(pchild0, recursive=F)
pparen1 = unlist(pparen0, recursive=F)
names(pchild1) = names(pparen1) = meta_cell$phenotype
pchild = pchild1[sapply(pchild1, function(x) unlist(x)[1]!="NA")]
pparen = pparen1[sapply(pparen1, function(x) x[1]!="NA")]



#save
# save(meta_cell_child, file=paste0(meta_cell_child_dir, ".Rdata"))
# save(meta_cell_child_ind, file=paste0(meta_cell_child_ind_dir, ".Rdata"))
# save(pchild, file=paste0(meta_cell_child_names_dir, ".Rdata")) #list of children for each named parent
# save(meta_cell_childpn, file=paste0(meta_cell_childpn_dir, ".Rdata"))
# save(meta_cell_childpn_ind, file=paste0(meta_cell_childpn_ind_dir, ".Rdata"))
save(pchild, file=paste0(meta_cell_childpn_names_dir, ".Rdata"))

# save(meta_cell_parent, file=paste0(meta_cell_parent_dir, ".Rdata"))
# save(meta_cell_parent_ind, file=paste0(meta_cell_parent_ind_dir, ".Rdata"))
save(pparen, file=paste0(meta_cell_parent_names_dir, ".Rdata"))

time_output(start)
}



