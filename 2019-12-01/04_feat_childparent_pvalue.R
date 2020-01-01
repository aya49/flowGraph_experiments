## input: meta_cell & base features e.g. countadj, props (each cell population's counts as a percentage of the total count), p values (for if we want to delete rows/columns with no significant values)
## output: features
## - parent_contrib (ratio of the difference between the avg control fcm file's and experiment fcm file's child population over that of the parent)
## - parent_effort (same as above except instead of difference, use log fold, and instead of ratio use different as values are in log scale)

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
  meta_cell_childpn_names_dir = paste(meta_dir, "/cell_childpn_names",sep="") #specifies a phenotypes children and splits them into +/- (only for when both -/+ exists)
  # meta_cell_childpn_ind_dir = paste(meta_dir, "/cell_childpn_ind",sep="")
  meta_cell_parent_names_dir = paste(meta_dir, "/cell_parent_names",sep="") #specifies a phenotypes parents
  
  
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
  
  m = as.matrix(Matrix(get(load(paste0(feat_file_cell_countAdj_dir,".Rdata")))))
  
  meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
  mp = Matrix(get(load(paste0(feat_file_cell_prop_dir,".Rdata"))))
  meta_cell_childpn_names = get(load(paste0(meta_cell_childpn_names_dir, ".Rdata")))
  meta_cell_parent_names = get(load(paste0(meta_cell_parent_names_dir, ".Rdata")))
  
  
  
  ## features that require p value script -----------------------
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
time_output(start)

