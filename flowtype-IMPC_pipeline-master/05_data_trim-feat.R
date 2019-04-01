# IMPC; Trim feature matrices using all pvalues
# aya43@sfu.ca 20180405

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

feat_dir = paste(result_dir, "/feat", sep=""); dir.create(feat_dir, showWarnings=F)
# feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")

## output directories

## libraries
source("code/_funcAlice.R")
libr("stringr")
libr("entropy")
libr("foreach")
libr("doMC")





## cores
no_cores = 15#detectCores() - 1
registerDoMC(no_cores)




## options
options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

writecsv = F

#matrices to not trim
notrimmatrix = "Max|pval"



writecsv = T




start = Sys.time()








start1 = Sys.time()

#get list of children for each non-leaf node & save
cat("\ncreating child matrix")

# m = get(load(paste0(feat_file_cell_count_dir,".Rdata")))
meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))

pvalTRIM_paths = list.files(feat_dir, full.names=T, pattern=".Rdata")
pvalTRIM_paths = pvalTRIM_paths[grepl("TRIM",pvalTRIM_paths) & grepl("pval",pvalTRIM_paths)]
pvalTRIM_labels = gsub(".Rdata","",fileNames(pvalTRIM_paths))
pvalTRIM_labels = sapply(str_split(pvalTRIM_labels,"[-]"), 
                         function(x) paste(x[3:length(x)],collapse="-") )

pvalTRIMs <- lapply(pvalTRIM_paths, function(x) get(load(x)))
names(pvalTRIMs) = pvalTRIM_labels




feat_paths = list.files(feat_dir, full.names=T, pattern=".Rdata")
feat_paths = feat_paths[!grepl("TRIM|FULL",feat_paths)]
feat_paths = feat_paths[!grepl(notrimmatrix,feat_paths)]

for (feat_path in feat_paths) {
  feat_m0 = get(load(feat_path))
  feat_mcol0 = sapply(str_split(colnames(feat_m0),"_"), function(x) x[length(x)])
  for (pvalTRIM_name in names(pvalTRIMs)) {
    feat_p = pvalTRIMs[[pvalTRIM_name]]
    feat_mcol_ind = feat_mcol0%in%colnames(feat_p)
    feat_mrow_ind = rownames(feat_m0)%in%rownames(feat_p)
    
    feat_m = feat_m0[feat_mrow_ind, feat_mcol_ind]
    
    feat_mcol = sapply(str_split(colnames(feat_m),"_"), function(x) x[length(x)])
    pis0 = which(as.matrix(feat_p)==0,arr.ind=T)
    for (pis0c in unique(pis0[,2])) {
      rowind = rownames(feat_m) %in% rownames(feat_p)[ pis0[pis0[,2]%in%pis0c,1] ]
      colind = feat_mcol %in% colnames(feat_p)[pis0c]
      feat_m[rowind, colind] = 0
    }
    
    feat_m = Matrix(feat_m, sparse=T)
    save(feat_m,file=gsub(".Rdata",paste0(".",pvalTRIM_name,".Rdata"),feat_path))
    if (writecsv) write.csv(feat_m,file=gsub(".Rdata",paste0(".",pvalTRIM_name,".csv"),feat_path))
  }
}


TimeOutput(start)





