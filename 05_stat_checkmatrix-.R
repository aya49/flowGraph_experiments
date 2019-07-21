## Input: original features --> Output: prints whether there are irregular values in features
# aya43@sfu.ca 20170316

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr("Matrix")

for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  # result_dir = paste0(root, "/result/impc_panel1_sanger-spleen") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  cat("\n\n",result_dir)
  
  ## input directories
  feat_dir = paste(result_dir, "/feat", sep="")
  
  
  ## feat_types
  feat_types = gsub(".Rdata","",list.files(path=feat_dir,pattern=".Rdata"))
  
  ## loop through feat and check for irregular values
  start = Sys.time()
  
  for (feat_type in sort(feat_types)) {
    m = get(load(paste0(feat_dir, "/", feat_type,".Rdata")))
    
    options(sep="")
    pl = paste0(
    "\n", feat_type,": ",
    ", dim ", paste0(dim(m), collapse="x"),
    ", NA ", sum(is.na(m)),
    # cat(", NaN ", sum(is.nan(m)))
    ", Inf ", sum(m==Inf),
    ", neg ", sum(m<0),
    ", 0/+/- ", sum(m==0),"/",sum(m>0),"/",sum(m<0),
    ", neg ", sum(m<0),
    ", max ", max(m[is.finite(m)]))
    # cat(", rowname ", rownames(m)[1])
    # cat(", colname ", colnames(m)[10])
    # cat(", ",checkm(m," "))
    # if(feat_typeneg) cat(" neg ")
    # cat("rownames; ")
    
    cat(pl)
    
    png(paste0(feat_dir, "/", feat_type,".png"), width=800)
    hist(m[is.finite(m)], main=pl)
    graphics.off()
  }
}