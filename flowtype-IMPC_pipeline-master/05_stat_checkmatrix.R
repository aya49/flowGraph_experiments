## Input: original features --> Output: prints whether there are irregular values in features
# aya43@sfu.ca 20170316

## root directory
root = "~/projects/IMPC"
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#,"Sanger_MLN","CIPHE","TCP","H")
# controlL = c("+_+|+_Y","+_+|+_Y","WildType","WildType","WildType") #control value in target_col column
ci = 1; panel = panelL[ci]; centre = centreL[ci]

result_dir = paste0("result/", panelL, "/", centreL); suppressWarnings(dir.create (result_dir))

## input directories
feat_dir = paste(result_dir, "/feat", sep="")

## libraries
source("~/projects/IMPC/code/_funcAlice.R")
libr("Matrix")

## feat_types
feat_types = gsub(".Rdata","",list.files(path=feat_dir,pattern=".Rdata"))

## loop through feat and check for irregular values
start = Sys.time()

for (feat_type in sort(feat_types)) {
  start2 = Sys.time()
  m = get(load(paste0(feat_dir, "/", feat_type,".Rdata")))
  
  options(sep="")
  cat("\n\n", feat_type,": ")
  cat("\nNA ", sum(is.na(m)))
  # cat(", NaN ", sum(is.nan(m)))
  cat(", Inf ", sum(m==Inf))
  cat(", dim ", dim(m))
  cat(", neg ", sum(m<0))
  cat(", 0/+/- ", sum(m==0),"/",sum(m>0),"/",sum(m<0))
  cat(", rowname ", rownames(m)[1])
  cat(", colname ", colnames(m)[10])
  # cat(", ",checkm(m," "))
  # if(feat_typeneg) cat(" neg ")
  # cat("rownames; ")
}
TimeOutput(start)
