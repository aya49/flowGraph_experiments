# aya43@sfu.ca 20170316
# checkss for irregular values in matrices going into calculating distance

root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/matrix", sep="")
# matrix_count = c("CountAdj")
matrix_type = c("Child_entropyTRIM_CountAdj", "Child_entropyTRIM_Prop", "Parent_entropyTRIM_CountAdj", "Parent_entropyTRIM_Prop", 
                "LogFoldTRIM_CountAdj", "LogFoldTRIM_Prop", "PvalTRIM_CountAdj", "PvalTRIM_Prop",
                "CountAdj", "Parent_entropy", "Child_entropy",
                "LogFold_CountAdj", "LogFold_Prop", "Pval_CountAdj", "Pval_Prop",
                "Parent_effortTRIM_CountAdj", "Parent_effortTRIM_Prop", "Parent_contribTRIM_CountAdj", "Parent_contribTRIM_Prop", 
                "Child_pnratioTRIM_CountAdj", "Child_pnratioTRIM_Prop", "Child_propTRIM_CountAdj", "Child_propTRIM_Prop",
                "Parent_effort_CountAdj", "Parent_effort_Prop", "Parent_contrib_CountAdj", "Parent_contrib_Prop",
                "Child_pnratio", "Child_prop")

libr(Matrix)
libr(foreach)
libr(doMC)
source("~/projects/IMPC/code/_funcAlice.R")

no_cores = detectCores()-1
registerDoMC(no_cores)








start = Sys.time()


#for (ci in 4:1) {

# m0 = get(load(paste0(matrix_dir, matrix_count,".Rdata")))

#load different matrices
# for (mcp in matrix_type_) {
for (mcp in sort(matrix_type)) {
  cat("\n", mcp, " ",sep="")
  start2 = Sys.time()
  if (!file.exists(paste0(matrix_dir, mcp,".Rdata"))) {cat("doesn't exist"); next}
  mm = get(load(paste0(matrix_dir, mcp,".Rdata")))
  
  #check for irregular values in matrix
  m = mm
  naed = F
  Infed = F
  mcpneg = F
  if (!is.null(dim(m))) {
    if (sum(m[!is.na(m)]<0)>1) mcpneg = T
    cat(dim(m))
    cat(" ",checkm(m," "))
    if(mcpneg) cat(" neg; ")
  } else {
    cat(dim(m[[length(m)]])," ")
    nainfneg = foreach (mind = 1:length(m),.combine="rbind") %dopar% {
      nai = F
      infi = F
      negi = F
      checkmed = checkm(m[[mind]],paste(" ",mind))
      if (grepl("_NA",checkmed,ignore.case=T)) {nai=T}
      if (grepl("_Inf",checkmed,ignore.case=T)) {infi=T}
      if (sum(m[[mind]][!is.na(m[[mind]])]<0)>1) { negi=T }
      return(c(nai,infi,negi))
    }
    if (sum(nainfneg[,1])>0) cat (min(which(nainfneg[,1]==T)),".na ",sep="")
    if (sum(nainfneg[,2])>0) cat (min(which(nainfneg[,2]==T)),".inf ",sep="")
    if (sum(nainfneg[,3])>0) cat (min(which(nainfneg[,3]==T)),".neg ",sep="")
  }
}


TimeOutput(start)

# for(i in 1:nrow(m)) {if (sum(m[i,]==Inf)>0) cat("\n",i," ",which(m[i,]==Inf))}
