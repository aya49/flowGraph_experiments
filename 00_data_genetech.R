## input: flowtype by daniel yokosawa on the genetech data (bone marrow + blood for 3 patients)
## output: feat_file_cell_count, meta_file, meta_cell


## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)


## input
input_dir = "/mnt/f/Brinkman group/current/Daniel/Genentech/TAP/comp_mix/flowType/Tube_003"


## ouput
result_dir = paste0(root, "/result/genetech"); dir.create(result_dir, showWarnings=F, recursive=T)
feat_dir = paste0(result_dir,"/feat"); dir.create(feat_dir, showWarnings=F, recursive=T)
feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")
feat_file_cell_prop_dir = paste(feat_dir, "/file-cell-prop", sep="")
meta_dir = paste(result_dir, "/meta", sep=""); dir.create(meta_dir, showWarnings=F, recursive=T)
meta_cell_dir = paste(meta_dir, "/cell", sep="")
meta_file_dir = paste(meta_dir, "/file", sep="")


## libraries
source("source/_func.R")
libr(c("flowCore", "flowDensity", "flowType",
       "foreach", "doMC", "plyr", "stringr"))
# BiocManager::install(version="3.9")

## cores
no_cores = detectCores()-1
registerDoMC(no_cores)

writecsv = F


ftl = get(load(paste0(input_dir,"/flowType_myeloid.Rdata")))
markers = get(load(paste0(input_dir,"/MarkerNames_myeloid.Rdata")))

ft = as.matrix(ldply(ftl, function(ft) ft@CellFreqs)[,-1])
ftcell = laply(ftl[[1]]@PhenoCodes, function(x) decodePhenotype(x, markers, ftl[[1]]@PartitionsPerMarker) )
ftp = foreach(xi=1:ncol(ft), .combine='cbind') %dopar% { return(ft[,xi]/ft[,1]) }
colnames(ft) = colnames(ftp) = ftcell
rownames(ft) = rownames(ftp) = names(ftl)

ft = as.matrix(ft)
save(ft, file=paste0(feat_file_cell_count_dir,".Rdata"))
if (writecsv) write.csv(ft, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)

ftp = as.matrix(ftp)
save(ftp, file=paste0(feat_file_cell_prop_dir,".Rdata"))
if (writecsv) write.csv(ftp, file=paste0(feat_file_cell_prop_dir,".csv"), row.names=T)


## meta ------------------------------

temp_ = Reduce("rbind", str_split(rownames(ft),"_"))
meta_file_ = data.frame(id=rownames(ft), class=temp_[,3], patient=as.numeric(gsub(".fcs","",temp_[,4])))
meta_file_$class[meta_file_$class=="100%WB"] = "control"

save(meta_file_, file=paste0(meta_file_dir,".Rdata"))
if (writecsv) write.csv(meta_file_, file=paste0(meta_file_dir,".csv"), row.names=T)

meta_cell = getPhen(ftcell)
save(meta_cell, file=paste0(meta_cell_dir,".Rdata"))
if (writecsv) write.csv(meta_cell, file=paste0(meta_cell_dir,".csv"), row.names=T)


time_output(start)

