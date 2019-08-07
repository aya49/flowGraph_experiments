## input: flowtype by daniel yokosawa on the genetech data (bone marrow + blood for 3 patients)
## output: feat_file_cell_count, meta_file, meta_cell


## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)


## input
input_dir = "/mnt/f/Brinkman group/current/Alice/gating_projects/HDCytoData_Bodenmiller"
data_dir = paste0(input_dir,"/data")
gate_dir = paste0(input_dir,"/gates")
fcs_dir = paste0(input_dir,"/fcs")
gate_dir = paste0(input_dir,"/gates")


## ouput
result_dir = paste0(root, "/result/bodenmiller"); dir.create(result_dir, showWarnings=F, recursive=T)
feat_dir = paste0(result_dir,"/feat"); dir.create(feat_dir, showWarnings=F, recursive=T)
feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")
feat_file_cell_prop_dir = paste(feat_dir, "/file-cell-prop", sep="")
meta_dir = paste(result_dir, "/meta", sep=""); dir.create(meta_dir, showWarnings=F, recursive=T)
meta_cell_dir = paste(meta_dir, "/cell", sep="")
meta_file_dir = paste(meta_dir, "/file", sep="")


## libraries
source("source/_func.R")
libr(c("HDCytoData", # requires bioconductor 3.9 which requires R 3.6
       "flowCore", "flowDensity", "diffcyt", "flowType",
       "MASS", "RColorBrewer",
       "foreach", "doMC", "plyr", "stringr"))
# BiocManager::install(version="3.9")

## cores
no_cores = detectCores()-1
registerDoMC(no_cores)

writecsv = F


ftl = get(load(load("/mnt/f/Brinkman group/current/Daniel/Genentech/TAP/comp_mix/flowType/Tube_003/flowType_myeloid.Rdata")))

ft = ldply(ftl, function(ft) ft@CellFreqs)
ftcell = unlist(lapply(ftl[[1]]@PhenoCodes, function(x){return(decodePhenotype(x, markers, ftl[[1]]@PartitionsPerMarker) )}))
ftp = foreach(xi=1:ncol(ft), .combine='cbind') %dopar% { return(ft[,xi]/ft[,1]) }
colnames(ft) = colnames(ftp) = ftcell
rownames(ft) = rownames(ftp) = names(fslist)

ft = as.matrix(ft)
save(ft, file=paste0(feat_file_cell_count_dir,".Rdata"))
if (writecsv) write.csv(ft, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)

ftp = as.matrix(ftp)
save(ftp, file=paste0(feat_file_cell_prop_dir,".Rdata"))
if (writecsv) write.csv(ftp, file=paste0(feat_file_cell_prop_dir,".csv"), row.names=T)


## meta ------------------------------

meta_file_ = meta_file[!duplicated(meta_file$sample_id),-4]
for (ci in 1:ncol(meta_file_)) 
  meta_file_[,ci] = as.character(meta_file_[,ci])
meta_file_ = as.data.frame(meta_file_)
colnames(meta_file_) = c("class","patient","id")
meta_file_$class = gsub("reference","control",meta_file_$class,ignore.case=T)

save(meta_file_, file=paste0(meta_file_dir,".Rdata"))
if (writecsv) write.csv(meta_file_, file=paste0(meta_file_dir,".csv"), row.names=T)

meta_cell = getPhen(ftcell)
save(meta_cell, file=paste0(meta_cell_dir,".Rdata"))
if (writecsv) write.csv(meta_cell, file=paste0(meta_cell_dir,".csv"), row.names=T)


time_output(start)

