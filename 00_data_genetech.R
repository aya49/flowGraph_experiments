## input: genetech by daniel yokosawa on the genetech data (bone marrow + blood for 3 patients x 5 samples; computationally mixed)
## output: feat_file_cell_count, meta_file


## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)


## input directories
input_dir = "/mnt/f/Brinkman group/current/Daniel/Genentech/TAP/comp_mix/flowType/Tube_003"


## ouput directories
result_dir = paste0(root, "/result/genetech"); dir.create(result_dir, showWarnings=F, recursive=T)
feat_dir = paste0(result_dir,"/feat"); dir.create(feat_dir, showWarnings=F, recursive=T)
feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")
meta_dir = paste(result_dir, "/meta", sep=""); dir.create(meta_dir, showWarnings=F, recursive=T)
meta_file_dir = paste(meta_dir, "/file", sep="")


## libraries
source("source/_func.R")
libr(c("flowCore", "flowType",
       "foreach", "doMC", "plyr", 
       "stringr"))


## cores
no_cores = detectCores()-1
registerDoMC(no_cores)


## options
writecsv = F



start = Sys.time()


## feat/file-cell-count: load and compile flowtype count files
ftl = get(load(paste0(input_dir,"/flowType_myeloid.Rdata")))
markers = get(load(paste0(input_dir,"/MarkerNames_myeloid.Rdata")))

m0 = as.matrix(ldply(ftl, function(ft) ft@CellFreqs)[,-1])
ftcell = laply(ftl[[1]]@PhenoCodes, function(x) decodePhenotype(x, markers, ftl[[1]]@PartitionsPerMarker) )
colnames(m0) = ftcell
rownames(m0) = names(ftl)

save(m0, file=paste0(feat_file_cell_count_dir,".Rdata"))
if (writecsv) write.csv(m0, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)


## meta/file
temp_ = Reduce("rbind", str_split(rownames(m0),"_"))
meta_file = data.frame(id=rownames(m0), class=temp_[,3], patient=as.numeric(gsub(".fcs","",temp_[,4])))
meta_file$class[meta_file$class=="100%WB"] = "control"
meta_file$class = gsub("%","",meta_file$class)

save(meta_file, file=paste0(meta_file_dir,".Rdata"))
if (writecsv) write.csv(meta_file, file=paste0(meta_file_dir,".csv"), row.names=T)



time_output(start, "data_genetech")

