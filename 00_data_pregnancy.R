## input: gates + fcm files,  meta data paths
## output: feat_file_cell_count, meta_file
## process: 
## - takes gates + fcm files, outputs flowtype vectors 
## - compiles flowtype vectors together to create cell count matrix
## - reformats meta file (meta info for fcm files)
## - creates cell meta file (meta info for fcm cell populations)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## input directories
input_dir = "/mnt/f/Brinkman group/current/Alice/gating_projects/pregnancy"
meta_file_dir = paste0(input_dir, "/meta_file.Rdata")
flowtype_dir = paste0(input_dir, "/flowtype")

## ouput directories
result_dir0 = paste0(root, "/result/pregnancy"); dir.create(result_dir0, showWarnings=F, recursive=T)


## libraries
source("source/_func.R")
libr(c("flowCore", "flowType", "CytoML", "flowWorkspace",
       "foreach", "doMC", "plyr"))

## cores
no_cores = detectCores()-1
registerDoMC(no_cores)


## options
writecsv = F



start = Sys.time()


## prepare flowtype directories & load meta data
ft_dirs = list.files(flowtype_dir, full.names=T)
ft_names = gsub(".Rdata|.fcs|Gates_|_Unstim|_Repeat","",fileNames(ft_dirs))
ft_names = gsub("BL","4",ft_names)

meta_file0 = get(load(meta_file_dir))
colnames(meta_file0)[1] = "subject"

fto = match(meta_file0$id,ft_names)


## feat/file-cell-count: load and compile flowtype count files
m00 = llply(loopInd(ft_dirs[fto],no_cores), function(ii) {
  llply(ii, function(i) get(load(i))@CellFreqs)
}, .parallel=T)
m00 = as.matrix(Reduce('rbind',llply(m00,function(x)Reduce(rbind,x))))
ft = get(load(ft_dirs[1]))
ftcell = unlist(lapply(ft@PhenoCodes, function(x)
  decodePhenotype(x, ft@MarkerNames, rep(2,length(ft@MarkerNames))) ))
colnames(m00) = ftcell
rownames(m00) = meta_file0$id


meta_file0$class[meta_file0$class==4] = "control"
meta_file0$train = ifelse(meta_file0$type=="train",T,F)
meta_file0 = meta_file0[,-4]


## save

# output directories
result_dir = result_dir0
meta_dir = paste(result_dir, "/meta", sep=""); dir.create(meta_dir, showWarnings=F, recursive=T)
feat_dir = gsub("meta","feat",meta_dir); dir.create(feat_dir, showWarnings=F, recursive=T)

meta_file = meta_file0
save(meta_file, file=paste0(meta_dir,"/file.Rdata"))
if (writecsv) write.csv(meta_file, file=paste0(meta_dir,"/file.csv"), row.names=F)

feat_dir = gsub("meta","feat",meta_dir); dir.create(feat_dir, showWarnings=F, recursive=T)
feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")

m0 = m00
save(m00, file=paste0(feat_file_cell_count_dir,".Rdata"))
if (writecsv) write.csv(m00, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)


time_output(start, "data_pregnancy")

