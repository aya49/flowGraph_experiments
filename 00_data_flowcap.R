## input: flowtype file,  meta data paths
## output: feat_file_cell_count, meta_file
## process: 
## - takes flowtype output files and compiles them together to create cell count matrix
## - reformats meta file (meta info for fcm files)
## - creates cell meta file (meta info for fcm cell populations)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)


## input directories
data_dir = "/mnt/f/Brinkman group/current/Alice/flowCAP-II/data" #main data directory
fcs_dir = paste0(data_dir,"/FCS") #fcs file directory
ft_dir = paste0(data_dir,"/FT") #flowtype file directory
csv_dir = paste0(data_dir,"/AML.csv") #meta file directory
csv_dir2 = paste0(data_dir,"/AMLTraining.csv") #meta file directory


## output directories (NOT last section, output split by tube/panel)
result_dir0 = paste0(root, "/result/flowcap") #; dir.create(result_dir, showWarnings=F, recursive=T)


## libraries
source("source/_func.R")
libr(c("flowCore", "flowType",
       "doMC", "foreach", "plyr","stringr"))


## cores
no_cores = 5#detectCores()-1
registerDoMC(no_cores)

## options
writecsv = F
options(stringsAsFactors=F)



start = Sys.time()


## prepare flowtype directories & load meta data
ft_dirs = sort( dir(ft_dir, pattern=".Rda", all.files=T, full.names=T, recursive=T) )
ftnames = fileNames(ft_dirs, "Rda")

meta_filetemp = data.frame(read.csv(csv_dir)) # meta file
meta_file_trt = read.csv(csv_dir2) # meta file  with with training/testing


## feat/file-cell-count: load and compile flowtype count files
m00 = llply(loopInd(1:length(ft_dirs),no_cores),function(ii){
  llply(ii, function(i) get(load(ft_dirs[i]))@CellFreqs)
}, .parallel=T)
m00 = as.matrix(Reduce('rbind',llply(m00,function(x)Reduce(rbind,x))))
rownames(m00) = ftnames
colnames(m00) = rownames(get(load(ft_dirs[1]))@MFIs)
markers_ = c("FS","SS",paste0("FL",1:5)) #unique(unlist(str_split(colnames(m00),"[+-]")))

## meta/file
tube_subject = ldply(1:length(ft_dirs), function(i) {
  fn = strsplit(ftnames[i], "S")[[1]]
  c(substr(fn[1],2,nchar(fn[1])), gsub("FT","",fn[2]))
})
meta_filetemp$train = ifelse(is.na(meta_file_trt$Label),F,T)
meta_file0 = meta_filetemp[match(
  apply(tube_subject,1,function(x) paste0(x,collapse=" ")), 
  apply(meta_filetemp[,c(2,3)],1,function(x) paste0(x,collapse=" "))),]
meta_file0$id = ftnames
meta_file0 = meta_file0[,c("Label","id","train")]
meta_file0$subject = as.numeric(gsub("T[0-9]S|FT","",meta_file0$id))
colnames(meta_file0)[1] = "class"
meta_file0$class[meta_file0$class=="normal"] = "control"
meta_file0$tube = as.numeric(substr(meta_file0$id,2,2))
meta_file0 = meta_file0[match(rownames(m00),meta_file0$id),]


## save: split data by tube/panel
# controli = NULL
for (tube in unique(meta_file0$tube)) {
  # if (tube!=6) next # save only tube 6 for now
  
  # output directories
  result_dir = paste0(result_dir0,"_",tube)
  meta_dir = paste(result_dir, "/meta", sep=""); dir.create(meta_dir, showWarnings=F, recursive=T)
  feat_dir = gsub("meta","feat",meta_dir); dir.create(feat_dir, showWarnings=F, recursive=T)
  
  # prepare tube sample indices
  tubei = meta_file0$tube==tube
  
  meta_file = meta_file0[tubei,]
  meta_file$id = as.numeric(gsub("T[0-9]S|FT","",meta_file$id))
  
  # # randomly pick normal patients as controls
  # if (is.null(controli)) {
  #   normn = sum(meta_file_$class=="normal")
  #   amln = sum(meta_file_$class=="aml")
  #   controli = grep("normal", meta_file_$class)[sample(normn, normn-amln)]
  # }
  # meta_file_$class[controli] = "control"
  
  meta_file_dir = paste(meta_dir, "/file", sep="")
  save(meta_file, file=paste0(meta_file_dir,".Rdata"))
  if (writecsv) write.csv(meta_file, file=paste0(meta_file_dir,".csv"), row.names=F)
  
  f = read.FCS(paste0(fcs_dir,"/",str_pad(meta_file$id[1], 4, pad="0"),".fcs"))
  markers = str_extract(f@parameters@data$desc,"[A-Za-z0-9]+")
  m0 = m00[tubei,]
  rownames(m0) = meta_file$id
  #colnames(m0) = laply(ft@PhenoCodes, function(x) decodePhenotype(x, markers, ft@PartitionsPerMarker) )
  for (i in 1:length(markers_))
    colnames(m0) = gsub(markers_[i],markers[i],colnames(m0))
  

  feat_file_cell_count_dir = paste(feat_dir,"/file-cell-count",sep="")
  save(m0, file=paste0(feat_file_cell_count_dir,".Rdata"))
  if (writecsv) write.csv(m0, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)
}


time_output(start, "data_flowcap")


