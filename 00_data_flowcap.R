## input: flowtype file,  meta data paths
## output: feat_file_cell_count, meta_file, meta_cell
## process: 
## - takes flowtype output files and compiles them together to create cell count matrix
## - reformats meta file (meta info for fcm files)
## - creates cell meta file (meta info for fcm cell populations)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metrics"
setwd(root)

result_dir = paste0(root, "/result/flowcap") #; dir.create(result_dir, showWarnings=F, recursive=T)
# result_dir = "results"
# suppressWarnings(dir.create (result_dir, recursive=T))
data_dir = "/mnt/f/Brinkman group/current/Alice/flowCAP-II/data" #main data directory



## input directories
ft_dir = paste0(data_dir,"/FT") #flowtype file directory
csv_dir = paste0(data_dir,"/AML.csv") #meta file directory
csv_dir2 = paste0(data_dir,"/AMLTraining.csv") #meta file directory

## output directories (see last section, output split by tube/panel)


## libraries
source("source/_funcAlice.R")
libr(c("flowCore", "flowType",
       "foreach", "doMC"))

## cores
no_cores = detectCores()-1
registerDoMC(no_cores)

## options
writecsv = F
options(stringsAsFactors=F)
options(device="cairo")
countThres = 0 #delete columns/rows where all values equal or below countThres
levelThres = 8 #Inf if delete no layers; >levelcutoff are deleted



start = Sys.time()

## prepare directories
ftFile_dir = sort( dir(ft_dir, pattern=".Rda", all.files=T, full.names=T, recursive=T) )
ftGT = folderNames(ftFile_dir)
ftFileNames = fileNames(ftFile_dir, "Rda")

## get meta data on files for meta_file
meta_filetemp = data.frame(read.csv(csv_dir))
colnames(meta_filetemp) = c("id", "tube", "specimen", "class") #rename columns

## get meta data on phenotypes (cell populations) for meta_cell
ft = get(load(ftFile_dir[1]))
markers = ft@MarkerNames
# save(markers, file=paste0(markers_dir,".Rdata"))

meta_cell = getPhen(rownames(ft@MFIs))

## load flowtype files for feat_file_cell_count & fill in meta_file
feat_file_cell_count = NULL
meta_file = NULL

result = foreach (i=1:length(ftFile_dir), .combine="rbind") %dopar% {
  ft = get(load(ftFile_dir[i]))
  print(ft@MarkerNames)
  npheno = length(ft@CellFreqs)
  # feat_file_cell_count = rbind(feat_file_cell_count, ft@CellFreqs)
  cat("\n", i, "Loaded file: ", ftGT[i], "/", ftFileNames[i], "; # of phenotypes x samples: ", npheno, " x ", length(ft), "", sep="")
  
  fn = ftFileNames[i]
  fn = strsplit(fn, "\\.")[[1]][1]
  fn = strsplit(fn, "S")[[1]]
  tube = substr(fn[1],2,nchar(fn[1]))
  specimen = gsub("FT","",fn[2])
  # meta_file = rbind(meta_file, meta_filetemp[which(as.numeric(meta_filetemp$specimen)==as.numeric(specimen) & as.numeric(meta_filetemp$tube)==as.numeric(tube)),])
  return(c(meta_filetemp[which(as.numeric(meta_filetemp$specimen)==as.numeric(specimen) & as.numeric(meta_filetemp$tube)==as.numeric(tube)),], ft@CellFreqs))
  
  rm(ft) #save memory
}

#collate feat_file_cell_count
feat_file_cell_count = as.matrix(apply(result[,(ncol(meta_filetemp)+1):ncol(result)], 2, as.numeric))
sm = result[,1:(ncol(meta_filetemp))]

#order meta_cell and feat_file_cell_count
pheno_order = order(meta_cell$phenocode)
meta_cell = meta_cell[pheno_order,]
feat_file_cell_count = feat_file_cell_count[,pheno_order]

#make meta_file
rownames(feat_file_cell_count) = ftFileNames
meta_file = as.data.frame(ftFileNames)
meta_file$tube = unlist(sm[,2])
meta_file$specimen = unlist(sm[,3])
meta_file$class = unlist(sm[,4])
colnames(meta_file) = colnames(meta_filetemp)
meta_file_trt = read.csv(csv_dir2) # meta file with training/testing, not used.

rownames(feat_file_cell_count) = meta_file[,1] = ftFileNames
colnames(feat_file_cell_count) = meta_cell$phenotype
feat_file_cell_prop = feat_file_cell_count/feat_file_cell_count[,1]
dimnames(feat_file_cell_prop) = dimnames(feat_file_cell_count)

#rename classes, so there is a control group
# meta_file$class[meta_file$class=="normal"] = "control"



## trim matrix/meta_cell ------------------------------
rowIndex = apply(feat_file_cell_count, 1, function(x) any(x > countThres)) #delete rows of all 0 or too little count
colIndex1 = apply(feat_file_cell_count, 2, function(x) any(x > countThres)) #delete cols of all 0 or too little count
colIndex2 = meta_cell$phenolevel <= levelThres #delete cols of too high level
colIndex = colIndex1 & colIndex2

feat_file_cell_count <- feat_file_cell_count[rowIndex,colIndex]
feat_file_cell_prop <- feat_file_cell_prop[rowIndex,colIndex]
meta_cell <- meta_cell[colIndex,]
meta_file <- meta_file[rowIndex,]


## split data by tube/panel and save; also randomly pick same controls to be a non-control class so to make features-------------------------
controli = NULL
for (tube in unique(meta_file$tube)) {
  tubei = meta_file$tube==tube
  
  meta_file_ = meta_file[tubei,]
  tubeorder = order(meta_file_$specimen)
  meta_file_ = meta_file_[tubeorder,c("specimen","class"), drop=F]
  colnames(meta_file_)[1] = "id"
  
  # # randomly pick normal patients as controls
  # if (is.null(controli)) {
  #   normn = sum(meta_file_$class=="normal")
  #   amln = sum(meta_file_$class=="aml")
  #   controli = grep("normal", meta_file_$class)[sample(normn, normn-amln)]
  # }
  # meta_file_$class[controli] = "control"
  
  meta_dir = paste(result_dir, "_panel", tube, "/meta", sep=""); dir.create(meta_dir, showWarnings=F, recursive=T)
  meta_cell_dir = paste(meta_dir, "/cell", sep="")
  meta_file_dir = paste(meta_dir, "/file", sep="")
  save(meta_cell, file=paste0(meta_cell_dir,".Rdata"))
  if (writecsv) write.csv(meta_cell, file=paste0(meta_cell_dir,".csv"), row.names=F)
  save(meta_file_, file=paste0(meta_file_dir,".Rdata"))
  if (writecsv) write.csv(meta_file_, file=paste0(meta_file_dir,".csv"), row.names=F)
  
  feat_file_cell_count_ = feat_file_cell_count[tubei,]
  feat_file_cell_prop_ = feat_file_cell_prop[tubei,]
  
  feat_file_cell_count_ = feat_file_cell_count_[tubeorder,]
  feat_file_cell_prop_ = feat_file_cell_prop_[tubeorder,]
  rownames(feat_file_cell_prop_) = rownames(feat_file_cell_count_) = meta_file_$id
  
  feat_dir = gsub("meta","feat",meta_dir); dir.create(feat_dir, showWarnings=F, recursive=T)
  feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")
  feat_file_cell_prop_dir = paste(feat_dir, "/file-cell-prop", sep="")
  
  feat_file_cell_count_ = as.matrix(feat_file_cell_count_)
  save(feat_file_cell_count_, file=paste0(feat_file_cell_count_dir,".Rdata"))
  if (writecsv) write.csv(feat_file_cell_count_, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)
  
  feat_file_cell_prop_ = as.matrix(feat_file_cell_prop_)
  save(feat_file_cell_prop_, file=paste0(feat_file_cell_prop_dir,".Rdata"))
  if (writecsv) write.csv(feat_file_cell_prop_, file=paste0(feat_file_cell_prop_dir,".csv"), row.names=T)
}

time_output(start)


