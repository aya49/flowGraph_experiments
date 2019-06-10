## input: gates + fcm files,  meta data paths
## output: feat_file_cell_count, meta_file, meta_cell
## process: 
## - takes gates + fcm files, outputs flowtype vectors 
## - compiles flowtype vectors together to create cell count matrix
## - reformats meta file (meta info for fcm files)
## - creates cell meta file (meta info for fcm cell populations)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## input
input_dir = "/mnt/f/Brinkman group/current/Alice/gating_projects/pregnancy"
meta_file_dir = paste0(input_dir, "/meta_file.Rdata")
gs_dir = paste0(input_dir, "/gs")
channelsind_dir = paste0(input_dir, "/channels_ind.Rdata")
# gthres_dir = paste0(input_dir, "/gthres.Rdata")
gthres_dir = paste0(input_dir, "/gates_flowLearn.Rdata")
filters_dir = paste0(input_dir, "/filters.Rdata")

## ouput
result_dir = paste0(root, "/result/pregnancy"); dir.create(result_dir, showWarnings=F, recursive=T)




## libraries
source("source/_func.R")
libr(c("flowCore", "flowType", "flowDensity", "flowViz",
       "CytoML", "flowWorkspace",
       "pracma", "tools", "MASS", "KernSmooth",
       "colorRamps",
       "foreach", "doMC", "plyr", "stringr"))

## cores
no_cores = 5#detectCores()-1
registerDoMC(no_cores)


writecsv = F

meta_file = get(load(meta_file_dir))
gthres = get(load(gthres_dir))
filters = get(load(filters_dir))
channels.ind = get(load(channelsind_dir))


## flowtype -----------------------------------------
start = Sys.time()

gs = load_gs(gs_dir)


fslist = as(getData(gs),Class="list")
markers = c(# "CD11b", "CD11c", 
  "CD123", "CD14", "CD16", 
  # "CD19", 
  # "CD25",
  "CD3", "CD4","CD45", "CD45RA", "CD56", 
  "CD66", 
  "CD7", "CD8", 
  # "FoxP3", 
  # "HLADR", 
  "Tbet", "TCRgd")
gthresm = c(# "cd11b", "cd11c", 
  "cd123", "cd14.low", "cd16.gate.mid", 
  # "cd19", 
  # "cd25...",# slanted
  "cd3", "cd4","cd45", "cd45ra.cd8", "cd56.gate.mid", 
  "cd66", #.low
  "cd7", "cd8a", 
  # "foxp3...", # slanted
  # "hladr.MMDSCs", 
  "tbet.cd8t", "tcrd")
ftl = llply(1:length(fslist), function(i) {
  flowType(Frame=fslist[[i]], 
           PropMarkers=channels.ind[markers], MarkerNames=markers, 
           MaxMarkersPerPop=6, PartitionsPerMarker=2, Methods='Thresholds', 
           Thresholds=as.list(gthres[[i]][gthresm]), 
           verbose=F, MemLimit=60)
}, .parallel=F)
ft = ldply(ftl, function(ft) ft@CellFreqs)
ftcell = unlist(lapply(ftl[[1]]@PhenoCodes, function(x){return( decodePhenotype(x, markers, ftl[[1]]@PartitionsPerMarker) )}))
meta_cell = getPhen(ftcell)
ftp = foreach(xi=1:ncol(ft), .combine='cbind') %dopar% { return(ft[,xi]/ft[,1]) }
colnames(ft) = colnames(ftp) = ftcell
rownames(ft) = rownames(ftp) = meta_file$id


# compile ---------------------------
# make a control class (1st week) for meta_file class

meta_file$class[meta_file$class==1] = "control"

# for (typed in unique(meta_file$type)) {
#   typei = meta_file$type==typed
#   
#   meta_dir = paste(result_dir, "_",typed,"/meta", sep=""); dir.create(meta_dir, showWarnings=F, recursive=T)
meta_dir = paste(result_dir, "/meta", sep=""); dir.create(meta_dir, showWarnings=F, recursive=T)
# meta_file_ = meta_file[typei,!colnames(meta_file)%in%"type"]
meta_file_ = meta_file
save(meta_file_, file=paste0(meta_dir,"/file.Rdata"))
if (writecsv) write.csv(meta_file_, file=paste0(meta_dir,"/file.csv"), row.names=F)
save(meta_cell, file=paste0(meta_dir,"/cell.Rdata"))
if (writecsv) write.csv(meta_cell, file=paste0(meta_dir,"/cell.csv"), row.names=F)

feat_dir = gsub("meta","feat",meta_dir); dir.create(feat_dir, showWarnings=F, recursive=T)
feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")
feat_file_cell_prop_dir = paste(feat_dir, "/file-cell-prop", sep="")
# ft_ = ft[typei,]
# ftp_ = ftp[typei,]
ft_ = ft
ftp_ = ftp

ft_ = as.matrix(ft_)
save(ft_, file=paste0(feat_file_cell_count_dir,".Rdata"))
if (writecsv) write.csv(ft_, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)

ftp_ = as.matrix(ftp_)
save(ftp_, file=paste0(feat_file_cell_prop_dir,".Rdata"))
if (writecsv) write.csv(ftp_, file=paste0(feat_file_cell_prop_dir,".csv"), row.names=T)
# }

time_output(start)

