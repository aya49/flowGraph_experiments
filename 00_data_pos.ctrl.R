## input: gates + fcm files,  meta data paths
## output: feat_file_cell_count, meta_file, meta_cell
## process: 
## - takes gates + fcm files (values randomly generated via exponential distribution; for class 3,4 FCS files, we increase gate values of the first 2 markers by 25% for artificial positive experiment files), outputs flowtype vectors 
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
# see for loop below



## libraries
source("source/_func.R")
libr(c("flowCore", "flowType", "flowDensity", "flowViz",
       "CytoML", "flowWorkspace",
       "pracma", "tools", "MASS", "KernSmooth", "fitdistrplus",
       "colorRamps",
       "foreach", "doMC", "plyr", "stringr"))

## cores
no_cores = 8#detectCores()-1
registerDoMC(no_cores)


writecsv = F

meta_file = get(load(meta_file_dir))
gthres = get(load(gthres_dir))
filters = get(load(filters_dir))
channels.ind = get(load(channelsind_dir))

# make a control class and experiment class
meta_file$class[meta_file$class%in%c(1,2)] = "control"
meta_file$class[meta_file$class%in%c(3,4)] = "exp"



## flowtype -----------------------------------------
start = Sys.time()

if (!exists("gs")) gs = load_gs(gs_dir)
fslist = as(getData(gs),Class="list")
f = fslist[[1]]

markers = c(# "CD11b", "CD11c", 
  "CD123", "CD14", "CD16", 
  # "CD19", 
  # "CD25",
  # "CD3", "CD4","CD45", "CD45RA", "CD56", 
  "CD66", 
  "CD7", "CD8", 
  # "FoxP3", 
  # "HLADR", 
  "Tbet", "TCRgd")
gthresm = c(# "cd11b", "cd11c", 
  "cd123", "cd14.low", "cd16.gate.mid", 
  # "cd19", 
  # "cd25...",# slanted
  # "cd3", "cd4","cd45", "cd45ra.cd8", "cd56.gate.mid", 
  "cd66", #.low
  "cd7", "cd8a", 
  # "foxp3...", # slanted
  # "hladr.MMDSCs", 
  "tbet.cd8t", "tcrd")
thress = as.list(gthres[[1]][gthresm])
f = fslist[[1]]
for (jj in 1:length(markers)) {
  j = channels.ind[markers[jj]]
  # a = 
  # a = a-min(a)
  # a = a*(max(f@exprs[,j])-max(a))
  f@exprs[,j] = rnorm(nrow(f@exprs),2,1)
  thress[[gthresm[jj]]] = mean(f@exprs[,j])
}
thress1 = thress2 = thress
for (jj in 1:2) {
  j = channels.ind[markers[jj]]
  thress2[[gthresm[jj]]] = quantile(f@exprs[,j], .75)
}

# start = Sys.time()
for (ds in c("ctrl","pos")) {
  ## ouput
  result_dir = paste0(root, "/result/",ds)
  meta_dir = paste(result_dir, "/meta", sep=""); dir.create(meta_dir, showWarnings=F, recursive=T)
  feat_dir = paste(result_dir, "/feat", sep=""); dir.create(feat_dir, showWarnings=F)
  feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")
  feat_file_cell_prop_dir = paste(feat_dir, "/file-cell-prop", sep="")
  
  
  if (ds=="pos") {
    ftl = llply(1:length(fslist), function(i) {
      ii = names(fslist)[i]
      f = fslist[[i]]
      if (grepl("_1_|_2_",ii)) {
        thress = thress1
      } else {
        thress = thress2
      }
      for (jj in 1:length(markers)) {
        j = channels.ind[markers[jj]]
        f@exprs[,j] = rexp(nrow(f@exprs))
      }
      flowType(Frame=f, 
               PropMarkers=channels.ind[markers], MarkerNames=markers, 
               MaxMarkersPerPop=6, PartitionsPerMarker=2, Methods='Thresholds', 
               Thresholds=thress, 
               verbose=F, MemLimit=60)
    }, .parallel=T)
  } else if (ds=="ctrl") {
    ftl = llply(1:length(fslist), function(i) {
      f = fslist[[i]]
      for (jj in 1:length(markers)) {
        j = channels.ind[markers[jj]]
        f@exprs[,j] = rexp(nrow(f@exprs))
      }
      flowType(Frame=f, 
               PropMarkers=channels.ind[markers], MarkerNames=markers, 
               MaxMarkersPerPop=6, PartitionsPerMarker=2, Methods='Thresholds', 
               Thresholds=thress, 
               verbose=F, MemLimit=60)
    }, .parallel=T)
  }
  ft = ldply(ftl, function(ft) ft@CellFreqs)
  ftcell = unlist(lapply(ftl[[1]]@PhenoCodes, function(x){return( decodePhenotype(x, markers, ftl[[1]]@PartitionsPerMarker) )}))
  meta_cell = getPhen(ftcell)
  ftp = foreach(xi=1:ncol(ft), .combine='cbind') %dopar% { return(ft[,xi]/ft[,1]) }
  colnames(ft) = colnames(ftp) = ftcell
  rownames(ft) = rownames(ftp) = meta_file$id
  
  save(meta_file, file=paste0(meta_dir,"/file.Rdata"))
  if (writecsv) write.csv(meta_file, file=paste0(meta_dir,"/file.csv"), row.names=F)
  save(meta_cell, file=paste0(meta_dir,"/cell.Rdata"))
  if (writecsv) write.csv(meta_cell, file=paste0(meta_dir,"/cell.csv"), row.names=F)
  
  # compile ---------------------------
  ft_ = as.matrix(ft)
  save(ft_, file=paste0(feat_file_cell_count_dir,".Rdata"))
  if (writecsv) write.csv(ft_, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)
  
  ftp_ = as.matrix(ftp)
  save(ftp_, file=paste0(feat_file_cell_prop_dir,".Rdata"))
  if (writecsv) write.csv(ftp_, file=paste0(feat_file_cell_prop_dir,".csv"), row.names=T)
}
time_output(start)

