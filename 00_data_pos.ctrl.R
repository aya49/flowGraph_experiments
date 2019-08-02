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
# meta_file_dir = paste0(input_dir, "/meta_file.Rdata")
gs_dir = paste0(input_dir, "/gs")
# channelsind_dir = paste0(input_dir, "/channels_ind.Rdata")
# gthres_dir = paste0(input_dir, "/gthres.Rdata")
gthres_dir = paste0(input_dir, "/gates_flowLearn.Rdata")
# filters_dir = paste0(input_dir, "/filters.Rdata")

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
no_cores = 10#detectCores()-1
registerDoMC(no_cores)


writecsv = F
nsample = 1000
nctrl = .5

# meta_file = get(load(meta_file_dir))
gthres = get(load(gthres_dir))
# filters = get(load(filters_dir))
# channels.ind = get(load(channelsind_dir))



## flowtype -----------------------------------------
start = Sys.time()

if (!exists("gs")) gs = load_gs(gs_dir)
fslist = as(getData(gs),Class="list")
f = fslist[[1]]

# number of cells in each sample
ncell0 = laply(fslist, function(x) nrow(x@exprs))
ncells = floor(rnorm(nsample,mean(ncell0),sd(ncell0))) # number of cells for each sample
meta_file = data.frame(id=paste0("a",1:nsample),class="exp")
meta_file$class[1:(nctrl*nsample)] = "control"
rm(fslist,gs); gc()

# markers
markers = c(# "CD11b", "CD11c", 
  "CD123", "CD14", "CD16", 
  # "CD19", 
  # "CD25",
  # "CD3", "CD4","CD45", "CD45RA", "CD56", 
  "CD66", 
  "CD7", "CD8"#, 
  # "FoxP3", 
  # "HLADR", 
  #"Tbet", "TCRgd"
  )
gthresm = c(# "cd11b", "cd11c", 
  "cd123", "cd14.low", "cd16.gate.mid", 
  # "cd19", 
  # "cd25...",# slanted
  # "cd3", "cd4","cd45", "cd45ra.cd8", "cd56.gate.mid", 
  "cd66", #.low
  "cd7", "cd8a"#, 
  # "foxp3...", # slanted
  # "hladr.MMDSCs", 
  #"tbet.cd8t", "tcrd"
  )

# marker thresholds
cvd = rnorm(nrow(f@exprs),2,1)
thress = as.list(gthres[[1]][gthresm])
for (jj in 1:length(markers)) 
  thress[[gthresm[jj]]] = mean(cvd)
thress0 = thress1 = thress2 = thress4 = thress
thress1[[gthresm[1]]] = thress2[gthresm[1:2]] = quantile(cvd, .75)
thress4[gthresm[1:4]] = quantile(cvd, .52)

# marker indices in f@exprs
ci = c(1:length(markers)); names(ci) = markers 

save.image(paste0(root,"/temp.Rdata"))

# start = Sys.time()
for (ds in c("ctrl","pos1","pos2","pos4")) {
  start2 = Sys.time()
  # clear/load memory
  a = ls(all=T); a = a[!a%in%c("ds","root","start")]
  rm(list=a); gc()
  load(paste0(root,"/temp.Rdata"))
  setwd(root)
  registerDoMC(no_cores)
  
  ## ouput
  result_dir = paste0(root, "/result/",ds)
  meta_dir = paste(result_dir, "/meta", sep=""); dir.create(meta_dir, showWarnings=F, recursive=T)
  feat_dir = paste(result_dir, "/feat", sep=""); dir.create(feat_dir, showWarnings=F)
  feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")
  feat_file_cell_prop_dir = paste(feat_dir, "/file-cell-prop", sep="")
  
  ftl = llply(loopInd(1:nsample,no_cores), function(ii) {
    llply(ii, function(i) {
      f@exprs = matrix(rnorm(ncells[i]*length(markers),2,1), nrow=ncells[i])
      thress = thress0
      if (ds!="ctrl" & i>(nsample*nctrl)) {
        thress = switch(ds, 
                        pos1 = thress1,
                        pos2 = thress2,
                        pos4 = thress4)
      }
      flowType(Frame=f, PropMarkers=ci, MarkerNames=markers, 
               MaxMarkersPerPop=6, PartitionsPerMarker=2, 
               Thresholds=thress, 
               Methods='Thresholds', verbose=F, MemLimit=60)
    })
  }, .parallel=T)
  ftl = unlist(ftl,recursive=F)
  ft = ldply(ftl, function(ft) ft@CellFreqs)[,-1]
  ftcell = unlist(lapply(ftl[[1]]@PhenoCodes, function(x){return( decodePhenotype(x, LETTERS[1:length(markers)], ftl[[1]]@PartitionsPerMarker) )}))
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
  
  time_output(start2, ds)
}
time_output(start)

