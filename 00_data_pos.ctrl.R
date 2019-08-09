## input: gates + fcm files,  meta data paths
## output: feat_file_cell_count, meta_file
## process: 
## - takes gates + fcm files (values randomly generated via exponential distribution; for class 3,4 FCS files, we increase gate values of the first 2 markers by 25% for artificial positive experiment files), outputs flowtype vectors 
## - compiles flowtype vectors together to create cell count matrix
## - reformats meta file (meta info for fcm files)
## - creates cell meta file (meta info for fcm cell populations)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## input directories
input_dir = "/mnt/f/Brinkman group/current/Alice/gating_projects/pregnancy" # make controls based off of pregnancy data set

## ouput directories
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


## options
writecsv = F
nsample = 1000 # # of samples to create per data set
nctrl = .5 # % of control samples
markern = 6 # # of markers

normean = 300000 # 301234.7 for pregnancy
normsd = 0 # 64734.05 for pregnancy


start = Sys.time()


# define number of cells in each fcs file
# for pregnancy data set, mean=301234.7, sd=64734.05
# ncells = floor(rnorm(nsample,300000,65000)) # number of cells for each sample
ncells = rnorm(nsample,normean,normsd) # number of cells for each sample

## meta/file
meta_file = data.frame(id=paste0("a",1:nsample), class="exp")
meta_file$class[1:(nctrl*nsample)] = "control"


## feat/file-count-count: flowtype

# load sample fcs file
f = read.FCS(paste0(input_dir,"/samplefcs.fcs"))

# markers
gthresm = markers = LETTERS[1:markern]

# marker indices in f@exprs
ci = c(1:markern); names(ci) = markers 

# marker thresholds
cvd = rnorm(ncells[1],2,1)
p50 = quantile(cvd, .5)
p75 = quantile(cvd, .75)
p25 = quantile(cvd, .25)
thress0 = llply(markers, function(x) p50); names(thress0) = markers
thress1 = thress2 = thress4 = thress0
thress1[[gthresm[1]]] = thress2[gthresm[1:2]] = p75
thress4[gthresm[1:4]] = quantile(cvd, .501)

save.image(paste0(root,"/temp.Rdata"))

# start = Sys.time()
for (ds in c("ctrl","pos1","pos2","pos3","pos4")) {
  # clear/load memory
  a = ls(all=T); a = a[!a%in%c("ds","root")]
  rm(list=a); gc()
  
  # start
  start2 = Sys.time()
  load(paste0(root,"/temp.Rdata"))
  setwd(root)
  registerDoMC(no_cores)
  

  # ouput directories
  result_dir = paste0(root, "/result/",ds)
  meta_dir = paste(result_dir, "/meta", sep=""); dir.create(meta_dir, showWarnings=F, recursive=T)
  feat_dir = paste(result_dir, "/feat", sep=""); dir.create(feat_dir, showWarnings=F)
  feat_file_cell_count_dir = paste(feat_dir,"/file-cell-count",sep="")

  # flowtype
  ftl = llply(loopInd(1:nsample,no_cores), function(ii) {
    llply(ii, function(i) {
      f@exprs = matrix(rnorm(ncells[i]*length(markers),2,1), nrow=ncells[i])
      thress = thress0
      if (i>(nsample*nctrl)) {
        if (ds=="pos3") {
          # .125 -> .19 a+b+c+
          f@exprs = matrix(rnorm(ncells[i]*length(markers),2,1), nrow=ncells[i])
          ap = f@exprs[,1]>thress[[1]]
          bp = f@exprs[,2]>thress[[2]]
          cp = f@exprs[,3]>thress[[3]]
          ab = (ap & bp)
          ac = (ap & cp)
          bc = (bp & cp)
          triple = (ap & bp & cp)
          
          all = ap | bp | cp
          non = !all
          doublea = bc & cp & !ap
          doubleb = ap & cp & !bp
          doublec = ap & bp & !cp
          singlea = ap & !bp & !cp
          singleb = bp & !ap & !cp
          singlec = cp & !ab & !bp
          
          # single = singlea | singleb | singlec
          # double = (ab | bc | ac) & !triple
          tm = sum(triple)/2
          f@exprs[sample(which(doublea),tm),1] = p75
          f@exprs[sample(which(doubleb),tm),1] = p25
          f@exprs[sample(which(doublec),tm),1] = p25
          f@exprs[sample(which(non),tm),1] = p75
          f@exprs[sample(which(doublea),tm),2] = p25
          f@exprs[sample(which(doubleb),tm),2] = p75
          f@exprs[sample(which(doublec),tm),2] = p25
          f@exprs[sample(which(non),tm),2] = p75
          f@exprs[sample(which(doublea),tm),3] = p25
          f@exprs[sample(which(doubleb),tm),3] = p25
          f@exprs[sample(which(doublec),tm),3] = p75
          f@exprs[sample(which(non),tm),3] = p75
          
         ft = flowType(Frame=f, PropMarkers=ci, MarkerNames=markers, 
                   MaxMarkersPerPop=4, PartitionsPerMarker=2, 
                   Thresholds=thress, 
                   Methods='Thresholds', verbose=F, MemLimit=60)
         ftcell = unlist(lapply(ft@PhenoCodes, function(x)
           decodePhenotype(x, ft@MarkerNames, rep(2,length(ft@MarkerNames))) ))
         ftv = ft@CellFreqs
         ftv = round(ftv/ftv[1],3)
         names(ftv) = ftcell
         a = getPhenCP(cp=ftcell,no_cores=no_cores)
         al = layout_gr(a$gr$e,a$gr$v)
         
        } else if (ds!="ctrl") {
          thress = switch(ds, 
                          pos1 = thress1,
                          pos2 = thress2,
                          pos4 = thress4)
        }
      }
      flowType(Frame=f, PropMarkers=ci, MarkerNames=markers, 
               MaxMarkersPerPop=6, PartitionsPerMarker=2, 
               Thresholds=thress, 
               Methods='Thresholds', verbose=F, MemLimit=60)
    })
  }, .parallel=T)
  ftl = unlist(ftl,recursive=F)
  ft = ldply(ftl, function(ft) ft@CellFreqs)[,-1]
  ftcell = unlist(lapply(ftl[[1]]@PhenoCodes, function(x){return( decodePhenotype(x, markers, ftl[[1]]@PartitionsPerMarker) )}))
  colnames(ft) = ftcell
  rownames(ft) = meta_file$id
  
  save(meta_cell, file=paste0(meta_dir,"/cell.Rdata"))
  if (writecsv) write.csv(meta_cell, file=paste0(meta_dir,"/cell.csv"), row.names=F)
  
  # compile ---------------------------
  ft_ = as.matrix(ft)
  save(ft_, file=paste0(feat_file_cell_count_dir,".Rdata"))
  if (writecsv) write.csv(ft_, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)

  time_output(start2, ds)
}
time_output(start)

