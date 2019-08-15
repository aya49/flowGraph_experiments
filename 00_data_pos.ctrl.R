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
no_cores = detectCores()-1
registerDoMC(no_cores)


## options
writecsv = F
nsample = 1000 # # of samples to create per data set
nctrl = .5 # % of control samples
markern = 6 # # of markers

normean = 300000 # 301234.7 for pregnancy
normsd = 0 # 64734.05 for pregnancy; if >0, remember to save as count not countAdj & run 01_feat_normalizeadj


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
for (ds in c("ctrl1","ctrl2","ctrl3","ctrl4","ctrl5","ctrl6","ctrl7","ctrl8","ctrl9","ctrl10","pos1","pos2","pos3","pos4")) {
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
  feat_file_cell_count_dir = paste(feat_dir,"/file-cell-countAdj",sep="")
  
  # cell names
  f@exprs = matrix(rnorm(ncells[1]*length(markers),2,1),nrow=ncells[1])
  a = flowType(Frame=f, PropMarkers=ci, MarkerNames=markers, 
               MaxMarkersPerPop=6, PartitionsPerMarker=2, 
               Thresholds=thress0, 
               Methods='Thresholds', verbose=F, MemLimit=60)
  ftcell = unlist(lapply(a@PhenoCodes, function(x){return( decodePhenotype(x, markers, a@PartitionsPerMarker) )}))

  # flowtype
  ftl = llply(loopInd(1:nsample,no_cores), function(ii) {
    llply(ii, function(i) {
      f@exprs = matrix(rnorm(ncells[i]*length(markers),2,1), nrow=ncells[i])
      thress = thress0
      if (i>(nsample*nctrl) & !grepl("ctrl",ds)) {
        # make base graph for plot
        if (i == nsample*nctrl+1 & grepl("pos",ds)) {
          ft = flowType(Frame=f, PropMarkers=ci, MarkerNames=markers,
                        MaxMarkersPerPop=4, PartitionsPerMarker=2,
                        Thresholds=thress,
                        Methods='Thresholds', verbose=F, MemLimit=60)
          ftcell = unlist(lapply(ft@PhenoCodes, function(x)
            decodePhenotype(x, ft@MarkerNames, rep(2,length(ft@MarkerNames))) ))
          ftv0 = ft@CellFreqs
          ftv0 = round(ftv0/ftv0[1],3)
        }
        # change f values
        if (ds=="pos3") {
          # .125 -> .19 a+b+c+
          
          ap = f@exprs[,1]>thress[[1]]
          bp = f@exprs[,2]>thress[[2]]
          cp = f@exprs[,3]>thress[[3]]
          dp = f@exprs[,4]>thress[[4]]
          # ep = f@exprs[,5]>thress[[5]]
          triple = (ap & bp & cp)
          # quad = ap & bp & cp & dp
          # quint = ap & bp & cp & dp & ep
          
          
          tm = sum(triple)/2
          f@exprs[sample(which(bp & cp & !ap),tm),1] = p75 # bc
          f@exprs[sample(which(ap & cp & !bp),tm),1] = p25 # ac
          f@exprs[sample(which(ap & bp & !cp),tm),1] = p25 # ab
          f@exprs[sample(which(!ap & !bp & !cp),tm),1] = p75 # a
          
          # tn = sum(quint)/2
          # f@exprs[sample(which(!ap & bp & cp & dp & ep),tn),1] = p75
          # f@exprs[sample(which(ap & !bp & cp & dp & ep),tn),1] = p25
          # f@exprs[sample(which(ap & bp & !cp & dp & ep),tn),1] = p25
          # f@exprs[sample(which(ap & bp & cp & !dp & ep),tn),1] = p25
          # f@exprs[sample(which(ap & bp & cp & dp & !ep),tn),1] = p25
          # 
          # f@exprs[sample(which(!ap & bp & cp & !dp & !ep),tn),1] = p75
          # f@exprs[sample(which(!ap & bp & !cp & dp & !ep),tn),1] = p75
          # f@exprs[sample(which(!ap & bp & !cp & !dp & ep),tn),1] = p75
          # f@exprs[sample(which(!ap & !bp & cp & dp & !ep),tn),1] = p75
          # f@exprs[sample(which(!ap & !bp & cp & !dp & ep),tn),1] = p75
          # f@exprs[sample(which(!ap & !bp & !cp & dp & ep),tn),1] = p75
          # 
          # f@exprs[sample(which(ap & bp & !cp & !dp & !ep),tn),1] = p25
          # f@exprs[sample(which(ap & !bp & cp & !dp & !ep),tn),1] = p25
          # f@exprs[sample(which(ap & !bp & !cp & dp & !ep),tn),1] = p25
          # f@exprs[sample(which(ap & !bp & !cp & !dp & ep),tn),1] = p25
          # 
          # f@exprs[sample(which(!ap & !bp & !cp & !dp & !ep),tn),1] = p75
          
        } 
        else if (ds=="pos4") {
          
          ap = f@exprs[,1]>thress[[1]]
          bp = f@exprs[,2]>thress[[2]]
          cp = f@exprs[,3]>thress[[3]]
          dp = f@exprs[,4]>thress[[4]]
          # ep = f@exprs[,5]>thress[[5]]
          # triple = (ap & bp & cp)
          quad = ap & bp & cp & dp
          # quint = ap & bp & cp & dp & ep
          
          tn = sum(quad)/2
          f@exprs[sample(which(!ap & bp & cp & dp),tn),1] = p75 # bcd
          f@exprs[sample(which(ap & !bp & cp & dp),tn),1] = p25 # acd
          f@exprs[sample(which(ap & bp & !cp & dp),tn),1] = p25 # abd
          f@exprs[sample(which(ap & bp & cp & !dp),tn),1] = p25 # abc
          f@exprs[sample(which(!ap & bp & !cp & !dp),tn),1] = p75 # ab
          f@exprs[sample(which(!ap & !bp & cp & !dp),tn),1] = p75 # ac
          f@exprs[sample(which(!ap & !bp & !cp & dp),tn),1] = p75 # ad
          f@exprs[sample(which(ap & !bp & !cp & !dp),tn),1] = p25 # a
        } 
        else if (ds=="pos1") {
          thress = thress1
        } 
        else if (ds=="pos2") {
          thress = thress2
        }
        if (i == nsample*nctrl+1 & grepl("pos",ds)) {
          if (ds%in%c("pos1","pos2")) la=1
          if (ds%in%c("pos3")) la=3
          if (ds%in%c("pos4")) la=4
          
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
          alp = layout_gr(a$grp$e,a$grp$v)
          al = gpdf(al)
          al$v$label = paste0(al$v$name,":",ftv)
          al$v$size=1
          al$v$sizeb=1
          
          vind = abs(ftv-ftv0)/ftv0 >.05
          gp = gggraph(al, v_ind=vind, vb_ind = rep(F,nrow(al$v)),
                        e_ind=al$e[,1]%in%al$v$name[vind] & al$e[,2]%in%al$v$name[vind], 
                        label_ind=str_count(al$v$name,"[-+]")==la & vind)
          ggsave(paste0(meta_dir,"/all_sig.png"), plot=gp, scale = 1, width =11, height =7.5, units = "in", dpi = 300, limitsize = TRUE)
          
          vind = !grepl("-",al$v$name) & abs(ftv-ftv0)/ftv0 >.05
          gp = gggraph(al, v_ind=vind, vb_ind = rep(F,nrow(al$v)),
                        e_ind=al$e[,1]%in%al$v$name[vind] & al$e[,2]%in%al$v$name[vind], 
                        label_ind=str_count(al$v$name,"[-+]")==la & vind)
          ggsave(paste0(meta_dir,"/all_sigpos.png"), plot=gp, scale = 1, width =11, height =7.5, units = "in", dpi = 300, limitsize = TRUE)
        }
      }
      flowType(Frame=f, PropMarkers=ci, MarkerNames=markers, 
               MaxMarkersPerPop=6, PartitionsPerMarker=2, 
               Thresholds=thress, 
               Methods='Thresholds', verbose=F, MemLimit=60)@CellFreqs
    })
  }, .parallel=T)
  ft = Reduce(rbind,unlist(ftl,recursive=F))
  colnames(ft) = ftcell
  rownames(ft) = meta_file$id
  
  save(meta_file, file=paste0(meta_dir,"/file.Rdata"))
  if (writecsv) write.csv(meta_file, file=paste0(meta_dir,"/file.csv"), row.names=F)
  
  # compile ---------------------------
  m0 = as.matrix(ft)
  save(m0, file=paste0(feat_file_cell_count_dir,".Rdata"))
  if (writecsv) write.csv(m0, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)
  
  time_output(start2, ds)
}
time_output(start)

