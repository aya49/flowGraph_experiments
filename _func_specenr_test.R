
## start -------------------------

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
# root = "/home/ayue/projects/flowtype_metric"
setwd(root)

## libraries
source("source/_func_specenr.R")
libr = function(pkgs) {
  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) 
    install.packages(setdiff(pkgs, rownames(installed.packages())), verbose=F)
  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
    if (!requireNamespace("BiocManager", quietly=T)) install.packages("BiocManager")
    BiocManager::install(setdiff(pkgs, rownames(installed.packages())), ask=F)
  }
  sapply(pkgs, require, character.only=T)
}
libr(c("flowCore", "flowType", 
       "doMC", "foreach", 
       "plyr","stringr", "gsubfn"))

## cores
no_cores = detectCores()-1
registerDoMC(no_cores)

## options
options(stringsAsFactors=F)


## bodenmiller ------------------------------------

## input: gates + fcm files,  meta data paths
## output: feat/file-cell-count, meta_file
## HDCytoDate: e.g. Levine_32dim_SE(metadata = FALSE) Levine_32dim_flowSet(metadata = FALSE) 
## remember to transform values! for cytof, usually use asinh with cofactor = 5 (cofactor = 150 for flow cytometry)
## - Clustering:
##   - Levine_32dim: 
##     * cytof "Data-driven phenotypic dissection of AMLreveals progenitor-like cells that correlate with prognosis" 2015; 
##     * human bone marrow cells from 2 healthy subjects H1, H2; 
##     * (265627 (104184 manually gated, 161443 ungated) x cytof 32 surface markers) 
##     * manually gated 14 cell populations
##   - Levine_13dim: 
##     * cytof 
##     * human bone marrow cells from 1 healthy subject 
##     * (167044 (81747 manually gated, 85297 ungated) x 13 surface markers) 
##     * manually gated 24 cell populations
##   - Samusik_01: 
##     * (86864 (53173 manually gated, 33691 ungated) x 39 + ungated surface markers)
##   - Samusik_all:  
##     * cytof "Automated mapping of phenotype space with single-cell data" 2016; 
##     * mouse bone marrow from 10 C57BL/6J mice clones 
##     * (841644 (514386 manually gated, 327258 ungated) x 39 surface markers) 
##     * manually gated 24 + ungated cell populations
##   - Nilsson_rare: 
##     * flow "Frequency determination of rare populations by flow cytometry: A hematopoietic stem cell perspective" 2013; 
##     * human bone marrow cells from 1 healthy subject
##     * (44140 (358 manually gated hematopoietic stem cells) x 13 surface markers)
##   - Mosmann_rare:
##     * flow "SWIFT - Scalable clustering for automated identification of rare cell populations in large, high-dimensional flow cytometry datasets, Part 2: Biological evaluation" 2014
##     * human peripheral blood cells exposed to influenza agents from 1 healthy subject
##     * (296460 (109 manually gated rare activated cytokine producing memory CD4 T cells) x 14 (7 surface + 7 signalling) markers)
## - Differential analysis:
##   - Krieg_Anti_PD_1: strong batch effect, 2 different days ('batch23' and 'batch29')
##     * cytof "High-dimensional single-cell analysispredicts response to anti-PD-1 immunotherapy" 2018
##     * human peripheral blood from 20 melanoma skin cancer patients treated with anti-PD-1 immunotherapy x 2 days pre/post treatment (9/11 non/responders)
##     * CD14+CD16-HLA-DRhi monocytes (a small subpopulation of CD14+CD33+HLA-DRhiICAM-1+CD64+CD141+CD86+CD11c+CD38+PD-L1+CD11b+ monocytes) prior to treatment is strong predictor of survival status following immunotherapy treatment.
##     * (85715 cells x 24 cell type markers (exclude CD45 b/c all cells show high so CD45=none))
##   - Bodenmiller_BCR_XL 
##     * cytof "Multiplexed masscytometry profiling of cellular states perturbed by small-molecule regulators" 2012
##     * human peripheral blood from 8 x 2 BCR-XL (b cell receptor / Fc receptor cross linker) un/stimulated healthy subjects; 
##     * differentially expressed signalling markers in many cell populations e.g. phosphorylated S6 (pS6) in B cells
##     * (172791 x 10 surface markers (cell type) + 14 intracellular signalling functional markers (cell state))

## input
input_dir = "/mnt/f/Brinkman group/current/Alice/gating_projects/HDCytoData_Bodenmiller"
data_dir = paste0(input_dir,"/data")
gate_dir = paste0(input_dir,"/gates")
fcs_dir = paste0(input_dir,"/fcs")
gate_dir = paste0(input_dir,"/gates")


## ouput
result_dir = paste0(root, "/result/bodenmiller"); dir.create(result_dir, showWarnings=F, recursive=T)


start = Sys.time()

## load data
load(paste0(fcs_dir,".Rdata")) # fslist; fcs files
meta_file0 = get(load(paste0(input_dir,"/meta_file.Rdata")))
load(paste0(input_dir,"/meta_mark.Rdata")) # meta_mark; markers in fcs
# load(paste0(data_dir,"/se.Rdata")) # se; data

## feat/file-count-count: flowtype
markers = c("CD3","CD4","CD20","CD33","CD14","IgM","CD7") # HLA-DR
gthresm = c("cd3.gate","cd4.gate","cd20.gate","cd33.gate","cd14.gate","igm.gate","cd7.gate")
gatesfd = gates = get(load(paste0(gate_dir,"/gates.Rdata"))) #gates

for (i in 1:length(fslist)) {
  f = fslist[[i]]
  print(f@exprs[1:10,1])
}

ftl = llply(1:length(fslist), function(i) {
  flowType(Frame=fslist[[i]], 
           PropMarkers=match(markers, meta_mark$marker_name), 
           MarkerNames=markers, 
           MaxMarkersPerPop=6, PartitionsPerMarker=2, Methods='Thresholds', 
           Thresholds=as.list(gates[i,gthresm]), 
           verbose=F, MemLimit=60)
}, .parallel=T)

## prepare meta
meta_file = meta_file0[!duplicated(meta_file0$sample_id),-4]
for (ci in 1:ncol(meta_file)) 
  meta_file[,ci] = as.character(meta_file[,ci])
meta_file = as.data.frame(meta_file)
colnames(meta_file) = c("class","subject","id")
meta_file$class = gsub("reference","control",meta_file$class,ignore.case=T)
meta_file$class[meta_file$class!="control"] = "exp"
meta_file = meta_file[match(names(fslist),meta_file$id),]
for (uc in unique(meta_file$class)) {
  uci = meta_file$class==uc
  meta_file$train[uci] = ifelse(which(uci)%in%sample(which(uci),sum(uci)/2),T,F)
}

fg = flowgraph(ftl, no_cores=no_cores, meta=meta_file, norm_path=paste0(result_dir,"/count_norm"))
fg = flowgraph_mean_class(fg, class="subject", no_cores=no_cores)
save(fg, file=paste0(result_dir,"/fg.Rdata"))

time_output(start, "data_bodenmiller")




## flowcap ---------------------------------------------------------

## input directories
fc_dir = "/mnt/f/Brinkman group/current/Alice/flowCAP-II"
ft_dir = paste0(fc_dir,"/data/FT") # flowtype file directory
csv_dir = paste0(fc_dir,"/meta_file.csv") # meta file directory
markers_dir = paste0(fc_dir,"/markers.Rdata") # fcs file directory

## output directories (NOT last section, output split by tube/panel)
result_dir0 = paste0(root, "/result/flowcap")

## options
matchsamples = 5 #pos: number of samples from each normal and aml to mix
amlprop = 1/2 #ctrl: proportion of control sample to make as "aml"


start = Sys.time()

## prepare flowtype directories
ft_dirs = sort( dir(ft_dir, pattern=".Rda", all.files=T, full.names=T, recursive=T) )
ft_names = sapply(strsplit(ft_dirs,"/"), function(a) gsub(".Rda","",a[length(a)]) )

## prepare meta
meta_file0 = read.csv(csv_dir)[,-1]
for (i in 1:ncol(meta_file0)) 
  if (is.factor(meta_file0[,i])) 
    meta_file0[,i] = as.character(meta_file0[,i])

markers = get(load(markers_dir))

# ## this or
# m00 = llply(loop_ind_f(1:length(ft_dirs),no_cores),function(ii)
#   llply(ii, function(i) get(load(ft_dirs[i]))@CellFreqs), .parallel=T)
# m00 = as.matrix(Reduce('rbind',llply(m00,function(x)Reduce(rbind,x))))
# rownames(m00) = ft_names
# colnames(m00) = rownames(get(load(ft_dirs[1]))@MFIs)
# fg0 = flowgraph(m00, no_cores=no_cores, meta=meta_file0, normalize=F, specenr=F)

## that
fg0 = flowgraph(ft_dirs, no_cores=no_cores, meta=meta_file0, normalize=F, specenr=T)

randomindp = randomindc = NULL
for (tube in unique(meta_file0$tube)) {
  if (!tube==6) next
  
  ## split data by tube/panel
  fg = extract_samples(fg0, fg0@meta$id[fg0@meta$tube==tube])
  fg = gsub_ids(fg, as.numeric(gsub("T[0-9]S|FT","",fg@meta$id)) )
  fg = flowgraph_normalize(fg, no_cores=no_cores, norm_path=paste0(result_dir0, "_", tube, "/count_norm")) # normalize count
  marker = markers[[tube]]
  fg = gsub_markers(fg, marker)
  
  mc = fg@feat$node$count_norm
  
  ## extract controls
  fg_c = extract_samples(fg, fg@meta$id[fg@meta$class=="control"])
  fg_c@meta$class[1:(nrow(fg_c@meta)/2)] = "experiment"
  dir.create(paste0(result_dir0, "_", tube, "_ctrl"), showWarnings=F)
  save(fg_c, file=paste0(result_dir0, "_", tube, "_ctrl/fg.Rdata"))
  
  
  ## make a new class: mixed
  normali = which(fg@meta$class=="control")
  amli = which(fg@meta$class=="aml")
  if (is.null(randomindp))
    for (i in 1:min(length(normali),length(amli))) {
      randomindp[[i]] = list()
      randomindp[[i]]$normal = sample(normali, matchsamples)
      randomindp[[i]]$aml = sample(amli, matchsamples)
    }
  weight = 1/(2*matchsamples)
  
  fg2 = fg
  for (nfeat in names(fg2@feat)) {
    for (nne in names(fg2@feat[[nfeat]])) {
      m = as.matrix(fg2@feat[[nfeat]][[nne]])
      mc2 = as.matrix(do.call(rbind,llply(randomindp,function(ri)
        weight*colSums(m[append(ri$normal,ri$aml),,drop=F]) )))
      rownames(mc2) = c((1+nrow(m)):(nrow(mc2)+nrow(m)))
      fg2@feat[[nfeat]][[nne]] = mc2
    }
  }
  
  # meta
  fg2@meta = meta2 = data.frame(
    class=rep("mix", length(randomindp)), 
    id=rownames(mc2), 
    train=append(rep(F,floor(length(randomindp)/2)), 
                 rep(T,ceiling(length(randomindp)/2))),
    subject=0,
    tube=tube
  )
  fg1 = fg
  fg = merge_samples(fg1, fg2)
  dir.create(paste0(result_dir0, "_", tube), showWarnings=F)
  save(fg, file=paste0(result_dir0, "_", tube, "/fg.Rdata"))
}

time_output(start, "data_flowcap")


## genentech -------------------------------------

## input: genetech by daniel yokosawa on the genetech data (bone marrow + blood for 3 patients x 5 samples; computationally mixed)
## output: feat_file_cell_count, meta_file

## input directories
input_dir0 = "/mnt/f/Brinkman group/current/Alice/gating_projects/genentech"

for (tube in 2:4) {
  input_dir = paste0(input_dir0,"/Tube_", str_pad(tube,3,"left",0),"/comp_mix")


## ouput directories
result_dir = paste0(root, "/result/genentech_",tube); dir.create(result_dir, showWarnings=F, recursive=T)


start = Sys.time()

## prepare flowtype directories
# ftl = get(load(paste0(input_dir,"/flowType_myeloid.Rdata")))
ftl = get(load(paste0(input_dir,"/flowType_sing.Rdata")))
names(ftl) = gsub("[%]","",names(ftl))
markers = get(load(paste0(input_dir,"/MarkerNames_sing.Rdata")))

## prepare meta
temp_ = Reduce("rbind", str_split(gsub(".fcs","",names(ftl)),"_"))
meta_file = data.frame(id=names(ftl), class=temp_[,3], patient=as.numeric(gsub(".fcs","",temp_[,4])))
meta_file$class[meta_file$class=="100WB"] = "control"
for (uc in unique(meta_file$class)) {
  uci = meta_file$class==uc
  meta_file$train[uci] = ifelse(which(uci)%in%sample(which(uci),sum(uci)/2),T,F)
}

fg = flowgraph(ftl, markers=markers, no_cores=no_cores, meta=meta_file, norm_path=paste0(result_dir,"/count_norm"))
save(fg, file=paste0(result_dir,"/fg.Rdata"))
}
time_output(start, "data_genetech")


## pregnancy ----------------------------------

## input: gates + fcm files,  meta data paths
## output: feat_file_cell_count, meta_file
## process: 
## - takes gates + fcm files, outputs flowtype vectors 
## - compiles flowtype vectors together to create cell count matrix
## - reformats meta file (meta info for fcm files)
## - creates cell meta file (meta info for fcm cell populations)

## input directories
input_dir = "/mnt/f/Brinkman group/current/Alice/gating_projects/pregnancy"
meta_file_dir = paste0(input_dir, "/meta_file.Rdata")
flowtype_dir = paste0(input_dir, "/flowtype")

## ouput directories
result_dir0 = paste0(root, "/result/pregnancy"); dir.create(result_dir0, showWarnings=F, recursive=T)


start = Sys.time()

## prepare flowtype directories
ft_dirs = list.files(flowtype_dir, full.names=T)
ft_names = sapply(strsplit(ft_dirs,"/"), function(a) gsub(".Rdata","",a[length(a)]) )
ft_names = gsub(".Rdata|.fcs|Gates_|_Unstim|_Repeat","",ft_names)
ft_names = gsub("BL","4",ft_names)

## prepare meta
meta_file0 = get(load(meta_file_dir))
colnames(meta_file0)[1] = "subject"
meta_file0 = meta_file0[match(ft_names, meta_file0$id),,drop=F]
meta_file0$class[meta_file0$class==4] = "control"
meta_file0$train = ifelse(meta_file0$type=="train",T,F)
meta_file0 = meta_file0[,-4]

## feat/file-cell-count: load and compile flowtype count files
fg = flowgraph(ft_dirs, meta=meta_file0, no_cores=no_cores, norm_path=paste0(result_dir0,"/count_norm"))
fg = flowgraph_mean_class(fg, class="subject", no_cores=no_cores)

save(fg, file=paste0(result_dir0,"/fg.Rdata"))

time_output(start, "data_pregnancy")



## positive/negative control -----------------------

## input: gates + fcm files,  meta data paths
## output: feat_file_cell_count, meta_file
## process: 
## - takes gates + fcm files (values randomly generated via normal distribution; see different pos_ for how positive controls are changed
## - compiles flowtype vectors together to create cell count matrix
## - reformats meta file (meta info for fcm files)
## - creates cell meta file (meta info for fcm cell populations)

## options
nsample = 1000 # # of samples to create per data set; don't change this, i hard coded in cytodx and save fcm below on which files to use
nctrl = .5 # % of control samples
markern = 4 # # of markers
maxmarker = 6

normean = 300000 # 301234.7 for pregnancy
normsd = 0 # 64734.05 for pregnancy; if >0, remember to save as count not countAdj & run 01_feat_normalizeadj

lastlsdp = .1 # when generating only last layer cell proportions, use lastlsdp*normean/(2^markern) as sd

## cores
no_cores = 6
registerDoMC(no_cores)


start = Sys.time()

# define number of cells in each fcs file
# for pregnancy data set, mean=301234.7, sd=64734.05
# ncells = floor(rnorm(nsample,300000,65000)) # number of cells for each sample
ncells = rnorm(nsample,normean,normsd) # number of cells for each sample

## meta/file
meta_file = data.frame(id=paste0("a",1:nsample), class="exp")
meta_file$class[1:(nctrl*nsample)] = "control"
meta_file$train = rep(F,nrow(meta_file))
meta_file$train[(nctrl*nsample+1):(nctrl*nsample+(1-nctrl)*nsample/2)] = meta_file$train[1:(nctrl*nsample/2)] = T

## prepare flowtype files
# load sample fcs file
# f = read.FCS("/mnt/f/Brinkman group/current/Alice/gating_projects/pregnancy/samplefcs.fcs")
f = new("flowFrame")

markers = LETTERS[1:markern] # markers

# marker thresholds
cvd = rnorm(ncells[1],2,1)
p50 = quantile(cvd, .5)
p60 = quantile(cvd, .6)
p75 = quantile(cvd, .75)
p25 = quantile(cvd, .25)
thress0 = llply(markers, function(x) p50); names(thress0) = markers
# thress1 = thress2 = thress4 = 
  thress5 = thress0
# thress1[[markers[1]]] = thress2[markers[1:2]] = p25
# thress4[markers[1:4]] = quantile(cvd, .501)
thress5[[markers[1]]] = c(p25,p50)
thress5[[markers[2]]] = c(p25,p50,p60)

#paste0("ctrl",c(0:9)), 
#paste0("pos",c(1:30))
for (ds in c(paste0("pos",c(1:30)),paste0("ctrl",c(0:9)))) {
  start2 = Sys.time()
  
  # ouput directories
  result_dir = paste0(root, "/result/",ds)
  fcs_dir = paste0(result_dir,"/fcs"); dir.create(fcs_dir, showWarnings=F, recursive=T)

  # make cell names
  f@exprs = matrix(rnorm(ncells[1]*length(markers),2,1),nrow=ncells[1])
  # a = flowType(Frame=f, PropMarkers=ci, MarkerNames=markers, 
  #              MaxMarkersPerPop=min(markern,maxmarker), PartitionsPerMarker=2, 
  #              Thresholds=thress0, 
  #              Methods='Thresholds', verbose=F, MemLimit=60)
  # ftcell = unlist(lapply(a@PhenoCodes, function(x){return( decodePhenotype(x, markers, a@PartitionsPerMarker) )}))
  # ftcell_ = str_count(ftcell,"[+|-]")
  # lastlcp = ftcell[ftcell_==length(markers)]
  # lastlcpm = llply(str_extract_all(lastlcp,"[A-Z][+|-]"), function(x) 
  #   grepl("[+]",x) )
  # lastlallposi = which(sapply(lastlcpm, function(x) all(x)))
  # lastlallnegi = which(sapply(lastlcpm, function(x) all(!x)))
  
  # flowtype
  loop_ind = loop_ind_f(1:nsample,no_cores)
  ftl = llply(loop_ind, function(ii) {
    llply(ii, function(i) {
      # v1 randomized matrix
      f@exprs = matrix(rnorm(ncells[i]*length(markers),2,1), nrow=ncells[i])
      colnames(f@exprs) = markers
      ci = c(1:ncol(f@exprs)); names(ci) = colnames(f@exprs) # marker indices in f@exprs
      
      thress = thress0
      if (i>(nsample*nctrl) & grepl("pos",ds)) {
        # make base graph for plot
        # if (i == nsample*nctrl+1 & grepl("pos",ds)) {
        #   ft = flowType(Frame=f, PropMarkers=ci, MarkerNames=markers,
        #                 MaxMarkersPerPop=min(markern,maxmarker), PartitionsPerMarker=2,
        #                 Thresholds=thress,
        #                 Methods='Thresholds', verbose=F, MemLimit=60)
        #   ftcell = unlist(lapply(ft@PhenoCodes, function(x)
        #     decodePhenotype(x, ft@MarkerNames, rep(2,length(ft@MarkerNames))) ))
        #   ftv0 = ftv0_ = ft@CellFreqs
        #   ftv0 = round(ftv0/ftv0[1],3)
        # }
          dp = f@exprs[,4]>thress[[4]]
        ap = f@exprs[,1]>thress[[1]]
        bp = f@exprs[,2]>thress[[2]]
        cp = f@exprs[,3]>thress[[3]]
        # ep = f@exprs[,5]>thress[[5]]
        double = ap & bp
        triple = ap & bp & cp
          quad   = ap & bp & cp & dp
        # quint = ap & bp & cp & dp & ep
        
        # change f values
        if (ds=="pos1") { # A+ > .75; A- > .25
          tm = sum(ap)/2
          f@exprs[sample(which(!ap),tm),1] = p75 # 
        } 
        else if (ds=="pos2") { # A+, B+ > .75; A-, B- > .25
          tm = sum(ap)/2
          f@exprs[sample(which(!ap),tm),1] = p75 # 
          f@exprs[sample(which(!bp),tm),2] = p75 # 
        } 
        else if (ds=="pos3") { # A-B+ > A+B+ x1.5
          tm = sum(double)/2
          f@exprs[sample(which(bp & !ap),tm),1] = p75 # 
          f@exprs[sample(which(ap & !bp),tm),1] = p25 # 
          f@exprs[sample(which(!ap & !bp),tm),1] = p25 # 
        }
        else if (ds=="pos4") { # A-B+ > A+B+ x1.5; D+c- > D+c+ x1.5
          tm = sum(double)/2
          f@exprs[sample(which(bp & !ap),tm),1] = p75 # 
          f@exprs[sample(which(ap & !bp),tm),1] = p25 # 
          f@exprs[sample(which(!ap & !bp),tm),1] = p25 # 
          
          f@exprs[sample(which(dp & !cp),tm),3] = p75 # 
          f@exprs[sample(which(cp & !dp),tm),3] = p25 # 
          f@exprs[sample(which(!dp & !cp),tm),3] = p25 # 
        }
        else if (ds=="pos5") { # A-B+C+ > A+B+C+ x1.5
          tm = sum(triple)/2
          f@exprs[sample(which(bp & cp & !ap),tm),1] = p75 # bc
          f@exprs[sample(which(ap & cp & !bp),tm),1] = p25 # ac
          f@exprs[sample(which(ap & bp & !cp),tm),1] = p25 # ab
          f@exprs[sample(which(!ap & !bp & !cp),tm),1] = p75 # a
        }
        else if (ds=="pos6") { # A-B+C+ > A+B+C+ x1.5; b+D+c- > b+D+c+ x1.5
          tm = sum(triple)/2
          f@exprs[sample(which(bp & cp & !ap),tm),1] = p75 # bc
          f@exprs[sample(which(ap & cp & !bp),tm),1] = p25 # ac
          f@exprs[sample(which(ap & bp & !cp),tm),1] = p25 # ab
          f@exprs[sample(which(!ap & !bp & !cp),tm),1] = p75 # a
          
          f@exprs[sample(which(bp & cp & !dp),tm),3] = p75 # bc
          f@exprs[sample(which(bp & dp & !cp),tm),3] = p25 # ac
          f@exprs[sample(which(dp & cp & !bp),tm),3] = p25 # ab
          f@exprs[sample(which(!bp & !cp & !dp),tm),3] = p75 # a
        }
        else if (ds=="pos7") { # A+B+C+D+ > x2
          f@exprs = rbind(f@exprs, f@exprs[ap & bp & cp & dp,])
        } 
        else if (ds=="pos8") { # A-B-C-D-E- > x2
          f@exprs = rbind(f@exprs, f@exprs[!ap & !bp & !cp & !dp,])
        } 
        else if (ds=="pos9") { # same as above but both
          f@exprs = rbind(f@exprs, f@exprs[ap & bp & cp & dp,])
          f@exprs = rbind(f@exprs, f@exprs[!ap & !bp & !cp & !dp,])
        } 
        else if (ds=="pos10") { #1.5x A+, B+C+D+
          tm = sum(ap)/2
          tripleind = which(ap)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
          tm = sum(triple)/2
          tripleind = which(cp & dp & bp)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
        }
        else if (ds=="pos11") { #1.5x A+
          tm = sum(double)
          f@exprs = rbind(f@exprs,f@exprs[sample(which(ap),tm),])
        } 
        else if (ds=="pos12") { #1.5x A+, B+
          tm = sum(double)
          f@exprs = rbind(f@exprs,f@exprs[sample(which(ap),tm),])
          f@exprs = rbind(f@exprs,f@exprs[sample(which(bp),tm),])
        } 
        else if (ds=="pos13") { #1.5x A+B+
          tm = sum(double)/2
          f@exprs = rbind(f@exprs,f@exprs[sample(which(ap & bp),tm),])
        } 
        else if (ds=="pos14") { #1.5x A+B+, C+D+
          tm = sum(double)/2
          f@exprs = rbind(f@exprs,f@exprs[sample(which(ap & bp),tm),])
          f@exprs = rbind(f@exprs,f@exprs[sample(which(dp & cp),tm),])
        } 
        else if (ds=="pos15") {
          tm = sum(triple)/2
          f@exprs = rbind(f@exprs,f@exprs[sample(which(ap & bp & cp),tm),])
        }
        else if (ds=="pos16") { #1.5x A+B+C+, B+C+D+
          tm = sum(triple)/2
          f@exprs = rbind(f@exprs,f@exprs[sample(which(ap & bp & cp),tm),])
          f@exprs = rbind(f@exprs,f@exprs[sample(which(cp & dp & bp),tm),])
        } 
        else if (ds=="pos17") { #1.5x A+B+; decrease other ones accordingly
          tm = .33*nrow(f@exprs) - sum(ap & bp) #sum(double)/3
          f@exprs[c(sample(which(!ap & !bp),tm/3) , 
                    sample(which(ap & !bp),tm/3) , 
                    sample(which(!ap & bp),tm/3)),c(1,2)] = p75
        } 
        else if (ds=="pos18") { #1.5x A+B+, C+D+, like above
          tm = sum(double)/2
          f@exprs[c(sample(which(!ap & !bp),tm/3) , 
                    sample(which(ap & !bp),tm/3) , 
                    sample(which(!ap & bp),tm/3)),c(1,2)] = p75
        }
        else if (ds=="pos19") { #1.5x A+, A+B+C+
          tm = sum(triple)/2
          f@exprs = rbind(f@exprs,f@exprs[sample(which(ap & bp & cp),tm),])
          tm = sum(ap)*1.5-sum(ap)
          f@exprs = rbind(f@exprs,f@exprs[sample(
            which(f@exprs[,1]>thress[[1]] & f@exprs[,2]>thress[[2]]), tm),])
        }
        else if (ds=="pos20") { #1.5x A+, A+B+C+
          tm = sum(ap)/2
          f@exprs = rbind(f@exprs,f@exprs[sample(which(ap),tm),])
          tm = sum(triple)/2
          f@exprs = rbind(f@exprs,f@exprs[sample(which(cp & ap & bp),tm),])
        }
        else if (ds=="pos21") { # 1.5x A+
          tm = sum(double)
          f@exprs = f@exprs[-sample(which(ap),tm),]
        } 
        else if (ds=="pos22") { # 1.5x A+, B+
          tm = sum(double)
          f@exprs = f@exprs[-unique(c(sample(which(bp),tm), sample(which(ap),tm))),]
        } 
        else if (ds=="pos23") { #1.5x A+B+
          tm = sum(double)/2
          f@exprs = f@exprs[-sample(which(ap & bp),tm),]
        } 
        else if (ds=="pos24") { #1.5x A+B+, C+D+
          tm = sum(double)/2
          tind = sample(which(dp & cp),tm)
          tind2 = sample(which(ap & bp),tm)
          f@exprs = f@exprs[-c(unique(c(tind,tind2))),]
        } 
        else if (ds=="pos25") {
          tm = sum(triple)/2
          f@exprs = f@exprs[-sample(which(ap & bp & cp),tm),]
        }
        else if (ds=="pos26") { #1.5x A+B+C+, B+C+D+
          tm = sum(triple)/2
          f@exprs = f@exprs[unique(c(which(ap & bp & cp),sample(which(cp & dp & bp),tm))),]
        } 
        else if (ds=="pos27") { #1.5x A+, A+B+C+
          tm = sum(triple)/2
          tripleind = which(ap & bp & cp)
          f@exprs = rbind(f@exprs,f@exprs[-unique(c(sample(which(ap & bp & cp),sum(triple)/2),sample(which(ap),sum(ap)/2))),])
        }
        else if (ds=="pos28") { #1.5x A+, A+B+C+
          tm = sum(triple)/2
          tripleind = which(cp & ap & bp)
          tm1 = sum(ap)/2
          tripleind1 = which(ap)
          tripleind1 = tripleind1[!tripleind1%in%tripleind]
          f@exprs = f@exprs[-c(sample(tripleind1,tm1),sample(tripleind,tm)),]
        }
        
        else if (ds=="pos29") { #1.5x A+, B+C+D+
          tm = sum(ap)/2
          tripleind = which(ap)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
          tm = sum(triple)/2
          tripleind = which(cp & dp & bp)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
        }
        else if (ds=="pos30") { #same as 25; A+B+C+ + 50%
          tm = sum(triple)/2
          tripleind = which(ap & bp & cp)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
        }

        # if (i == nsample*nctrl+1 & grepl("pos",ds)) {
        # 
        #   ft = flowType(Frame=f, PropMarkers=ci, MarkerNames=markers,
        #                 MaxMarkersPerPop=min(markern,maxmarker), PartitionsPerMarker=2,
        #                 Thresholds=thress,
        #                 Methods='Thresholds', verbose=F, MemLimit=60)
        #   ftcell = unlist(lapply(ft@PhenoCodes, function(x)
        #     decodePhenotype(x, ft@MarkerNames, rep(2,length(ft@MarkerNames))) ))
        #   ftv = ftv_ = ft@CellFreqs
        #   ftv = round(ftv/ftv[1],3)
        #   names(ftv) = ftcell
        #   a = getPhenCP(cp=ftcell,no_cores=no_cores)
        # }
      }
      fe = f@exprs
      colnames(fe) = LETTERS[1:ncol(fe)]
      if (i%in%c(1:5,251:255,501:505,751:755) & grepl("pos",ds)) save(fe,file=paste0(fcs_dir,"/a",i,".Rdata"))
      if (ds%in%c("pos30")) {
        thress = thress5
        ppm = rep(2,markern)
        ppm[1] = 3
        ppm[2] = 4
        ft = flowType(Frame=f, PropMarkers=ci, MarkerNames=markers, 
                      MaxMarkersPerPop=min(markern,maxmarker), PartitionsPerMarker=ppm, 
                      Thresholds=thress, 
                      Methods='Thresholds', verbose=F, MemLimit=60)
        ftcell = unlist(lapply(ft@PhenoCodes, function(x)
          decodePhenotype(x, ft@MarkerNames, ppm) ))
        # ft = ft@CellFreqs
        names(ft@CellFreqs) = ftcell
      } else {
        thress = thress[colnames(f@exprs)]
        ft = flowType(Frame=f, PropMarkers=ci, MarkerNames=colnames(f@exprs), 
                      MaxMarkersPerPop=min(markern,maxmarker), PartitionsPerMarker=2, 
                      Thresholds=thress, 
                      Methods='Thresholds', verbose=F, MemLimit=60)#@CellFreqs
      }
      return(ft)
    })
  }, .parallel=T)
  ftl = unlist(ftl,recursive=F)

  if (ds=="pos30") {
    fg0 = flowgraph(ftl, no_cores=no_cores, meta=meta_file, prop=F, specenr=F, normalize=F)
    fg1 = flowgraph_cumsum(fg0, no_cores=no_cores)
    
    fg0 = flowgraph_prop(fg0)
    fg0 = flowgraph_prop_edge(fg0, no_cores=no_cores)
    fg0 = flowgraph_normalize(fg0, norm_path=paste0(result_dir,"/count_norm"), no_cores=no_cores)
    fg0 = flowgraph_specenr(fg0, no_cores=no_cores)
    
    fg1 = flowgraph_prop(fg1)
    fg1 = flowgraph_prop_edge(fg1, no_cores=no_cores)
    fg1 = flowgraph_normalize(fg1, norm_path=paste0(result_dir,"_cumsum/count_norm"), no_cores=no_cores)
    fg1 = flowgraph_specenr(fg1, no_cores=no_cores)
    
    save(fg0, file=paste0(result_dir,"/fg.Rdata"))
    save(fg1, file=paste0(result_dir,"_cumsum/fg.Rdata"))
    
  } else {
    fg = flowgraph(ftl, no_cores=no_cores, meta=meta_file, norm_path=paste0(result_dir,"/count_norm"))
    save(fg, file=paste0(result_dir,"/fg.Rdata"))
  }
  
  time_output(start2, ds)
  rm(list=c("ftl","fg")); gc()
}
time_output(start)


## p value ------------------------

## cores
no_cores = detectCores()-1
registerDoMC(no_cores)

result_dirs = list.dirs(paste0(root,"/result"), recursive=F)
for (result_dir in result_dirs) {
  # if (!grepl("genentech",result_dir)) next
  try ({
  start1 = Sys.time()
  cat(result_dir)
  fg = get(load(paste0(result_dir,"/fg.Rdata")))
  fg = flowgraph_clear_p(fg)
  fg = flowgraph_p(
    fg, no_cores=no_cores, class="class", control="control",
    overwrite=F, 
    test_name="t_BY",
    diminish=F,
    p_thres=.05, p_rate=2 # only used if diminish=T
    )
  fg = flowgraph_p(
    fg, no_cores=no_cores, class="class", control="control",
    overwrite=F, 
    test_name="t_BY_diminish",
    diminish=T,
    p_thres=.05, p_rate=2 # only used if diminish=T
    )
  
  save(fg, file=paste0(result_dir,"/fg.Rdata"))
  time_output(start1)
  
  })
  
  
  
  # ## test for pos15 data set
  # for (cpop in fg@graph$v$phenotype[-1]) {
  #   a1 = fg@feat$node$expect_prop[1:500,cpop]
  #   a2 = fg@feat$node$prop[1:500,cpop]
  #   a = a1/a2
  #   b1 = fg@feat$node$expect_prop[501:1000,cpop]
  #   b2 = fg@feat$node$prop[501:1000,cpop]
  #   b = b1/b2
  #   
  #   dir.create(paste0(result_dir,"/plots/ratios"), recursive=T, showWarnings=F)
  #   try ({
  #   gp = ggplot() + ggtitle(paste0("population ",cpop,"\nred=control, blue=experiment","\nexpected/actual proportion", "\n t test p value: ", t.test(a,b)$p.value,"\nmeans: ", round(mean(a1),3),"/",round(mean(a2),3),", ",round(mean(b1),3),"/",round(mean(b2),3))) + 
  #     geom_density(aes(x=x), colour="red",data=data.frame(x=a)) + 
  #     geom_density(aes(x=x), colour="blue",data=data.frame(x=b))
  #   
  #   ggsave(paste0(result_dir,"/plots/ratios/",cpop,".png"), plot=gp, scale=1, width=5, height=5, units="in", dpi=500, limitsize=T)
  #   
  #   })
  #   
  #   try({
  #     
  #   gp = ggplot() + ggtitle(paste0("population ",cpop,"\nred=control, blue=experiment","\nexpected/actual proportion", "\n t test p value: ", t.test(log(a),log(b))$p.value,"\nmeans: ", round(mean(log(a1)),3),"/",round(mean(log(a2)),3),", ",round(mean(log(b1)),3),"/",round(mean(log(b2)),3))) + 
  #     geom_density(aes(x=x), colour="red",data=data.frame(x=log(a))) + 
  #     geom_density(aes(x=x), colour="blue",data=data.frame(x=log(b)))
  #   
  #   ggsave(paste0(result_dir,"/plots/ratios/",cpop,"_log.png"), plot=gp, scale=1, width=5, height=5, units="in", dpi=500, limitsize=T)
  #   })
  #   
  # }
  
}


## plots --------------------------

result_dirs = list.dirs(paste0(root,"/result"), recursive=F)
for (result_dir in result_dirs) {
  if (!grepl("genentech",result_dir)) next
  start1 = Sys.time()
  cat(result_dir)
  
  fg = get(load(paste0(result_dir,"/fg.Rdata")))
  for (cls in unique(fg@meta$class)) {
    try ({
      if (cls=="control") next
      flowgraph_summary_plot(
        fg, sumplot=F,
        idname1="control", idname2=cls, class="class", 
        nodeft="specenr", show_label=rep(T,nrow(fg@graph$v)),
        label1="expect_prop", label2="prop", # node only
        show_bgedges=T, width=20, height=9,
        path=paste0(result_dir,"/plots"))
    })
  }
  time_output(start1)
  
}

for (result_dir in result_dirs) {
  if (!grepl("genentech",result_dir)) next
  start1 = Sys.time()
  cat(result_dir)
  
  fg = get(load(paste0(result_dir,"/fg.Rdata")))
  
  for (cls in unique(fg@meta$class)) {
    try ({
      if (cls=="control") next
      flowgraph_summary_plot(
        fg,
        method="t_BY", class="class", idname1="control", idname2=cls,
        nodeft="specenr", edgeft="prop", 
        label1="expect_prop", label2="prop", # node only
        p_thres=.05, show_bgedges=T,
        path=paste0(result_dir,"/plots"))
    })
  }
  time_output(start1)
  
}


