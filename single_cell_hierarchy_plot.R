rstudioapi::openProject("/mnt/f/Brinkman group/current/Alice/flowGraph") # CHANGE THIS
# install.packages("/mnt/f/Brinkman group/current/Alice/flowGraph",repos=NULL, type="source")
library(flowGraph)
library(flowType)

root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric" # CHANGE THIS

## options
nctrl = .5 # % of control samples
markern = maxmarker = 3 # of markers
nsample = 1

normean = 300000 # 301234.7 for pregnancy
normsd = 0 # 64734.05 for pregnancy; if >0, remember to save as count not countAdj & run 01_feat_normalizeadj

lastlsdp = .1 # when generating only last layer cell proportions, use lastlsdp*normean/(2^markern) as sd


# define number of cells in each fcs file
# for pregnancy data set, mean=301234.7, sd=64734.05
# ncells = floor(rnorm(nsample,300000,65000)) # number of cells for each sample
ncells = rnorm(nsample,normean,normsd) # number of cells for each sample

## prepare flowtype files
# load sample fcs file
# f = read.FCS("/mnt/f/Brinkman group/current/Alice/gating_projects/pregnancy/samplefcs.fcs")
f = new("flowFrame")

markers = LETTERS[1:markern] # markers

# marker thresholds
cvd = rnorm(ncells[1], 2, 1)
p50 = quantile(cvd, .5)

thress0 = as.list(rep(p50,10))
names(thress0) = markers
thress_ = purrr::map(seq_len(nsample), function(x) thress0)

# make flow files + flowtype files
ftl = purrr::map2(ncells, thress_, function(ncellsi, thress) {
f@exprs = matrix(rnorm(ncellsi*markern,2,1), nrow=ncellsi)
colnames(f@exprs) = markers

ci = c(1:ncol(f@exprs))
names(ci) = colnames(f@exprs) # marker indices in f@exprs

thress = thress[colnames(f@exprs)]
ft = flowType(Frame=f, PropMarkers=ci, MarkerNames=colnames(f@exprs),
              MaxMarkersPerPop=min(markern,maxmarker), PartitionsPerMarker=2,
              Thresholds=thress,
              Methods='Thresholds', verbose=FALSE, MemLimit=60)#@CellFreqs
})

# flowgraph them
fg=flowGraph(ftl)

# plot
gr=ggdf(fg@graph)
gr$v$colour = -1; gr$v$colour[!grepl("[-]",gr$v$phenotype)] = 1
gr$v$label = gr$v$phenotype
gr$v$label_ind = gr$v$v_ind = T

gp = plot_gr(gr)

ggplot2::ggsave(paste0(root,"/temp.png"), plot=gp,
                scale=1, width=9, height=9, units="in", dpi=600, limitsize=TRUE)

