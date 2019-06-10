## input: gates + fcm files,  meta data paths
## output: feat_file_cell_count, meta_file, meta_cell
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


## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)


## input
input_dir = "/mnt/f/Brinkman group/current/Alice/gating_projects/HDCytoData_Bodenmiller"
data_dir = paste0(input_dir,"/data")
gate_dir = paste0(input_dir,"/gates")
fcs_dir = paste0(input_dir,"/fcs")
gate_dir = paste0(input_dir,"/gates")


## ouput
result_dir = paste0(root, "/result/HDCytoData_Bodenmiller"); dir.create(result_dir, showWarnings=F, recursive=T)
feat_dir = paste0(result_dir,"/feat"); dir.create(feat_dir, showWarnings=F, recursive=T)
feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")
feat_file_cell_prop_dir = paste(feat_dir, "/file-cell-prop", sep="")
meta_dir = paste(result_dir, "/meta", sep=""); dir.create(meta_dir, showWarnings=F, recursive=T)
meta_cell_dir = paste(meta_dir, "/cell", sep="")
meta_file_dir = paste(meta_dir, "/file", sep="")


## libraries
source("source/_func.R")
libr(c("HDCytoData", # requires bioconductor 3.9 which requires R 3.6
       "flowCore", "flowDensity", "diffcyt", "flowType",
       "MASS", "RColorBrewer",
       "foreach", "doMC", "plyr", "stringr"))
# BiocManager::install(version="3.9")

## cores
no_cores = 3 #detectCores()-1
registerDoMC(no_cores)

writecsv = F



load(paste0(input_dir,"/meta_file.Rdata"))
load(paste0(input_dir,"/meta_mark.Rdata"))

load(paste0(data_dir,"/se.Rdata")) #se
# Transform data
se <- transformData(se, cofactor=5)
selist = llply(unique(meta_file$sample_id), function(x) se@assays$data$exprs[meta_file$sample_id==x,])
names(selist) = unique(meta_file$sample_id)





## flowtype ------------------------

start = Sys.time()
markers = c("CD3","CD4","CD20","CD33","CD14","IgM","CD7") # HLA-DR
gthresm = c("cd3.gate","cd4.gate","cd20.gate","cd33.gate","cd14.gate","igm.gate","cd7.gate")
gatesfd = gates = get(load(paste0(gate_dir,"/gates.Rdata"))) #gates
load(paste0(fcs_dir,".Rdata")) #fslist

# # alternative gates
# for (ei in names(selist)) {
#   pop = meta_file$population_id[meta_file$sample_id==ei]
#   exprs = selist[[ei]]
#   
#   gates$cd3.gate[ei] = min(exprs[grepl("T-cells",pop),"CD3"]) #low
#   gates$cd4.gate[ei] = min(exprs[grepl("CD4 T-cells",pop),"CD4"])
#   gates$cd20.gate[ei] = min(exprs[grepl("B-cells",pop),"CD20"]) #low
#   gates$cd33.gate[ei] = max(exprs[grepl("NK cells",pop),"CD33"]) #high
#   # gates$cd14.gate[ei] = max(exprs[grepl("",pop),"CD14"])
#   gates$igm.gate[ei] = min(exprs[grepl("B-cells IgM+",pop),"IgM"]) #low
#   gates$cd7.gate[ei] = min(exprs[grepl("T-cells",pop) | grepl("NK cells",pop), "CD7"]) #low
# }


ftl = llply(1:length(fslist), function(i) {
  flowType(Frame=fslist[[i]], 
           PropMarkers=match(markers, meta_mark$marker_name), 
           MarkerNames=markers, 
           MaxMarkersPerPop=6, PartitionsPerMarker=2, Methods='Thresholds', 
           Thresholds=as.list(gates[i,gthresm]), 
           verbose=F, MemLimit=60)
}, .parallel=F)
ft = ldply(ftl, function(ft) ft@CellFreqs)
ftcell = unlist(lapply(ftl[[1]]@PhenoCodes, function(x){return(decodePhenotype(x, markers, ftl[[1]]@PartitionsPerMarker) )}))
ftp = foreach(xi=1:ncol(ft), .combine='cbind') %dopar% { return(ft[,xi]/ft[,1]) }
colnames(ft) = colnames(ftp) = ftcell
rownames(ft) = rownames(ftp) = names(fslist)

ft = as.matrix(ft)
save(ft, file=paste0(feat_file_cell_count_dir,".Rdata"))
if (writecsv) write.csv(ft, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)

ftp = as.matrix(ftp)
save(ftp, file=paste0(feat_file_cell_prop_dir,".Rdata"))
if (writecsv) write.csv(ftp, file=paste0(feat_file_cell_prop_dir,".csv"), row.names=T)


## meta ------------------------------

meta_file_ = meta_file[!duplicated(meta_file$sample_id),-4]
for (ci in 1:ncol(meta_file_)) 
  meta_file_[,ci] = as.character(meta_file_[,ci])
meta_file_ = as.data.frame(meta_file_)
colnames(meta_file_) = c("class","patient","id")
meta_file_$class = gsub("reference","control",meta_file_$class,ignore.case=T)

save(meta_file_, file=paste0(meta_file_dir,".Rdata"))
if (writecsv) write.csv(meta_file_, file=paste0(meta_file_dir,".csv"), row.names=T)

meta_cell = getPhen(ftcell)
save(meta_cell, file=paste0(meta_cell_dir,".Rdata"))
if (writecsv) write.csv(meta_cell, file=paste0(meta_cell_dir,".csv"), row.names=T)


time_output(start)

