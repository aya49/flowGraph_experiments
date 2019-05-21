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
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metrics"
setwd(root)


## ouput
result_dir = paste0(root, "/result/diffcytof"); dir.create(result_dir, showWarnings=F, recursive=T)


## libraries
source("source/_funcAlice.R")
libr(c("HDCytoData",
       "foreach", "doMC", "plyr", "stringr"))

## cores
no_cores = 8#detectCores()-1
registerDoMC(no_cores)


writecsv = F




start = Sys.time()


time_output(start)

