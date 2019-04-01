## Input: raw cell counts --> Output: tables and plots
# aya43@sfu.ca 20170419

#Directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)
options(include.rownames=F)

#Input
meta_dir = paste0(result_dir,"/meta")
meta_cell_dir = paste(meta_dir, "/cell", sep="") #metadata for cell populations (phenotype)
meta_file_dir = paste(meta_dir, "/file", sep="") #metadata for FCM files
feat_dir = paste(result_dir, "/feat", sep=""); dir.create(feat_dir, showWarnings=F)
feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")
fcs_dir = paste0("fcs/0006.FCS")

#Output
plot_dir = paste(result_dir, "/plots", sep=""); dir.create(plot_dir, showWarnings=F)
plot_vis_dir = paste(plot_dir, "/vis.png", sep=""); dir.create(plot_vis_dir, showWarnings=F)

#Libraries/Functions
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
libr("stringr")
libr("foreach")
libr("doMC")
libr("flowDensity")
libr("flowCore")

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)






meta_file = get(load(paste0(meta_file_dir,".Rdata")))
meta_cell_dir = get(load(paste0(meta_cell_dir,".Rdata")))
feat_file_cell_count = get(load(paste0(feat_file_cell_count_dir,".Rdata")))








f = read.FCS(fcs_dir)
dups <- duplicated(f@exprs[,2:7])


#Tsne
libr(Rtsne)
libr(ggplot2)
# tsne = Rtsne(as.matrix(f2@exprs[,markercols]), perplexity=1, max_iter=10000)
# plot(tsne$Y,col=rhcl[as.numeric(factor(la))],pch=16,cex=.5)
out_Rtsne <- Rtsne(f@exprs[!dups,1:7], pca = FALSE, verbose = TRUE)

data_plot <- as.data.frame(out_Rtsne$Y)
colnames(data_plot) <- c("tSNE_1", "tSNE_2")

plot(data_plot,col=gg_color_hue(length(unique(la)))[as.numeric(factor(la))], pch=16,cex=1, main="TSNE")


#FlowSOM
libr(FlowSOM)
ff = ReadInput(f2,compensate = F,transform = F, scale = F)
fSOM <- BuildSOM(ff,colsToUse=markercols)
fSOM <- BuildMST(fSOM,tSNE=T)
labels_pre <- fSOM$map$mapping[, 1]
metaClustering <- metaClustering_consensus(fSOM$map$codes,k=8)

res <- data.frame(cluster = metaClustering[labels_pre]) #cluster assignments

PlotPies(fSOM,cellTypes=as.factor(la), backgroundValues = as.factor(metaClustering), colorPalette=grDevices::colorRampPalette(gg_color_hue(length(unique(la)))))








