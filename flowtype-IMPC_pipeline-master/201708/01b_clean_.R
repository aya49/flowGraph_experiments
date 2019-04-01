## Input: Normalized Count and Prop matrices --> Output: Clean out surrogate variable effects
# aya43@sfu.ca 20151228

## root directory
root = "~/projects/IMPC"
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#,"Sanger_MLN","CIPHE","TCP","H")
controlL = c("+_+|+_Y","+_+|+_Y","WildType","WildType","WildType") #control value in target_col column
ci = 1; panel = panelL[ci]; centre = centreL[ci]

result_dir = paste0("result/", panelL, "/", centreL); suppressWarnings(dir.create (result_dir))



## input directories
meta_cell_dir = paste(result_dir, "/meta-cell.Rdata", sep="")
meta_file_dir = paste(result_dir, "/meta-file.Rdata", sep="")
feat_dir = paste(result_dir, "/features", sep=""); dir.create(feat_dir, showWarnings=F)
feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")
feat_file_cell_countAdj_dir = paste(feat_dir, "/file-cell-countAdj", sep="")
feat_file_cell_prop_dir = paste(feat_dir, "/file-cell-prop", sep="")
phenoChild_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChild.Rdata",sep="")
phenoChild_ind_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChild_ind.Rdata",sep="")




interested = c("gene", "gender", "date", "colony", "strain", "fur", "birth_date", "specimen", "sample") #sampleMeta columns to plot
nosample = 40 #number of samples to plot >2
countThres = 500
dpp = c("leafonly","nonleaf")

libr(stringr)
libr(devtools)
libr(Biobase)
libr(preprocessCore)
libr(foreach)
libr(doMC)
source("~/projects/IMPC/code/_funcAlice.R")


#Output
plot_dir = paste(result_dir, "/plots", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dir[i])) }
matrixCountAdj_qn_dir = paste(result_dir, "/matrixCountAdj_qn_all_",dpp,countThres,".Rdata", sep="")
matrixCountAdj_new_dir = paste(result_dir, "/matrixCountAdj_all_",dpp,countThres,".Rdata", sep="")



start = Sys.time()

no_cores = detectCores() - 1
registerDoMC(no_cores)


control = controlL[ci]

sampleMeta = get(load(meta_file_dir))
phenoMeta = get(load(meta_cell_dir))
mm = get(load(feat_file_cell_countAdj_dir))
pc = get(load(phenoChild_dir))
pc_ind = get(load(phenoChild_ind_dir))


interestedCols = which(colnames(sampleMeta)%in%interested)

## Quantile normalization

#features to delete
dsi = which(!sampleMeta$gene%in%control)
lowCounts = which(apply(mm,2,function(x) all(x<=countThres))) #deletes about 20k
leaves = unlist(getleaves(phenoMeta,no_cores))

for (dp in 1:length(dpp)) { #leaves, noleaves
  if (dp==1) dpi = intersect(setdiff(1:nrow(phenoMeta),lowCounts),leaves)
  if (dp==2) dpi = union(lowCounts,setdiff(1:nrow(phenoMeta),pc_ind))
  # m = mm[-dsi,-dpi]
  # sm = sampleMeta[-dsi,]
  m = mm[,-dpi]
  save(m, file=matrixCountAdj_new_dir[dp])
  sm = sampleMeta
  pm = phenoMeta[-dpi,]
  
  
  
  #sample original distributions
  #quantile normalization
  m_qn = t(normalize.quantiles(t(m)))
  dimnames(m_qn) = dimnames(m)
  save(m_qn,file=matrixCountAdj_qn_dir[dp])
  
  pngname = paste0(plot_dir, "/qn_all_", dpp[dp], ".png")
  png(filename=pngname, width=2*700, height=700)
  par(mfrow=c(1,2), mar=rep(2,4))
  
  colramp = colorRampPalette(c(3,"white",2))(nosample)
  s = sample((1:nrow(sm)),nosample)
  plot(density(log(m[s[1],]),na.rm=T),col=colramp[1],lwd=3, main=paste0("Quantile Normalization ln() BEFORE ",dpp[dp],"; countThres ", countThres, "; "))#,ylim=c(0,.30))
  for(i in 2:nosample){lines(density(log(m[s[i],]),na.rm=T),lwd=3,col=colramp[i])}
  
  plot(density(log(m_qn[s[1],]), na.rm=T),col=colramp[1],lwd=3, main=paste0("AFTER : ",dpp[dp],"; countThres ", countThres, "; "))#,ylim=c(0,.20))
  for(i in 2:nosample){lines(density(log(m_qn[s[i],]),na.rm=T),lwd=3,col=colramp[i])}
  
  graphics.off()
  
}

TimeOutput(start)






