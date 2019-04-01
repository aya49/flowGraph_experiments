# aya43@sfu.ca 20151228
# Normalizes cell count matrix

#root = "~/projects/IMPC"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN") #,"Sanger_MLN","CIPHE","TCP","H")

options(stringsAsFactors=FALSE)
#options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
matrixCountAdj_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixCountAdj.Rdata", sep="")
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
plot_dir = paste(result_dir, "/", panelL, "/", centreL, "/plots", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dir[i])) }
matrixCountAdj_qn_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixCountAdj_qn_all_",dpp,countThres,".Rdata", sep="")
matrixCountAdj_new_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixCountAdj_all_",dpp,countThres,".Rdata", sep="")



start = Sys.time()

no_cores = detectCores() - 1
registerDoMC(no_cores)

for (ci in 1:length(paste0(panelL,centreL))) {
  start1 = Sys.time()
  
  sampleMeta = get(load(sampleMeta_dir[ci]))
  phenoMeta = get(load(phenoMeta_dir[ci]))
  mm = get(load(matrixCountAdj_dir[ci]))
  pc = get(load(phenoChild_dir[ci]))
  pc_ind = get(load(phenoChild_ind_dir[ci]))

  control = c("+_+", "+_Y")
  if (!grepl("Sanger",paste0(panelL,centreL)[ci])) control = "WildType"
  
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

    pngname = paste0(plot_dir[ci], "/qn_all_", dpp[dp], ".png")
    png(filename=pngname, width=2*700, height=700)
    par(mfrow=c(1,2), mar=rep(2,4))
    
    colramp = colorRampPalette(c(3,"white",2))(nosample)
    s = sample((1:nrow(sm)),nosample)
    plot(density(log(m[s[1],]),na.rm=T),col=colramp[1],lwd=3, main=paste0("Quantile Normalization ln() BEFORE ",dpp[dp],"; countThres ", countThres, "; "))#,ylim=c(0,.30))
    for(i in 2:nosample){lines(density(log(m[s[i],]),na.rm=T),lwd=3,col=colramp[i])}
    
    plot(density(log(m_qn[s[1],]), na.rm=T),col=colramp[1],lwd=3, main=paste0("AFTER : ",dpp[dp],"; countThres ", countThres, "; "))#,ylim=c(0,.20))
    for(i in 2:nosample){lines(density(log(m_qn[s[i],]),na.rm=T),lwd=3,col=colramp[i])}

    graphics.off()
    TimeOutput(start1)
  }
  
}

TimeOutput(start)






