# Plot Distance + its NCA
# aya43@sfu.ca 20161220

#Directory
root = "~/projects/IMPC"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#,"Sanger_MLN","CIPHE","TCP","H")

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrix", sep="")
matrix_type = c("CountAdj", "Count", "Prop")
dist_dir = paste(result_dir, "/", panelL, "/", centreL, "/dist", sep="")
dist_score_dir = paste(dist_dir, "/dist_score", sep=""); for (i in 1:length(dist_score_dir)) { suppressWarnings(dir.create(dist_score_dir[i])) }
dist_clusterf_dir = paste(dist_score_dir, "/dist_cluster_f", sep=""); for (i in 1:length(dist_clusterf_dir)) { suppressWarnings(dir.create(dist_clusterf_dir[i])) }


sampleCountThres = 3 #only compare if >=3 samples available

savefit = F


source("~/projects/IMPC/code/_funcAlice.R")
libr("stringr")
libr("colorspace")
libr("fastcluster")
libr("dendextend")
libr("circlize")
libr("Rtsne")
libr("MASS")
libr("RDRToolbox")
libr("foreach")
libr("doMC")
libr("lubridate") #if there are date variables

no_cores = detectCores() - 1
registerDoMC(no_cores)

doHC = F
mds = c("tsne") #"mds","iso",
theta=.5

link = c("ward.D", "ward.D2", "mcquitty") #, "single", "complete") # , "median", "average", "centroid"

interested = c("gene", "gender", "date", "colony", "strain", "birth_date", "specimen", "sample") # "fur", ; #sampleMeta columns to plot

controlCol = NULL #if there are plots where you want only to plot controls of a specific column #
control = NULL #name of control factor in control column


#Output
plot_dir = paste(result_dir, "/", panelL, "/", centreL, "/plots", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dir[i])) }
plot_dist_dir = paste(plot_dir, "/dist", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dist_dir[i])) }
plot_dist_source_dir = paste0(plot_dist_dir, "/source"); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dist_source_dir[i])) }




start = Sys.time()




for (ci in 1:length(paste0(panelL,centreL))) {
  start2 = Sys.time()
  centre = paste0(panelL,centreL)[ci]
  cat("\n",paste0(panelL," ",centreL)[ci])
  
  sampleMeta0 = get(load(sampleMeta_dir[ci]))
  distmfile = list.files(dist_dir[ci], recursive=F, full.names=T, pattern=".Rdata")
  delpvalind = grep(".csv",distmfile,ignore.case=T)
  if (length(delpvalind)>0) distmfile = distmfile[-delpvalind]
  distmfilenames = fileNames(distmfile)
  
  # No of factors in each sampleMeta column
  # sm = sapply(sampleMeta, function(x) as.numeric(factor(x, ordered=T)))
  
  #each distance file will plot out all interested columns, set dimensions for the plot
  width = 400; height = 400 # of each plot in png
  hcwidth = 700; hcheight = 400
  legperplot = 12
  perplot = height/legperplot # legend items per plot
  
  
  
  a = foreach(i=1:length(distmfile)) %dopar% {
    tryCatch({
      fm = get(load(paste0(dist_clusterf_dir[ci],"/",fileNames(distmfile[i]))))
      a = fm[[1]][[1]][[1]][[1]][c("NCA","silmed","pearsongamma")]
      scoress = paste0("\n",paste0(names(a),"-",signif(a,digits=4),collapse="   "))
      
      d = get(load(distmfile[i]))
      dcol = rep(0,ncol(sampleMeta0))
      for (j in 1:ncol(sampleMeta0)) {
        l = sum(!is.na(match(colnames(as.matrix(d)),sampleMeta0[,j])))
        if (l>0) dcol[j] = l
      }
      dcol = which.max(dcol)
      if (length(unique(sampleMeta0[,dcol]))==nrow(sampleMeta0)) {
        sampleMeta <- sampleMeta0[match(colnames(as.matrix(d)),sampleMeta0[,dcol]),]
        interestedCols = which(colnames(sampleMeta)%in%interested)
      } else {
        interestedCols = 1
        sampleMeta = as.data.frame(colnames(as.matrix(d)))
      }
      
      numplot = NULL
      ic = 0
      for (col in interestedCols) {
        ic = ic+1
        numplot[ic] = 2*ceiling( length(unique(sampleMeta0[,col]))/perplot ) # each legend row takes up 18 pixels; prevent legend from overspilling
        if (numplot[ic]>2)numplot[ic] = 2+numplot[ic]
        if (grepl("date", colnames(sampleMeta0)[col], ignore.case=T)) numplot[ic] = 2
      }
      numplots = sum(numplot)
      rowplot = 3
      colplot = ceiling(numplots/rowplot) # all genes + WT only, >3 samples
      
      for (mdsi in 1:length(mds)) {
        dsname = paste0(plot_dist_source_dir[ci], "/", mds[mdsi], "_", fileNames(distmfile[i]))
        if (T) {
          #if (!file.exists(dsname)) {
          if (mds[mdsi]=="mds") { fit = cmdscale(d,eig=TRUE, k=2)
          } else if (mds[mdsi]=="iso") { fit = isoMDS(d, k=2)
          } else if (mds[mdsi]=="tsne") { fit = Rtsne(d, is_distance=T, theta=theta) }
        } else {
          fit = get(load(dsname))
        }
        if (sum(c("mds","iso")==mds[mdsi])>0) { X = fit$points[,c(1,2)] #Metric MDS fit$points[,c(1,2)]
        } else if (mds[mdsi]=="tsne") { X = fit$Y }
        rownames(X) = rownames(as.matrix(d))
        if (savefit) save(X, file=dsname)
        
        pngname = paste0(plot_dist_dir[ci], "/", mds[mdsi], "_", fileNames(distmfile[i]), ".png")
        png(pngname, width=width*colplot, height=height*rowplot)
        par(mfrow=c(rowplot,colplot),mar=c(1,1,10,1))
        
        #sampleNo = 3 #only consider variables with more than sampleNo samples per factor
        ic = 0
        for (col in interestedCols) {
          ic = ic+1
          varno = table(sampleMeta[,col])
          var = names(varno)
          rownames(X) = sampleMeta[,col]
          if (grepl("date", colnames(sampleMeta)[col], ignore.case=T)) {
            if (!is.null(control)) {
              ftWTGT0 = unique(sampleMeta[,controlCol,grep(control, sampleMeta[,controlCol])])
              plotTsne(X[which(sampleMeta[,controlCol]%in%ftWTGT0),], continuous=T, main=paste(mds[mdsi], " (",fileNames(distmfile[i]),"); \ntheta=",theta,"; \ndays since ", min(sampleMeta[,col]),"; control only",scoress, sep=""))
            }
            plotTsne(X, continuous=T, main=paste(mds[mdsi], " (",fileNames(distmfile[i]),"); \ntheta=",theta,"; \ndays since ", min(sampleMeta[,col]), scoress, sep=""))
          } else {
            #plotTsne(X[sampleMeta$gene%in%c(control,as.character(ftKOGT)[which(varno>sampleNo)]),], leg=leg, main=paste(mds[mdsi], " (",fileNames(distmfile[i]),"); \ntheta=",theta,"; ", colnames(sampleMeta)[col], " with ",sampleNo,"< samples", sep=""))
            plotTsne(X, leg=T, main=paste(mds[mdsi], " (",fileNames(distmfile[i]),"); \ntheta=",theta,"; ", colnames(sampleMeta)[col], scoress, sep=""))
            if (length(var)>legperplot/2) plotTsne(X, leg=F, main=paste(mds[mdsi], " (",fileNames(distmfile[i]),"); \ntheta=",theta, scoress, sep=""))
            #if too many factors, plot one without legend, and plot rest of them with legend
            if (numplot[ic]>2) {
              for (j in 1:((numplot[ic]-2)/2)) {
                e = min(j*perplot, length(unique(sampleMeta[,col])))
                gt = unique(sampleMeta[,col])[c(((j-1)*perplot+1):e)]
                plotTsne(X[sampleMeta[,col]%in%gt,], leg=T, main=paste(mds[mdsi], " (",fileNames(distmfile[i]),"); \ntheta=",theta, scoress, sep=""))
                plotTsne(X[sampleMeta[,col]%in%gt,], leg=F, main=paste(mds[mdsi], " (",fileNames(distmfile[i]),"); \ntheta=",theta, scoress, sep=""))
              }
            }
          }
        }
        graphics.off()
      }
    
    
    ## hclust ------------------------------------------------
    if (doHC) { cat(" HClust")
      for (j in 1:length(link)) { cat(link(j)," ",sep="") #all links
        #Hclust
        try({
          hc = hclust(d, method=link[j])
          dend = as.dendrogram(hc)
          dend = rotate(dend, 1:nrow(sampleMeta))
        })
        
        #plot
        pngname = paste0(plot_dist_dir[ci], "/hclust_", fileNames(distmfile[i]), "_", link[j], ".png")
        png(file=pngname , width=hcwidth*length(interestedCols), height=hcheight)
        par(mar=(c(5,5,5,40) + 0.1), mfrow=c(1,length(interestedCols)))
        
        dend0 = dend
        for (col in interestedCols) {
          dend = dend0
          var = factor(sampleMeta[,col])
          suppressWarnings({ labels_colors(dend) = rainbow_hcl(length(rev(levels(var))))[ sort_levels_values( as.numeric(var)[order.dendrogram(dend)] ) ]  }) 
          labels(dend) = paste(as.character(var)[order.dendrogram(dend)], sep = "")
          dendh = hang.dendrogram(dend,hang_height=0.1)
          plot(dendh, main = paste("aml samples (labels = ", colnames(sampleMeta)[col],"); \ndist=",dis[i],"; link=",link[j],sep=""), horiz=TRUE,  nodePar=list(cex = .007))
          # dend = color_branches(dend, k=length(gene))
          # legend("topleft", legend = gene, fill = rainbow_hcl(length(gene)))
        }
      }
      graphics.off()
    }
    }, error = function(e) {
      cat(paste("ERROR:  ",e));
    })
    
  }
  
  cat("\n centre ", centre, " ",TimeOutput(start2)," \n",sep="")
}

TimeOutput(start)













## AMI evaluation sci.kit.learn ==FAST on python, not good on R NMI

## libr (fpc)

#set.seed(4634)
#face <- rFace(300,dMoNo=2,dNoEy=0)
#grface <- as.integer(attr(face,"grouping"))
#plotcluster(face,grface==1)
#plotcluster(face,grface, clnum=1, method="vbc")








