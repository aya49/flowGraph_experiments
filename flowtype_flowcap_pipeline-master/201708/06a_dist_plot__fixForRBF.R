# plot dist samples
# aya43@sfu.ca 20161220

root = "~/projects/flowCAP-II"
result_dir = "result"
setwd(root)


options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/matrix", sep="")
matrix_type = c("CountAdj", "Count", "Prop")
dist_dir = paste(result_dir, "/dist", sep="")
sampleCountThres = 3 #only compare if >=3 samples available

savefit = T # save the md reduction
ignoredist = ".csv|simmatrix"


libr(stringr)
libr(colorspace)
libr(fastcluster)
libr(dendextend)
libr(circlize)
libr(Rtsne)
libr(MASS)
libr(RDRToolbox)
libr(foreach)
libr(doMC)
libr(lubridate) #if there are date variables
source("~/projects/IMPC/code/_funcAlice.R")

no_cores = detectCores() - 1
registerDoMC(no_cores)

doHC = F
mds = c("tsne") #"mds","iso",
theta=.5

link = c("ward.D", "ward.D2", "mcquitty") #, "single", "complete") # , "median", "average", "centroid"

interested = c("aml","tube") #,"specimen",) #sampleMeta columns to plot

controlCol = NULL #if there are plots where you want only to plot controls of a specific column #
control = NULL #name of control factor in control column


#Output
plot_dir = paste(result_dir, "/plots", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dir[i])) }
plot_dist_dir = paste(plot_dir, "/dist", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dist_dir[i])) }
plot_dist_source_dir = paste0(plot_dist_dir, "/source"); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dist_source_dir[i])) }




start = Sys.time()






sampleMeta0 = get(load(sampleMeta_dir))
distmfile = list.files(dist_dir, recursive=F, full.names=T, pattern=".Rdata")
delpvalind = grepl(ignoredist,distmfile,ignore.case=T)
if (sum(delpvalind)>0) distmfile = distmfile[!delpvalind]
distmfile = distmfile[grepl("Freqp0.5_orig",distmfile)]
distmfilenames = fileNames(distmfile)
distmfilelayer = str_extract(distmfile, "layer[0-9]+")

# No of factors in each sampleMeta column
# sm = sapply(sampleMeta, function(x) as.numeric(factor(x, ordered=T)))

#each distance file will plot out all interested columns, set dimensions for the plot
width = 700; height = 700 # of each plot in png
hcwidth = 1000; hcheight = 7000
legperplot = 12
perplot = height/legperplot # legend items per plot



a = foreach(i=order(distmfilelayer,decreasing=T)) %dopar% {
  tryCatch({
    
    d = get(load(distmfile[i])); if (is.null(dim(as.matrix(d)))) return(NULL)
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
      dsname = paste0(plot_dist_source_dir, "/", mds[mdsi], "_", fileNames(distmfile[i]))
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
      
      pngname = paste0(plot_dist_dir, "/", mds[mdsi], "_", fileNames(distmfile[i]), ".png")
      png(pngname, width=width*colplot, height=height*rowplot)
      par(mfrow=c(rowplot,colplot))
      
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
            plotTsne(X[which(sampleMeta[,controlCol]%in%ftWTGT0),], continuous=T, main=paste(mds[mdsi], " (",distmfile[i],") distance matrix; \ntheta=",theta,"; \ndays since ", min(sampleMeta[,col]),"; control only", sep=""))
          }
          plotTsne(X, continuous=T, main=paste(mds[mdsi], " (",distmfile[i],") distance matrix; \ntheta=",theta,"; \ndays since ", min(sampleMeta[,col]), sep=""))
        } else {
          #plotTsne(X[sampleMeta$gene%in%c(control,as.character(ftKOGT)[which(varno>sampleNo)]),], leg=leg, main=paste(mds[mdsi], " (",distmfile[i],") distance matrix; \ntheta=",theta,"; ", colnames(sampleMeta)[col], " with ",sampleNo,"< samples", sep=""))
          plotTsne(X, leg=T, main=paste(mds[mdsi], " (",distmfile[i],") distance matrix; \ntheta=",theta,"; ", colnames(sampleMeta)[col], sep=""))
          if (length(var)>legperplot/2) plotTsne(X, leg=F, main=paste(mds[mdsi], " (",distmfile[i],") distance matrix; \ntheta=",theta, sep=""))
          #if too many factors, plot one without legend, and plot rest of them with legend
          if (numplot[ic]>2) {
            for (j in 1:((numplot[ic]-2)/2)) {
              e = min(j*perplot, length(unique(sampleMeta[,col])))
              gt = unique(sampleMeta[,col])[c(((j-1)*perplot+1):e)]
              plotTsne(X[sampleMeta[,col]%in%gt,], leg=T, main=paste(mds[mdsi], " (",distmfile[i],") distance matrix; theta=",theta, sep=""))
              plotTsne(X[sampleMeta[,col]%in%gt,], leg=F, main=paste(mds[mdsi], " (",distmfile[i],") distance matrix; theta=",theta, sep=""))
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
        pngname = paste0(plot_dist_dir, "/hclust_", fileNames(distmfile[i]), "_", link[j], ".png")
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
  }, error = function(err) {
    cat(paste("ERROR:  ",err))
  })
  
}



TimeOutput(start)













## AMI evaluation sci.kit.learn ==FAST on python, not good on R NMI

## libr (fpc)

#set.seed(4634)
#face <- rFace(300,dMoNo=2,dNoEy=0)
#grface <- as.integer(attr(face,"grouping"))
#plotcluster(face,grface==1)
#plotcluster(face,grface, clnum=1, method="vbc")








