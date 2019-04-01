# aya43@sfu.ca 20161220
# Uses different distance measures to calculate distance & plot samples (for pvalue matrices)

root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrix", sep="")
matrix_type = list.files(result_dir, pattern="Pval", full.names=T)
matrix_type = matrix_type[-grep("FULL",matrix_type)]
matrix_type = gsub("matrix","",fileNames(matrix_type, ext="Rdata")) 


#Output
dist_dir = paste(result_dir, "/dist", sep=""); suppressWarnings(dir.create (dist_dir))

libr(stringr)
libr(colorspace)
libr(vegan) # libr(proxy)
libr(prabclus)
libr(fastcluster)
libr(dendextend)
libr(circlize)
libr(Rtsne)
libr(foreach)
libr(doMC)
libr(lubridate) #if there are date variables
source("code/_funcAlice.R")



dodist = F
doHC = F
doTsne = T
pvalThres = .05

dis = c("canberra","bray","kulczynski", "binomial","cao", "jaccardInd") # "jaccard","mahalanobis", "gower","altGower","morisita","horn", "euclidean", ""manhattan", "maximum", "binary", "minkowski", "raup", "mountford" (all at same level),)

start = Sys.time()

no_cores = detectCores() - 2
registerDoMC(no_cores)

load(sampleMeta_dir)
load(phenoMeta_dir)

k0 = c(1,max(phenoMeta[,3])+1) # how many markers to consider i.e. k=max(phenolevel) only

for (mcp in matrix_type) { # Load & fix pvalue matrices
  loop.ind = 1:length(dis)
  result = foreach(i = loop.ind) %dopar% { #for each phenotype
    
    for (k in k0) {
      start1 = Sys.time()
      
      dname = paste(dist_dir[ci], "/distance_", mcp, "_layer-", str_pad(k, 2, pad = "0"), "_", dis[i], ".Rdata", sep = "" )
      
      if (file.exists(dname) & !dodist) { cat("\nLoading Distance Object #", length(dis)-i+1, " ", dis[i]," ",sep="")
        d <- get(load(dname))
      } else { cat("\nCalculating Distance Object ", dis[i], sep="")
        mm = get(load(paste0(matrix_dir[ci], mcp,".Rdata")))
        phenoLevel = str_count(colnames(mm), "[+-]")
        m = mm[,which(phenoLevel<=k)]
        sigInd = which(m<pvalThres)
        insigInd = which(m>=pvalThres)
        
        start2 = Sys.time()
        #Dist matrix
        if (dis[i]=="jaccardInd") {
          m[sigInd] = 1
          m[insigInd] = 0
          d = jaccard(t(m))
        } else {
          m[sigInd] = log(m[sigInd])
          m[insigInd] = -1
          m = -m
          d = vegdist(m, method=dis[i])
        }
        save(d,file=dname)
        cat(" ",TimeOutput(start2), sep="")
      }
      
      if (doTsne) {
        cat("; Tsne theta: ")
        for (theta in c(0,.5)) { cat(theta," ", sep="")
          tryCatch({
            
            tsne = Rtsne(d, is_distance=T, theta=theta)
            colnames(tsne$Y) = c("x","y")
            # tsnem = Rtsne(m)
            # rownames(tsnem$Y) = sampleMeta$gene
            # palette = choose_palette(pal=rainbow_hcl, n=length(unique(sampleMeta$gene)))
            
            width = 700; height = 700 # of each plot in png
            perplot = height/12 # legend items per plot
            
            numplots = NULL
            numplots = 2*ceiling( length(unique(rownames(m)))/perplot ) # each legend row takes up 18 pixels; prevent legend from overspilling
            if (numplots>2) {
              numplots = 2+numplots
            }
            rowplot = 3
            colplot = ceiling(numplots/rowplot) # all genes + WT only, >3 samples
            
            pngname = paste0(dist_dir[ci], "/tsne_", mcp, "_layer-", str_pad(k, 2, pad = "0"), "_", dis[i], "_tsne",theta, ".png")
            
            png(pngname, width=width*colplot, height=height*rowplot)
            par(mfrow=c(rowplot,colplot))
            
            var = unique(rownames(m))
            names(var) = var
            varno = unlist(lapply(var, function(x) length(which(rownames(m)==x)))) # of samples per factor
            rownames(tsne$Y) = rownames(m)
            plotTsne(tsne$Y, leg=T, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta,"; speciman", sep=""), text=T)
            plotTsne(tsne$Y, leg=F, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta, sep=""))
            #if too many factors, plot one without legend, and plot rest of them with legend
            if (numplots>2) {
              for (j in 1:((numplots-2)/2)) {
                e = min(j*perplot, length(unique(rownames(m))))
                gt = unique(rownames(m))[c(((j-1)*perplot+1):e)]
                plotTsne(tsne$Y[rownames(m)%in%gt,], leg=T, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta, sep=""))
                plotTsne(tsne$Y[rownames(m)%in%gt,], leg=F, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta, sep=""))
              }
            }
            graphics.off()
          }, error = function(err) { graphics.off()})
        }
      }
      
    }
  }
}

TimeOutput(start)
