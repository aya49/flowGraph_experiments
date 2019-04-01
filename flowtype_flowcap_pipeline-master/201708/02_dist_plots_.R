# aya43@sfu.ca 20161220
# Uses different distance measures to calculate distance & plot samples

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
matrix_type = c("Count", "CountAdj", "Prop")

#Output
dist_dir = paste(result_dir, "/dist", sep=""); suppressWarnings(dir.create (dist_dir))

libr(stringr)
libr(colorspace)
libr(vegan) # libr(proxy)
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

dis = c("canberra","bray","kulczynski", "binomial","cao", "canberra2") # "jaccard","mahalanobis", "gower","altGower","morisita","horn", "euclidean", ""manhattan", "maximum", "binary", "minkowski", "raup", "mountford" (all at same level),)
link = c("ward.D", "ward.D2", "mcquitty") #, "single", "complete") # , "median", "average", "centroid"

interestedCols = c(3,4) #from sampleMeta, mainly for plotting colours in RTsne


start = Sys.time()

no_cores = detectCores() - 2
registerDoMC(no_cores)

load(sampleMeta_dir)
load(phenoMeta_dir)

k0 = c(1,max(phenoMeta[,3])+1) # how many markers to consider i.e. k=max(phenolevel) only

for (mcp in matrix_type) { # Load & fix cell count/countAdj/proportion matrix; 1:3 are regular matrices, 4... are pvalues
  
  loop.ind = 1:length(dis)
  result = foreach(i = loop.ind) %dopar% { #for each phenotype
    #for (i in 1:length(dis)) {
    
    start1 = Sys.time()
    
    
    for (k in k0) {
      
      dname = paste(dist_dir[ci], "/distance_", mcp, "_layer-", str_pad(k, 2, pad = "0"), "_", dis[i], ".Rdata", sep = "" )
      
      if (file.exists(dname) & !dodist) { cat("\nLoading Distance Object #", length(dis)-i+1, " ", dis[i]," ",sep="")
        d <- get(load(dname))
      } else { cat("\nCalculating Distance Object #", length(dis)-i+1, " ", dis[i], sep="")
        mm = get(load(paste0(matrix_dir[ci], mcp,".Rdata")))
        if (nrow(mm)!=nrow(sampleMeta)) { mm = t(mm) } #phenotype on cols
        if (grepl("Count",mcp)) {
          delPhenoInd = colIndBelow(mm, countThres)
          mm = mm[,-delPhenoInd]
          m = mm[,which(phenoMeta[-delPhenoInd,3]<=k)]
        } else {
          m = mm[,which(phenoMeta[,3]<=k)]
        }
        
        start2 = Sys.time()
        #Dist matrix
        if (dis[i]=="canberra2") {
          m[which(m<1)] = 1
          d = vegdist(m, method="canberra")
        } else {
          d = vegdist(m, method=dis[i])
        }
        save(d,file=dname)
        cat(" ",TimeOutput(start2), sep="")
      }
      
      start2 = Sys.time()
      
      ## create hclust for distance object i
      if (doHC) {
        cat("; HClust: ")
        for (j in 1:length(link)) { cat(length(link)-j+1," ",sep="") #all links
          #Hclust
          try({
            hc = hclust(d, method=link[j])
            dend = as.dendrogram(hc)
            dend = rotate(dend, 1:nrow(sampleMeta))
          })
          
          #plot
          numplots = length(interestedCols)
          
          pngname = paste0(dist_dir[ci], "/hclust_", mcp, "_layer-", str_pad(k, 2, pad = "0"), "_", dis[i], "_", link[j], ".png")
          
          png (file=pngname , width=1000*numplots, height=8000)
          par(mar=(c(5,5,5,40) + 0.1), mfrow=c(1,numplots))
          
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
      
      
      ## rtsne plot for distance object i
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
            
            numplot = NULL
            ic = 0
            for (col in interestedCols) {
              ic = ic+1
              numplot[ic] = 2*ceiling( length(unique(sampleMeta[,col]))/perplot ) # each legend row takes up 18 pixels; prevent legend from overspilling
              if (numplot[ic]>2) {
                numplot[ic] = 2+numplot[ic]
              }
              # if (grepl("date", colnames(sampleMeta)[col], ignore.case=T)) {
              #   numplot[ic] = 2
              # }
            }
            numplots = sum(numplot)
            
            rowplot = 3
            colplot = ceiling(numplots/rowplot) # all genes + WT only, >3 samples
            
            pngname = paste0(dist_dir[ci], "/tsne_", mcp, "_layer-", str_pad(k, 2, pad = "0"), "_", dis[i], "_tsne",theta, ".png")
            
            png(pngname, width=width*colplot, height=height*rowplot)
            par(mfrow=c(rowplot,colplot))
            
            #sampleNo = 3 #only consider variables with more than sampleNo samples per factor
            ic = 0
            for (col in interestedCols) {
              ic = ic+1
              var = unique(sampleMeta[,col])
              names(var) = var
              varno = unlist(lapply(var, function(x) length(which(sampleMeta[,col]==x)))) # of samples per factor
              rownames(tsne$Y) = sampleMeta[,col]
              # if (grepl("date", colnames(sampleMeta)[col], ignore.case=T)) {
              #   ftWTGT0 = unique(sampleMeta$gene[grep(ftWTGT, sampleMeta$gene)])
              #   plotTsne(tsne$Y[which(sampleMeta$gene%in%ftWTGT0),], continuous=T, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta,"; \ndays since ", min(sampleMeta[,col]),"; WT only", sep=""))
              #   plotTsne(tsne$Y, continuous=T, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta,"; \ndays since ", min(sampleMeta[,col]), sep=""))
              # } else {
              #plotTsne(tsne$Y[sampleMeta$gene%in%c(ftWTGT,as.character(ftKOGT)[which(var>3)]),], leg=leg, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta,"; ", colnames(sampleMeta)[col], " with ",sampleNo,"< samples", sep=""))
              plotTsne(tsne$Y, leg=T, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta,"; ", colnames(sampleMeta)[col], sep=""))
              plotTsne(tsne$Y, leg=F, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta, sep=""))
              #if too many factors, plot one without legend, and plot rest of them with legend
              if (numplot[ic]>2) {
                for (j in 1:((numplot[ic]-2)/2)) {
                  e = min(j*perplot, length(unique(sampleMeta[,col])))
                  gt = unique(sampleMeta[,col])[c(((j-1)*perplot+1):e)]
                  plotTsne(tsne$Y[sampleMeta[,col]%in%gt,], leg=T, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta, sep=""))
                  plotTsne(tsne$Y[sampleMeta[,col]%in%gt,], leg=F, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta, sep=""))
                }
                # }
              }
            }
            
            # rownames(tsne$Y) = ymd(sampleMeta$date)
            # plotTsne(tsne$Y[which(sampleMeta$gene%in%ftWTGT),], continuous=T, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta,"; \ndays since ", min(sampleMeta$date),"; WT only", sep=""))
            # plotTsne(tsne$Y[sampleMeta$gene%in%c(ftWTGT,as.character(ftKOGT)[which(ftKOIndexL>3)]),], continuous=T, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta,"; \ndays since ", min(sampleMeta$date),"; gene with ",sampleNo,"< samples", sep=""))
            # plotTsne(tsne$Y, continuous=T, main=paste("TSNE ",dis[i]," distance matrix; theta=",theta, "\n days since ", min(sampleMeta$date), sep=""))
            # 
            # rownames(tsne$Y) = as.character(sampleMeta[,3])
            # plotTsne(tsne$Y, main=paste("TSNE ",dis[i]," distance matrix; gender", sep=""))
            
            graphics.off()
            
            # require(cluster)
            # require(geometry)
            
            #           cNo = 7
            #           km = kmeans(tsne$Y,cNo)
            #           for (i in 1:cNo) {
            #             line(convhulln(tsne$Y[which(km$cluster==i),]))
            #           }
            #           convhulln
            #           
            #           require(RDRToolbox)
            #           isomap = Isomap(data=m, dims=2, k=10)
            #           ecb(isomap$dim2, unique(sampleMeta$gene), paste("Isomap ",dis[i]," distance matrix; theta=0", sep="") )
            #           ecb(isomap$dim2, ftWTGT, paste("Isomap ",dis[i]," distance matrix; theta=0; WildType only", sep=""))
            
          }, error = function(err) { graphics.off()})
        }
      }
      
      cat(" viz ",TimeOutput(start2)," ",sep="")
    }
    cat(" dis ",TimeOutput(start1)," ",sep="")
    
    
    
  }
  
}


TimeOutput(start)



#require(RDRToolbox)
#libr(rgl)
## Isomap ##



# libr(dbscan)
# oc = optics(d, eps=10, minPts=3) # no results (gravitate towards each other)
