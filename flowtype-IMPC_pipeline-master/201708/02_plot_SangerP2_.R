# Plot: plot samples & divide dates up for p value calculation
# aya43@sfu.ca 20161220

root = "~/projects/IMPC/SangerP2"
result_dir = "result"
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN") #"Sanger_SPLEEN","Sanger_MLN","CIPHE","TCP",

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrix", sep="")
matrix_type = c("CountAdj") #, "Count", "Prop", "CountAdj_qn_all_leafonly500", "CountAdj_qn_all_nonleaf500", "CountAdj_all_leafonly500", "CountAdj_all_nonleaf500"
dist_dir = paste(result_dir, "/", panelL, "/", centreL, "/dist", sep="")
sampleCountThres = 3 #only compare if >=3 samples available

#Output
plot_dir = paste(result_dir, "/", panelL, "/", centreL, "/plots", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dir[i])) }
cp_dir = paste(plot_dir, "/changepoint", sep=""); for(i in 1:length(cp_dir)) { suppressWarnings(dir.create(cp_dir[i])) }
single_dir = paste(plot_dir, "/singlephen", sep=""); for(i in 1:length(single_dir)) { suppressWarnings(dir.create(single_dir[i])) }


libr(stringr)
libr(colorspace)
libr(changepoint) # libr(proxy)
libr(FKF)
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
source("~/projects/IMPC/code/_funcdist.R")


countThres = 20
plotpc = 5 #number of pca pc to plot
doHC = F
doISO = T
doTsne = T
theta=.5 #for Tsne
dofullPCA = F #do PCA for all cell popoulations not just first level

dis = c("canberra", "canberra2","bray","kulczynski", "binomial","cao") # "jaccard","mahalanobis", "gower","altGower","morisita","horn", "euclidean", ""manhattan", "maximum", "binary", "minkowski", "raup", "mountford" (all at same level),)
methods = c("BinSeg","AMOC","PELT") #Changepoint analysis; AMOC at most one change, PELT is poisson (quick exact), BinSeg (quick, approx), SegNeigh (slow, exact)
usemethod="AMOC" #use this method to determine changepoints for pvalue calculation
mds_type = c("iso", "mds")
link = c("ward.D", "ward.D2", "mcquitty") #, "single", "complete") # , "median", "average", "centroid"

interested = c("gene", "gender", "date", "colony", "strain", "fur", "birth_date", "specimen", "sample") #sampleMeta columns to plot



start = Sys.time()

no_cores = 9#detectCores() - 1
registerDoMC(no_cores)

#for (ci in 4:1) {
for (ci in 1:length(paste0(panelL,centreL))) {
  start2 = Sys.time()
  centre = paste0(panelL,centreL)[ci]
  cat("\n",paste0(panelL," ",centreL)[ci])
  
  load(sampleMeta_dir[ci])
  #interestedCols0 = which(colnames(sampleMeta)%in%interested)
  
  #order samples by date, exclude genotypes with less than 3 samples
  if (length(grep("barcode", colnames(sampleMeta)))>0) sampleMeta2 <- sampleMeta[order(sampleMeta$barcode),] #sanger centre
  if (length(grep("sample", colnames(sampleMeta)))>0) sampleMeta2 <- sampleMeta[order(sampleMeta$sample),] #other centres
  if (length(grep("specimen", colnames(sampleMeta)))>0) sampleMeta2 <- sampleMeta[order(sampleMeta$specimen),] #other centres
  sampleMeta2 <- sampleMeta2[order(sampleMeta2$date),]
  
  ftWTGT = c("+_+","+_Y"); if (!grepl("Sanger",centre)) ftWTGT = "WildType"
  g = getGTindex(sampleMeta2$gene, ftWTGT, sampleCountThres)
  goodind = as.vector(c(unlist(g$expIndex[g$goodexpIndex]), g$controlIndex))
  allind = c(1:nrow(sampleMeta2))
  # #only keep samples with >3 samples
  # if (length(goodind)<length(allind)) {
  #   delind = sort(setdiff(allind, goodind))
  #   sampleMeta2 = sampleMeta2[-delind,]
  #   g = getGTindex(sampleMeta2$gene, ftWTGT, sampleCountThres)
  # }
  morder = match(sampleMeta2$fileName, sampleMeta$fileName)
  #if (sum(is.na(morder)>0)) morder = morder[!is.na(morder)]
  ko = g$expIndex
  wt = g$controlIndex #wildtype
  
  #prepare plot colours by date
  ts = as.numeric(factor(sampleMeta2$date))
  tswt = ts[wt]
  tsc = c(1,1+which(diff(ts)!=0))
  tscwt = c(1,1+which(diff(tswt)!=0))
  colour = heat.colors(length(unique(ts))+25)[1:length(unique(ts))]
  
  for (mcp in matrix_type) { cat("\n  ", mcp, ": ")
    start2 = Sys.time()
    
    sm0 = sapply(sampleMeta2, function(x) as.numeric(factor(x, ordered=T)))
    uniquecols = apply(sm0, 2, function(x) nrow(sm0)==length(unique(x)))
    if (grepl("Sanger",centre)) { sm = sm0[,!uniquecols]
    } else { sm = sm0[,!uniquecols] }
    
    ## get interested columns
    interestedCols = which(colnames(sm)%in%interested)
    
    mm = get(load(paste0(matrix_dir[ci], mcp,".Rdata")))
    
    ## get only 1st layer phenotypes
    phenoLevel = str_count(colnames(mm), "[+-]")
    if (length(which(phenoLevel<=1))>0) {
      m = mm[,which(phenoLevel<=1)]
      #colnames(m)[1] = "all"
      cols = seq(from=1,to=ncol(m),by=2)
      
      ## create line plots for WT only (Kalman filtering), for use in plots -----------------------------------------
      cat("Kalman filtering for WT; ")
      fkffit = NULL
      statsfit = NULL
      for (i in cols) {
        ywt = as.numeric(m[morder,i])[wt]
        if (length(unique(ywt))==1) {
          fkffit[[i]] <- statsfit[[i]] <- ywt
        } else {
          ## Set constant parameters:
          dt <- ct <- matrix(0) 
          Zt <- Tt <- matrix(1)
          a0 <- ywt[1]           # Estimation of the first sample count
          P0 <- matrix(100)     # Variance of 'a0'
          ## Estimate parameters 23min if TS
          fit.fkf <- optim(c(HHt = var(ywt, na.rm = TRUE) * .5,
                             GGt = var(ywt, na.rm = TRUE) * .5),
                           fn = function(par, ...) -fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik,
                           yt = rbind(ywt), a0 = a0, P0 = P0, dt = dt, ct = ct,
                           Zt = Zt, Tt = Tt, check.input = FALSE)
          ## Filter Nile data with estimated parameters:
          fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, HHt = matrix(fit.fkf$par[1]), GGt = matrix(fit.fkf$par[2]), yt = rbind(ywt))
          fkffit[[i]] <- fkf.obj$att[1,]
          ## Compare with the stats' structural time series implementation: 5min if TS
          statsfit[[i]] <- fitted(StructTS(ywt, type = "level"))
        }
      }
      
      cat("Kalman filtering for all; ")
      fkffitall = NULL
      statsfitall = NULL
      for (i in cols) {
        yall = as.numeric(m[morder,i])
        if (length(unique(yall))==1) {
          fkffitall[[i]] <- statsfitall[[i]] <- yall
        } else {
          ## Set constant parameters:
          dt <- ct <- matrix(0) 
          Zt <- Tt <- matrix(1)
          a0 <- yall[1]           # Estimation of the first sample count
          P0 <- matrix(100)     # Variance of 'a0'
          ## Estimate parameters 23min if TS
          fit.fkf <- optim(c(HHt = var(yall, na.rm = TRUE) * .5,
                             GGt = var(yall, na.rm = TRUE) * .5),
                           fn = function(par, ...) -fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik,
                           yt = rbind(yall), a0 = a0, P0 = P0, dt = dt, ct = ct,
                           Zt = Zt, Tt = Tt, check.input = FALSE)
          ## Filter Nile data with estimated parameters:
          fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, HHt = matrix(fit.fkf$par[1]), GGt = matrix(fit.fkf$par[2]), yt = rbind(yall))
          fkffitall[[i]] <- fkf.obj$att[1,]
          ## Compare with the stats' structural time series implementation: 5min if TS
          statsfitall[[i]] <- fitted(StructTS(yall, type = "level"))
        }
      }
      
      
      
      ## plot all single phenotypes onto an image ------------------------------------------------
      cat("single phen;")
      pngname <- paste0(single_dir[ci], "/singlePhen_", mcp, ".png")
      png(filename=pngname, width=2*1000, height=length(cols)*400)
      par(mfrow=c(length(cols),2), mar=rep(2,4))
      for (i in cols) {
        y <- as.numeric(m[morder,i])
        ylim <- c(min(m[morder,]),max(m[morder,]))
        for (j in 1:2) {
          if (j==1) { #everything to same scale
            plot(y, main=paste0(colnames(m)[i], "; heat colours = days since ", min(sampleMeta$date), "; vertical lines = first sample on day"), col=colour[ts], pch=19,cex=.4, ylim=ylim)
          } else {
            plot(y, main=paste0(colnames(m)[i], "; heat colours = days since ", min(sampleMeta$date), "; vertical lines = first sample on day"), col=colour[ts], pch=19,cex=.4)
          }
          abline(v=tsc, col="#DCDCDC")
          points(wt ,y[wt], col="black")
          lines(wt, statsfit[[i]], col = "green")
          lines(wt, fkffit[[i]], col = "blue")
          legend("top", c("Actual datapoints (clack=WT)", "Local level (StructTS)", "Local level (fkf)"), col = c("red", "green", "blue"), lty = 1)
          legend.col(col=colour, lev=ts)
        }
      }
      graphics.off()
      
      
      
      ## plot one single phenotypes and its changepoints ------------------------------------------------
      cat(" changepoint; ")
      for (i in cols) {
        y <- as.numeric(m[morder,i])
        ywt <- y[wt]
        
        
        pngname <- paste0(cp_dir[ci], "/cp_", mcp, "_", colnames(m)[i], ".png")
        png(filename=pngname, width=length(methods)*800, height=(1+3)*400)
        layout(matrix(c(1,1,1, 2:(3*length(methods)+1)),ncol=length(methods),byrow=T))
        par(mar=rep(2,4))
        
        mvalueall <- cpt.mean(fkffitall[[i]],method=usemethod, Q=20, penalty="MBIC", minseglen=5)
        plot(mvalueall, main=paste0(colnames(m)[i], "; kalman filtered on WT only; heat colours = days since ", min(sampleMeta$date), "; vertical lines = first sample on day"))
        abline(v=tsc, col="#DCDCDC")
        points(y, col=colour[ts], cex=.4)
        points(wt ,ywt, cex=.4, col="blue")#pch=19, 
        lines(wt, statsfit[[i]], col = "green")
        lines(wt, fkffit[[i]], col = "blue")
        legend("top", c("Actual datapoints (blue=WT)", "Local level (StructTS)", "Local level (fkf)"), col = c("red", "green", "blue"), lty = 1)
        legend.col(col=colour, lev=ts)
        
        ylim <- c(min(ywt), max(ywt))
        for (j in 1:length(methods)) {
          mvalue <- cpt.mean(fkffit[[i]],method=methods[j], Q=20, penalty="MBIC", minseglen=5)
          if (methods[j]==usemethod & i==1 & mcp==matrix_type[1]) { mvaluewt = mvalue }
          plot(mvalue, main=paste("WT only; mean change: ",methods[j], "; Penalty MBIC; ",colnames(m)[i], sep=""), ylim=ylim)
          points(ywt, col=colour[tswt], pch=19,cex=.4)
          lines(statsfit[[i]], col = "green")
          lines(fkffit[[i]], col = "blue")
          legend.col(col=colour, lev=ts)
        }
        for (j in 1:length(methods)) {
          vvalue <- cpt.var(diff(fkffit[[i]]), method=methods[j], penalty="MBIC")
          plot(vvalue, main=paste("variance change: ",methods[j],sep=""))
          #points(wt ,ywt, col=colour[tswt], pch=19,cex=.4)
          legend.col(col=colour, lev=ts)
        }
        for (j in 1:length(methods)) {
          mvvalue <- cpt.meanvar(diff(fkffit[[i]]), method=methods[j], penalty="MBIC")
          plot(mvvalue, main=paste("variance/mean change: ",methods[j],sep=""))
          #points(wt ,ywt, col=colour[tswt], pch=19,cex=.4)
          legend.col(col=colour, lev=ts)
        }
        graphics.off()
      }
      
      
      ##Group days by change in mean; incomplete, this only works for AMOC, large amounts of dates grouped into one.
      if (mcp==matrix_type[1]) {
        group = rep(0,nrow(sampleMeta))
        enddate=0
        for (i in 1:length(mvaluewt@cpts)) {
          startdate = min(sampleMeta2$date[sampleMeta2$date>enddate])
          enddate = sampleMeta2$date[wt[mvaluewt@cpts[i]]]
          group[sampleMeta$date>=startdate & sampleMeta$date<=enddate] <- i
        }
        sampleMeta$group = group
        save(sampleMeta, file=sampleMeta_dir[ci])
      }      
      
    }
    
    
    
    
    
    ## pca analysis ------------------------------------------------
    iso <- 0
    if (doISO) { cat("iso; ")
      iso <- 1
      fit <- Isomap(m[morder,],k=2)
    }
    cat("pca; ")
    pc <- prcomp(m[morder,])
    pngname <- paste0(plot_dir[ci], "/pca_iso_", mcp, "_layer-", str_pad(1, 2, pad = "0"), ".png")
    png(filename=pngname, width=length(interestedCols)*400, height=(plotpc+1+iso)*400)
    layout(matrix(c(rep(1,length(interestedCols)), 2:((plotpc+iso)*length(interestedCols)+1)),ncol=length(interestedCols),byrow=T))
    
    plot(pc$sdev^2/sum(pc$sdev^2), xlim = c(0, 15), type = "b", pch = 16, xlab = "principal components", ylab = "variance explained")
    for (i in 1:plotpc) {
      for (col in interestedCols) {
        coloursm <- sm[,col]
        if (grepl("date",colnames(sm)[col])) coloursm = colour[as.numeric(factor(sm[,col]))]
        plot(pc$x[,i], pc$x[,i+1], col = coloursm, main = paste0("PCA ", colnames(sm)[col]), xlab = paste0("PC_",i), ylab = paste0("PC_",i+1))
        points(0, 0, pch = 3, cex = 4, lwd = 4)
      }
    }
    if (doISO) {
      for (col in interestedCols) {
        coloursm <- sm[,col]
        if (grepl("date",colnames(sm)[col])) coloursm <- colour[as.numeric(factor(sm[,col]))]
        plot(fit$dim2, col = coloursm, main = paste0("ISO ", colnames(sm)[col]))
        points(0, 0, pch = 3, cex = 4, lwd = 4)
      }
    }
    graphics.off()
    
    #same as SVD on centred data
    # cx <- sweep(cbind(sm,m[morder,]), 2, colMeans(x), "-")
    # sv <- svd(cx)
    
    
  }
  
  cat("\n centre ", centre, " ",TimeOutput(start2)," \n",sep="")
}

TimeOutput(start)





## See dates without WT
#Harwell
ftWTGT <- "FR-FCM-ZYCB-WildType_01_All"
KOdays <- table(sampleMeta[!sampleMeta$gene%in%ftWTGT,]$date)
colnames(sampleCountThres) <- table(sampleMeta[sampleMeta$gene%in%ftWTGT,]$date)
noWTdays <- KOdays [!(names(KOdays) %in% names(WTdays))]



noWTdaysInd <- which(as.character(sampleMeta2$date)%in%names(noWTdays))
plot()








## AMI evaluation sci.kit.learn ==FAST on python, not good on R NMI

## libr (fpc)

#set.seed(4634)
#face <- rFace(300,dMoNo=2,dNoEy=0)
#grface <- as.integer(attr(face,"grouping"))
#plotcluster(face,grface==1)
#plotcluster(face,grface, clnum=1, method="vbc")








