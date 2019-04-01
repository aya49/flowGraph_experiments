# PEER coundounding factor analysis
# aya43@sfu.ca 20170924

root = "~/projects/IMPC"
result_dir = "result"
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN") #"Sanger_SPLEEN","Sanger_MLN","CIPHE","TCP",

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrix", sep="")
dist_dir = paste(result_dir, "/", panelL, "/", centreL, "/dist", sep="")

#Output
peer_dir = paste(result_dir, "/", panelL, "/", centreL, "/PEER", sep=""); for(i in 1:length(peer_dir)) { suppressWarnings(dir.create(peer_dir[i])) }
#peer_dir = paste(result_dir, "/", panelL, "/", centreL, "/PEER_withoutDate", sep=""); for(i in 1:length(peer_dir)) { suppressWarnings(dir.create(peer_dir[i])) }
resid_dir = paste(peer_dir, "/peer_residual_analysis.txt", sep="")


#Libraries/Functions
libr(stringr)
libr(colorspace)
libr(lubridate) #if there are date variables
libr(peer)
libr(qtl)
libr(arules)
libr(foreach)
libr(doMC)
source("~/projects/IMPC/code/helpers.R")
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")









#Options for script
matrix_type = c("CountAdj") #, "Count", "Prop", "CountAdj_qn_all_leafonly500", "CountAdj_qn_all_nonleaf500", "CountAdj_all_leafonly500", "CountAdj_all_nonleaf500"
sampleCountThres = 3 #only compare if >=3 samples available

layers = c(1,3,5)
factors = c(2,5,10,15,20,25)

target_col = "gene" #column with control/experiment
controlL = c("+_+|+_Y","+_+|+_Y","WildType","WildType","WildType") #control value in target_col column
date_col = "date"

interested = c("gene", "date", "birth_date", "gender", "colony", "strain", "fur") #sampleMeta columns to plot against generated covriates
interestedCovs = list(NULL,c("gene", "date", "gender"),c("gene","gender"),c("gene", "date", "birth_date", "gender", "colony", "strain", "fur")) #sampleMeta columns to input into model as known covariates
interestedCont = c("date","birth_date") #continuous variables
wtonly = "" #"_WTonly" if only analyzing wildtypes, else ""
doNoCov = T # if you want to do the model without known covariates




#Setup Cores
no_cores = min(length(layers),detectCores()-1)
registerDoMC(no_cores)










start = Sys.time()

#for (ci in 4:1) {
for (ci in 1:length(paste0(panelL,centreL))) {
  start2 = Sys.time()
  centre = paste0(panelL,centreL)[ci]
  cat("\n",paste0(panelL," ",centreL)[ci])
  
  
  #Prepare data
  sampleMeta = get(load(sampleMeta_dir[ci]))
  #interestedCols0 = which(colnames(sampleMeta)%in%interested)
  
  # #order samples by date, exclude genotypes with less than 3 samples
  # if (length(grep("barcode", colnames(sampleMeta)))>0) sampleMeta2 <- sampleMeta[order(sampleMeta$barcode),] #sanger centre
  # if (length(grep("sample", colnames(sampleMeta)))>0) sampleMeta2 <- sampleMeta[order(sampleMeta$sample),] #other centres
  # if (length(grep("specimen", colnames(sampleMeta)))>0) sampleMeta2 <- sampleMeta[order(sampleMeta$specimen),] #other centres
  # sampleMeta2 <- sampleMeta2[order(sampleMeta2[,date_col]),]
  
  for (mcp in matrix_type) { cat("\n  ", mcp, ": ")
    start2 = Sys.time()
    
    # sm0 = sapply(sampleMeta2, function(x) as.numeric(factor(x, ordered=T)))
    # uniquecols = apply(sm0, 2, function(x) nrow(sm0)==length(unique(x)))
    # if (grepl("Sanger",centre)) { sm = sm0[,!uniquecols]
    # } else { sm = sm0[,!uniquecols] }
    
    ## get interested columns
    interestedCols = which(colnames(sampleMeta)%in%interested)
    
    mm = get(load(paste0(matrix_dir[ci], mcp,".Rdata")))
    
    ## get only 1st layer phenotypes
    phenoLevel = str_count(colnames(mm), "[+-]")
    
    ## peer analysis ------------------------------------------------
    a = foreach (layer = layers) %dopar% {
      Y = mm[,which(phenoLevel<=layer)]
      sampleMeta2 = sampleMeta[match(sampleMeta$fileName,rownames(mm)),]
      
      
      model = PEER()
      PEER_setPhenoMean(model, as.matrix(Y)) # data for inference - note the as.matrix() !
      # set priors (these are the default settings of PEER)
      PEER_setPriorAlpha(model,0.1,0.1);
      PEER_setPriorEps(model,0.1,10.);
      PEER_setNmax_iterations(model,1000)
      PEER_setAdd_mean(model, T)
      
      colour = rainbow(length(factors))
      
      ## no covariates included (sampleMeta not used)
      for (interestedCov in interestedCovs) {
        
        noCov = T
        columnssm = match(interested,colnames(sampleMeta2))
        columnssm = columnssm[!is.na(columnssm)]
        covs0 = sampleMeta2[,columnssm]
        
        #create genotype matrix
        sampleMeta2[,columnssm]
        
        if (!is.null(interestedCov)) {
          covs1 = sampleMeta2[,columnssm]
          if (wtonly!="") {
            wtind = grepl(controlL[ci],covs1[,target_col])
            Y = Y[wtind,]
            covs1 = covs1[wtind,]
            covs1 = covs1[,!colnames(covs1)%in%target_col]
          }
          covs = sapply(covs1, function(x) as.numeric(factor(x, ordered=T)))
          
          noCov = F
        }
        
        
        
        
        fname = paste0(peer_dir[ci], "/peer_", mcp,wtonly, "_layer-", str_pad(layer, 2, pad = "0"),ifelse(noCov,),"",paste("_knownCov-",paste(interested,collapse=".")))
        png(filename=paste0(fname, ".png"), width=700, height=700)
        
        for (factori in 1:length(factors)) {
          kfactor = factors[factori]
          
          # set data and parameters
          PEER_setNk(model, kfactor) #number of factor for learning
          if (!noCov) PEER_setCovariates(model, as.matrix(covs)) # covariates (e.g batch, RNA quality, ...) - not the as.matrix()!
          
          
          # perform inference
          PEER_update(model)
          
          peer = list()
          #investigate results
          #factors:
          peer$factors = X = PEER_getX(model)
          #weights:
          peer$weights = W = PEER_getW(model)
          #ARD parameters
          peer$ARD = Alpha = PEER_getAlpha(model)
          #get corrected dataset:
          peer$corrected = Yc = PEER_getResiduals(model)
          
          # plot variance of factors - in this case, we expect a natural elbow where there are 5 active factors, as 5 were simulated
          
          tryCatch({
            lines(1.0 / Alpha, type="l", col=colour[factori])
          }, error = function(err) {
            plot(1.0 / Alpha,xlab="Factors", ylab="Factor relevance", main="PEER inference; no covariance", type="l")#,xlim=nrow(Alpha),ylim=c(0,max(Alpha[is.finite(Alpha)])))       
          })
          
        }
        graphics.off()
        
        Alpha[!is.finite(Alpha)] = 0
        AlphaPlot = which(Alpha!=0)
        if (length(AlphaPlot)>0) {
          
          save(peer,file=paste0(fname,".Rdata"))
          
          png(filename=paste0(fname, "_compareWithKnown.png"), width=length(ncol(covs))*900, height=(length(factors))*400)
          par(mfrow=c(length(AlphaPlot),length(ncol(covs))),mar=c(.5,.5,.5,.5))
          
          colourH = heat.colors(max(sapply(which(colnames(covs)%in%interestedCont), function(x)length(unique(covs[,x])))))
          for (ai in AlphaPlot) {
            for (col in ncol(covs)) {
              cont = colnames(covs)[col]%in%interestedCont
              
              x_value = coloursm = covs[,col]
              coloursm[is.na(coloursm)] = "gray" #grey for non-finite values
              x_split = as.integer(factor(x_value))
              
              cor = cor(x_split, X[,ai])
              
              if (cont) {
                coloursm = colourH[coloursm]
                x_splitnames = sort(unique(covs0[,col]))
                if (grepl("date",colnames(covs0)[col])) {
                  x_value = as_date(covs0[,col])
                  x_splitnames = sort(unique(x_value))
                }
                plot(x_value, X[,ai], col = coloursm, 
                     main = paste0("Factor ",ai," and ", colnames(covs)[col]," Pearson Cor = ", cor), 
                     #ylim=c(min(X[is.finite(X[,ai]),ai]),max(X[X[,ai],ai]))+1, 
                     xlab = colnames(covs)[col], ylab = paste0("factor ",ai))
              } else {
                xy_boxplot = lapply(sort(unique(x_split)), function(x) X[x_split==x,ai])
                boxplot(xy_boxplot, lwd = 1, outline=F, #ylim=c(min(X[is.finite(X[,ai]),ai]),max(X[X[,ai],ai])),
                        main = paste0("Factor ",ai," and ", colnames(covs)[col]," Pearson Cor (ok if 2 var values) = ", cor),
                        xaxt = 'n', xlab = colnames(covs)[col], ylab = paste0("factor",ai)) #,xaxt ='n', xlab=testcol
                axis(1, at=1:length(x_splitnames), labels=x_splitnames)
                points(jitter(x_split, factor=1), X[,ai], col = coloursm)
              }
              
            }
            
          }
          graphics.off()
        }
        
      }
      
      
      ## covariates used
      fname = paste0(peer_dir[ci], "/peer_Cov_", mcp, wtonly, "_layer-", str_pad(layer, 2, pad = "0"),"_knownCov-",paste(interested,collapse="."))
      png(filename=paste0(fname, ".png"), width=700, height=700)
      
      for (factori in 1:length(factors)) {
        kfactor = factors[factori]
        
        # set data and parameters
        PEER_setNk(model, kfactor) #number of factor for learning
        
        PEER_update(model)
        
        #investigate results
        #factors:
        X = PEER_getX(model)
        #weights:
        W = PEER_getW(model)
        #ARD parameters
        Alpha = PEER_getAlpha(model)
        #get corrected dataset:
        Yc = PEER_getResiduals(model)
        
        # plot variance of factors - in this case, we expect a natural elbow where there are 5 active factors, as 5 were simulated
        
        
        # plot variance of factors - in this case, we expect a natural elbow where there are 5 active factors, as 5 were simulated
        tryCatch({
          lines(1.0 / Alpha, type="l", col=colour[factori])
        }, error = function(err) {
          plot(1.0 / Alpha,xlab="Factors", ylab="Factor relevance", main=paste0("PEER inference; with covariances:\n",paste(interested,collapse=".")), type="l")
        })
      }
      graphics.off()
      
      
      
      
      
      
      png(filename=paste0(fname, "_compareWithKnown.png"), width=length(ncol(covs))*400, height=(length(factors))*400)
      par(mfrow=c(length(AlphaPlot),length(ncol(covs))))
      
      colourH = heat.colors(max(sapply(which(colnames(covs)%in%interestedCont), function(x)length(unique(covs[,x])))))
      for (ai in AlphaPlot) {
        for (col in ncol(covs))
          cont = colnames(covs)[col]%in%interestedCont
        
        x_value = coloursm = covs[,col]
        x_split = as.integer(factor(x_value))
        
        cor = cor(x_split, X[,ai])
        
        if (cont) {
          coloursm = colourH[coloursm]
          x_splitnames = sort(unique(covs0[,col]))
          if (grepl("date",colnames(covs0)[col])) {
            x_value = as_date(covs0[,col])
            x_splitnames = sort(unique(x_value))
          }
          plot(x_value, X[,ai], col = coloursm, main = paste0("Factor ",ai," and ", colnames(covs)[col]," Pearson Cor = ", cor), xlab = colnames(covs)[col], ylab = paste0("factor ",ai), ylim=c(min(X[is.finite(X[,ai]),ai]),max(X[X[,ai],ai])))
        } else {
          xy_boxplot = lapply(sort(unique(x_split)), function(x) X[x_split==x,ai])
          boxplot(xy_boxplot, lwd = 1, outline=F, ylim=c(min(X[is.finite(X[,ai]),ai]),max(X[X[,ai],ai])),
                  main = paste0("Factor ",ai," and ", colnames(covs)[col]," Pearson Cor (ok if 2 var values) = ", cor),
                  xaxt = 'n', xlab = colnames(covs)[col], ylab = paste0("factor",ai)) #,xaxt ='n', xlab=testcol
          axis(1, at=1:length(x_splitnames), labels=x_splitnames)
          points(jitter(x_split, factor=1), X[,ai], col = coloursm)
        }
        
      }
      graphics.off()
    }
    
    
    
    #try a whole bunch of residual parameters
    sink(file=resid_dir,append=T)
    print("\n\n\n======================================================\nLayer",layer,": plot factor weight variances for a large set of maximum number K of factors, and see if K has an effect\n")
    # plot factor weight variances for a large set of maximum number K of factors, and see if K has an effect
    for(pa in c(0.0001, 0.1, 1000)){
      for(pb in c(0.1,10,1000)){
        model = get_simple_PEER_object() # see simple_unsupervised_demo for how it is constructed
        PEER_setPriorEps(model,0.1, pb);
        sink()
        PEER_update(model)
        sink(file=resid_dir,append=T)
        print(paste("Eps pa=", pa, "pb=", pb, "mean(residuals^2)=",mean((PEER_getResiduals(model))**2)))
      }
    }
    sink()
    
    
    
    
  }
  
  cat("\n centre ", centre, " ",TimeOutput(start2)," \n",sep="")
}

TimeOutput(start)





# ## See dates without WT
# #Harwell
# control <- "FR-FCM-ZYCB-WildType_01_All"
# KOdays <- table(sampleMeta[!sampleMeta$gene%in%control,]$date)
# colnames(sampleCountThres) <- table(sampleMeta[sampleMeta$gene%in%control,]$date)
# noWTdays <- KOdays [!(names(KOdays) %in% names(WTdays))]
# 
# 
# 
# noWTdaysInd <- which(as.character(sampleMeta2$date)%in%names(noWTdays))
# plot()








## AMI evaluation sci.kit.learn ==FAST on python, not good on R NMI

## libr (fpc)

#set.seed(4634)
#face <- rFace(300,dMoNo=2,dNoEy=0)
#grface <- as.integer(attr(face,"grouping"))
#plotcluster(face,grface==1)
#plotcluster(face,grface, clnum=1, method="vbc")








