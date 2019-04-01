# HCP coundounding factor analysis
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
hcp_dir = paste(result_dir, "/", panelL, "/", centreL, "/hcp", sep=""); for(i in 1:length(hcp_dir)) { suppressWarnings(dir.create(hcp_dir[i])) }


#Libraries/Functions
libr(stringr)
libr(colorspace)
libr(lubridate) #if there are date variables
libr(Rhcpp)
libr(robustHD)
libr(foreach)
libr(doMC)
source("~/projects/IMPC/code/helpers.R")
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")









#Options for script
matrix_type = c("CountAdj","Prop") #, "Count", "Prop", "CountAdj_qn_all_leafonly500", "CountAdj_qn_all_nonleaf500", "CountAdj_all_leafonly500", "CountAdj_all_nonleaf500"
sampleCountThres = 3 #only compare if >=3 samples available

layers = c(1,3,5)
factors = c(15,5,10)


target_col = "gene" #column with control/experiment
controlL = c("+_+|+_Y","+_+|+_Y","WildType","WildType","WildType") #control value in target_col column
date_col = "date"

interested = c("gene", "date", "birth_date", "gender", "colony", "strain", "fur") #sampleMeta columns to plot against generated covriates
interestedCovs = list(c("gene","gender"),c("gene", "date", "gender"),c("gene", "date"),c("gene", "date", "birth_date", "gender", "colony", "strain", "fur")) #sampleMeta columns to input into model as known covariates
# interestedCovs = list(NULL,c("gene","gender")) #sampleMeta columns to input into model as known covariates
interestedCont = c("date","birth_date") #continuous variables
wtonly = "_WTonly" #"_WTonly" if only analyzing wildtypes, else ""
layerequal = c(T,F) # only include cell populations in the layer (T) or it and all above (F)

replaceModel = F #if true, run the whole thing, else use existing predictions
dateRestriction = c("_dateRestrict","2015-08-25","2016-03-03") #if "", don't consider, if "_dateRestrict", cut samples into date interval
niter = 10000




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
  
  # sm0 = sapply(sampleMeta, function(x) as.numeric(factor(x, ordered=T))) #convert all values to integers
  # uniquecols = apply(sm0, 2, function(x) nrow(sm0)==length(unique(x))) #get rid of unique columns such as file name
  # sm = sm0[,!uniquecols]
  
  # #order samples by date, exclude genotypes with less than 3 samples
  # if (length(grep("barcode", colnames(sampleMeta)))>0) sampleMeta2 <- sampleMeta[order(sampleMeta$barcode),] #sanger centre
  # if (length(grep("sample", colnames(sampleMeta)))>0) sampleMeta2 <- sampleMeta[order(sampleMeta$sample),] #other centres
  # if (length(grep("specimen", colnames(sampleMeta)))>0) sampleMeta2 <- sampleMeta[order(sampleMeta$specimen),] #other centres
  # sampleMeta2 <- sampleMeta2[order(sampleMeta2[,date_col]),]
  
  for (mcp in matrix_type) { cat("\n  ", mcp, ": ") #for each type of feature
    start2 = Sys.time()
    
    # ## get interested columns
    # interestedCols = which(colnames(sm)%in%interested)
    
    mm = get(load(paste0(matrix_dir[ci], mcp,".Rdata")))
    sampleMeta2 = sampleMeta[match(sampleMeta$fileName,rownames(mm)),]
    
    #if wildtype only, clip samplemeta and feature matrix
    if (wtonly!="") {
      wtind = grepl(controlL[ci],sampleMeta2[,target_col])
      mm = mm[wtind,]
      sampleMeta2 = sampleMeta2[wtind,!colnames(sampleMeta2)%in%target_col]
    }
    if (dateRestriction[1]!="") {
      dateind = ymd(sampleMeta2$date)>ymd(dateRestriction[2]) & ymd(sampleMeta2$date)<ymd(dateRestriction[3]) & mm[,5]<250000 & mm[,5]>75000
      sampleMeta2 = sampleMeta2[dateind,]
      mm = mm[dateind,]
    }
    
    
    ## get only 1st layer phenotypes
    phenoLevel = str_count(colnames(mm), "[+-]")
    
    ## hcp analysis ------------------------------------------------
    # a = foreach (layer = layers) %dopar% {
    for (layer in layers) { cat("layer ",layer,", ")
      for (le in layerequal) { cat("... ")
        
        # trim feature matrix
        if (le) { Y = mm[,which(phenoLevel==layer)]; lei = "_only" }
        if (!le) { Y = mm[,which(phenoLevel<=layer)]; lei = "" }
        
        # keep only needed columns in sampleMeta3 (interested)
        columnssm = match(interested,colnames(sampleMeta2))
        columnssm = columnssm[!is.na(columnssm)]
        sampleMeta3 = sampleMeta2[,columnssm]
        
        colour = rainbow(length(factors)) #compare predictions given different number of factors to find
        
        ## no covariates included (sampleMeta not used)
        for (interestedCov in interestedCovs) { cat(", cov")
          
          #create known covariate matrix
          noCov = T
          covs = NULL
          if (!is.null(interestedCov)) {
            columnssm = match(interestedCov,colnames(sampleMeta3))
            columnssm = columnssm[!is.na(columnssm)]
            covs1 = sampleMeta3[,columnssm]
            if (is.null(dim(covs1))) { covs1 = matrix(covs1,ncol=1) }
            covs = sapply(1:ncol(covs1), function(x) as.numeric(factor(covs1[,x], ordered=T)))
            if (is.null(dim(covs))) { covs = matrix(covs,ncol=1) }
            colnames(covs) = colnames(sampleMeta3)[columnssm] 
            noCov = F
          }
          
          
          
          fnamepartt = ""
          if (!noCov) fnamepartt = paste0("_knownCov-",paste(colnames(covs),collapse="."))
          fname = paste0(hcp_dir[ci], "/hcp_", mcp,wtonly,lei,"_layer-", str_pad(layer, 2, pad = "0"),"_n.",ncol(Y),fnamepartt)
          if (dateRestriction[1]!="") fname = paste0(hcp_dir[ci], "/hcp_", mcp,wtonly,dateRestriction[1], lei,"_layer-", str_pad(layer, 2, pad = "0"),"_n.",ncol(Y),fnamepartt)
          
          
          for (factori in 1:length(factors)) {
            kfactor = factors[factori]
            
            if (F & !replaceModel & file.exists(paste0(fname,".Rdata"))) {
              Rres = get(load(paste0(fname,".Rdata")))
            } else {
              # build model
              #do HCP
              covs0 = standardize(covs)
              # covs0 = covs0-min(covs0)+min(abs(covs0))
              covs0 = sapply(data.frame(covs0), function(x) x-min(x)+min(abs(x)))
              covs0 = covs-min(covs)+.1
              covs0 = log(covs0)
              Y0 = standardize(Y)
              # Y0 = Y-min(Y)+.1
              Y0 = sapply(data.frame(Y0), function(x) x-min(x)+min(abs(x)))
              Y0 = log(Y0)
              Rres <- hcp(Z=covs0, Y=Y0, k=kfactor, lambda1=100, lambda2=1, lambda3=1, iter=niter, stand=F, fast=F, log=F)
              
              X = Rres$W
              AlphaPlot = sapply(1:nrow(Rres$B), function(x) !all(Rres$B[x,]==0)) #check which factors to plot; just do top 10 for now
              # Alpha = #add inverse variance weighted mean!
              AlphaPlot = which(AlphaPlot)
            }
            
            ## does learned factors matter
            if (length(AlphaPlot)>0) {
              save(Rres,file=paste0(fname,"_factor.",kfactor,".Rdata"))
              sampleMeta3.1 = sapply(sampleMeta3, function(x) as.numeric(factor(x, ordered=T))) #numeric
              
              ## plot factors against variables if they do
              png(filename=paste0(fname, "_factor.", kfactor,"_compareWithKnown.png"), width=ncol(sampleMeta3)*900, height=(length(AlphaPlot))*400)
              par(mfrow=c(length(AlphaPlot),ncol(sampleMeta3)),mar=c(.5,.5,.5,.5))
              
              #heat colours for continuous values
              colourH = heat.colors(max(sapply(which(colnames(sampleMeta3)%in%interestedCont), function(x)length(unique(sampleMeta3[,x])))))
              
              for (ai in AlphaPlot) { #for each known/learned factor
                for (col in 1:ncol(sampleMeta3)) { #for each interested known factor
                  cont = colnames(sampleMeta3)[col]%in%interestedCont #is interested factor continuous
                  
                  # interested factor values to plot
                  x_value = sampleMeta3[,col]
                  coloursm = x_split = as.integer(factor(x_value))
                  coloursm[is.na(coloursm)] = "gray" #grey for non-finite values
                  x_splitnames = sort(unique(sampleMeta3[,col]))
                  
                  cor = cor(x_split, X[,ai])
                  
                  try({
                    if (cont) {
                      coloursm = colourH[as.integer(coloursm)]
                      if (grepl("date",colnames(sampleMeta3)[col])) {
                        x_value = as_date(sampleMeta3[,col])
                        x_splitnames = sort(unique(x_value))
                      }
                      plot(x_value, X[,ai], col = coloursm,
                           main = paste0("Factor ",ai," and ", colnames(sampleMeta3)[col]," Pearson Cor = ", cor),
                           #ylim=c(min(X[is.finite(X[,ai]),ai]),max(X[X[,ai],ai]))+1,
                           xlab = colnames(sampleMeta3)[col], ylab = paste0("factor ",ai))
                    } else {
                      xy_boxplot = lapply(sort(unique(x_split)), function(x) X[x_split==x,ai])
                      boxplot(xy_boxplot, lwd = 1, outline=F, #ylim=c(min(X[is.finite(X[,ai]),ai]),max(X[X[,ai],ai])),
                              main = paste0("Factor ",ai," and ", colnames(sampleMeta3)[col]," Pearson Cor (ok if 2 var values) = ", cor),
                              xaxt = 'n', xlab = colnames(sampleMeta3)[col], ylab = paste0("factor",ai), ylab = "Covariate value") #,xaxt ='n', xlab=testcol
                      axis(1, at=1:length(x_splitnames), labels=x_splitnames)
                      points(jitter(x_split, factor=1), X[,ai], col = coloursm, cex=.5, pch=16)
                    }
                  })
                  
                }
                
              }
              graphics.off()
            }
            
          }
          
          
        }
        
        
        
      }
    }
    
    
    
    # #try a whole bunch of residual parameters
    # sink(file=resid_dir,append=T)
    # cat("\n\n\n======================================================\nLayer",layer,": plot factor weight variances for a large set of maximum number K of factors, and see if K has an effect\n")
    # # plot factor weight variances for a large set of maximum number K of factors, and see if K has an effect
    # for(pa in c(0.0001, 0.1, 1000)){
    #   for(pb in c(0.1,10,1000)){
    #     model = get_simple_PEER_object() # see simple_unsupervised_demo for how it is constructed
    #     PEER_setPriorEps(model,0.1, pb);
    #     PEER_update(model)
    #     cat(paste("\nEps pa=", pa, "pb=", pb, "mean(residuals^2)=",mean((PEER_getResiduals(model))**2)))
    #   }
    # }
    # sink()
    
    
    
    
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








