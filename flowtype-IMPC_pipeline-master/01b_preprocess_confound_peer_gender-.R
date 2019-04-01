# PEER coundounding factor analysis (gender comparison)
# aya43@sfu.ca 20170924

## root directory
root = "~/projects/IMPC"
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#,"Sanger_MLN","CIPHE","TCP","H")
controlL = c("+_+|+_Y","+_+|+_Y","WildType","WildType","WildType") #control value in target_col column
ci = 1; panel = panelL; centre = centreL

result_dir = paste0("result/", panelL, "/", centreL); suppressWarnings(dir.create (result_dir))

## input directories
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste(meta_dir, "/file", sep="")
feat_dir = paste(result_dir, "/feat", sep="")

## output directories
scale = T
peer_dir = paste(result_dir, "/PEER", sep="")
if (scale) peer_dir = paste0(peer_dir, "_scale", sep="")
dir.create(peer_dir, showWarnings=F)
resid_dir = paste(peer_dir, "/peer_residual_analysis.txt", sep="")


#Libraries/Functions
source("~/projects/IMPC/code/helpers.R")
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
libr("stringr")
libr("colorspace")
libr("lubridate") #if there are date variables
libr("peer")
libr("qtl")
libr("DMwR")
libr("arules")
libr("foreach")
libr("doMC")









#Options for script
feat_types = c("CountAdj")#,"Prop") #, "Count", "CountAdj_qn_all_leafonly500", "CountAdj_qn_all_nonleaf500", "CountAdj_all_leafonly500", "CountAdj_all_nonleaf500"
# sampleCountThres = 3 #only compare if >=3 samples available

layers = c(1:6)
factors = c(15,10,5)


target_col = "gene" #column with control/experiment
controlL = c("[+]_[+]|[+]_Y","[+]_[+]|[+]_Y","WildType","WildType","WildType") #control value in target_col column
date_col = "date"

genderCompare = c(0,1) #female only=0, male only=1, in this order, compare at male
interested = c("gene", "date", "birth_date", "gender", "colony", "strain", "fur") #mta_file columns to plot against generated covriates
interestedCovs = list(NULL,c("gene"),c("gender"),c("gene","gender"),c("gene", "date", "gender"),c("gene", "date", "birth_date", "gender", "colony", "strain", "fur")) #mta_file columns to input into model as known covariates
# interestedCovs = list(NULL,c("gene","gender")) #mta_file columns to input into model as known covariates
interestedCont = c("date","birth_date") #continuous variables
wtonlyL = c("","_WTonly") #"_WTonly" if only analyzing wildtypes, else ""; code modified, just leave this! do both!
layerequal = c(T,F) # only include cell populations in the layer (T) or it and all above (F)


replaceModel = T #if true, run the whole thing, else use existing predictions
dateRestriction = c("_dateRestrict","2015-08-25","2016-03-03") #if "", don't consider, if "_dateRestrict", cut samples into date interval
niter = 10000

scale_center = T
scale_scale = T

plotCountThres = 500 #how high must mean cell population count be ot get plotted for before and after plot
plotCountNo = 4 #how many of those cell population counts to plot




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
  mta_file = get(load(mta_file_dir[ci]))
  mta_file$gene[grepl("[+]_Y",mta_file$gene)] = "+_+"
  
  #interestedCols0 = which(colnames(mta_file)%in%interested)
  
  # sm0 = sapply(mta_file, function(x) as.numeric(factor(x, ordered=T))) #convert all values to integers
  # uniquecols = apply(sm0, 2, function(x) nrow(sm0)==length(unique(x))) #get rid of unique columns such as file name
  # sm = sm0[,!uniquecols]
  
  # #order samples by date, exclude genotypes with less than 3 samples
  # if (length(grep("barcode", colnames(mta_file)))>0) mta_file2 <- mta_file[order(mta_file$barcode),] #sanger centre
  # if (length(grep("sample", colnames(mta_file)))>0) mta_file2 <- mta_file[order(mta_file$sample),] #other centres
  # if (length(grep("specimen", colnames(mta_file)))>0) mta_file2 <- mta_file[order(mta_file$specimen),] #other centres
  # mta_file2 <- mta_file2[order(mta_file2[,date_col]),]
  
  for (mcp in matrix_type) { cat("\n  ", mcp, ": ") #for each type of feature
    start2 = Sys.time()
    
    # ## get interested columns
    # interestedCols = which(colnames(sm)%in%interested)
    
    mm = get(load(paste0(matrix_dir[ci], mcp,".Rdata")))
    mta_file2 = mta_file20 = mta_file[match(mta_file$fileName,rownames(mm)),]
    
    if (dateRestriction[1]!="") {
      dateind = ymd(mta_file2$date)>ymd(dateRestriction[2]) & ymd(mta_file2$date)<ymd(dateRestriction[3]) & mm[,5]<250000 & mm[,5]>75000
      mta_file2 = mta_file20 = mta_file2[dateind,]
      mm = mm[dateind,]
    }
    
    for (wtonly in wtonlyL) {
      #if wildtype only, clip mta_file and feature matrix
      mta_file200 = mta_file20
      if (wtonly!="") {
        wtind = grepl(controlL[ci],mta_file2[,target_col])
        mm = mm[wtind,]
        mta_file2 = mta_file2[wtind,!colnames(mta_file2)%in%target_col]
        mta_file200 = mta_file20[wtind,]
      }
      
      
      
      ## get only 1st layer phenotypes
      phenoLevel = str_count(colnames(mm), "[+-]")
      
      ## peer analysis ------------------------------------------------
      a = foreach (layer = layers) %dopar% {
      #for (layer in layers) { cat("layer ",layer,", ")
        for (le in layerequal) { cat("... ")
          if (le==F & layer==1) next()
          
          
          # trim feature matrix
          if (le) { Y = mm[,which(phenoLevel==layer)]; lei = "_only" }
          if (!le) { Y = mm[,which(phenoLevel<=layer)]; lei = "" }
          
          # keep only needed columns in mta_file3 (interested)
          columnssm0 = match(interested,colnames(mta_file2))
          columnssm0 = columnssm0[!is.na(columnssm0)]
          mta_file3 = mta_file2[,columnssm0]
          
          colour = rainbow(length(factors)) #compare predictions given different number of factors to find
          
          ## no covariates included (mta_file not used)
          for (interestedCov in interestedCovs) { cat(", cov")
            
            
            #create known covariate matrix
            noCov = T
            covs = NULL
            if (!is.null(interestedCov)) {
              columnssm = match(interestedCov,colnames(mta_file3))
              columnssm = columnssm[!is.na(columnssm)]
              if (length(columnssm)>0) {
                covs1 = mta_file3[,columnssm]
                if (is.null(dim(covs1))) { covs1 = matrix(covs1,ncol=1) }
                covs = sapply(1:ncol(covs1), function(x) as.numeric(factor(covs1[,x], ordered=T)))
                if (is.null(dim(covs))) { covs = matrix(covs,ncol=1) }
                colnames(covs) = colnames(mta_file3)[columnssm] 
                noCov = F
              }
            }
            
            for (gender in genderCompare) {
              genderind = mta_file200$gender==gender
              Ytemp = Y
              Y = Y[genderind,]
              
              
              fnamepartt = ""
              if (!noCov) fnamepartt = paste0("_knownCov-",paste(interestedCov,collapse="."))
              fname0 = paste0(peer_dir[ci], "/peer_", mcp,wtonly,lei,"_layer-", str_pad(layer, 2, pad = "0"),"_n.",ncol(Y))
              if (dateRestriction[1]!="") fname0 = paste0(peer_dir[ci], "/peer_", mcp,wtonly,dateRestriction[1], lei,"_layer-", str_pad(layer, 2, pad = "0"),"_n.",ncol(Y))
              fname = paste0(fname0,"_gender-",gender,fnamepartt)
              
              
              #prepare list of models ((list of the predicted factors from PEER))
              modelexists = T
              if (!replaceModel & file.exists(paste0(fname,".Rdata"))) {
                models = get(load(paste0(fname,".Rdata")))
                if (length(models)==0) next()
                modelsfactors = as.numeric(names(models))
                #if models was already made, just replot it, likely plotting's changed that's why i'm coming back; and it doesn't take long
              } else {
                models = list()
                modelsfactors = factors
                modelexists = F
              }
              if (gender==genderCompare[2]) {
                modelsF = get(load(paste0(gsub(paste0("_gender-",gender),paste0("_gender-",genderCompare[1]),fname),".Rdata")))
                modelsF = Filter(Negate(is.null), modelsF)
                modelsfactors = as.numeric(intersect(modelsfactors,names(modelsF)))
              }
              modelsfactors = modelsfactors[!is.na(modelsfactors)]
              
              
              
              for (factori in 1:length(modelsfactors)) {
                kfactor = modelsfactors[factori]
                
                if (gender==genderCompare[2]) {
                  AlphaF = modelsF[[factori]]$Alpha
                  AlphaF = as.vector(AlphaF)
                  
                  #factors:
                  XF = modelsF[[factori]]$X
                  #weights:
                  WF = modelsF[[factori]]$W
                }
                if (modelexists) {
                  Alpha = models[[factori]]$Alpha
                  AlphaPlot = which(Alpha!=0)
                  #factors:
                  X = models[[factori]]$X
                  #weights:
                  W = models[[factori]]$W
                } else {
                  # build model
                  # set data and parameters
                  model = PEER()
                  if (scaleL) {
                    PEER_setPhenoMean(model, scale(as.matrix(Y),center=scale_center,scale=scale_scale)) # data for inference - note the as.matrix() !
                  } else {
                    PEER_setPhenoMean(model, as.matrix(Y)) # data for inference - note the as.matrix() !
                  }
                  # set priors (these are the default settings of PEER)
                  PEER_setPriorAlpha(model,0.1,0.1)
                  PEER_setPriorEps(model,0.1,10)
                  PEER_setNmax_iterations(model,niter)
                  PEER_setNk(model, kfactor) #number of factor for learning
                  PEER_setAdd_mean(model, TRUE) #doesn't work otherwise; adds one factor that is all 1 with weights as mean count?
                  if (!noCov) PEER_setCovariates(model, as.matrix(covs)) # covariates (e.g batch, RNA quality, ...) - not the as.matrix()!
                  
                  # perform inference
                  PEER_update(model)
                  
                  #investigate results
                  #ARD parameters
                  Alpha = PEER_getAlpha(model)
                  Alpha[!is.finite(Alpha)] = 0
                  #factors:
                  X = PEER_getX(model)
                  #delete factors and fix
                  AlphaPlot = which(Alpha!=0 & sapply(1:ncol(X), function(x) length(unique(X[,x]))>1))
                  if (!length(AlphaPlot)>0) next() #skip if no valid factors found
                  AlphaPlot = AlphaPlot[order(Alpha[AlphaPlot,],decreasing=F)]
                  Alpha = Alpha[AlphaPlot,]
                  if (is.null(dim(Alpha))) Alpha = matrix(Alpha,ncol=1)
                  X = X[,AlphaPlot]
                  if (is.null(dim(X))) X = matrix(X,ncol=1)
                  #weights:
                  W = PEER_getW(model)[,AlphaPlot]
                  if (is.null(dim(W))) W = matrix(W,ncol=1)
                  #get corrected dataset:
                  Yc = PEER_getResiduals(model)
                  dimnames(Yc) = dimnames(Y)
                  if (scaleL) {
                    YcUnscaled = unscale(Yc,scale(as.matrix(Y)))
                    dimnames(YcUnscaled) = dimnames(as.matrix(Y))
                  } else {
                    YcUnscaled = NULL
                  }
                  
                  models[[as.character(kfactor)]] = list(X=X, Alpha=Alpha, W=W)
                }
                if (gender==genderCompare[1]) next()
                
                
                
                #compare PEER on Male data vs Female only
                
                png(filename=paste0(fname, "_factor-",kfactor, "_compareWTandALL.png"), width=length(AlphaPlot)*500, height=ncol(WF)*300)
                layout(matrix(c(1,rep(2,length(AlphaPlot)-1), 3:(ncol(WF)*length(AlphaPlot)+2)),ncol=length(AlphaPlot),byrow=T))
                
                par(mar=c(3,3,3,3)) #mfrow=c(length(AlphaPlot),ncol(WF))
                
                # WF.1 = sapply(WF, function(x) as.numeric(factor(x, ordered=T))) #numeric
                
                
                
                if (max(1/Alpha)<max(1/as.numeric(AlphaF))) {
                  ## plot inverse variance of factors - in this case, we expect a natural elbow where there are 5 active factors, as 5 were simulated
                  plot(1.0 /  as.numeric(AlphaF),xlab="Factors", ylab="Factor relevance (factor inverse variance over cell pops)", main=paste0("PEER inference; dashed=",gender), type="l",col=colour[factori])#,xlim=nrow(Alpha),ylim=c(0,max(Alpha[is.finite(Alpha)])))       
                  lines(1.0 / Alpha,col="black",lty=2)
                  legend("topright",legend=factors,text.col=colour,title="no. factors to predict")
                  
                  plot(1.0 /  as.numeric(AlphaF),xlab="Factors", ylab="Factor relevance (factor inverse variance over cell pops) (dashed line = aLL files)", main=paste0("PEER inference; dashed=",gender), type="l",col=colour[factori])#,xlim=nrow(Alpha),ylim=c(0,max(Alpha[is.finite(Alpha)]))) 
                  lines(1.0 /Alpha,col="black",lty=2)
                  legend("topright",legend=factors,text.col=colour,title="no. factors to predict")
                } else {
                  ## plot inverse variance of factors - in this case, we expect a natural elbow where there are 5 active factors, as 5 were simulated
                  plot(1.0 / Alpha,xlab="Factors", ylab="Factor relevance (factor inverse variance over cell pops)", main=paste0("PEER inference; dashed=",gender), type="l",col=colour[factori])#,xlim=nrow(Alpha),ylim=c(0,max(Alpha[is.finite(Alpha)])))       
                  lines(1.0 / as.numeric(AlphaF),col="black",lty=2)
                  legend("topright",legend=factors,text.col=colour,title="no. factors to predict")
                  
                  plot(1.0 / Alpha,xlab="Factors", ylab="Factor relevance (factor inverse variance over cell pops) (dashed line = aLL files)", main=paste0("PEER inference; dashed=",gender), type="l",col=colour[factori])#,xlim=nrow(Alpha),ylim=c(0,max(Alpha[is.finite(Alpha)]))) 
                  lines(1.0 / as.numeric(AlphaF),col="black",lty=2)
                  legend("topright",legend=factors,text.col=colour,title="no. factors to predict")
                }
                
                
                
                ## plot factors against variables if they do
                
                #heat colours for continuous values
                colourH = heat.colors(max(sapply(1:ncol(WF), function(x)length(unique(WF[,x])))))
                
                for (col in 1:ncol(WF)) { #for each interested known factor
                  for (ai in 1:length(AlphaPlot)) { #for each known/learned factor
                    
                    cont = AlphaF[col]%in%interestedCont #is interested factor continuous
                    
                    # interested factor values to plot
                    w_value = WF[,col]
                    coloursm = w_split = as.integer(factor(w_value))
                    coloursm[is.na(coloursm)] = "gray" #grey for non-finite values
                    w_splitnames = sort(unique(WF[,col]))
                    
                    if (length(unique(W[,ai]))==1) { #should only happen in ploti=2 when we compare gene factor
                      cort = "NA"
                    } else {
                      cort = specify_decimal(cor.test(w_split, W[,ai])$p.value,4)
                    }
                    cor = specify_decimal(cor(w_split, W[,ai]),4)
                    
                    if (cont) {
                      
                      coloursm = colourH[as.integer(coloursm)]
                      if (grepl("date",colnames(WF)[col])) {
                        w_value = as_date(WF[,col])
                        w_splitnames = sort(unique(w_value))
                      }
                      mainn = paste0("WEIGHTS: Factor ",ai," and gender=",gender,"'s factor ", col," Pearson Cor = ", cor,"; p value=",cort)
                      plot(w_value, W[,ai], col = coloursm,
                           main = mainn,
                           #ylim=c(min(W[is.finite(W[,ai]),ai]),max(W[W[,ai],ai]))+1,
                           xlab = colnames(WF)[col], ylab = paste0("factor ",ai))
                    } else {
                      if (length(unique(w_split))>3 & length(unique(W[,ai]))>1) {
                        if (length(unique(W[,ai]))>1.5*length(unique(w_split))) {
                          Xai = as.numeric(factor(discretize(W[,ai],categories=2*length(unique(w_split)))))
                        } else {
                          Xai = W[,ai]
                        }
                        cort = specify_decimal(chisq.test(w_split, Xai)$p.value,4)
                        cort = paste0(cort," (Chi2)")
                      } 
                      mainn = paste0("WEIGHTS: Factor ",ai," and gender=",gender,"'s factor ", col," Pearson Cor (ok if 2 var values) = ", cor,"; p value=",cort)
                      
                      xy_boxplot = lapply(sort(unique(w_split)), function(x) W[w_split==x,ai])
                      boxplot(xy_boxplot, lwd = 1, outline=F, #ylim=c(min(W[is.finite(W[,ai]),ai]),max(W[W[,ai],ai])),
                              main = mainn,
                              xaxt = 'n', xlab = colnames(WF)[col], ylab = paste0("factor",ai), ylab = "Covariate value") #,xaxt ='n', xlab=testcol
                      axis(1, at=1:length(w_splitnames), labels=w_splitnames)
                      points(jitter(w_split, factor=1), W[,ai], col = coloursm, cex=.5, pch=16)
                    }
                    
                  }
                }
                
                
                

                graphics.off()
                
                
                
                
                
                
              }
              Y = Ytemp
              if (!modelexists & length(models)>0) {
                save(models,file=paste0(fname,".Rdata"))
              }
              
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
  
}

cat("\n centre ", centre, " ",TimeOutput(start2)," \n",sep="")
}

TimeOutput(start)





# ## See dates without WT
# #Harwell
# control <- "FR-FCM-ZYCB-WildType_01_All"
# KOdays <- table(mta_file[!mta_file$gene%in%control,]$date)
# colnames(sampleCountThres) <- table(mta_file[mta_file$gene%in%control,]$date)
# noWTdays <- KOdays [!(names(KOdays) %in% names(WTdays))]
# 
# 
# 
# noWTdaysInd <- which(as.character(mta_file2$date)%in%names(noWTdays))
# plot()








## AMI evaluation sci.kit.learn ==FAST on python, not good on R NMI

## libr (fpc)

#set.seed(4634)
#face <- rFace(300,dMoNo=2,dNoEy=0)
#grface <- as.integer(attr(face,"grouping"))
#plotcluster(face,grface==1)
#plotcluster(face,grface, clnum=1, method="vbc")








