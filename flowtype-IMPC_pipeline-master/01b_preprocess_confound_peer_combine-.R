# PEER coundounding factor analysis; combine layers
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
scaleL = T
if (scaleL) {
  peer_dir = paste(result_dir, "/", panelL, "/", centreL, "/PEER_scale", sep=""); for(i in 1:length(peer_dir)) { suppressWarnings(dir.create(peer_dir[i])) }
} else {
  peer_dir = paste(result_dir, "/", panelL, "/", centreL, "/PEER", sep=""); for(i in 1:length(peer_dir)) { suppressWarnings(dir.create(peer_dir[i])) }
}
resid_dir = paste(peer_dir, "/peer_residual_analysis.txt", sep="")


#Libraries/Functions
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
libr("stringr")
libr("lubridate") #if there are date variables
libr("foreach")
libr("doMC")









#Options for script
matrix_type = c("CountAdj")#,"Prop") #, "Count", "Prop", "CountAdj_qn_all_leafonly500", "CountAdj_qn_all_nonleaf500", "CountAdj_all_leafonly500", "CountAdj_all_nonleaf500"
# sampleCountThres = 3 #only compare if >=3 samples available

factors = c(5) #c(15,5,10)


target_col = "gene" #column with control/experiment
controlL = c("[+]_[+]|[+]_Y","[+]_[+]|[+]_Y","WildType","WildType","WildType") #control value in target_col column
date_col = "date"

interested = c("gene", "date", "birth_date", "gender", "colony", "strain", "fur") #sampleMeta columns to plot against generated covriates
interestedCovs = list(NULL)#,c("gene"),c("gender"),c("gene","gender"),c("gene", "date", "gender"),c("gene", "date", "birth_date", "gender", "colony", "strain", "fur")) #sampleMeta columns to input into model as known covariates
# interestedCovs = list(NULL,c("gene","gender")) #sampleMeta columns to input into model as known covariates
interestedCont = c("date","birth_date") #continuous variables
wtonlyL = c("") #this script is only for c(""); if only analyzing wildtypes, else ""; code modified, just leave this! do both!
# layerequal = c(T,F) # T!!! only include cell populations in the layer (T) or it and all above (F)


# replaceModel = F #if true, run the whole thing, else use existing predictions
dateRestriction = c("","2015-08-25","2016-03-03") #if "", don't consider, if "_dateRestrict", cut samples into date interval

scale_center = T
scale_scale = T





# #Setup Cores
# no_cores = min(length(layers),detectCores()-1)
# registerDoMC(no_cores)










start = Sys.time()

#for (ci in 4:1) {
for (ci in 1:length(paste0(panelL,centreL))) {
  start2 = Sys.time()
  centre = paste0(panelL,centreL)[ci]
  cat("\n",paste0(panelL," ",centreL)[ci])
  
  for (mcp in matrix_type) { cat("\n  ", mcp, ": ") #for each type of feature
    start2 = Sys.time()
    
    
    mm = get(load(paste0(matrix_dir[ci], mcp,".Rdata")))
    phenoLevel = str_count(colnames(mm), "[+-]")
    
    wtonly=""
    
    
    # ## get only 1st layer phenotypes
    # phenoLevel = str_count(colnames(mm), "[+-]")
    # 
    # ## peer analysis ------------------------------------------------
    # a = foreach (layer = layers) %dopar% {
    le=T
    lei = "_only"
    
    ## no covariates included (sampleMeta not used)
    for (interestedCov in interestedCovs) { cat(", cov")
      
      layer="00" #placeholder
      
      fnamepartt = ""
      if (!is.null(interestedCov)) fnamepartt = paste0("_knownCov-",paste(interestedCov,collapse="."))
      fname0 = paste0(peer_dir[ci], "/peer_", mcp,wtonly,lei,"_layer-", str_pad(layer, 2, pad = "0"),"_n.","[0-9]+")
      if (dateRestriction[1]!="") fname0 = paste0(peer_dir[ci], "/peer_", mcp,wtonly,dateRestriction[1], lei,"_layer-", str_pad(layer, 2, pad = "0"),"_n.","[0-9]+")
      fname = paste0(fname0,fnamepartt)
      
      fpaths = paste(peer_dir[ci],"/",list.files(peer_dir[ci],pattern=".Rdata"),sep="")
      fpaths = fpaths[!grepl("_gender-",fpaths)]
      fpaths = fpaths[!grepl("Prop",fpaths)]
      fpaths
      layers = unique(as.numeric(gsub("layer","",str_extract(fpaths,"layer[0-9]+"))))
      
      models_combine = list()
      for (layer in sort(layers)) {
        fpath = fpaths[grepl(paste0(gsub("00",str_pad(layer, 2, pad = "0"),fname),".Rdata"), fpaths)]
        model = get(load(fpath))
        model[[as.character(factors)]]$YcUnscaled = unscale(model[[as.character(factors)]]$Yc, scale(as.matrix(mm),center=scale_center,scale=scale_scale))
        Yc = model[[as.character(factors)]]$YcUnscaled
        if (layer==1) { # add on layer 0 node
          markers = gsub("[+]|[-]","",colnames(mm)[phenoLevel==1])
          candidate_layer0 = lapply(markers[!duplicated(markers)], function(marker) if (sum(markers==marker)>1) { rowSums(Yc[,markers==marker]) } else { NULL } )
          candidate_layer0_mean = rowMeans(Reduce("cbind",candidate_layer0))
          Yc = cbind(candidate_layer0_mean,Yc)
        }
        models_combine[[as.character(layer)]] = Yc
      }
      Y = Reduce("cbind", models_combine)
      dimnames(Y) = dimnames(mm)
      save(Y,file=paste0(matrix_dir[ci], mcp,"PEER.Rdata"))
      
      
      
      
      
    }
    
  }
  
  cat("\n centre ", centre, " ",TimeOutput(start2)," \n",sep="")
}

TimeOutput(start)







