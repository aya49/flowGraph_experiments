# aya43@sfu.ca 20170116
# Evaluate all distance metrics

#root = "~/projects/IMPC"
result_dir = "result"
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#"Sanger_MLN","CIPHE","TCP","H")

options(stringsAsFactors=FALSE)
#options(device="cairo")
options(na.rm=T)

#Input
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
dist_dir = paste(result_dir, "/", panelL, "/", centreL, "/dist", sep="")
doUnderflow = T #if numbers used to calculate scores --> 0, do work-around (slow)

KOonly = T #delete all WT samples
adjustD = T #divide by a number if the numbers are too big e.g. manhattan

#Output
dist_score_dir = paste(dist_dir, "/dist_score", sep=""); for (i in 1:length(dist_score_dir)) { suppressWarnings(dir.create(dist_score_dir[i])) }

source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
libr("stringr")
libr("foreach")
libr("doMC")
libr("lubridate") #if there are date variables

interested = c("gene", "gender", "date", "colony", "strain", "fur", "birth_date", "specimen", "sample") #sampleMeta columns to plot

start = Sys.time()

no_cores = 6#detectCores()-1
registerDoMC(no_cores)

for (ci in 1:length(paste0(panelL,centreL))) {
  start2 = Sys.time()
  cat("\n", paste0(panelL," ",centreL)[ci], sep="")
  centre = paste0(panelL," ",centreL)[ci]
  ftWTGT = c("+_+","+_Y"); if (!grepl("Sanger",centre)) ftWTGT = "WildType"
  
  sampleMeta0 = get(load(sampleMeta_dir[ci]))
  distmfile = list.files(dist_dir[ci], recursive=F, full.names=T, pattern=".Rdata")
  #delpvalind = grep("pval",distmfile,ignore.case=T)
  #if (length(delpvalind)>0) distmfile = distmfile[-delpvalind]
  distmfilenames = fileNames(distmfile)
  interestedCols0 = which(colnames(sampleMeta0)%in%interested)
  
  #calculate score_nac for each distance matrix & interested classes
  #result = list()
  #for(i in 1:length(distmfile)) {
  result = foreach(i=1:length(distmfile)) %dopar% {
    #loop.ind = 1:length(distmfile)
    #for(i in 1:length(distmfile)) {
    #result = foreach(i=loop.ind) %dopar% { 
    start1 = Sys.time()
    cat("\n", distmfile[i], ", loading dist matrix", sep="")
    #result[[i]] = list()
    d = get(load(distmfile[i]))
    if (sum(is.na(as.matrix(d)))>0) return(NULL)
    drows = colnames(as.matrix(d))
    
    #adjust sampleMeta
    dcol = rep(0,ncol(sampleMeta0))
    for (j in 1:ncol(sampleMeta0)) {
      l = sum(!is.na(match(drows,sampleMeta0[,j])))
      if (l>0) dcol[j] = l
    }
    dcol = which.max(dcol)
    if (length(unique(sampleMeta0[,dcol]))==nrow(sampleMeta0)) {
      sampleMeta <- sampleMeta0[match(drows,sampleMeta0[,dcol]),]
      interestedCols = interestedCols0
      #distmfile[i] = paste0("noScore_",distmfile[i])
    } else {
      interestedCols = 1
      sampleMeta = as.data.frame(colnames(as.matrix(d)))
    }
    
    if (KOonly) {
      delind = which(sampleMeta$gene%in%ftWTGT)
      if (length(delind)>0) {
        sampleMeta = sampleMeta[-delind,]
        d = as.dist(as.matrix(d)[-delind,-delind])
      }
    }
    
    if (adjustD) {
      d = as.matrix(d)
      if (min(d[which(d>0)])>1 & mean(d)>10000) d = d/mean(d)
      d = as.dist(d)
    }
    
    #   pp = list()
    #   for (j in 1:length(interestedCols)) { #14.5min each; total 8.5hr
    #     cat(", ", colnames(sampleMeta)[interestedCols[j]], " ", sep="")
    #     #result[[i]][[colnames(sampleMeta)[interestedCols[j]]]] = NCA_score(d, sampleMeta[,interestedCols[j]], doUnderflow=doUnderflow) # 15s for 2500*2500 d
    #     pp[[colnames(sampleMeta)[interestedCols[j]]]] = NCA_score(d, sampleMeta[,interestedCols[j]], doUnderflow=doUnderflow) #2min
    #     TimeOutput(start1)
    #   }
    #   cat("; Done ", sep="")
    #   TimeOutput(start1)
    #   return(pp)
    # }
    pp = NULL
    cat(", calculating score for feature ", sep="")
    
    tryCatch({
      for (j in 1:length(interestedCols)) {
        cat(", ", colnames(sampleMeta)[interestedCols[j]], " ", sep="")
        #result[[i]][[colnames(sampleMeta)[interestedCols[j]]]] = NCA_score(d, sampleMeta[,interestedCols[j]], doUnderflow=doUnderflow) # 15s for 2500*2500 d
        pp[[colnames(sampleMeta)[interestedCols[j]]]] = NCA_score(d, sampleMeta[,interestedCols[j]], doUnderflow=doUnderflow) #2min, 21s
        TimeOutput(start1)
      }
      
      cat("; Done ", sep="")
      TimeOutput(start1)
      return(list(pp=pp,drows=drows))
    }, error = function(e) {
      cat(paste("ERROR:  ",e));
    })
  }
  
  names(result) = distmfile
  distmfile_working = which(vapply(result, Negate(is.null), NA))
  result = result[distmfile_working]
  distmfile = distmfile[distmfile_working]
  distmfilenames = distmfilenames[distmfile_working]
  if (length(which(vapply(result, Negate(is.null), NA)))>0) result = result[which(vapply(result, Negate(is.null), NA))]
  save(result,file=paste0(dist_score_dir[ci],"/score_nac_list.Rdata"))
  
  # split scores from distance matrices on differenct sets of objects
  drows = NULL
  for (i in 1:length(result)) {
    drows[[i]] = result[[i]]$drows
    result[[i]] = result[[i]]$pp
  }
  names(result) = distmfile
  
  drows_ind = list()
  drows_ind[[1]] = 1
  drows_factor = list()
  drows_factor[[1]] = drows[[1]]
  if (length(drows)>1) {
    for (i in 2:length(drows)) {
      matched = Position(function(x) identical(x, drows[[i]]), drows_factor, nomatch = 0)
      if (matched>0) {
        a = drows_ind[[matched]]
        drows_ind[[matched]] = append(a,i)
      } else {
        drows_factor[[length(drows_factor)+1]] = drows[[i]]
        drows_ind[[length(drows_ind)+1]] = c(i)
      }
    }
  }
  
  #compile score_nac results & save
  cat("; Compiling scores ")
  for (k in 1:length(drows_ind)) {
    resultinds = drows_ind[[k]]
    p_nac = matrix(0, nrow=length(resultinds), ncol=length(result[[resultinds[1]]]), dimnames=list(distmfilenames[resultinds], names(result[[resultinds[1]]])))
    filenames = distmfilenames[resultinds]
    if (!grepl("linear",filenames[1])) {
      filnamestable = matrix("",nrow=length(filenames),ncol=length(unlist(strsplit(filenames[1],"_"))))
      for (fnamesi in 1:length(filenames)) {
        a = unlist(strsplit(gsub(".Rdata","",filenames[fnamesi]),"_"))
        filnamestable[fnamesi,min(length(a),ncol(filnamestable))] = a[1:min(length(a),ncol(filnamestable))]
      }
    }
    for (j in 1:length(result[[resultinds[1]]])) {
      pyt_nac = matrix(0, nrow=length(resultinds), ncol=1+length(result[[resultinds[1]]][[j]]$pyt), dimnames=list(filenames,as.vector(c("score",paste(result[[resultinds[1]]][[j]]$yt,names(result[[resultinds[1]]][[j]]$yt))))))
      for (i in 1:length(resultinds)) {
        tryCatch({
          pyt_nac[i,] = c(result[[resultinds[i]]][[j]]$p, result[[resultinds[i]]][[j]]$pyt)
          p_nac[i,j] = result[[resultinds[i]]][[j]]$p
        }, error = function(err) { print(paste("MY_ERROR: ", i, j,err)) })
      }
      if (!grepl("linear",filenames[1])) { write.csv(cbind(filnamestable,pyt_nac), file=paste0(dist_score_dir[ci],"/score_nac_",k,"_",names(result[[resultinds[1]]])[j],".csv"), row.names=F)
      } else {write.csv(pyt_nac, file=paste0(dist_score_dir[ci],"/score_nac_",k,"_",names(result[[resultinds[1]]])[j],".csv")) }
    }
    write.csv(cbind(filnamestable,p_nac), file=paste0(dist_score_dir[ci],"/score_nac_",k,".csv"))
  }
  TimeOutput(start2)
  
}

TimeOutput(start)









