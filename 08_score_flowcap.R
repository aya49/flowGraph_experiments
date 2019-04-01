# Evaluate all distance metrics
# aya43@sfu.ca 20170116

root = "~/projects/flowCAP-II"
result_dir = "result"
setwd(root)

options(stringsAsFactors=FALSE)
#options(device="cairo")
options(na.rm=T)

#Input
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
dist_dir = paste(result_dir, "/dist", sep="")
doUnderflow = T #if numbers used to calculate scores --> 0, do work-around (slow)

KOonly = F #delete all WT samples
adjustD = T #divide by a number if the numbers are too big e.g. manhattan
overwrite = F #redo and overwrite all past scores

ignoredist = ".csv"


#Output
dist_score_dir = paste(dist_dir, "/dist_score", sep=""); for (i in 1:length(dist_score_dir)) { suppressWarnings(dir.create(dist_score_dir[i])) }
dist_score_result_dir = paste0(dist_score_dir,"/score_nac_list.Rdata")

source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
libr("stringr")
libr("foreach")
libr("doMC")
libr("Brobdingnag")
libr("lubridate") #if there are date variables

interested = c("aml","specimen","tube") #sampleMeta columns to plot
splitby = c("tube","tube","none")
dis = c("euclidean", "manhattan", "canberra")#, "mahalanobis") #assume after splitting dist filename by "_", distance is second element

start = Sys.time()

no_cores = detectCores()-1
registerDoMC(no_cores)

ftWTGT = "normal"

sampleMeta0 = get(load(sampleMeta_dir))
distmfile = list.files(dist_dir, recursive=F, full.names=T, pattern=".Rdata")
delpvalind = grepl(ignoredist,distmfile,ignore.case=T)
if (sum(delpvalind)>0) distmfile = distmfile[!delpvalind]
if (!overwrite) {
  result0 = get(load(dist_score_result_dir))
  doneind = which(!is.na(match(names(result0),distmfile)))
  if (length(doneind)>0) distmfile = distmfile[-doneind]
}
#distmfile = sample(distmfile,15)
distmfilenames = fileNames(distmfile)
interestedCols0 = which(colnames(sampleMeta0)%in%interested)


distMeta = distMetafun(distmfile,dis)


#calculate score_nac for each distance matrix & interested classes

result = foreach(i=1:length(distmfile)) %dopar% {
  tryCatch({
    start1 = Sys.time()
    cat("\n", distmfile[i], ": loading dist matrix", sep="")
    d0 = as.matrix(get(load(distmfile[i]))); if (sum(is.na(as.matrix(d0)))>0 | sum(as.matrix(d0)==Inf)>0) return(NULL)
    drows = colnames(d0)
    
    #adjust distance values
    if (adjustD) if (min(d0[which(d0>0)])>1 & mean(d0)>10000) d0 = d0/mean(d0)
    
    #adjust sampleMeta
    dcol = rep(0,ncol(sampleMeta0))
    for (j in 1:ncol(sampleMeta0)) {
      l = sum(!is.na(match(colnames(d0),sampleMeta0[,j])))
      if (l>0) dcol[j] = l
    }
    dcol = which.max(dcol)
    if (length(unique(sampleMeta0[,dcol]))==nrow(sampleMeta0)) {
      sampleMeta = sampleMeta0[match(colnames(d0),sampleMeta0[,dcol]),]
    } else { return(NULL) }
    
    #upload original matrix when needed
    
    #for each interested column
    pp = list()
    for (col in 1:length(interested)) { 
      colnam = paste0("interested-",interested[col],"_splitby-",splitby[col])
      
      #split by another column first, assuming d is by filename
      split=T
      if (splitby[col]=="none") split=F
      d = splitmatrix(as.matrix(d0),sampleMeta[,colnames(sampleMeta)==splitby[col]],split=split)
      sm = splitmatrix(sampleMeta,sampleMeta[,colnames(sampleMeta)==splitby[col]],split=split,missym=F)
      
      
      #for each split
      for (dind in 1:length(d)) {
        dindname0 = paste0(colnam,"_",names(d)[dind])
        dindname = names(d)[dind]
        cat("\n",interested[col],": ",names(d)[dind],".", sep="")
        
        #result[[i]][[colnames(sampleMeta)[interestedCols[j]]]] = NCA_score(d, sampleMeta[,interestedCols[j]], doUnderflow=doUnderflow) # 15s for 2500*2500 d
        pp[[dindname0]] = NCA_score(d[[dindname]], sm[[dindname]][,interested[col]], doUnderflow=doUnderflow) #2min, 21s
        TimeOutput(start1)
      }
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
if (!length(distmfile_working)>0) { result = NULL; break }
result = result[distmfile_working]
if (!overwrite) result = append(result0,result)
distmfile = names(result)
distmfilenames = names(result) = fileNames(distmfile)
save(result,file=dist_score_result_dir)

# split scores from distance matrices on differenct sets of objects
  drows = lapply(result, function(x) x$drows)
  result = lapply(result, function(x) x$pp)

names(result) = names(drows) = distmfilenames

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
  distMeta1 = distMeta[resultinds,]
  filenames = distmfilenames[resultinds]
  
  metatable = list()
  metawidth = c()
  for (fnamesi in 1:length(filenames)) {
    a = unlist(strsplit(gsub(".Rdata","",filenames[fnamesi]),"_"))
    astart = which(a%in%dis)
    if (!length(astart)>0) astart = grep("layer",a) #for linear, because all manhattan
    if (astart>2) a = c(paste(a[1:(astart-1)],collapse=""), a[astart:length(a)])
    # filnamestable[fnamesi,] = a
    metatable[[fnamesi]] = a
    metawidth = append(metawidth,length(a))
  }
  filnamestable = foreach(r=1:length(metatable), .combine='rbind') %dopar% {
    if (metawidth[r]!=max(metawidth)) return(append(metatable[[r]],rep("",max(metawidth)-metawidth[r])))
    return(metatable[[r]])
  }
  if (is.null(dim(filnamestable))) filnamestable = matrix(filnamestable,nrow=1)
  for (j in 1:length(result[[resultinds[1]]])) {
    pyt_nac = matrix(0, nrow=length(resultinds), ncol=1+length(result[[resultinds[1]]][[j]]$pyt), dimnames=list(filenames,as.vector(c("score",paste(result[[resultinds[1]]][[j]]$yt,names(result[[resultinds[1]]][[j]]$yt))))))
    for (i in 1:length(resultinds)) {
      tryCatch({
        pyt_nac[i,] = c(result[[resultinds[i]]][[j]]$p, result[[resultinds[i]]][[j]]$pyt)
        p_nac[i,j] = result[[resultinds[i]]][[j]]$p
      }, error = function(err) { print(paste("ERROR: ", i, j,err)) })
    }
    if (!grepl("linear",filenames[1])) { write.csv(cbind(filnamestable,pyt_nac), file=paste0(dist_score_dir,"/score_nac_",k,"_",names(result[[resultinds[1]]])[j],".csv"), row.names=F)
    } else {write.csv(pyt_nac, file=paste0(dist_score_dir,"/score_nac_",k,"_",names(result[[resultinds[1]]])[j],".csv")) }
  }
  write.csv(cbind(filnamestable,p_nac), file=paste0(dist_score_dir,"/score_nac_",k,".csv"))
}


TimeOutput(start)








