## input: a specified feature matrix with control and experiment fcm files (must run 01.5_meta_timegroup.R before this script for impc)
## output: p value feature and their adjustments for experiment files (we split up the control in flowcap, and label some of them as simply healthy files and treat them as experiment as well)


## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr(c("stringr", "lubridate", "Matrix",
       "foreach", "doMC"))


## cores
no_cores = detectCores() - 2
registerDoMC(no_cores)



## options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

writecsv = F

adjust = c("","BY","BH","bonferroni") #pvalue adjustment
#test = "wilcox" #pvalue test
cellCountThres = 1000 #insignificant if count under
pval_thres = .025 #delete phenotypes/rows without any significant changes from the pVal matrix
good_sample = 3 #only compare if >=3 samples available
good_sample_wt = 15 #min 70 wt used to compare with KO

control = "control"

id_col = "id"
target_col = "class"
split_col = "gender"
split_date = NULL
# split_col = "group"
# split_date = "date" # null if no date col associated with split_col groups
max_control = 70 # set to null if don't require a max amount of controls samples to compare other samples to, only needed when split_date!=NULL

feat_types = c("file-cell-countAdj") #, "file-cell-countAdj.PEER-layerbylayer","file-cell-countAdj.PEER-all")


for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)[-16]) {
  # result_dir = paste0(root, "/result/impc_panel1_sanger-spleen") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta")
  meta_file_dir = paste(meta_dir, "/file", sep="")
  feat_dir = paste(result_dir, "/feat", sep="")
  
  ## output directories
  feat_file_cell_pval_dir = paste(feat_dir, "/file-cell-pval",sep="")
  feat_file_cell_logfold_dir = paste(feat_dir, "/file-cell-logfold",sep="")
  feat_file_cell_countAdjMax_dir = paste(feat_dir, "/file-cell-countAdjMax",sep="")
  feat_file_cell_countAdjKO_dir = paste(feat_dir, "/file-cell-countAdjKO",sep="")
  
  
  
  start = Sys.time()
  
  
  meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
  if (length(unique(meta_file0[,target_col]))<3) next
  
  for (feat_type in feat_types) {
    start1 = Sys.time()
    
    #load matrices
    m = as.matrix(get(load(paste0(feat_dir,"/", feat_type,".Rdata"))))
    morder = match(rownames(m),meta_file0[,id_col]); morder = morder[!is.na(morder)]
    meta_file = meta_file0[morder,]
    
    #wildtypes
    g = getGTindex(meta_file[,target_col], control, good_sample, meta_file[,id_col])
    ftGT = g$attlist; controli = g$controlIndex; ftKOGT = g$exp; ftKOIndex = g$expIndex; ftKOgoodGTIndex = g$goodexpIndex; rm(g)  #Index of unique KO genotypes that has 3+ samples available ftKOIndex(ftKOGTIndexIndex)
    
    #split days up based on mean of wt; pvalues will only be calculated with WT samples that were made on similar days as the ko in question
    controli = grep(control, meta_file[,target_col])
    rowcombos = NULL
    for (i in nrow(m):1) {
      rowcombos[[i]] = list()
      rowcombos[[i]][[1]] = controli
      rowcombos[[i]][[2]] = i
    } 
    if (!is.null(split_col) & split_col%in%colnames(meta_file)) {
      for (i in nrow(m):1) {
        i_group = meta_file[i,split_col]
        control_group = meta_file[controli,split_col]
        rowcombos[[i]][[1]] = wtind = controli[control_group==i_group]
        if (!is.null(split_date)) {
          control_date = meta_file[controli, split_date]
          datediff = abs(ymd(meta_file[wtind,split_date])-ymd(meta_file[i,split_date]))
          rowcombos[[i]][[1]] = wtind[order(datediff)[1:ifelse(!is.null(max_control),max_control,length(datediff))]]
        } 
      }
    }
    # #split days up based on mean of wt; pvalues will only be calculated with WT samples that were made on similar days as the ko in question
    # rowcombos = NULL
    # for (i in nrow(meta_file):1) {
    #   kodate_group = meta_file$group[i]
    #   wtdates = meta_file[ftWTIndex,time_col]
    #   wtdates_group = meta_file$group[ftWTIndex]
    #   wtind = ftWTIndex[wtdates_group==kodate_group]
    #   # wtind = ftKOIndex[[i]] #use above if care about date groups, else, comment out above and use this line
    #   datediff = abs(ymd(meta_file[wtind,time_col])-ymd(meta_file[i,time_col]))
    #   wtindchosen = NULL
    #   for (j in 1:length(datediff)) {
    #     wtindchosen = append(wtindchosen,wtind[datediff==min(datediff)])
    #     datediff[datediff==min(datediff)] = Inf
    #     if (length(wtindchosen)>=70) break
    #   }
    #   #wtfn = meta_file$fileName[wtind[order(datediff)]]
    #   #rowcombos[[i]][[1]] = match(wtfn[1:(min(good_sample_wt,length(wtdates)))], meta_file$fileName) #should be max(good_sample_wt,length(wtdates)) assuming that wtdates should have >70 samples
    #   wtindchoseng = wtindchosen[meta_file[wtindchosen,split_col]==meta_file[i,split_col]]
    #   if (length(wtindchoseng)<min(good_sample_wt,length(wtindchosen))) wtindchoseng = wtindchosen
    #   rowcombos[[i]][[1]] = wtindchoseng
    #   rowcombos[[i]][[2]] = i
    # }
    
    cat(paste("getting pValues of ", ncol(m), " possible cell populations for ", length(rowcombos), " genotypes ", sep="")) #3iTCell specific
    
    loop.ind = 1:ncol(m)
    loop.inds = split(loop.ind, sort(loop.ind%%(no_cores)))
    result = foreach(loop.ind = loop.inds) %dopar% { #for each phenotype
      #for (k in 1:ncol(m)){ cat(paste(" ", j, sep="")) {
      
      stuff = list()
      for (k in loop.ind) {
        
        pvalcol = rep(0,length(rowcombos)) 
        logfold = rep(0,length(rowcombos))
        maxcount = rep(0,length(rowcombos))
        kocount = rep(0,length(rowcombos))
        for (j in 1:length(rowcombos)) { #for each KO gene
          compare1 = as.numeric(m[ rowcombos[[j]][[1]],k ])
          compare2 = as.numeric(m[ rowcombos[[j]][[2]],k ])
          if ((!median(compare1)<cellCountThres & !mean(compare2)<cellCountThres) | !grepl("prop",feat_type,ignore.case=T)) { #TRIM: if both WT and KO medians (to avoid outler influence) is < cell count threshold then change is not significant (set as 1)
            # if (test=="wilcox") {
            #   pvalcol[j] = wilcox.test(compare1, compare2)$p.value
            # } 
            # else if (test=="ttest") {
            #   try({ pvalcol[j] = t.test(compare1, compare2)$p.value })
            # } 
            all0 = all(append(compare1,compare2)==0) | all(compare2==compare1)
            pvalcol[j] = ifelse(all0,1, t.test.single(compare1, compare2))
            logfold[j] = ifelse(all0,1, log(compare2/exp(mean(log(compare1)))))
            maxcount[j] = max(compare2,exp(mean(log(compare1))))
            kocount[j] = compare2
          }
        }
        stuff[[k]] = list(pvalcol=pvalcol,logfold=logfold,maxcount=maxcount, kocount=kocount)
      }
      pvalcols = foreach(l=1:length(stuff),.combine="cbind") %dopar% { return(stuff[[l]]$pvalcol) }
      logfolds = foreach(l=1:length(stuff),.combine="cbind") %dopar% { return(stuff[[l]]$logfold) }
      maxcounts = foreach(l=1:length(stuff),.combine="cbind") %dopar% { return(stuff[[l]]$maxcount) }
      kocounts = foreach(l=1:length(stuff),.combine="cbind") %dopar% { return(stuff[[l]]$kocount) }
      
      return(list(pvalcol=pvalcols,logfold=logfolds,maxcount=maxcounts, kocount=kocounts))
    }
    
    # feat_file_cell_pval = lapply(1:length(result), function(i) return(result[[i]]$pvalcol))
    # feat_file_cell_logfold = lapply(1:length(result), function(i) return(result[[i]]$logfold))
    # feat_file_cell_countAdjMax = lapply(1:length(result), function(i) return(result[[i]]$maxcount))
    # feat_file_cell_countAdjKO = lapply(1:length(result), function(i) return(result[[i]]$kocount))
    # time_output(start1)
    
    # feat_file_cell_pvalFULL1 = Reduce("cbind",feat_file_cell_pval)
    feat_file_cell_pvalFULLori = foreach(k=1:length(result),.combine="cbind") %dopar% { return(result[[k]]$pvalcol) }
    feat_file_cell_logfoldFULL = foreach(k=1:length(result),.combine="cbind") %dopar% { return(result[[k]]$logfold) }
    feat_file_cell_countAdjMaxFULL = foreach(k=1:length(result),.combine="cbind") %dopar% { return(result[[k]]$maxcount) }
    feat_file_cell_countAdjKOFULL = foreach(k=1:length(result),.combine="cbind") %dopar% { return(result[[k]]$kocount) }
    time_output(start1)
    
    colnames(feat_file_cell_pvalFULLori) = colnames(feat_file_cell_logfoldFULL) = colnames(feat_file_cell_countAdjMaxFULL) = colnames(feat_file_cell_countAdjKOFULL) = colnames(m)
    rownames(feat_file_cell_pvalFULLori) = rownames(feat_file_cell_logfoldFULL) = rownames(feat_file_cell_countAdjMaxFULL) = rownames(feat_file_cell_countAdjKOFULL) = rownames(m)
    
    save(feat_file_cell_countAdjMaxFULL, file=paste0(feat_file_cell_countAdjMax_dir, "FULL.", feat_type,".Rdata"))
    save(feat_file_cell_countAdjKOFULL, file=paste0(feat_file_cell_countAdjKO_dir,"FULL.",feat_type,".Rdata"))
    
    #trim/mod pvalues
    for (adj in adjust) {
      feat_file_cell_pvalFULL0 = feat_file_cell_pvalFULLori
      if (adj!="") {
        feat_file_cell_pvalFULL0 = foreach(i=1:ncol(feat_file_cell_pvalFULLori), .combine = 'cbind') %dopar% {
          return(p.adjust(feat_file_cell_pvalFULLori[,i], method=adj))
        }
        colnames(feat_file_cell_pvalFULL0) = colnames(feat_file_cell_pvalFULLori)
        rownames(feat_file_cell_pvalFULL0) = rownames(feat_file_cell_pvalFULLori)
      }
      
      feat_file_cell_pvalFULL = -log(feat_file_cell_pvalFULL0)
      feat_file_cell_pvalFULL[is.na(feat_file_cell_pvalFULL)] = 0
      feat_file_cell_pvalFULL[feat_file_cell_pvalFULL==Inf] = 10^(ceiling(log(max(feat_file_cell_pvalFULL[feat_file_cell_pvalFULL!=Inf]),10)))
      feat_file_cell_pvalFULL = Matrix(feat_file_cell_pvalFULL, sparse=T)
      save(feat_file_cell_pvalFULL, file=paste0(feat_file_cell_pval_dir, adj, "FULL.", feat_type,".Rdata"))
      
      trimRowIndex <- apply(feat_file_cell_pvalFULL0[,-1], 1, function(x) all(x<=(pval_thres)))
      trimColIndex <- apply(feat_file_cell_pvalFULL0[-1,], 2, function(x) all(x<=(pval_thres)))
      
      feat_file_cell_pval = feat_file_cell_pvalTRIM = feat_file_cell_pvalFULL[!trimRowIndex,!trimColIndex]
      save(feat_file_cell_pval, file=paste0(feat_file_cell_pval_dir, adj, ".", feat_type,".Rdata"))
      if (writecsv) write.csv(feat_file_cell_pval, file=paste0(feat_file_cell_pval_dir, adj, ".", feat_type,".csv"), row.names=T)
      
      trimIndex = feat_file_cell_pval <= -log(pval_thres)
      
      feat_file_cell_pvalTRIM[trimIndex] = 0
      if (writecsv) write.csv(feat_file_cell_pvalTRIM, file=paste0(feat_file_cell_pval_dir, adj, ".", feat_type,".TRIM.csv"), row.names=T)
      feat_file_cell_pvalTRIM = Matrix(feat_file_cell_pvalTRIM, sparse=T, dimnames=dimnames(feat_file_cell_pvalTRIM))
      save(feat_file_cell_pvalTRIM, file=paste0(feat_file_cell_pval_dir, adj, ".", feat_type,".TRIM.Rdata"))
      
      mt = m
      mt[as.matrix(trimIndex)] = 0
      if (writecsv) write.csv(mt, file=paste0(feat_dir,"/",feat_type, adj, ".", feat_type,".TRIM.csv"), row.names=T)
      mt = Matrix(mt, sparse=T, dimnames=dimnames(m))
      save(mt, file=paste0(feat_dir,"/",feat_type, adj, ".", feat_type,".TRIM.Rdata"))
      
      if (adj=="") { #don't need to trim others with adjusted p values, too much space taken up lol
        save(feat_file_cell_logfoldFULL, file=paste0(feat_file_cell_logfold_dir, "FULL.", feat_type,".Rdata"))
        
        feat_file_cell_logfold = feat_file_cell_logfoldTRIM = feat_file_cell_logfoldFULL[!trimRowIndex,!trimColIndex]
        save(feat_file_cell_logfold, file=paste0(feat_file_cell_logfold_dir, ".", feat_type,".Rdata"))
        if (writecsv) write.csv(feat_file_cell_logfold, file=paste0(feat_file_cell_logfold_dir, ".", feat_type,".csv"))
        
        feat_file_cell_logfoldTRIM[as.matrix(trimIndex)] = 0
        if (writecsv) write.csv(feat_file_cell_logfoldTRIM, file=paste0(feat_file_cell_logfold_dir, ".", feat_type,".TRIM.csv"), row.names=T)
        feat_file_cell_logfoldTRIM = Matrix(feat_file_cell_logfoldTRIM, sparse=T, dimnames=dimnames(feat_file_cell_logfoldTRIM))
        save(feat_file_cell_logfoldTRIM, file=paste0(feat_file_cell_logfold_dir, ".", feat_type,".TRIM.Rdata"))
        
        feat_file_cell_countAdjMax = feat_file_cell_countAdjMaxFULL[!trimRowIndex,!trimColIndex]
        save(feat_file_cell_countAdjMax, file=paste0(feat_file_cell_countAdjMax_dir, ".", feat_type,".Rdata"))
        if (writecsv) write.csv(feat_file_cell_countAdjMax, file=paste0(feat_file_cell_countAdjMax_dir, ".", feat_type,".csv"))
      }
      
    }
    
    cat("\n feat_type", feat_type,": ",time_output(start1), "\n", sep="") #IMPC Sanger P1 ~3h
  }
  
  cat("\nTime taken to calculate p values & barcode matrices is: ",time_output(start), "\n", sep="") #3iTcell ~40min
  
  # change "normal" to "control", normal is just used to do p value in flowcap
  if (grepl("flowcap", result_dir)) {
    meta_file$class[meta_file$class=="normal"] = "control"
    save(meta_file,file=paste0(meta_file_dir,".Rdata"))
  }
}




