# aya43@sfu.ca 20151228
# Calculate PValues; must run 02_plot.R before this script

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
meta_cell_parent_dir = paste(result_dir, "/phenoParent.Rdata",sep="")
meta_cell_parent_ind_dir = paste(result_dir, "/phenoParent_ind.Rdata",sep="")
feat_dir = paste(result_dir, "/feat", sep="")

## output directories
feat_file_cell_pval_dir = paste(result_dir, "/file-cell-pval",sep="")
feat_file_cell_pvalFULL_dir = paste(result_dir, "/file-cell-pvalFULL",sep="")
feat_file_cell_pvalTRIM_dir = paste(result_dir, "/file-cell-pvalTRIM",sep="")
feat_file_cell_logfold_dir = paste(result_dir, "/file-cell-logfold",sep="")
feat_file_cell_logfoldFULL_dir = paste(result_dir, "/file-cell-logfoldFULL",sep="")
feat_file_cell_logfoldTRIM_dir = paste(result_dir, "/file-cell-logfoldTRIM",sep="")
feat_file_cell_countAdjMax_dir = paste(result_dir, "/file-cell-countAdjMax",sep="")
feat_file_cell_countAdjMaxFULL_dir = paste(result_dir, "/file-cell-countAdjMaxFULL",sep="")
file_cell_countAdjKOFULL_dir = paste(result_dir, "/file-cell-countAdjKOFULL",sep="")

## libraries
libr(stringr)
libr(foreach)
libr(doMC)
libr(lubridate)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")


## cores
no_cores = detectCores() - 1
registerDoMC(no_cores)



## options
#adjust = "BY" #"BH" #pvalue adjustment
#test = "wilcox" #pvalue test
cellCountThres = 200 #insignificant if count under
pval_thres = .025 #delete phenotypes/rows without any significant changes from the pVal matrix
good_sample = 3 #only compare if >=3 samples available
good_sample_wt = 15 #min 70 wt used to compare with KO

control = str_split(controlL,"[|]")[[1]]



start = Sys.time()


meta_file = get(load(paste0(meta_file_dir,".Rdata")))
                meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
                
                for (feat_type in feat_types) {
                  start1 = Sys.time()
                  
                  #load matrices
                  m = get(load(paste0(feat_dir,"/", feat_type,".Rdata")))
                  
                  #wildtypes
                  ftWTGT = c("+_+","+_Y"); if (!grepl("Sanger",centre)) ftWTGT = "WildType"
                  g = getGTindex(meta_file$gene, ftWTGT, good_sample, meta_file$fileName)
                  ftGT = g$attlist; ftWTIndex = g$controlIndex; ftKOGT = g$exp; ftKOIndex = g$expIndex; ftKOgoodGTIndex = g$goodexpIndex; rm(g)  #Index of unique KO genotypes that has 3+ samples available ftKOIndex(ftKOGTIndexIndex)
                  
                  #split days up based on mean of wt; pvalues will only be calculated with WT samples that were made on similar days as the ko in question
                  rowcombos = NULL
                  for (i in length(ftKOIndex):1) {
                    kodate_group = meta_file$group[ftKOIndex[[i]]]
                    wtdates = meta_file$date[ftWTIndex]
                    wtdates_group = meta_file$group[ftWTIndex]
                    wtind = ftWTIndex[wtdates_group==kodate_group]
                    # wtind = ftKOIndex[[i]] #use above if care about date groups, else, comment out above and use this line
                    datediff = abs(ymd(meta_file$date[wtind])-ymd(meta_file$date[ftKOIndex[[i]]]))
                    wtindchosen = NULL
                    for (j in 1:length(datediff)) {
                      wtindchosen = append(wtindchosen,wtind[datediff==min(datediff)])
                      datediff[datediff==min(datediff)] = Inf
                      if (length(wtindchosen)>=70) break
                    }
                    #wtfn = meta_file$fileName[wtind[order(datediff)]]
                    #rowcombos[[i]][[1]] = match(wtfn[1:(min(good_sample_wt,length(wtdates)))], meta_file$fileName) #should be max(good_sample_wt,length(wtdates)) assuming that wtdates should have >70 samples
                    wtindchoseng = wtindchosen[meta_file$gender[wtindchosen]==meta_file$gender[ftKOIndex[[i]]]]
                    if (length(wtindchoseng)<min(good_sample_wt,length(wtindchosen))) wtindchoseng = wtindchosen
                    rowcombos[[i]][[1]] = wtindchoseng
                    rowcombos[[i]][[2]] = ftKOIndex[[i]] 
                  }
                  
                  cat(paste("getting pValues of ", ncol(m), " possible cell populations for ", length(rowcombos), " genotypes ", sep="")) #3iTCell specific
                  
                  loop.ind = 1:ncol(m)
                  result = foreach(k = loop.ind) %dopar% { #for each phenotype
                    #for (k in 1:ncol(m)){ cat(paste(" ", j, sep="")) {
                    
                    pvalcol = rep(0,length(rowcombos))
                    logfold = rep(0,length(rowcombos))
                    maxcount = rep(0,length(rowcombos))
                    kocount = rep(0,length(rowcombos))
                    for (j in 1:length(rowcombos)) { #for each KO gene
                      compare1 = as.numeric(m[ rowcombos[[j]][[1]],k ])
                      compare2 = as.numeric(m[ rowcombos[[j]][[2]],k ])
                      if (!exp(mean(log(compare1)))<cellCountThres & !exp(mean(log(compare2)))<cellCountThres & !grepl("Prop",feat_type)) { #TRIM: if both WT and KO medians (to avoid outler influence) is < cell count threshold then change is not significant (set as 1)
                        # if (test=="wilcox") {
                        #   pvalcol[j] = wilcox.test(compare1, compare2)$p.value
                        # } 
                        # else if (test=="ttest") {
                        #   try({ pvalcol[j] = t.test(compare1, compare2)$p.value })
                        # } 
                        pvalcol[j] = -log(t.test.single(compare1, compare2))
                        logfold[j] = log(compare2/exp(mean(log(compare1))))
                        maxcount[j] = max(compare2,exp(mean(log(compare1))))
                        kocount[j] = compare2
                      }
                    }
                    return(list(pvalcol=pvalcol,logfold=logfold,maxcount=maxcount, kocount=kocount))
                  }
                  
                  feat_file_cell_pval = lapply(1:length(result), function(i) return(result[[i]]$pvalcol))
                  feat_file_cell_logfold = lapply(1:length(result), function(i) return(result[[i]]$logfold))
                  feat_file_cell_countAdjMax = lapply(1:length(result), function(i) return(result[[i]]$maxcount))
                  file_cell_countAdjKO = lapply(1:length(result), function(i) return(result[[i]]$kocount))
                  TimeOutput(start1)
                  
                  feat_file_cell_pvalFULL = foreach(k=1:length(feat_file_cell_pval),.combine="cbind") %dopar% { return(feat_file_cell_pval[[k]]) }
                  feat_file_cell_logfoldFULL = foreach(k=1:length(feat_file_cell_logfold),.combine="cbind") %dopar% { return(feat_file_cell_logfold[[k]]) }
                  feat_file_cell_countAdjMaxFULL = foreach(k=1:length(feat_file_cell_countAdjMax),.combine="cbind") %dopar% { return(feat_file_cell_countAdjMax[[k]]) }
                  file_cell_countAdjKOFULL = foreach(k=1:length(file_cell_countAdjKO),.combine="cbind") %dopar% { return(file_cell_countAdjKO[[k]]) }
                  TimeOutput(start1)
                  
                  feat_file_cell_pvalFULL[is.nan(feat_file_cell_pvalFULL)] = 0
                  feat_file_cell_pvalFULL[feat_file_cell_pvalFULL==Inf] = 10^(ceiling(log(max(feat_file_cell_pvalFULL[feat_file_cell_pvalFULL!=Inf]),10)))
                  colnames(feat_file_cell_pvalFULL) = colnames(feat_file_cell_logfoldFULL) = colnames(feat_file_cell_countAdjMaxFULL) = colnames(file_cell_countAdjKOFULL) = colnames(m)
                  rownames(feat_file_cell_pvalFULL) = rownames(feat_file_cell_logfoldFULL) = rownames(feat_file_cell_countAdjMaxFULL) = rownames(file_cell_countAdjKOFULL) = ftKOGT
                  
                  # adjust when doing Wilcox
                  # feat_file_cell_pvalAdj = foreach(i=1:ncol(feat_file_cell_pval), .combine='cbind') %dopar% {
                  #   return(p.adjust(feat_file_cell_pval[,i], method=adjust))
                  # }
                  # feat_file_cell_pvalFULL = feat_file_cell_pvalAdj
                  
                  #trim pvalues
                  trimRowIndex <- which(apply(feat_file_cell_pvalFULL[,-1], 1, function(x) all(x<=-log(pval_thres)))==T)
                  trimColIndex <- which(apply(feat_file_cell_pvalFULL[-1,], 2, function(x) all(x<=-log(pval_thres)))==T)
                  
                  feat_file_cell_pvalTRIM = feat_file_cell_pval = feat_file_cell_pvalFULL[-trimRowIndex,-trimColIndex]
                  feat_file_cell_logfoldTRIM = feat_file_cell_logfold = feat_file_cell_logfoldFULL[-trimRowIndex,-trimColIndex]
                  feat_file_cell_countAdjMax = feat_file_cell_countAdjMaxFULL[-trimRowIndex,-trimColIndex]
                  
                  feat_file_cell_pvalTRIM[feat_file_cell_pvalTRIM<=-log(pval_thres)] = 0
                  feat_file_cell_logfoldTRIM[feat_file_cell_pvalTRIM<=-log(pval_thres)] = 0
                  TimeOutput(start1)
                  
                  save(feat_file_cell_pval, file=paste0(feat_file_cell_pval_dir, "_", feat_type,".Rdata"))
                  write.table(feat_file_cell_pval, file=paste0(feat_file_cell_pval_dir, "_", feat_type,".csv"))
                  save(feat_file_cell_pvalFULL, file=paste0(feat_file_cell_pvalFULL_dir, "_", feat_type,".Rdata"))
                  save(feat_file_cell_pvalTRIM, file=paste0(feat_file_cell_pvalTRIM_dir, "_", feat_type,".Rdata"))
                  a = feat_file_cell_pvalTRIM; rownames(a) = meta_file$gene[match(rownames(a),meta_file$fileName)]
                  write.table(a, file=paste0(feat_file_cell_pvalTRIM_dir, "_", feat_type,".csv"), sep=",", col.names=F)
                  save(feat_file_cell_logfold, file=paste0(feat_file_cell_logfold_dir, "_", feat_type,".Rdata"))
                  write.table(feat_file_cell_logfold, file=paste0(feat_file_cell_logfold_dir, "_", feat_type,".csv"))
                  save(feat_file_cell_logfoldFULL, file=paste0(feat_file_cell_logfoldFULL_dir, "_", feat_type,".Rdata"))
                  save(feat_file_cell_logfoldTRIM, file=paste0(feat_file_cell_logfoldTRIM_dir, "_", feat_type,".Rdata"))
                  b = feat_file_cell_logfoldTRIM; rownames(b) = meta_file$gene[match(rownames(b),meta_file$fileName)]
                  write.table(b, file=paste0(feat_file_cell_logfoldTRIM_dir, "_", feat_type,".csv"), sep=",", col.names=F)
                  save(feat_file_cell_countAdjMax, file=paste0(feat_file_cell_countAdjMax_dir, "_", feat_type,".Rdata"))
                  write.table(feat_file_cell_countAdjMax, file=paste0(feat_file_cell_countAdjMax_dir, "_", feat_type,".csv"))
                  save(feat_file_cell_countAdjMaxFULL, file=paste0(feat_file_cell_countAdjMaxFULL_dir, "_", feat_type,".Rdata"))
                  
                  save(file_cell_countAdjKOFULL, file=paste0(file_cell_countAdjKOFULL_dir,"_",feat_type,".Rdata"))
                  
                  meta_file_trim = meta_file[-trimRowIndex,]
                  meta_cell_trim = meta_cell[-trimColIndex,]
                  save(meta_file_trim, file=paste0(meta_file_trim_dir, "_", feat_type,".Rdata"))
                  save(meta_cell_trim, file=paste0(meta_cell_trim_dir, "_", feat_type,".Rdata"))
                  
                  cat("\n centre", centre,": ",TimeOutput(start1), "\n", sep="") #3iTcell ~40min
                }
                }

cat("\nTime taken to calculate p values & barcode matrices is: ",TimeOutput(start), "\n", sep="") #3iTcell ~40min
