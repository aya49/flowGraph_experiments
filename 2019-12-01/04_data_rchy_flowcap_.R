## Input: original/trimmed features --> Output: trimmed features + rchyoptimyx additional nodes
# aya43@sfu.ca 20161220

#Directory
root = "~/projects/flowtype_metrics"
setwd(root)

start = Sys.time()
for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  
  #Options
  options(stringsAsFactors=FALSE)
  # options(device="cairo")
  options(na.rm=T)
  
  #Input
  meta_dir = paste0(result_dir,"/meta") # meta files directory
  meta_file_dir = paste(meta_dir, "/file", sep="") #meta for rows (sample)
  meta_cell_dir = paste(meta_dir, "/cell", sep="") #meta for rows (sample)
  feat_dir = paste(result_dir, "/feat", sep="") #feature files directory
  
  #Output
  rchy_dir = paste(result_dir, "/rchy", sep=""); dir.create(rchy_dir, showWarnings=F)
  
  
  #Libraries/Functions
  source("source/_funcAlice.R")
  libr(c("Matrix","stringr",
         "foreach","doMC",
         "RchyOptimyx"))
  
  #Setup Cores
  no_cores = detectCores()-1
  registerDoMC(no_cores)
  
  
  
  
  
  
  
  
  #Options for script
  overwrite = T #overwrite biclust?
  # writecsv = F
  
  good_count = 3 #trim matrix; only keep col/rows that meet criteria for more than 3 elements
  good_sample = 3 #trim matrix; only keep rows that are a part of a class with more than 3 samples
  
  k = 4 #1,4,max(meta_cell$phenolevel)) # how many markers to consider i.e. k=max(phenolevel) only
  countThres = 1000 #insignificant if count under
  
  target_col = "class" #the interested column in meta_file
  control = "control" #control value in target_col column
  id_col = "id" #the column in meta_file matching rownames in feature matrices
  order_cols = NULL #if matrix rows should be ordered by a certain column
  split_col = NULL # if certain rows in matrices should be analyzed in isolation, split matrix by this column in meta_file
  
  kpaths = 5
  
  # caluclating score; weigh different factors
  score_layer = .2
  score_count = .2
  score_score = .6 #of feat_score
  score_discretize = F #discretize score e.g. p value
  
  notrimmatrix = "Max"
  
  feat_count = c("file-cell-countAdj")
  feat_score = paste0(c("file-cell-pval."), feat_count) #matrix with 0s for insignificant nodes
  feat_value = paste0(c("file-cell-pval."), feat_count) #matrix as input into rchy
  
  feat_weight = paste0(c("file-cell-countAdjMax."), feat_count) #assume by cell pop
  
  # matrix_type = c("file-cell-countAdj.pval.file-cell-countAdj.TRIM")
  # matrix_all_type = c("file-cell-pval.file-cell-countAdj")
  # matrix_edge_type = c("Parent_contrib_CountAdj") #plot only
  
  feat_types = list.files(path=feat_dir,pattern=".Rdata")
  feat_types = feat_types[!grepl("TRIM|FULL",feat_types) & !grepl(notrimmatrix,feat_types)]
  feat_types = gsub(".Rdata","",feat_types)
  
  
  #Prepare data
  mc = get(load(paste0(feat_dir,"/", feat_count,".Rdata")))
  meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
  meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
  
  
  
  
  
  
  
  
  
  
  
  
  start = Sys.time()
  
  for (feat_type in feat_types) {
    tryCatch({
      # map = matrix_all_type[feat_type]
      # mep = matrix_edge_type[feat_type]
      # cat("\n", feat_type, " ",sep="")
      # start2 = Sys.time()
      # 
      # #start a list of phenotypes included in each distance matrix calculation such that none are repeated
      # leavePhenotype = list()
      # doneAll = F
      
      #load feature matrix
      m0 = as.matrix(get(load(paste0(feat_dir,"/", feat_type,".Rdata"))))
      ms0 = as.matrix(get(load(paste0(feat_dir,"/", feat_score,".Rdata"))))
      # mresult = Loadintermatrices(c(paste0(matrix_dir, feat_type,".Rdata"),paste0(matrix_dir, map,".Rdata"),paste0(matrix_dir, feat_type,".Rdata")))
      # mml0 = mresult$mml
      # mmlname = names(mml0)
      # pt = mpi = mresult$pt
      # gt = mresult$gt
      
      
      #trim / order feature & score matrix
      mm = trimMatrix(m0,TRIM=T, mc=mc, sampleMeta=meta_file0, sampleMeta_to_m1_col=id_col, target_col=target_col, control=control, order_cols=order_cols, colsplitlen=NULL, k=k, countThres=countThres, goodcount=good_count, good_sample=good_sample)
      if (is.null(mm)) next
      m = mm$m
      meta_file = mm$sm
      if (all(meta_file[,target_col]==meta_file[1,target_col])) next
      
      msm = trimMatrix(ms0,TRIM=T, mc=mc, sampleMeta=meta_file0, sampleMeta_to_m1_col=id_col, target_col=target_col, control=control, order_cols=order_cols, colsplitlen=NULL, k=k, countThres=countThres, goodcount=good_count, good_sample=good_sample)
      if (is.null(mm)) next
      ms = msm$m
      metas_file = msm$sm
      if (all(metas_file[,target_col]==metas_file[1,target_col])) next
      
      colorder = sort(intersect(colnames(m),colnames(ms)))
      roworder = sort(intersect(rownames(m),rownames(ms)))
      m = m[roworder,match(colorder,colnames(m))]
      ms = ms[roworder,match(colorder,colnames(ms))]
      meta_file = meta_file[roworder,]
      metas_file = metas_file[roworder,]
      
      phenocodes = meta_cell$phenocode[match(colnames(m),meta_cell$phenotype)]
      
      fname = paste0(feat_dir,"/",feat_type,".RCHY.dir.k",str_pad(kpaths, 2, pad = "0"),"_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres)
      fname0 = paste0(feat_dir,"/",feat_type,".RCHY.all.k",str_pad(kpaths, 2, pad = "0"),"_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres)
      
      a = foreach(ii=1:nrow(m)) %dopar% {
        sample = rownames(m)[ii]
        sample_feat = m[ii,]
        sample_score = ms[ii,]
        
        # get cell names that either change +/- in pvalues
        rchy0 = NULL
        if (! sum(sample_feat!=0) > good_count) next
        score = abs(sample_score)
        
        startpheno = meta_cell$phenocode[match(names(sample_feat)[sample_feat!=0], meta_cell$phenotype)]
        startpheno = setdiff(startpheno, paste0(rep(0,nchar(meta_cell$phenocode[1]))))
        
        ### TRIM STARTPHENO!!!!
        
        rchy0 = lapply(startpheno, function(sp) {
          rr = RchyOptimyx(pheno.codes=phenocodes, phenotypeScores=sample_score, startPhenotype=sp, pathCount=kpaths, trimPaths=F)
          rr@nodes[2,] = sample_score[match(rr@nodes[1,], phenocodes)]
          rownames(rr@nodes) = c("phenocode",feat_score,"colour")
          return(rr)
        })
        names(rchy0) = startpheno
        # rchy = rchy1[[1]]; if (length(rchy1)>1) for (ri in 2:length(rchy1)) { rchy = merge(rchy,rchy1[[ri]]) }
        
        
        rchy0_m = NULL
        if (!is.null(rchy0)) {
          rchy0_m = rchy0[[1]]
          if (length(rchy0)>1) {
            for (ri in 2:length(rchy0)) {
              rchy0_m = merge(rchy0_m,rchy0[[ri]])
            }
          }
        }
        
        
        
        rchy = NULL
        # positive change pvalues
        rchy1 = NULL
        if (sum(m[ii,]>0)>0) {
          phenoscore = ma[ii,]
          if (min(phenoscore)<0) { phenoscore = phenoscore+min(-ma[ii,]);  }
          startpheno = meta_cell$phenocode[match(colnames(m)[m[i,]>0],meta_cell$phenotype)]
          startpheno = setdiff(startpheno,paste0(rep(0,nchar(meta_cell$phenocode[1]))))
          rchy1 = lapply(startpheno, function(sp) {
            rr = RchyOptimyx(pheno.codes=phenocodes, phenotypeScores=phenoscore, startPhenotype=sp, pathCount=kpaths, trimPaths=F)
            rr@nodes[2,] = ma[ii,match(rr@nodes[1,],meta_cell$phenocode[match(colnames(ma),meta_cell$phenotype)])]
            rownames(rr@nodes) = c("phenocode",matrix_type[feat_type],"colour")
            return(rr)
          })
          names(rchy1) = startpheno
          # rchy = rchy1[[1]]; if (length(rchy1)>1) for (ri in 2:length(rchy1)) { rchy = merge(rchy,rchy1[[ri]]) }
        }
        # negative change pvalues
        rchy2 = NULL
        if (sum(m[ii,]<0)>0) {
          phenoscore = -ma[ii,]
          if (min(phenoscore)<0) phenoscore = phenoscore+min(-ma[ii,])
          startpheno = meta_cell$phenocode[match(colnames(m)[m[i,]<0],meta_cell$phenotype)]
          startpheno = setdiff(startpheno,paste0(rep(0,nchar(meta_cell$phenocode[1])),collapse=""))
          rchy2 = lapply(startpheno, function(sp) {
            rr = RchyOptimyx(pheno.codes=phenocodes, phenotypeScores=phenoscore, startPhenotype=sp, pathCount=kpaths, trimPaths=F)
            rr@nodes[2,] = ma[ii,match(rr@nodes[1,],meta_cell$phenocode[match(colnames(ma),meta_cell$phenotype)])]
            rownames(rr@nodes) = c("phenocode",matrix_type[feat_type],"colour")
            return(rr)
          })
          names(rchy2) = startpheno
          # if (!sum(m[i,]>0)>0) rchy = rchy2[[1]]
          # rchy = merge(rchy,rchy2[[1]]); if (length(rchy2)>1) for (ri in 2:length(rchy2)) { rchy = merge(rchy,rchy2[[ri]]) }
        }
        rchy = append(rchy1,rchy2)
        
        rchy_m = NULL
        if (!is.null(rchy)) {
          rchy_m = rchy[[1]]
          if (length(rchy)>1) {
            for (ri in 2:length(rchy)) {
              rchy_m = merge(rchy_m,rchy[[ri]])
            }
          }
        }
        
        if (!is.null(rchy0)) save(rchy0,file=paste0(fname0,"/",i,".Rdata"))
        if (!is.null(rchy)) save(rchy,file=paste0(fname,"/",i,".Rdata"))
        if (!is.null(rchy0_m)) save(rchy0,file=paste0(fname0,"/",i,"_merged.Rdata"))
        if (!is.null(rchy_m)) save(rchy,file=paste0(fname,"/",i,"_merged.Rdata"))
        
      }
      
      TimeOutput(start2)
    }, error = function(err) { cat(paste("ERROR:  ",err)) })
  }
  
  TimeOutput(start)
  
  
  
  
  
  
  
  
  
  
  
  
  #delete distance matrices with all 0
  distmfile = list.files(dist_dir, recursive=F, full.names=T, pattern=".Rdata$")
  a = foreach (i = 1:length(distmfile),.combine="c") %dopar% {
    a = F
    d = get(load(distmfile[i]))
    if (any(is.na(as.matrix(d)))) { rm(d); return(T) }
    if (all(as.matrix(d)==as.matrix(d)[1,1])) { rm(d); return(T) }
    rm(d)
    return(a)
  }
  file.remove(distmfile[a])
  
  
}







