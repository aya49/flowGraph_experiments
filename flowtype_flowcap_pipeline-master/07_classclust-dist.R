## Input: dist --> Output: clusters/classifications
# aya43@sfu.ca 20180526

## root directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

## input directories
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste(meta_dir, "/file", sep="")
meta_cell_dir = paste(meta_dir, "/cell", sep="")
feat_dir = paste(result_dir, "/feat", sep="")
meta_train_dir = paste0("attachments/AMLTraining.csv") #which FCM files are testing/training files

## output directories
dist_dir = paste(result_dir, "/dist", sep=""); dir.create(dist_dir, showWarnings=F)
clust_dir = paste(result_dir, "/clust", sep=""); dir.create(clust_dir, showWarnings=F)
clust_source_dir = paste0(clust_dir,"/clust_source"); dir.create(clust_source_dir, showWarnings=F)

## libraries
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
libr("FastKNN")
libr("cluster")
libr("mclust")
libr("kernlab")
libr("densitycut") #devtools::install_bitbucket("jerry00/densitycut_dev")
libr("foreach")
libr("doMC")
libr("stringr")

#Setup Cores
no_cores = 10#detectCores()-3
registerDoMC(no_cores)






#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)


















## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

readcsv = F

cellCountThres = c(1200)

readcsv = F #read features as csv or Rdata
overwrite = T #overwrite clustering?

good_count = 3 #trim matrix; only keep col/rows that meet criteria for more than 3 elements
good_sample = 3 #trim matrix; only keep rows that are a part of a class with more than 3 samples
cellCountThres = c(200) #a cell is insignificant if count under cell CountThres so delete -- only for matrices that have cell populations as column names

target_col = "aml" #the interested column in meta_file
control = "normal" #control value in target_col column
id_col = "fileName" #the column in meta_file matching rownames in feature matrices
order_cols = NULL #if matrix rows should be ordered by a certain column
split_col = "tube" # if certain rows in matrices should be analyzed in isolation, split matrix by this column in meta_file

cmethods = c("distmatrix","knn","kmed","lv","spec","spec1","hc") #rw1 #clustering methods
cmethodsclass = c("knn") # classification
# clustering/classification parameters
cmethodspar = list(knn=c(1:6), kmed=c(1,3,5), lv=c(0,.05,.1,.2), spec=list(methods=c("rbf"),kpar="automatic",tries=1), spec1=1, hc=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))



feat_count = "file-cell-countAdj" # cell count features used to trim matrix




#dist matrix paths
dist_types = list.files(dist_dir, recursive=F, full.names=F, pattern=".Rdata")
dist_types = gsub(".Rdata","",dist_types)

dist_types = dist_types[grepl("rw",dist_types)]















start = Sys.time()

# read meta_file and count matrix (count matrix only to trim feature matrix of small cell populations -- only for matrices that have cell populations as column names)
if (readcsv) {
  # mc = read.csv(paste0(feat_dir,"/", feat_count,".csv"),row.names=1, check.names=F)
  meta_file = read.csv(paste0(meta_file_dir,".csv"),check.names=F)
} else {
  # mc = get(load(paste0(feat_dir,"/", feat_count,".Rdata")))
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
}
meta_train = read.csv(meta_train_dir)
ts = strsplit(meta_file$fileName,"[A-Z]")
tu = sapply(ts, function(x) x[2])
sa = sapply(ts, function(x) x[3])
meta_file$label = sapply(1:nrow(meta_file), function(i) 
  return(meta_train[which(meta_train[,"TubeNumber"]==tu[i] & meta_train[,"SampleNumber"]==sa[i]),"Label"]))

## for each feature
a = foreach(dist_type=dist_types) %dopar% {
  tryCatch({
    cat("\n", dist_type, " ",sep="")
    start2 = Sys.time()
    
    ## upload and prep feature matrix
    if (readcsv) {
      d0 = as.matrix(read.csv(paste0(dist_dir,"/", dist_type,".csv"),row.names=1, check.names=F))
    } else {
      d0 = as.matrix(get(load(paste0(dist_dir,"/", dist_type,".Rdata"))))
    }
    if (!rownames(d0)[1]%in%meta_file[,id_col]) {
      cat("\nskipped: ",dist_type,", matrix rownames must match fileName column in meta_file","\n", sep="")
      return(T)
    }
    
    meta_file_ordered = meta_file[match(rownames(d0),meta_file[,id_col]),]
    
    ## split up analysis of feature matrix rows by split_col
    if (is.null(split_col)) {
      split_ind = list(all = 1:nrow(meta_file_ordered))
    } else {
      split_ids = unique(meta_file_ordered[,split_col])
      split_ids = split_ids[!is.na(split_ids)]
      split_ind = lapply(split_ids, function(split_id) which(meta_file_ordered[,split_col]==split_id) )
      names(split_ind) = split_ids
    }
    
    for (tube in names(split_ind)) {
      d = d0[split_ind[[tube]],split_ind[[tube]]]
      if (!sum(d!=0)>0) next()
      sm = meta_file_ordered[split_ind[[tube]],]
      if (length(unique(sm[,target_col]))<=1) next()
      if (length(unique(sm[,target_col]))<2) next()
      if (min(sapply(unique(sm[,target_col]), function(x) length(sm[,target_col]==x)))<1) next()
      
      #list out class labels
      la0 = sm[,target_col]; la0n = length(unique(la0))
      las0 = sm$label
      
      ## for each clustering method
      for (cmethod in cmethods) {
        pars = cmethodspar[[cmethod]]
        if (is.null(pars)) pars = NA
        
        for (par in pars) {
          
          # where to save clustering
          cname = paste0(clust_source_dir, "/", dist_type, "_clust-",cmethod,".",par,"_splitby-",split_col,".", tube, "_class-", target_col,"_rowclust.csv")
          
          #to do or not to do
          if (file.exists(cname) & !overwrite) next
          # if (!cmethod%in%cmethodsorigm & i2==2) next
          # if (cmethod%in%cmethodsorigm) { #for other distance, these are clusterings made directly from original features
          #   if (i2!=2) next
          #   if (!"clt"%in%names(cm0[[colnam]][[dindname]][[cmethod]])) next
          #   clt = cm0[[colnam]][[dindname]][[cmethod]]$clt
          # }
          
          # #create class/cluster list
          # cm0[[colnam]][[dindname]][[cmethod]] = list()
          
          cat(" ",cmethod," ",sep="")
          start2 = Sys.time()
          
          ## knn for each k (classification)
          if (cmethod=="knn") {
            labelss = la0
            labelss[is.na(las0)] = NA
            clt1 = as.numeric(factor(knntable(as.matrix(d),par,labelss)[,1]))
            clt = rep(0,length(la0))
            clt[is.na(las0)] = clt1
          } #parameter k
          
          ## kmedoids
          if (cmethod=="kmed") { 
            clt = pamtable(d,la0n,par)[,1]
          } #number of tries
          
          ## louvain
          if (cmethod=="lv") { #input is a similarity matrix
            # if (grepl("dist.Rdata",distpath)) {
            # sim0 = get(load(paste0(dist_dir,"/",gsub("dist.Rdata","simmatrix.Rdata",distpath))))
            # simind = match(rownames(as.matrix(d[[dindname]])),rownames(sim0))
            # sim = sim0[simind,simind]
            # } else { 
            sim = get_graph(d)
            # }
            clt = lvtable(sim,par)[,1]
          }
          
          ## random walk refined distance matrices only
          if (cmethod=="rw1") {
            # if (grepl("dist.Rdata",distpath)) {
            #   sim0 = get(load(paste0(dist_dir,"/",gsub("dist.Rdata","simmatrix.Rdata",distpath))))
            #   simind = match(rownames(as.matrix(d[[dindname]])),rownames(sim0))
            #   sim = sim0[simind,simind]
            # } else if (!rwonly) {
            sim = get_graph(d) 
            # } else { next }
            clt = rw1table(sim,par)[,1]
          }
          
          ## spectral clustering via distance matrix
          if (cmethod=="spec1") {
            # if (grepl("dist.Rdata",distpath)) {
            #   sim0 = get(load(paste0(dist_dir,"/",gsub("dist.Rdata","simmatrix.Rdata",distpath))))
            #   simind = match(rownames(as.matrix(d)),rownames(sim0))
            #   sim = sim0[simind,simind]
            # } else { 
            sim = get_graph(d) 
            # }
            clt = spec1table(sim,la0n)[,1]
          }
          
          ## hierarchical clustering
          if (cmethod=="hc") { 
            clt = hctable(d,la0n,par)[,1]
          }
          
          ## densitycut via distance matrix (k=3)
          if (cmethod=="dc1") { 
            clt = dc1table(d)[,1]
          }
          
            write.csv(clt,file=cname)

          
          
          TimeOutput(start2)
          
          
          
        }
        
        
      }
    }
    return(F)
  }, error = function(e) {
    cat(paste("ERROR:  ",e)); return(T)
  })
}

TimeOutput(start)









# #check which features weren't finished
# aa = distMetafun(list.files(clust_dir,pattern="other"),dis=c(dis,"other"))
# mattype = sapply(1:nrow(aa), function(x) {
#   if (aa$rand[x]==0) { add = ""
#   } else { add = paste0("_",aa$rand[x]) }
#   return(paste0(aa$type[x],add))
# })
# mattype2 = sapply(unique(mattype), function(x) {
#   al = aa$layer[which(mattype==x)]
#   if (grepl("TRIM",x) & length(al)>1) return(T)
#   if (!grepl("TRIM",x) & length(al)>2) return(T)
#   return(F)
# })
# # unique(mattype)[mattype2]
# feat_types[!feat_types%in%unique(mattype)[mattype2]] #not done
# 


