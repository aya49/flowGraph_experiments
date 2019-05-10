## Input: feature --> Output: clusters/rbfdistance
# aya43@sfu.ca 20180526

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metrics"
setwd(root)
for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  # result_dir = paste0(root, "/result/impc_panel1_sanger-spleen") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  
## input directories
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste(meta_dir, "/file", sep="")
meta_cell_dir = paste(meta_dir, "/cell", sep="")
feat_dir = paste(result_dir, "/feat", sep="")

## output directories
dist_dir = paste(result_dir, "/dist", sep=""); dir.create(dist_dir, showWarnings=F)
clust_dir = paste(result_dir, "/clust", sep=""); dir.create(clust_dir, showWarnings=F)
clust_source_dir = paste0(clust_dir,"/clust_source"); dir.create(clust_source_dir, showWarnings=F)

## libraries
source("source/_funcAlice.R")
source("source/_funcdist.R")
libr("FastKNN")
libr("cluster")
libr("mclust")
libr("kernlab")
libr("densitycut") #devtools::install_bitbucket("jerry00/densitycut_dev")
libr("foreach")
libr("doMC")
libr("stringr")

#Setup Cores
no_cores = 14#detectCores()-3
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

cmethods = c("spec","dc") #clustering methods
# densitycut parameters
alpha=.85; nu=seq(0.0, 1.0, by=0.05)




#feature matrix paths
if (readcsv) {
  feat_types = list.files(path=feat_dir,pattern=".csv")
  feat_types = gsub(".csv","",feat_types)
} else {
  feat_types = list.files(path=feat_dir,pattern=".Rdata")
  feat_types = gsub(".Rdata","",feat_types)
}
feat_types = feat_types[!grepl("file-cell-count.",feat_types)]
feat_types = feat_types[!grepl("-fe-",feat_types)]
feat_types = feat_types[!grepl("KO|Max",feat_types)]

feat_count = "file-cell-countAdj" # cell count features used to trim matrix
#feat_types = list.files(path=result_dir,pattern=glob2rx("matrix*.Rdata"))














start = Sys.time()

# read meta_file and count matrix (count matrix only to trim feature matrix of small cell populations -- only for matrices that have cell populations as column names)
if (readcsv) {
  mc = read.csv(paste0(feat_dir,"/", feat_count,".csv"),row.names=1, check.names=F)
  meta_file = read.csv(paste0(meta_file_dir,".csv"),check.names=F)
} else {
  mc = get(load(paste0(feat_dir,"/", feat_count,".Rdata")))
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
}

## for each feature
a = foreach(feat_type=feat_types) %dopar% {
  tryCatch({
    cat("\n", feat_type, " ",sep="")
    start2 = Sys.time()
    
    ## upload and prep feature matrix
    if (readcsv) {
      m0 = as.matrix(read.csv(paste0(feat_dir,"/", feat_type,".csv"),row.names=1, check.names=F))
    } else {
      m0 = as.matrix(get(load(paste0(feat_dir,"/", feat_type,".Rdata"))))
    }
    if (!rownames(m0)[1]%in%meta_file[,id_col]) {
      cat("\nskipped: ",feat_type,", matrix rownames must match fileName column in meta_file","\n", sep="")
      return(T)
    }
    
    ## does feature matrix have cell populations on column names?
    layers = 0
    countThres = 0
    colhaslayer = grepl("_layer-",feat_type)
    colhascell = str_split(feat_type,"-")[[1]][2]=="cell"
    if (colhascell) {
      layers = c(1,2,4,max(unique(sapply(unlist(str_split(colnames(m0),"_")[[1]]), function(x) str_count(x,"[+-]")))))
      countThres = cellCountThres
    }
    
    ## for each layer, trim feature matrix accordingly
    for (k in layers) {
      #trim matrix
      mm = trimMatrix(m0,TRIM=T, mc=mc, sampleMeta=meta_file, sampleMeta_to_m1_col=id_col, target_col=target_col, control=control, order_cols=order_cols, colsplitlen=NULL, k=k, countThres=countThres, goodcount=good_count, good_sample=good_sample)
      if (is.null(mm)) next
      m_ordered = mm$m
      meta_file_ordered = mm$sm
      if (all(meta_file_ordered[,target_col]==meta_file_ordered[1,target_col])) next
      
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
        m = m_ordered[split_ind[[tube]],]
        if (!sum(m!=0)>0) next
        sm = meta_file_ordered_split = meta_file_ordered[split_ind[[tube]],]
        if (length(unique(sm[,target_col]))<=1) next
        
        #ensure every matrices have values and more than one class
        # for more than one column
        # if (length(Reduce('intersect',lapply(sm, function(x) unique(x[,target_col]))))<length(unique(sm[,target_col]))) next
        # if (length(Reduce('intersect',lapply(sm, function(x) x[duplicated(x[,target_col]),target_col])))<length(unique(sm[,target_col]))) next  #if not both aml and healthy in each tube
        
        if (length(unique(sm[,target_col]))<2) next()
        if (min(sapply(unique(sm[,target_col]), function(x) length(sm[,target_col]==x)))<1) next()
        
        #list out class labels
        la0 = sm[,target_col]
        
        ## for each clustering method
        for (cmethod in cmethods) {
          
          # where to save clustering
          if (colhaslayer) {
            dname0 = paste0("/",feat_type, "_dist-NADIST", sep = "")
          } else {
            dname0 = paste("/",feat_type, "_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_dist-NADIST")
          }
          dname = paste0(dist_dir,dname0)
          
          cname = paste0(clust_source_dir, dname0, "_clust-",cmethod,"_splitby-",split_col,".", tube, "_class-", target_col,"_rowclust.csv")
          if (file.exists(cname) & !overwrite) next()
          
          ## desnitycut clustering
          if (cmethod=="dc") {
            clt = DensityCut(X=m, alpha=alpha, nu=nu, show.plot=F)$cluster
          }
          
          ## spectral clustering and rbf kernel parameter learning
          if (cmethod=="spec") {
            dname1 = gsub("NADIST","rbf",dname)
            
            sp = NULL
            tr = try({ sp = specc(x=m, kernel="rbfdot",centers=length(unique(la0))) })
            if ("try-error" %in% class(tr))
              sp = specc(x=(m-min(m)) / (max(m)-min(m)), kernel="rbfdot",centers=length(unique(la0)))
            
            if (is.null(sp)) next()
            clt = sp@.Data

            ## calculate rbf distance
            parlist0 = unlist(sp@kernelf@kpar)
            simrbf = kernelMatrix(kernel=rbfdot(as.numeric(parlist0)), x=m)
            drbf = get_graphd(simrbf)
            save(drbf, file=paste0(dname1,".Rdata"))
            
          }
          names(clt) = rownames(m)
          write.csv(clt, file=cname)
          
        }
        
        
      }
    }
    return(T)
  }, error = function(e) {
    cat(paste("ERROR:  ",e)); graphics.off(); return(T)
  })
}

time_output(start)









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

}
