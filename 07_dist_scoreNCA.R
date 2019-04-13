# Evaluate all distance metrics
# aya43@sfu.ca 20170116

## root directory
root = "~/projects/flowtype_metrics"
setwd(root)
for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  # result_dir = paste0(root, "/result/impc_panel1_sanger-spleen") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  
  options(stringsAsFactors=FALSE)
  #options(device="cairo")
  options(na.rm=T)
  
  #Input
  meta_dir = paste0(result_dir,"/meta")
  meta_file_dir = paste(meta_dir, "/file", sep="")
  dist_dir = paste(result_dir, "/dist", sep="")
  
  #Output
  dist_score_dir = paste(dist_dir, "/dist_score", sep="")
  # dist_score_result_dir = paste0(dist_score_dir,"/score_nac_list.Rdata")
  
  source("source/_funcAlice.R")
  source("source/_funcdist.R")
  libr(c("stringr", "plyr",
         "foreach","doMC",
         "Brobdingnag","lubridate")) #if there are date variables
  
  doUnderflow = T #if numbers used to calculate scores --> 0, do work-around (slow)
  nocontrol = F #delete all WT samples
  adjustD = T #divide by a number if the numbers are too big e.g. manhattan
  overwrite = F #redo and overwrite all past scores
  
  id_col = "id"
  target_cols = c("class","gender","group") #meta_file columns to plot
  split_cols = c("gender", "group","none")
  
  
  
  
  
  
  
  start = Sys.time()
  
  
  #dist matrix paths
  dist_types = gsub(".Rdata","", list.files(dist_dir, recursive=F, full.names=F, pattern=".Rdata"))
  # dist_types = dist_types[grepl("rw",dist_types)]
  
  # read meta_file and count matrix (count matrix only to trim feature matrix of small cell populations -- only for matrices that have cell populations as column names)
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
  # meta_train = read.csv(meta_train_dir)
  
  ## for each feature
  ncasl = llply(dist_types, function(dist_type) { 
    tryCatch({ cat("\n", dist_type, " ",sep="")
      start2 = Sys.time()
      
      ## upload and prep feature matrix
      d0 = as.matrix(get(load(paste0(dist_dir,"/", dist_type,".Rdata"))))
      sm0 = meta_file[match(rownames(d0),meta_file[,id_col]),]
      
      ncas = list()
      for (split_col in split_cols) {
        if (split_col=="none") {
          dl = list(none=d0)
          sml = list(none=sm0)
        } else if (!split_col%in%colnames(meta_file)) { next
        } else {
          inds = llply(unique(meta_file[,split_col]), function(x) meta_file[,split_col]==x)
          dl = llply(inds, function(xi) d0[xi,xi])
          sml = llply(inds, function(xi) sm0[xi,xi])
        }
        
        ncas[[split_col]] = list()
        for (split_i in names(dl)) {
          d = dl[[split_i]]
          sm = sml[[split_i]]
          
          ncas[[split_col]][[paste0(split_i,".",nrow(sm))]] = list()
          for (target_col in target_cols) {
            if (!target_col%in%colnames(sm)) next
            
            #list out class labels
            class = sm[,target_col]; 
            class_unique = unique(class)
            # las0 = sm$label
            if (length(class_unique)<2 | any(table(class)<2)) next
            
            dirname = paste0(target_col, ".",
                             ifelse(length(class_unique)<10, 
                                    paste0(paste0(class_unique, "-",sapply(class_unique, function(x) sum(sm[,target_col]==x))), collapse="/"), length(class_unique)))
            # dist_score_dir1 = paste0(dist_score_dir,"/",dirname,"/splitby_",split_col,"-",split_i,".",nrow(sm)); dir.create(dist_score_dir1, showWarnings=F, recursive=T)
            
            # dcname1 = paste0(dist_score_dir1, "/", dist_type)
            # #to do or not to do
            # if (file.exists(dcname1) & !overwrite) next
            
            ncas[[split_col]][[paste0(split_i,"-",nrow(sm))]][[dirname]] = NCA_score(d, class, doUnderflow=doUnderflow)$p #2min, 21s
          }# target col
        }# split_i
      }# split_col
      return(ncas)
      # }, error = function(e) {
      #   cat(paste("ERROR:  ",e)); return(T)
    })
  }, .parallel=T)
  names(ncasl) = dist_types
  
  # save
  ncaslul = unlist(ncasl)
  ncasna = str_split(names(ncaslul),"[.]")
  dists = sapply(ncasna, function(x) paste0(x[1:(length(x)-4)],sep="."))
  splitcols = sapply(ncasna, function(x) x[length(x)-3])
  splitis = sapply(ncasna, function(x) x[length(x)-2])
  targcols = sapply(ncasna, function(x) paste0(x[(length(x)-1):length(x)],": "))
  scores = data.frame(distance=dists, split_by_col=splitcols, split_variable_nosamples=splitis, class=targcols, NCA_score=scores)
  write.csv(scores, paste0(dist_score_dir, ".csv"))
  
  time_output(start)
  
}
