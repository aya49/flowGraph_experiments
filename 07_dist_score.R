## input: distance matrices
## output: nca & dunn (didn't work well lol) table of scores + bar plots based on how well distance matrix aligns with known classes
## note: classes "normal" and "control" are treated the same with not phenodeviant features

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

options(stringsAsFactors=FALSE)
#options(device="cairo")
options(na.rm=T)

source("source/_func.R")
libr(c("stringr", "plyr",
       "foreach","doMC",
       "Brobdingnag","lubridate", "clValid", "lattice")) #if there are date variables

#Setup Cores
no_cores = 6#detectCores()-3
registerDoMC(no_cores)

doUnderflow = T #if numbers used to calculate scores --> 0, do work-around (slow)
nocontrol = F #delete all WT samples
adjustD = T #divide by a number if the numbers are too big e.g. manhattan
overwrite = F #redo and overwrite all past scores

id_col = "id"
target_cols = c("class","gender","group") #meta_file columns to plot
split_cols = c("gender", "group","none")



for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  # result_dir = paste0(root, "/result/impc_panel1_sanger-spleen") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  
  #Input
  meta_dir = paste0(result_dir,"/meta")
  meta_file_dir = paste(meta_dir, "/file", sep="")
  dist_dir = paste(result_dir, "/dist", sep="")
  
  #Output
  score_dir = paste(result_dir, "/score_dist", sep=""); dir.create(score_dir, showWarnings=F)
  # dist_score_result_dir = paste0(score_dir,"/score_nac_list.Rdata")
  
  
  
  
  
  
  start = Sys.time()
  
  
  #dist matrix paths
  dist_types = gsub(".Rdata","", list.files(dist_dir, recursive=F, full.names=F, pattern=".Rdata"))
  # dist_types = dist_types[grepl("rw",dist_types)]
  
  # read meta_file and count matrix (count matrix only to trim feature matrix of small cell populations -- only for matrices that have cell populations as column names)
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
  # meta_train = read.csv(meta_train_dir)
  
  ## for each feature
  scoresl = llply(dist_types, function(dist_type) {
    tryCatch({ cat("\n", dist_type, " ",sep="")
      start2 = Sys.time()
      
      ## upload and prep feature matrix
      d0 = as.matrix(get(load(paste0(dist_dir,"/", dist_type,".Rdata"))))
      sm0 = meta_file[match(rownames(d0),meta_file[,id_col]),]
      
      if (grepl("pregnancy",result_dir) & !grepl("TRIM",dist_type)) { try({
        mfri = ceiling(sqrt(length(unique(sm0$patient))))
        png(paste0(result_dir,"/", dist_type, "_intrapatient.png"), width=mfri*300, height=mfri*200)
        par(mfrow=c(mfri,mfri))
        cl = length(unique(sm0$class))
        for (pi in unique(sm0$patient)) {
          pind = sm0$patient==pi
          pind2 = pind & sm0$class=="control"
          if (sum(pind)<cl | sum(pind2)==0) next
          datas = d0[pind2,pind]
          plot(datas, xlab="class (weeks into pregnancy)", ylab="distance from control (1st week)", main=paste0("patient ", pi))
        }
        graphics.off()
      }) }
      
      scores = list()
      for (split_col in split_cols) {
        if (split_col=="none") {
          dl = list(none=d0)
          sml = list(none=sm0)
        } else if (!split_col%in%colnames(sm0)) { next
        } else {
          inds = llply(unique(sm0[,split_col]), function(x) sm0[,split_col]==x)
          dl = llply(inds, function(xi) d0[xi,xi])
          sml = llply(inds, function(xi) sm0[xi,])
        }
        
        # scores[[split_col]] = list()
        for (split_i in names(dl)) {
          d = dl[[split_i]]
          sm = sml[[split_i]]
          
          # scores[[split_col]][[paste0(split_i,".",nrow(sm))]] = list()
          for (target_col in target_cols) {
            if (!target_col%in%colnames(sm)) next
            
            #list out class labels
            class = sm[,target_col]; if (grepl("control",class) & grepl("normal",class)) class[grepl("normal",class)] = "control"
            class_unique = unique(class)
            # las0 = sm$label
            if (length(class_unique)<2 | any(table(class)<2)) next
            
            dirname = paste0(target_col, ".",
                             ifelse(length(class_unique)<10, 
                                    paste0(paste0(class_unique, "-",sapply(class_unique, function(x) sum(sm[,target_col]==x))), collapse="/"), length(class_unique)))
            # score_dir1 = paste0(score_dir,"/",dirname,"/splitby_",split_col,"-",split_i,".",nrow(sm)); dir.create(score_dir1, showWarnings=F, recursive=T)
            
            # dcname1 = paste0(score_dir1, "/", dist_type)
            # #to do or not to do
            # if (file.exists(dcname1) & !overwrite) next
            
            # NCA
            nca_score = NCA_score(d, class, doUnderflow=doUnderflow)$p
            
            # dunn
            dunn_score = dunn(as.dist(d), as.numeric(factor(class)))
            
            scores[[split_col]][[paste0(split_i,"-",nrow(sm))]][[dirname]] = list(nca=nca_score,dunn=dunn_score) #2min, 21s
            
            
            # plot -------------------
            
            
          }# target col
        }# split_i
      }# split_col
      return(scores)
      # }, error = function(e) {
      #   cat(paste("ERROR:  ",e)); return(T)
    })
  }, .parallel=T)
  names(scoresl) = paste0(dist_types,"\\")
  
  # save
  scoreslul = unlist(scoresl)
  temp = str_split(names(scoreslul),"\\\\")
  dists = sapply(temp, function(x) x[1])
  distss = sapply(str_split(dists,"-"), function(x) x[length(x)])
  other = sapply(temp, function(x) gsub("^[.]","",paste0(x[-1], collapse="/")))
  scoresna = str_split(other,"[.]")
  splitcols = sapply(scoresna, function(x) x[1])
  splitis = sapply(scoresna, function(x) x[2])
  targcols = sapply(scoresna, function(x) paste0(x[3:(length(x)-1)],": "))
  targcols = apply(targcols, 2, function(x) paste0(x, collapse=""))
  targcols_col = sapply(strsplit(targcols,"[:]"), function(x) x[1])
  score_types = sapply(scoresna, function(x) x[length(x)])
  score_table = data.frame(distance=distss, split_by_col=splitcols, split_variable_nosamples=splitis, NoOfClassesOrClasses=targcols, score_type=score_types, score=scoreslul)
  write.csv(score_table, paste0(score_dir, ".csv"))
  
  feat = sapply(str_split(dists, "_layer"), function(x) x[1])
  layer = as.numeric(gsub("layer-","",str_extract(dists, "layer-[0-9]+")))
  # dist = gsub("dist-","",str_extract(dists, "dist-[a-zA-Z]+"))
  
  # plot
  dir.create(score_dir, showWarnings=F)
  for (targ in unique(targcols_col)) {
    ti = targcols_col==targ
    for (split in unique(splitcols)) {
      si = splitcols==split
      for (lay in unique(layer)) {
        li = layer==lay
        # for (dis in unique(dist)) {
        #   di = dist==dis
        for (sco in unique(score_types)) {
          sci = score_types==sco
          
          score_temp = score_table[ti&si&li&sci,]
          score_temp$feat = feat[ti&si&li&sci]
          
          png(paste0(result_dir, "/splitby-", split, "_class-", targ, "_layer-", lay, "_scoretype-",sco, ".png"))
          pl = barchart(score~feat,data=score_temp[!grepl("paired",score_temp[,"feat"]),],groups=distance, 
                        scales=list(x=list(rot=90,cex=0.8)), main=paste0(sco," scores by feature type for class ", targ))
          print(pl)
          
          if (any(grepl("paired",score_temp[,"feat"]))) {
            graphics.off()
            png(paste0(result_dir, "/splitby-", split, "_class-", targ, "_layer-", lay, "_scoretype-",sco, "_paired.png"))
            pl = barchart(score~feat,data=score_temp[grepl("paired",score_temp[,"feat"]),],groups=distance, 
                          scales=list(x=list(rot=90,cex=0.8)), main=paste0(sco," PAIRED scores by feature type for class ", targ))
            print(pl)
          } 
          print(pl)
          graphics.off()
        }
        # }
      }
    }
  }
  
  time_output(start)
  
}
