## input: bi/clusters
## output: table of scores + bar plots based on how well clustering did according to known classes
## note: classes "normal" and "control" are treated the same with not phenodeviant features

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

options(stringsAsFactors=FALSE)
#options(device="cairo")
options(na.rm=T)

source("source/_funcAlice.R")
source("source/_funcdist.R")
libr(c("stringr", "plyr", "Matrix",
       "foreach","doMC",
       "clues", "PerfMeas", "cluster", "lattice")) #if there are date variables

#Setup Cores
no_cores = 5#detectCores()-3
registerDoMC(no_cores)

overwrite = F #redo and overwrite all past scores

id_col = "id"
target_cols = c("class","gender","group") #meta_file columns to plot
split_cols = c("gender", "group","none")


for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)[16:17]) {
  # result_dir = paste0(root, "/result/impc_panel1_sanger-spleen") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  ## input directories
  meta_dir = paste0(result_dir,"/meta")
  meta_file_dir = paste(meta_dir, "/file", sep="")
  clust_dir = paste(result_dir,  "/clust", sep="")
  dist_dir = paste(result_dir, "/dist", sep="")
  
  ## output directories
  score_dir = paste(result_dir, "/score_clust", sep=""); dir.create(score_dir, showWarnings=F)
  
  
  
  start = Sys.time()
  
  
  #dist matrix paths
  clust_types = gsub(".Rdata","", list.files(clust_dir, recursive=T, full.names=F, pattern=".Rdata"))
  # clust_types = clust_types[grepl("rw",clust_types)]
  
  # read meta_file and count matrix (count matrix only to trim feature matrix of small cell populations -- only for matrices that have cell populations as column names)
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
  # meta_train = read.csv(meta_train_dir)
  
  ## for each feature
  ncasl = llply(clust_types, function(clust_type) { 
    # ncasl = lapply(clust_types, function(clust_type) { 
    # for (clust_type in clust_types) {
    tryCatch({ cat("\n", clust_type, " ",sep="")
      clust_typedir0 = str_split(clust_type,"/")
      clust_typedir = paste0(clust_typedir0[[1]][1:(length(clust_typedir0[[1]])-1)], collapse="/")
      clust_typedir = paste0(score_dir,"/",clust_typedir)
      dir.create(clust_typedir, recursive=T, showWarnings=F)
      
      fname1 = str_split(clust_type, "[/]")[[1]]
      d = NULL
      if (grepl("dist",fname1[1]))
        d = as.matrix(get(load(paste0(dist_dir,"/",fname1[length(fname1)],".Rdata"))))
        
      ## upload and prep clust
      c00 = get(load(paste0(clust_dir,"/", clust_type,".Rdata")))
      
      # for each cluster result
      scores = list()
      for (i in 1:length(c00)) {
        clust_name = names(c00[i])
        c0 = c00[[i]]$x
        sm0 = meta_file[match(as.character(names(c0)), as.character(meta_file[,id_col])),]
        
        for (split_col in split_cols) {
          if (split_col=="none") {
            cl = list(none=c0)
            sml = list(none=sm0)
          } else if (!split_col%in%colnames(sm0)) { next
          } else {
            inds = llply(unique(sm0[,split_col]), function(x) sm0[,split_col]==x)
            cl = llply(inds, function(xi) c0[xi,xi])
            sml = llply(inds, function(xi) sm0[xi,xi])
          }
          
          # scores[[split_col]] = list()
          for (split_i in names(cl)) {
            c = cl[[split_i]]
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
              
              # dcname1 = paste0(score_dir1, "/", clust_type)
              # #to do or not to do
              # if (file.exists(dcname1) & !overwrite) next
              
              list_temp = list()
              
              
              
              
              
              
              
              # start -----------------------------------
              
              
              # ## make an adjusted version of clustering, such that they match with actual labels
              # clt = c
              # # for (j in 1:ncol(clt)) {
              # cl = clt
              # la = la0 = class #for all metrics
              # 
              # # cl1 s.t. each cluster is labeled by what majority of its real contents; if label is 
              # cl1 = rep(NA,length(la))
              # if (T | length(unique(la))>2) {
              #   tubes0 = unique(cl)[order(table(cl))]
              #   for (tubei in tubes0) {
              #     tci = which(cl==tubei) #index of cluster tubei
              #     tubej = Mode(la[tci]) #tubei is label taking up majority of cluster tubei
              #     cl1[tci] = tubej ## MAJORITY IN 2+ classes?
              #   }
              # } else if (F & length(unique(la))==2) {
              #   if (sum(la==unique(la)[1])>sum(la==unique(la)[2])) {
              #     tci = which(la==unique(la)[1])
              #     cl1[cl[tci]==Mode(cl[tci])] = unique(la)[1]
              #     cl1[cl[tci]!=Mode(cl[tci])] = unique(la)[2]
              #   } else {
              #     tci = which(la==unique(la)[2])
              #     cl1[cl[tci]==Mode(cl[tci])] = unique(la)[2]
              #     cl1[cl[tci]!=Mode(cl[tci])] = unique(la)[1]
              #   }
              # }
              # 
              # # cl such that class matches up with cl1; prioritize clusters to larger size labels
              # tubes = unique(cl1)[order(table(cl1),decreasing=F)] #don't switch these again!
              # for (k in 1:length(tubes)) {
              #   tubek = tubes[k]
              #   tck = which(cl1==tubek)
              #   tubei = Mode(cl[tck])
              #   stop = F
              #   if (k>1) {
              #     tck = tck[cl[tck]!=tubei & cl[tck]>k]
              #     while (length(tck)>0 & tubei%in%tubes[1:(k-1)] & !stop) {
              #       tubei = Mode(cl[tck])
              #       tck = tck[cl[tck]!=tubei]
              #       if (!length(tck)>0) stop = T
              #     }
              #   }
              #   if (tubei!=tubek & !stop) {
              #     clk = which(cl==tubek)
              #     cli = which(cl==tubei)
              #     if (length(clk)>0) cl[clk] = tubei
              #     cl[cli] = tubek
              #   }
              # }
              # 
              # c_a = cl1
              # # done!
              
              c_df = cluster_v2m(c)
              # c1_df = model.matrix(~ factor(c1) - 1); 
              # colnames(c1_df) = sort(unique(c1)); 
              # rownames(c1_df) = names(c1)
              
              class_df = cluster_v2m(class)
              # colnames(class_df) = sort(class_unique)
              rownames(class_df) = rownames(c_df)
              
              
              
              
              
              ## external validation F1 (classification & clustering)
              # if (cmethod%in%cmethodclass) {
              #   F1 = F.measure.single.over.classes(class_df,c_df)$average[-6]
              #   names(F1) = paste0(names(F1),"_1")
              #   r1 = adjustedRand(class,c); 
              #   names(r1) = paste0(names(r1), "_1")
              #   score = append(score, F1)
              #   score = append(score, r1)
              # } else {
              f11c = f.measure.comembership(class,c); 
              names(f11c) = paste0(names(f11c), "_co_1")
              list_temp = append(list_temp, as.list(f11c))
              # score = c(score, f11c, r1)
              # }
              
              
              ## external validation f1 (clustering)
              # if (!cmethod%in%cmethodclass) {
              #   f1c = f.measure.comembership(class_df,c_df); names(f1c) = paste0(names(f1c), "_co")
              #   # r = adjustedRand(class,c)
              #   # score = c(score,f1c,r)
              #   score = append(score,unlist(f1c))
              # }
              
              
              
              # ## internal validation (adjusted clustering)
              # if (length(unique(c1))==1) { sil = NA
              # } else { sil = median(silhouette(c1,d[[dindname]])[,3]) }
              # fm[[colnam]][[dindname]][[cltype]][[par]]["silmed_1"] = sil
              # 
              # if (length(unique(c1))==1) {
              #   score = rep(NA,9)
              #   names(score) = c("average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")
              # } else {
              #   score = unlist(cluster.stats(as.dist(d[[dindname]]),c1))[c("average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")] 
              # }
              # names(score) = paste0(names(score),"_1")
              # fm[[colnam]][[dindname]][[cltype]][[par]] = append(fm[[colnam]][[dindname]][[cltype]][[par]], score)
              
              
              
              
              
              
              
              if (!is.null(d)) {
              
              ## internal validation silmed (distance & clustering)
              # if (!cmethod%in%cmethodclass & dist_type!="NA") {
              # if (!"silmed"%in%names(fm[[colnam]][[dindname]][[cltype]][[par]]) | overwritef) {
              # if (length(unique(c))==1) { 
              #   sil = NA
              # } else { 
              sil = median(silhouette(c,d)[,3])
              names(sil) = "silmed"
              # }
              # fm[[colnam]][[dindname]][[cltype]][[par]]["silmed"] = sil
              list_temp = append(list_temp, sil)
              # }
              }
              
              # if (length(unique(cl))==1) {
              #   score = rep(NA,9)
              #   names(score) = c("average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")
              # } else {
              #   score = unlist(cluster.stats(as.dist(d[[dindname]]),cl))[c("average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")]
              # }
              # fm[[colnam]][[dindname]][[cltype]][[par]] = append(fm[[colnam]][[dindname]][[cltype]][[par]], score)
              
              # end ------------------------------------
              
              
              scores[[clust_name]][[split_col]][[paste0(split_i,"-",nrow(sm))]][[dirname]] = list_temp
              
              
              # plot ---------------------------
              cu = unique(class)
              stackedbar = sapply(unique(c), function(x) {
                sapply(cu, function(y) sum(class[c==x]==y)) })
              colnames(stackedbar) = unique(c)
              rownames(stackedbar) = cu
              png(paste0(clust_typedir, "/", clust_typedir0[[1]][length(clust_typedir0[[1]])], "_split-",split_col,"-",split_i,"_class-",target_col, ".png"))
              barplot(stackedbar, main="ground truth class distribution in each cluster",
                      xlab="cluster",
                      legend = rownames(stackedbar)) 
              graphics.off()
            }# target col
          }# split_i
        }# split_col
      }
      return(scores)
      # }, error = function(e) {
      #   cat(paste("ERROR:  ",e)); return(T)
    })
  }, .parallel=T)
  names(ncasl) = paste0(clust_types,"\\")
  
  # save ??? UNFINISHED
  ncaslul = unlist(ncasl)
  temp = str_split(names(ncaslul),"\\\\")
  clusts = sapply(temp, function(x) x[1])
  other = sapply(temp, function(x) gsub("^[.]","",paste0(x[-1], collapse="/")))
  ncasna0 = str_split(other,"[.]")
  splitcols_ind = sapply(ncasna0, function(x) {
    for (i in 1:(length(x)-1)) {
      if (grepl(paste0("^",x[i]),x[i+1]) & (x[i]%in%colnames(meta_file) | x[i]=="none")) return(i)
    }
  })
  clustmeth = sapply(1:length(ncasna0), function(xi) paste0(ncasna0[[xi]][1:(splitcols_ind[xi]-1)], collapse="."))
  ncasna = llply(1:length(ncasna0), function(xi) ncasna0[[xi]][splitcols_ind[xi]:length(ncasna0[[xi]])])

  splitcols = sapply(ncasna, function(x) x[1])
  splitis = sapply(ncasna, function(x) x[2])
  targcols = sapply(ncasna, function(x) paste0(x[3:(length(x)-1)], collapse="; "))
  targcols = sapply(targcols, function(x) paste0(x, collapse=""))
  targcols_col = sapply(strsplit(targcols,"[:|;]"), function(x) x[1])
  score_types = sapply(ncasna, function(x) x[5])
  score_table = data.frame(path=clusts, clustMethod.parameter=clustmeth, split_by_col=splitcols, split_variable_nosamples=splitis, NoOfClassesOrClasses=targcols, scoreType=score_types, score=ncaslul)
  if(file.exists(paste0(score_dir, ".csv"))) file.remove(paste0(score_dir, ".csv"))
  write.csv(score_table, paste0(score_dir, ".csv"), row.names=F)
  
  feat = sapply(str_split(clusts, "_layer"), function(x) x[1])
  feat = sapply(str_split(feat,"[/]"), function(x) x[length(x)])
  layer = as.numeric(gsub("layer-","",str_extract(clusts, "layer-[0-9]+")))
  dist = gsub("dist-","",str_extract(clusts, "dist-[a-zA-Z]+"))
  
  # plot
  dir.create(score_dir, showWarnings=F)
  targcols_cols = sapply(str_split(targcols_col,";"), function(x) x[1])
  for (targ in unique(targcols_cols)) {
    ti = targcols_cols==targ
    for (split in unique(splitcols)) {
      si = splitcols==split
      for (lay in unique(layer)) {
        li = layer==lay
        for (dis in unique(dist)) {
          di = dist==dis; di[is.na(di)] = F
          if (is.na(dis)) di = is.na(dist)
          # for (clus in unique(clustmeth)) {
          #   cli = clustmeth==clus
            for (sco in unique(score_types)) {
              sci = score_types==sco
              
              indsfull = ti&si&li&di&sci#&cli
              if (sum(indsfull)<2) next
              cat("\n", targ, ", ", split, ", ", lay, ", ", dis, ", ", sco)
              score_temp = score_table[indsfull,]
              score_temp$feat = feat[indsfull]
              
              png(paste0(score_dir, "/splitby-", split, "_class-", targ, "_layer-", lay, "_scoretype-",sco, ".png"), width=1200)
              pl = barchart(score~feat, groups=clustMethod.parameter, data=score_temp, auto.key = list(columns=2),
                       cex.axis=3, scales=list(x=list(rot=90,cex=0.8)), 
                       main=paste0(sco," scores by feature type for class ", targ, " grouped by cluster method (made from dist? ",dis,")"))
              print(pl)
              graphics.off()
            }
          # }
        }
      }
    }
  }
  
  time_output(start)
  
}
