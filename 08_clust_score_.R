## Input: bi/cluster --> Output: bi/cluster scores for samples
# aya43@sfu.ca 20180328

## root directory
root = "~/projects/flowtype_metrics"
setwd(root)
for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  # result_dir = paste0(root, "/result/impc_panel1_sanger-spleen") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  
## input directories
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste(meta_dir, "/file", sep="")
clust_dir = paste(result_dir,  "/clust", sep="")

## output directories
score_dir = paste(result_dir, "/score_clust", sep="")


## libraries
source("source/_funcAlice.R")
source("source/_funcdist.R")
libr(c("biclust","NMF","pheatmap",
       "clues",
       "foreach","doMC",
       "stringr"))

#Setup Cores
no_cores = 10#detectCores()-3
registerDoMC(no_cores)








## options for script
readcsv = F
overwrite = T
target_cols = c("aml")

control = "normal" #control value in target_col column for each centre
id_col = "fileName"
target_col = "aml"
split_col = "tube"

cmethodclass = "knn"

# bcmethods = c("plaid","CC","bimax","BB-binary","nmf-nsNMF","nmf-lee","nmf-brunet","CC")
# onlysigBB = T #only evaluate significant BB-binary clusters
# Kr = 20; Kc = 20 #number of row and column biclusters for binary bayesian biclustering
# pval_thres = .01 #pval_threshold for choosing biclusters out of all bayesian biclusters
# min_iter = 100 #min number of iterations for BB-binary (B2PS)
# nmf_thres = .05 # * max contribution: threshold at which a row/col can be considered a significant contribution to a factor in nmf

#data paths
# clust_paths = list.files(path=clust_source_dir,pattern="_rowclust.csv", full.names=T)
# list.files(path=clust_source_dir,pattern=".Rdata", full.names=T))


clust_paths = append(list.files(path=clust_source_dir,pattern="_rowclust.csv", full.names=T),
                     list.files(path=biclust_source_dir,pattern="_rowclust.csv", full.names=T) )


feat_count = "file-cell-countAdj"



#dist matrix paths
dist_types = list.files(dist_dir, recursive=F, full.names=F, pattern=".Rdata")
dist_types = gsub(".Rdata","",dist_types)


dist_paths = lapply(dist_types, function(dist_type) {
  ca = clust_paths[grepl(paste0("/",dist_type,"_"),clust_paths)]
  if (length(ca)>0) return(ca)
  return("NA")
})
names(dist_paths) = dist_types
clust_paths_rest = clust_paths[!clust_paths%in%unlist(dist_paths)]
if (length(clust_paths_rest)>0) dist_paths[["NA"]] = clust_paths_rest











start = Sys.time()


if (readcsv) {
  # mc = read.csv(paste0(feat_dir,"/", feat_count,".csv"),row.names=1, check.names=F)
  meta_file = read.csv(paste0(meta_file_dir,".csv"),check.names=F)
} else {
  # mc = get(load(paste0(feat_dir,"/", feat_count,".Rdata")))
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
}

score_list0 = foreach(dist_type=names(dist_paths)) %dopar% {
  tryCatch ({
    scores = list()
    
    ## upload and prep feature matrix
    d0 = NULL
    if (dist_type!="NA" & (overwrite | !file.exists(paste0(clust_score_dir, "/", dist_type, ".Rdata")))) {
      if (readcsv) {
        d0 = as.matrix(read.csv(paste0(dist_dir,"/", dist_type,".csv"),row.names=1, check.names=F))
      } else {
        d0 = as.matrix(get(load(paste0(dist_dir,"/", dist_type,".Rdata"))))
      }
      
      if (!rownames(d0)[1]%in%meta_file[,id_col]) {
        cat("\nskipped: ",dist_type,", matrix rownames must match fileName column in meta_file","\n", sep="")
      } else {
        score = NULL
        
        rowlabel = meta_file[match(rownames(d0),meta_file[,id_col]),target_col]
        names(rowlabel) = rownames(d0)
        
        rowlabel_df = model.matrix(~ factor(rowlabel) - 1)
        colnames(rowlabel_df) = laname = sort(unique(rowlabel))
        rownames(rowlabel_df) = names(rowlabel)
        
        
        ## internal validation NCA (distance)
        # if (dist_type!="NA") {
        # cl = la0
        # if (!"NCA"%in%names(fm[[colnam]][[dindname]][[cltype]][[par]]) | overwritef) {
        # if (length(unique(cl))==1) { nca = NA
        # } else { 
        nca = NCA_score(as.matrix(d0), rowlabel)$p
        names(nca) = "nca"
        score = append(score, nca)
        # }
        # fm[[colnam]][[dindname]][[cltype]][[par]]["NCA"] = nca
        # }
        # }
        score = unlist(score)
        scores[[dist_type]] = score
        
        save(score,file=paste0(clust_score_dir, "/", dist_type, ".Rdata"))
      }
      
      
    }
    
    if (dist_paths[[dist_type]]=="NA") return(scores)
    
    for (clust_path in dist_paths[[dist_type]]) {
      clust_path_ = gsub("_rowclust.csv","",fileNames(clust_path))
      if (!overwrite & file.exists(paste0(clust_score_dir, "/", clust_path_, ".Rdata"))) next()
      
      cat("\n", clust_path, " ",sep="")
      start2 = Sys.time()
      
      cmethod = sub(".*?clust-(.*?)_.*", "\\1", clust_path)
      score = NULL
      
      
      
      
      
      
      
      ## prep clusters / label
      rowclust0 = read.csv(clust_path, row.names=1)
      rowclust = rowclust0[,1]
      names(rowclust) = rownames(rowclust0)
      if (cmethod%in%cmethodclass) {
        rowclust = rowclust[rowclust!=0]
      }
      rowlabel = meta_file[match(names(rowclust),meta_file[,id_col]),target_col]
      names(rowlabel) = names(rowclust)
      
      
      # ## make an adjusted version of clustering, such that they match with actual labels
      # clt = rowclust
      # # for (j in 1:ncol(clt)) {
      # cl = clt
      # la = la0 = rowlabel #for all metrics
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
      # rowclust_a = cl1
      # # done!
      
      
      
      
      
      if (length(unique(rowclust[rowclust>0]))<2) next()
      
      rowclust_df = matrix(rep(1,length(rowclust)),ncol=1)
      rownames(rowclust_df) = names(rowclust)
      colnames(rowclust_df) = rowclust[1]
      
      # rowclust1_df = model.matrix(~ factor(rowclust1) - 1); 
      # colnames(rowclust1_df) = sort(unique(rowclust1)); 
      # rownames(rowclust1_df) = names(rowclust1)
      
      rowlabel_df = model.matrix(~rowlabel)
      colnames(rowlabel_df) = laname = sort(unique(rowlabel))
      rownames(rowlabel_df) = names(rowlabel)
      
      
      
      
      d = NULL
      if (dist_type!="NA") {
        rowind = match(names(rowclust),rownames(d0))
        d = d0[rowind,rowind]
      }
      
      
      
      
      
      
      
      
      
      
      ## external validation F1 (classification & clustering)
      if (cmethod%in%cmethodclass) {
        F1 = F.measure.single.over.classes(rowlabel_df,rowclust_df)$average[-6]
        names(F1) = paste0(names(F1),"_1")
        r1 = adjustedRand(rowlabel,rowclust); 
        names(r1) = paste0(names(r1), "_1")
        score = append(score, F1)
        score = append(score, r1)
      } else {
        f11c = f.measure.comembership(rowlabel,rowclust); 
        names(f11c) = paste0(names(f11c), "_co_1")
        score = append(score, f11c)
        # score = c(score, f11c, r1)
      }
      
      
      ## external validation f1 (clustering)
      # if (!cmethod%in%cmethodclass) {
      #   f1c = f.measure.comembership(rowlabel_df,rowclust_df); names(f1c) = paste0(names(f1c), "_co")
      #   # r = adjustedRand(rowlabel,rowclust)
      #   # score = c(score,f1c,r)
      #   score = append(score,unlist(f1c))
      # }
      
      
      
      # ## internal validation (adjusted clustering)
      # if (length(unique(rowclust1))==1) { sil = NA
      # } else { sil = median(silhouette(rowclust1,d[[dindname]])[,3]) }
      # fm[[colnam]][[dindname]][[cltype]][[par]]["silmed_1"] = sil
      # 
      # if (length(unique(rowclust1))==1) {
      #   score = rep(NA,9)
      #   names(score) = c("average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")
      # } else {
      #   score = unlist(cluster.stats(as.dist(d[[dindname]]),rowclust1))[c("average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")] 
      # }
      # names(score) = paste0(names(score),"_1")
      # fm[[colnam]][[dindname]][[cltype]][[par]] = append(fm[[colnam]][[dindname]][[cltype]][[par]], score)
      
      
      
      
      
      
      
      
      
      ## internal validation silmed (distance & clustering)
      if (!cmethod%in%cmethodclass & dist_type!="NA") {
        # if (!"silmed"%in%names(fm[[colnam]][[dindname]][[cltype]][[par]]) | overwritef) {
        # if (length(unique(rowclust))==1) { 
        #   sil = NA
        # } else { 
        sil = median(silhouette(rowclust,d)[,3])
        names(sil) = "silmed"
        # }
        # fm[[colnam]][[dindname]][[cltype]][[par]]["silmed"] = sil
        score = append(score, sil)
      }
      
      # if (length(unique(cl))==1) {
      #   score = rep(NA,9)
      #   names(score) = c("average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")
      # } else {
      #   score = unlist(cluster.stats(as.dist(d[[dindname]]),cl))[c("average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")]
      # }
      # fm[[colnam]][[dindname]][[cltype]][[par]] = append(fm[[colnam]][[dindname]][[cltype]][[par]], score)
      
      
      #get dist matrix between genes
      if (grepl("IMPC",root)) {
        rowgene = meta_file[match(names(rowclust),meta_file[,id_col]),gene_col]
        names(rowgene) = names(rowclust)
        for (netdisti in names(netdists)) {
          netdist = netdists[[netdisti]]
          netgenes = colnames(netdist)
          evalgenes = intersect(rowgene, netgenes)
          
          rowgene1 = rowgene[rowgene%in%evalgenes]
          rowclust1 = rowclust[match(names(rowgene1),names(rowgene))]
          rowclust1 = rowclust1[rowclust1>0]
          if (length(rowclust1)==0) next
          if (length(unique(rowclust1))==1) next
          rowgene1 = rowgene1[match(names(rowclust1),names(rowgene1))]
          
          rowclust1_df = model.matrix(~ factor(rowclust1) - 1); colnames(rowclust1_df) = sort(unique(rowclust1)); rownames(rowclust1_df) = names(rowclust1)
          netdist_rowcol = match((rowgene1),netgenes)
          netdist1 = netdist[netdist_rowcol,netdist_rowcol]
          netdist1 = apply(as.matrix(netdist1), 2, as.numeric)
          rownames(netdist1) = colnames(netdist1) = names(rowgene1)
          
          #metrics
          nca = NCA_score(netdist1, rowclust1)$p
          names(nca) = paste0("NCA_net-",netdisti,".g-",length(unique(rowgene1)))
          no_genes=length(evalgenes)
          score = append(score,no_genes,nca)
          
          silhouette = silhouette(rowclust1,netdist1)
          pngame = paste0(biclust_plot_dir,"/",fileNames(clust_path),"_silhouette.net-",netdisti,".g-",length(unique(rowgene1)),".png")
          png(pngame)
          plot(silhouette)
          graphics.off()
          sil = median(silhouette[,3])
          names(sil) = paste0("median-silhouette_net-",netdisti,".g-",length(unique(rowgene1)))
          score = append(score,sil)
        }
      }
      
      
      score = unlist(score)
      scores[[clust_path]] = score
      
      save(score,file=paste0(clust_score_dir, "/", clust_path_, ".Rdata"))
      time_output(start2)
      
    }
    return(scores)
  }, error = function(err) { cat(paste("error:  ",err)); return(NA) })
  
}

# score_list = unlist(score_list0, recursive=F)
# save(score_list, file = paste0(clust_score_dir,".Rdata"))

time_output(start)












score_files = list.files(clust_score_dir,full.names=T)
score_list = lapply(score_files, function(x) get(load(x)))
names(score_list) = score_files

## put scores into a table
error_ind = is.na(score_list)
score_names = unique(unlist(lapply(score_list, names)))
score_list = lapply(score_list[!error_ind], function(x) {
  y = x[match(score_names,names(x))]
  names(y) = score_names
  return(y)
})
score_table = Reduce("rbind",score_list)
colnames(score_table) = score_names

clust_files = fileNames(names(score_list))
clust_files = gsub("[.]csv|[.]Rdata","",clust_files)
clust_files = gsub("_rowclust","",clust_files)

clust_files_attr = str_split(clust_files,"_")
clust_files_list1 = lapply(clust_files_attr, function(x) {
  metav = append(x[1], sapply(x[2:(length(x)-1)], function(y) str_split(y,"-")[[1]][2]))
  names(metav) = append("feature", sapply(x[2:(length(x)-1)], function(y) str_split(y,"-")[[1]][1]))
  return(metav)
})
clust_files_names = Reduce("union",lapply(clust_files_list1, names))
clust_files_list = lapply(clust_files_list1, function(x) {
  y = x[match(clust_files_names,names(x))]
  names(y) = clust_files_names
  return(y)
})
clust_files_table = Reduce("rbind",clust_files_list)
colnames(clust_files_table) = clust_files_names

clust_files_table_final = cbind(clust_files_table, score_table)
write.csv(clust_files_table_final, file = paste0(clust_score_dir,".csv"))


time_output(start)






