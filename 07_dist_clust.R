## input: distance matrices
## output: clusters (if number of clusters/k required, the number of ground truth classes were provided)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)


## libraries
source("source/_funcAlice.R")
source("source/_funcdist.R")
libr(c("FastKNN","cluster","mclust","kernlab", "igraph", 
       "tsne",
       "densitycut", #devtools::install_bitbucket("jerry00/densitycut_dev")
       "foreach","doMC",
       "stringr", "plyr",
       "tcltk"))

#Setup Cores
no_cores = 10#detectCores()-3
setup_parallel(no_cores)




## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

readcsv = F


overwrite = T #overwrite clustering?

plotsize = 300

# countThres = 400 #a cell is insignificant if count under cell CountThres so delete -- only for matrices that have cell populations as column names
# good_count = 3 #trim matrix; only keep col/rows that meet criteria for more than 3 elements
# good_sample = 3 #trim matrix; only keep rows that are a part of a class with more than 3 samples
# feat_count = "file-cell-countAdj" # cell count features used to trim matrix, not used

target_cols = c("class","gender") #the interested column in meta_file
# control = "control" #control value in target_col column
id_col = "id" #the column in meta_file matching rownames in feature matrices
order_cols = NULL #if matrix rows should be ordered by a certain column
# split_col = NULL # if certain rows in matrices should be analyzed in isolation, split matrix by this column in meta_file

# cmethods = c("distmatrix","knn","kmed", "kmeans","lv","spec","spec1","hc") #rw1 #clustering methods
# cmethodsclass = c("knn") # classification
# # clustering/classification parameters
cmethodspar = list(#knn=c(1:6), 
  kmed=NA,
  kmeans=NA,
  lv=c(0,.05,.1,.2), 
  # spec=list(methods=c("rbf"), kpar="automatic",tries=1), 
  spec1=NA, 
  hc=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))


for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)[-15]) {
  # result_dir = paste0(root, "/result/impc_panel1_sanger-spleen") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta")
  meta_file_dir = paste(meta_dir, "/file", sep="")
  meta_cell_dir = paste(meta_dir, "/cell", sep="")
  feat_dir = paste(result_dir, "/feat", sep="")
  dist_dir = paste(result_dir, "/dist", sep="")
  # meta_train_dir = paste0("attachments/AMLTraining.csv") #which FCM files are testing/training files for flowCAP data set
  
  ## output directories
  clust_dir = paste0(result_dir,"/clust/dist")
  clust_plot_dir = paste(clust_dir, "/plot", sep=""); 
  dir.create(clust_plot_dir, showWarnings=F, recursive=T)
  
  
  start = Sys.time()
  
  
  #dist matrix paths
  dist_types = gsub(".Rdata","", list.files(dist_dir, recursive=F, full.names=F, pattern=".Rdata"))
  # dist_types = dist_types[grepl("rw",dist_types)]
  
  # read meta_file and count matrix (count matrix only to trim feature matrix of small cell populations -- only for matrices that have cell populations as column names)
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
  # meta_train = read.csv(meta_train_dir)
  
  a = llply(dist_types, function(dist_type) { 
    tryCatch({ cat("\n", dist_type, " ",sep="")
      start2 = Sys.time()
      
      ## upload and prep feature matrix
      d = as.matrix(get(load(paste0(dist_dir,"/", dist_type,".Rdata"))))
      sm = meta_file[match(rownames(d),meta_file[,id_col]),]
      
      # ## split up analysis of feature matrix rows by split_col
      # if (is.null(split_col)) {
      #   split_ind = list(all = 1:nrow(meta_file_ordered))
      # } else {
      #   split_ids = unique(meta_file_ordered[,split_col])
      #   split_ids = split_ids[!is.na(split_ids)]
      #   split_ind = lapply(split_ids, function(split_id) which(meta_file_ordered[,split_col]==split_id) )
      #   names(split_ind) = split_ids
      # }
      # 
      # for (tube in names(split_ind)) {
      #   d = d0[split_ind[[tube]],split_ind[[tube]]]
      #   if (!sum(d!=0)>0) next()
      # sm = meta_file_ordered[split_ind[[tube]],]
      # if (length(unique(sm[,target_col]))<=1) next()
      # if (length(unique(sm[,target_col]))<2) next()
      # if (min(sapply(unique(sm[,target_col]), function(x) length(sm[,target_col]==x)))<1) next()
      
      ## for each feature
      for (target_col in target_cols) {
        if (!target_col%in%colnames(sm)) next
        #list out class labels
        class = sm[,target_col]; 
        class_unique = unique(class)
        # las0 = sm$label
        
        dirname = paste0(target_col, "_",
                         ifelse(length(class_unique)<10, 
                                paste0(paste0(class_unique, ".",sapply(class_unique, function(x) sum(sm[,target_col]==x))), collapse="-"), length(class_unique)))
        clust_dir1 = paste0(clust_dir,"/",dirname); dir.create(clust_dir1, showWarnings=F)
        clust_plot_dir1 = paste0(clust_plot_dir,"/",dirname); dir.create(clust_plot_dir1, showWarnings=F)
        
        cname1 = paste0(clust_dir1, "/", dist_type)
        cname2 = paste0(clust_plot_dir1, "/", dist_type)
        #to do or not to do
        if ((file.exists(cname1) | file.exists(cname2)) & !overwrite) next
        
        # plot
        plotn = length(unlist(cmethodspar))
        plotn2 = ceiling(sqrt(plotn))
        tsned = tsne(dist(d), k=2)
        png(cname2, width=plotn2*plotsize, height=plotn2*plotsize)
        par(mfrow=c(plotn2,plotn2))
        plot(tsned, pch=16, cex=1, col=factor(class), main=paste0("tsne plot of dist"))
        legend("topleft", legend=levels(factor(class)), pch=16, col=unique(factor(class)))
        
        sim = get_graph(d) 
        clusts = list()
        ## for each clustering method
        for (cmethod in names(cmethodspar)) {
          for (par in cmethodspar[[cmethod]]) {
            cat(" ",cmethod," ",sep="")
            start2 = Sys.time()
            
            # ## knn for each k (classification)
            # if (cmethod=="knn") {
            #   labelss = class
            #   labelss[is.na(las0)] = NA
            #   clt1 = as.numeric(factor(knntable(as.matrix(d),par,labelss)[,1]))
            #   clt = rep(0,length(class))
            #   clt[is.na(las0)] = clt1
            # } #parameter k
            
            ## kmeans
            if (cmethod=="kmeans") {
              clt = kmeans(dist(d), centers=length(class_unique), nstart=10)$cluster
            } #number of tries
            
            ## kmedoids
            if (cmethod=="kmed") {
              # if (is.na(par)) {
              #   par = length(class_unique)
              #   if (length(class_unique)%in%pars) next
              # } 
              clt = pam(dist(d),length(class_unique))$clustering
            } #number of tries
            
            ## louvain
            if (cmethod=="lv") { #input is a similarity matrix
              sim1 = sim
              tops = quantile(as.vector(sim),par)
              sim1[sim1<tops] = 0
              gr = graph_from_adjacency_matrix(sim1, weighted=T, mode='undirected', diag=F)
              clt = cluster_louvain(gr)$membership
            }
            
            # ## random walk refined distance matrices only
            # if (cmethod=="rw1") {
            #   # } else { next }
            #   clt = rw1table(sim,par)[,1]
            # }
            
            ## spectral clustering via distance matrix
            if (cmethod=="spec1") {
              # matrix power operator: computes M^power (M must be diagonalizable)
              "%^%" <- function(M, power) 
                with(eigen(M), vectors %*% (values^power * solve(vectors)))
              
              D <- diag(apply(sim, 1, sum)) # sum rows
              # L <- D - sim # unnormalized laplacian
              # L <- diag(nrow(my.data)) - solve(D) %*% A  # simple Laplacian
              L <- (D %^% (-1/2)) %*% sim %*% (D %^% (-1/2))  # normalized Laplacian
              
              k   <- length(class_unique)
              evL <- eigen(L, symmetric=TRUE)
              Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
              # plot(Z, col=factor(class), pch=20)
              
              kmeans(Z, centers=length(class_unique), nstart=10)
              
              # pp0t = specc(as.kernelMatrix(d),centers=length(class_unique))@.Data
              # clt = spec1table(sim,class_unique)[,1]
            }
            
            ## hierarchical clustering
            if (cmethod=="hc") { 
              clt = cutree(hclust(as.dist(d),method=par),k=length(class_unique))
            }
            
            # ## densitycut via distance matrix (k=3)
            # if (cmethod=="dc1") { 
            #   clt = dc1table(d)[,1]
            # }
            
            plot(tsned, pch=16, cex=1, col=factor(class), main=paste0("method = ", cmethod,", parameter = ", par, " (NA if none or # of clusters); \n o = cluster, . = actual class"))
            points(tsned, cex=2, col=factor(clt))
            
            names(clt) = rownames(d)
            clusts[[paste0(cmethod,".",par)]] = list()
            clusts[[paste0(cmethod,".",par)]][["x"]] = clt
            
            time_output(start2)
            
          }
        }
        graphics.off()
        save(clusts, file=paste0(cname1,".Rdata"))
      }
      return(F)
    }, error = function(e) {
      cat(paste("ERROR:  ",e)); return(T)
    })
  }, .parallel=T)
  
  time_output(start)
  
  
  
  
  
  
  
  
  
  # #check which features weren't finished
  # aa = distMetafun(list.files(clust_plot_dir,pattern="other"),dis=c(dis,"other"))
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

