## input: features 
## output: distance matrices; see k to know how many layers of the cell hierarchy was used (4 for now)


## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr(c("stringr","colorspace", "Matrix", "plyr",
       "vegan", # libr(proxy)
       "foreach","doMC",
       "lubridate", #if there are date variables
       "kernlab"))



#Setup Cores
no_cores = 6 # detectCores()-5
registerDoMC(no_cores)



## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

overwrite = T #overwrite distances?
writecsv = F

readcsv = F

good_count = 3
good_sample = 3 # each class must have more thatn good_sample amount of samples or else prune
countThres = 1000 #insignificant if count under

id_col = "id"
target_col = "class"
order_cols = NULL

control = "control"

dis = c("manhattan", "euclidean", "parthasarathyOgihara", "focus") #distances metrics to use #, "binomial", "cao", "jaccard","mahalanobis", "gower","altGower","morisita","horn", "manhattan", "maximum", "binary", "minkowski", "raup", "mountford"
# disnoavg = c("euclidean", "manhattan", "canberra", "bray", "kulczynski", "morisita", "horn", "binomial")  #dis measures that don't average over number of features
disindist = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski") #use dist function when possible, instead of vegdist
disnoneg = c("canberra") #dis measures that can't handle negative values
# disinkernel = c("rbf","poly","vanilla","tanh","laplace","bessel","anova","spline") #kernels

# normalize = c("none","cellpop", "layer") # by none (for child matrices only), cell pop, layer

k = 4 # layers used


feat_count = c("file-cell-countAdj") #(needed if sample x cell pop matrices are used) count matrix, to get rid of low cell count cells


for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  # result_dir = paste0(root, "/result/impc_panel1_sanger-spleen") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta") # meta files directory
  meta_file_dir = paste(meta_dir, "/file", sep="") #meta for rows (sample)
  meta_cell_dir = paste(meta_dir, "/cell", sep="") #meta for rows (sample)
  meta_cell_childpn_names_dir = paste(meta_dir, "/cell_childpn_names",sep="") #specifies a phenotypes children and splits them into +/- (only for when both -/+ exists)
  meta_cell_parent_names_dir = paste(meta_dir, "/cell_parent_names",sep="") #specifies a phenotypes parents
  feat_dir = paste(result_dir, "/feat", sep="") #feature files directory
  
  
  ## output directories
  dist_dir = paste0(result_dir,"/dist"); suppressWarnings(dir.create (dist_dir))
  
  
  
  #data paths
  feat_types = list.files(path=feat_dir,pattern=".Rdata")
  feat_types = feat_types[!feat_types=="file-cell-count.Rdata"]
  feat_types = gsub(".Rdata","",feat_types)
  feat_types = feat_types[!grepl("KO|Max",feat_types)]
  # feat_types = feat_types[!grepl("Freqp_orig",feat_types)]
  
  
  
  # #convert to csv
  # m_paths = paste0(feat_dir,feat_types,".Rdata")
  # foreach(mp = m_paths) %dopar% {
  #   mm = get(load(mp))
  #   if (is.null(dim(mm))) mm = Reduce('cbind',mm)
  #   write.table(mm,sep=",",row.names=F,col.names=F,file=gsub(".Rdata",".csv",mp))
  # }
  # write.table(gsub("~/","/home/ayue/",gsub(".Rdata",".csv",m_paths)),file=paste0(result_dir,"/featlist.csv"),sep=",",row.names=F,col.names=F)
  
  
  mc = Matrix(get(load(paste0(feat_dir,"/", feat_count,".Rdata"))))
  meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
  
  ## made for li-ogihara-zhou, but takes too long
  # phenocodes = llply(meta_cell$phenocode, function(x) as.numeric(laply(seq(1, nchar(x), 1), function(i) substr(x, i, i))))
  # phenocodes = Reduce('rbind',phenocodes)
  # pd = Matrix(0, nrow=nrow(meta_cell), ncol=nrow(meta_cell), sparse=T)
  # colnames(pd) = rownames(pd) = meta_cell$phenotype
  # for (i in 2:nrow(pd)) 
  #   for (j in 1:(i-1)) 
  #     pd[i,j] = pd[j,i] = sum(phenocodes[i,]==phenocodes[j,])
  
  
  # layers = c(1,4,max(meta_cell$phenolevel)) # how many markers to consider i.e. k=max(phenolevel) only
  
  meta_cell_childpn_names = get(load(paste0(meta_cell_childpn_names_dir, ".Rdata")))
  meta_cell_parent_names = get(load(paste0(meta_cell_parent_names_dir, ".Rdata")))
  
  
  
  start = Sys.time()
  
  #load different features matrix and calculate distance matrix
  # for (feat_type in feat_types_) {
  a = llply(feat_types, function(feat_type) {
    tryCatch({
    cat("\n", feat_type, " ",sep="")
    start2 = Sys.time()
    
    #start a list of phenotypes included in each distance matrix calculation such that none are repeated
    # leavePhenotype = list()
    # doneAll = F
    
    ## upload and prep feature matrix
    m0 = Matrix(get(load(paste0(feat_dir,"/", feat_type,".Rdata"))))
    iscell = grepl("_",colnames(m0))[10]
    
    ## does feature matrix have cell populations on column names?
    # layers = 0
    # countThres = 0
    # colhascell = !grepl("_",colnames(m0)[1])
    # if (colhascell) {
    #   layers = c(1,2,4,max(unique(sapply(unlist(str_split(colnames(m0),"_")), function(x) str_count(x,"[+-]")))))
    #   countThres = cellCountThres
    # }
    
    ## for each layer, trim feature matrix accordingly
    ## trim matrix
    # mm = trimMatrix(m0,TRIM=T, mc=mc, sampleMeta=meta_file, sampleMeta_to_m1_col=id_col, target_col=target_col, control=control, order_cols=order_cols, colsplitlen=NULL, k=k, countThres=countThres, goodcount=good_count, good_sample=good_sample)
    sm = meta_file[match(rownames(m0),meta_file$id),]
    # good_sample
    tcl = table(sm$class); delrow = rownames(m0)%in%sm$id[sm$class%in%names(tcl)[tcl<=good_sample]]
    sm = sm[!delrow,]
    # layer k or less cell pops
    # pm = meta_cell[match(sapply(str_split(colnames(m0),"_"), function(x) x[1]),meta_cell$phenotype),]
    if (grepl("_",colnames(m0)[1])) {
      pm = getPhen(sapply(str_split(colnames(m0),"_"), function(x) x[1]))
    } else {
      pm = getPhen(colnames(m0))
    }
    m = m0[!delrow, pm$phenolevel<=k]
    if (all(sm$class==sm$class[1])) return(F)
    
    #for every distance type
    a = match(disnoneg,dis)
    loop.ind = 1:length(dis); if (sum(!is.na(a))>0) loop.ind = loop.ind[-a]
    # foreach(i=loop.ind) %dopar% { #for each phenotype
    for (i in loop.ind) {
      cat(", ", length(dis)-i+1, ":", dis[i], " ", sep="")
      
      #assume after splitting dist filename by "_", distance is second element
      dname = paste0(dist_dir, "/",feat_type, "_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_dist-",dis[i],".Rdata")
      
      ## calculate distances
      # if (is.null(dim(m))) { #edge feature
      #   
      #   if ("none"%in%normalize) { #none
      #     dname0 = paste0(dname, "none.Rdata", sep="")
      #     if (overwrite | !file.exists(dname0)) {
      #       if (dis[i]%in%disindist) { d = dist(a, method=dis[i])
      #       } else { d = vegdist(a, method=dis[i]) }
      #       
      #       save(d, file=dname0); if (writecsv) write.csv(as.matrix(d), file=gsub(".Rdata",".csv",checkm(d,dname0)))
      #     }
      #   } 
      #   if (("cellpop"%in%normalize | "layer"%in%normalize)) {
      #     if ("cellpop"%in%normalize) { #cell pop
      #       dc = matrix(0,nrow=nrow(m[[1]]),ncol=nrow(m[[1]]))
      #       dname1 = paste0(dname, "cellpop.Rdata", sep="")
      #     }
      #     if ("layer"%in%normalize) { #layer
      #       dl = lapply(unique(pm$phenolevel), function(y) return(matrix(0, nrow=nrow(m[[1]]), ncol=nrow(m[[1]]))) )
      #       names(dl) = as.character(unique(pm$phenolevel))
      #       dname2 = paste0(dname, "layer.Rdata", sep="")
      #     }
      #     if (overwrite | !file.exists(dname1) | !file.exists(dname2)) {
      #       for (x in 1:length(m)) {
      #         if (dis[i]%in%disindist) { dx = as.matrix(dist(m[[x]],method=dis[i]))
      #         } else { dx = as.matrix(vegdist(m[[x]],method=dis[i])) }
      #         if (dis[i]%in%disindist) dx = dx/ncol(m[[x]])
      #         if ("cellpop"%in%normalize) dc = dc + dx
      #         if ("layer"%in%normalize) dl[[as.character(pm$phenolevel[x])]] = dl[[as.character(pm$phenolevel[x])]] + dx #1 indexing (layer starts at 0)
      #       }
      #       if ("cellpop"%in%normalize) {
      #         save(dc, file=dname1); if (writecsv) write.csv(as.matrix(dc), file=gsub(".Rdata",".csv",checkm(d,dname1)))
      #       }
      #       if ("layer"%in%normalize) {
      #         dl = lapply(unique(pm$phenolevel), function(x) {
      #           a = dl[[as.character(x)]]
      #           a[a>0] = a[a>0]/sum(pm$phenolevel==x)
      #           return(a)
      #         })
      #         a = Reduce('+',dl); a[a>0] = a[a>0]/length(dl) #get mean of all matrices
      #         dl = as.dist(a)
      #         
      #         save(dl, file=dname2); if (writecsv) write.csv(as.matrix(dl), file=gsub(".Rdata",".csv",checkm(d,dname2)))
      #       }
      #     }
      #   }
      #   
      # } else { #node feature
      # if ("cellpop"%in%normalize) { #cell pop
      # dname1 = paste0(dname, "cellpop.Rdata")
      if (overwrite | !file.exists(dname)) {
        d = Matrix(0,nrow=nrow(m), ncol=nrow(m), sparse=T)
        colnames(d) = rownames(d) = rownames(m)
        m_ = m-min(m)
        m_ = m_/max(m_)
        
        if (dis[i]=="parthasarathyOgihara") { if (!grepl("_cell_",feat_type)) next
          ## - Parthasarathy-Ogihara for prop: 1-(sum(max(0,(1-theta*abs(a-b)))) / (length(union(a,b)))
          theta = 1
          minds = llply(1:nrow(m), function(i) m[i,]!=0)
          for (i in 1:(nrow(m)-1)) {
            for (j in (i+1):nrow(m)) {
              mij = m_[c(i,j), minds[[i]] & minds[[j]]]
              mdiff = 1-theta*abs(mij[1,]-mij[2,])
              mdiffsum = sum(sapply(a,function(x) max(x,0)))
              d[i,j] = d[j,i] = 1-(mdiffsum/sum(minds[[i]] | minds[[j]]))
            }
          }
        } else if (dis[i]=="focus") { if (!grepl("_cell_",feat_type)) next
          ## - FOCUS: sum(abs(a-b)) / (sum(a)+sum(b))
          for (i in 1:(nrow(m)-1)) {
            for (j in (i+1):nrow(m)) {
              mij = m_[c(i,j), minds[[i]] & minds[[j]]]
              mdiff = sum(abs(mij[1,]-mij[2,]))
              mdiff2 = sum(mij[1,]+mij[2,])
              d[i,j] = d[j,i] = mdiff / mdiff2
            }
          }
        } else if (dis[i]%in%disindist) {
          d = dist(m, method=dis[i])
        } else { 
          d = vegdist(m, method=dis[i]) 
        }
        
        ## - li-origihara-zhou: 1-(2*i3)/(i1+i2) -- normalization for mutual info
        #   - i1 = sum(sapply(every pair of cell pops in a, d(a1,a2)))
        #   - i2 = sum(sapply(every pair of cell pops in b, d(b1,b2)))
        #   - i3 = sum(sapply(every pair of cell pops in a,b, d(a1,b2))) -- mutual info between a,b
        #   - d(x,y) = c*log(1+c)*min(supp(x),supp(y))
        
        save(d, file=dname)
        if (writecsv) write.csv(as.matrix(d), file=gsub(".Rdata",".csv",checkm(d,dname)))
      }
      # }
      # if ("layer"%in%normalize) { #layer
      #   dname2 = paste0(dname, "layer.Rdata", sep="")
      #   if (overwrite | !file.exists(dname2)) {
      #     layers = unique(pm$phenolevel)
      #     #if (sum(layers==0)>0) layers = layers[-which(layers==0)]
      #     if (dis[i]%in%disindist) { dd = lapply(layers, function(x) return(as.matrix(dist(m[,which(pm$phenolevel==x)], method=dis[i]))) )
      #     } else { dd = lapply(layers, function(x) return(as.matrix(vegdist(m[,which(pm$phenolevel==x)], method=dis[i]))) ) }
      #     names(dd) = as.character(layers)
      #     
      #     #dis measures that don't average over number of features need extra processing
      #     if (dis[i]%in%disnoavg) dd = lapply(layers, function(x) {
      #       a = dd[[as.character(x)]]
      #       a[a>0] = a[a>0]/sum(pm$phenolevel==x)
      #       return(a)
      #     })
      #     a = Reduce('+',dd); a[a>0] = a[a>0]/length(dd) #get mean of all matrices
      #     dl = as.dist(a)
      #     
      #     save(dl,file=dname2) if (writecsv) write.csv(as.matrix(d), file=paste0(checkm(d,dname2),".csv"))
      #   }
      # }
      # }
      
    } #dis
    # } #countThres
    
    return(F)
    time_output(start2, feat_type)
    }, error = function(err) { cat(paste("ERROR:  ",err)); return(T) })
  }, .parallel=T)
  
  time_output(start)
  
  
  
  
  # if (convertcsv) {
  #   dist_paths = list.files(dist_dir,pattern=".Rdata",full.names=T,recursive=F)
  #   dist_paths = dist_paths[((grepl("effort|contrib|prop|pnratio",dist_paths) & grepl("normalize-none",dist_paths)) | (!grepl("effort|contrib|prop|pnratio",dist_paths) & grepl("normalize-cellpop",dist_paths))) & !grepl("simmatrix",dist_paths)]
  #   a = foreach(dp = dist_paths) %dopar% {
  #     dm = get(load(dp))
  #     write.table(as.matrix(dm),row.names=F,col.names=F,sep=",",file=gsub(".Rdata",".csv",dp))
  #   }
  #   
  # }
  
  
  
  
  
  
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


