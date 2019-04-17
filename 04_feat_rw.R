## Input: original features --> Output: random walk features (edge matrix + paths -- paths need to be saved one at a time, memory will overload otherwise)
# aya43@sfu.ca 20161220

## root directory
root = "~/projects/flowtype_metrics"
setwd(root)

## libraries
source("source/_funcAlice.R")
source("source/_funcdist.R")
libr(c("igraph","diffusr",
       "foreach","doMC",
       "stringr","plyr","Matrix"))

#Setup Cores
no_cores = 3 #detectCores()-6
registerDoMC(no_cores)

## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)


overwrite = T
writecsv = F

countThres = 1000 #(needed if sample x cell pop matrices are used) insignificant if count under
good_sample = c(3)
good_count = c(3)
id_col = "id"
target_col = "class"
# order_cols = c("barcode","sample","specimen","date")
control = "control"

min_steps = 100
times_steps = 50 #times number of edges; number of random walks to do
min_walks = 5 #each path should have at least 5 walks, if not, set to 0

feat_count = c("file-cell-countAdj") #(needed if sample x cell pop matrices are used) count matrix, to get rid of low cell count cells

for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  # result_dir = paste0(root, "/result/impc_panel1_sanger-spleen") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta") # meta files directory
  meta_file_dir = paste(meta_dir, "/file", sep="") #meta for rows (sample)
  meta_cell_dir = paste(meta_dir, "/cell", sep="") #meta for rows (sample)
  feat_dir = paste(result_dir, "/feat", sep="") #feature files directory
  meta_cell_child_names_dir = paste(meta_dir, "/cell_childpn_names",sep="") #specifies a phenotypes children
  
  ## output directories
  

  #data paths
  feat_types = gsub(".Rdata","",list.files(path=feat_dir,pattern=glob2rx("*.Rdata")))
  feat_types = lapply(feat_types[grepl("file-cell",feat_types) & !grepl("rw",feat_types)], function(x) c(x, F,F,T,F))
  # list(c("file-cell-pvalBH.file-cell-countAdj",F,F,T,F),c("file-cell-pvalBH.file-cell-countAdj.PEER-layerbylayer",F,F,T,F)),#c("file-cell-countAdj",T,F,F,T),c("file-cell-prop",F,F,F,F),c("file-cell-pvalbonferroni.file-cell-countAdj.PEER-all",T,F,F,T)) #node/edge, log or not, edge weight should be reversed or not; should edge weights be absolute, edge weight minimum at 0
  
  

  start = Sys.time()
  
  mc = get(load(paste0(feat_dir,"/", feat_count,".Rdata")))
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
  meta_cell0 = get(load(paste0(meta_cell_dir,".Rdata")))
  meta_cell_child_names0 = get(load(paste0(meta_cell_child_names_dir,".Rdata")))
  
  for (feat_type in feat_types) {
    # fe = foreach (feat_type = feat_types) %dopar% {
    # feat_typee = feat_type[[1]]
    tryCatch({
      cat("\n", feat_type, " ",sep="")
      start2 = Sys.time()
      print("preping matrix")
      
      ## upload and prep feature matrix
      m0 = Matrix(get(load(paste0(feat_dir,"/", feat_type[1],".Rdata"))))
      # if (!rownames(m0)[1]%in%meta_file[,id_col]) {
      #   cat("\nskipped: ",feat_type,", matrix rownames must match fileName column in meta_file","\n", sep="")
      # }
      
      ## does feature matrix have cell populations on column names?
      # layers = 0
      # countThres = 0
      colhascell = ifelse(!any(grepl("_",colnames(m0))),T,F)
      # if (colhascell) {
      #   layers = c(1,2,4,max(unique(sapply(unlist(str_split(colnames(m0),"_")[[1]]), function(x) str_count(x,"[+-]")))))
      #   countThres = cellCountThres
      # }
      
      #get matrix colnames/rownames (sample/celltype features)
      m0cn = colnames(m0)
      m0rn = rownames(m0)
      
      #get feature layers if feature names represent cell types
      if ( colhascell) meta_cell = meta_cell0[match(m0cn,meta_cell0$phenotype),]
      if (!colhascell) meta_cell = meta_cell0[match(sapply(m0cn, function(x) str_split(x,"_")[[1]][2]),meta_cell0$phenotype),]
      
      ## for each layer, trim feature matrix accordingly
      # for (k in layers) {
      #trim matrix
      m = m0
      # mm = trimMatrix(m0,TRIM=T, mc=mc, sampleMeta=meta_file, sampleMeta_to_m1_col=id_col, 
      #                 target_col=target_col, control=control, # order_cols=order_cols,
      #                 colsplitlen=meta_cell$phenolevel, k=max(meta_cell$phenolevel), 
      #                 countThres=countThres, goodcount=good_count, good_sample=good_sample)
      # if (is.null(mm)) next
      # m = m_ordered = mm$m
      # mf = meta_file_ordered = mm$sm
      mf = meta_file_ordered = meta_file
      if (all(meta_file_ordered[,target_col]==meta_file_ordered[1,target_col])) next
      
      # file name & check if overwrite
      ename = paste0(feat_dir, "/file-edge-rw.", feat_type[1]); # dir.create(ename, showWarnings=F)
      dname = paste0(feat_dir, "/file-cell-rw.", feat_type[1]); # dir.create(dname, showWarnings=F)
      if (overwrite | !file.exists(paste0(ename,".Rdata")) & !file.exists(paste0(dname,".Rdata")) ) {
        
        
        # define cell populations
        m0cnphen = sapply(str_split(m0cn,"_"), function(x) x[length(x)])
        pm = meta_cell[match(m0cnphen,meta_cell$phenotype),]
        
        time_output(start2)
        print("preping edgelist")
        
        # edge list
        if (colhascell) {
          meta_cell_child_names = meta_cell_child_names0[match(pm$phenotype, names(meta_cell_child_names0))]
          edgelist0 = sapply(1:length(meta_cell_child_names), function(xi) {
            xn = names(meta_cell_child_names)[xi]
            x = meta_cell_child_names[[xi]]
            el = sapply(unlist(x), function(y) c(xn,y))
          })
          edgelist0 = Reduce("cbind",edgelist0)
          edgelist0 = t(edgelist0)
          edgelist = edgelist0[edgelist0[,1]%in%m0cnphen & edgelist0[,2]%in%m0cnphen,]
        } else {
          edgelist = Reduce("rbind",str_split(colnames(m),"_"))
        }
        
        time_output(start2)
        print("preping random walks and tallies for each sample")
        
        
        # loop.ind = loopInd(1:nrow(mf), no_cores)
        # random_edges_all_files = list()
        # random_paths_all_files = list()
        loop_ind = 1:nrow(mf)
        tallies = foreach (i = loop_ind) %dopar% {
          # lname = paste0(ename,"/",mf[i,id_col],".Rdata")
          # pname = paste0(dname,"/",mf[i,id_col],".Rdata")
          
          # if (!overwrite & file.exists(lname) & file.exists(pname)) return(T)
          
          #prepare nodes and edges
          weight = as.numeric(m[i,])
          if (colhascell) weight = weight[match(edgelist[,2], meta_cell$phenotype)]
          edges = data.frame(from=unlist(edgelist[,1]), to=unlist(edgelist[,2]), weight)
          #log abs weight
          if (feat_type[2]!="FALSE") { 
            a = log(abs(edges$weight))
            negind = edges$weight<0
            a[negind] = -a[negind]
            a[a==-Inf] = min(a[a!=-Inf])
            edges$weight = a
          }
          #inverse weight by minus
          if (feat_type[3]!="FALSE") edges$weight = max(edges$weight) - edges$weight
          #abs weight
          if (feat_type[4]!="FALSE") edges$weight = abs(edges$weight)
          #min to 0 weight
          if (feat_type[5]!="FALSE") edges$weight = edges$weight - min(edges$weight)
          
          edges = edges[edges$from%in%meta_cell$phenotype & edges$to%in%meta_cell$phenotype,]
          from_layers = cell_type_layers(edges$from)
          e_order = order(from_layers)
          edges = edges[e_order,]
          from_layers = from_layers[e_order]
          from_layers_startind = from_layers_startind0 = sapply(sort(unique(from_layers)), function(x) min(which(from_layers==x)))
          
          #add leaf to root node edges + duplicate edges with higher weight as random_walk() doesn't recognize weights
          edges1 = edges
          
          q95 = quantile(edges1$weight, .95)
          edges1$weight = ceiling(100 * edges1$weight / q95)
          edges1$weight[edges1$weight > 100] = 100
          edges1$weight[is.nan(edges1$weight)] = min(edges1$weight[!is.nan(edges1$weight)])
          edges1$weight = edges1$weight - min(edges1$weight) + 1
          edges1 = edges1[rep(1:nrow(edges1), edges1$weight), 1:2]
          
          leaf_nodes = setdiff(edges1$to,edges1$from)
          leaf_edges = lapply(leaf_nodes, function(x) c(x,"")) #connect all leaf nodes to root node
          leaf_edges = Reduce("rbind",leaf_edges)
          
          edges1 = as.data.frame(rbind(as.matrix(edges1),as.matrix(leaf_edges))) #,weight=append(as.numeric(edges[,3]),as.numeric(leaf_edges[,3])))
          # edges1[edges1==""] = "0"
          # nodes = nodes[nodes$cell%in%unique_cells,]
          if (!""%in%edges1$from) {
            cly1 = cell_type_layers(edges1$from)
            edges1roottemp = unique(edges1$from[cly1==min(cly1)])
            edges1 = rbind(data.frame(from=rep("",length(edges1roottemp)), to=edges1roottemp), edges1)
          }
          
          #random walk
          g = graph_from_data_frame(edges1,directed=T,unique(unlist(edges1)))
          steps = max(min_steps,nrow(edges1)*times_steps)
          # walks0 = random_walk(g,start="",steps=max(min_steps,nrow(edges1)*1000),mode="out")
          walks0 = random_walk(g,start="",steps=steps,mode="out")
          walks = as_ids(walks0)
          
          
          #tally walks on each path from root to leaf by layer
          # if (overwrite | !file.exists(pname)) {
          # indof0 = which(walks=="")
          # #delete last path
          # walks = walks[-c(indof0[length(indof0)]:length(walks))]
          # indof0 = indof0[-length(indof0)]
          # 
          # paths = sapply(1:(length(indof0)-1), function(ii0) {
          #   starti = indof0[ii0]
          #   endi = indof0[ii0+1]-1
          #   return( paste(walks[starti:endi],collapse="_"))
          # })
          # paths = append(paths, paste(walks[indof0[length(indof0)]:length(walks)],collapse="_"))
          # path_tally = table(paths)
          #     save(path_tally, file=pname); rm(path_tally)
          # } else {
          #   path_tally = get(load(pname))
          # }
          
          # if (overwrite | !file.exists(lname)) {
          
          # alternative code for getting an edge tally
          # edge_tally_table = list()
          # edgeletss = str_split(names(path_tally),"_")
          # for (ei in 1:length(edgeletss)) {
          #   edgelets = edgeletss[[ei]]
          #   if (length(edgelets)<2) next()
          #   einedgelets = sapply(2:length(edgelets), function(x) paste(edgelets[(x-1):x], collapse="_") )
          #   edge_tally_table = append(edge_tally_table, 
          #                             data.frame(edge=einedgelets,
          #                                        freq=rep(path_tally[ei],length(einedgelets))))
          # }
          # edge_tally_table0 = Reduce("rbind",edge_tally_table)
          # edge_tally_table0 = as.data.frame(edge_tally_table0)
          # edge_tally_table1 = ddply(edge_tally_table0,"edge",numcolwise(sum))
          # edge_tally = edge_tally_table1$freq
          # names(edge_tally) = edge_tally_table1$edge
          
          
          edge_tally = rep(0, nrow(edges)) #without leaf to root edges
          names(edge_tally) = paste0(edges[,1], "_", edges[,2])
          for (ei in 1:nrow(edges)) {
            fromi = which(walks==edges$from[ei])
            if (length(fromi)==0) { edge_tally[ei] = 0; next }
            toi = which(walks==edges$to[ei])
            edge_tally[ei] = sum(toi%in%(fromi+1))
          }
          # names(edge_tally) = apply(edges[,1:2],1,function(x) paste(x,collapse="_"))
          # w0ind = which(walks=="")
          # to = 1 #to layer 1
          # while (length(from_layers_startind)>0) {
          #   startind = from_layers_startind[1]
          #   endind = from_layers_startind[2]-1
          #   if (length(from_layers_startind)==1) endind = nrow(edges)
          #   
          #   edges_temp = names(edge_tally)[startind:endind]
          #   wto = w0ind + to; wtodel = wto>length(walks); wto = wto[!wtodel]
          #   edges_walk0 = paste(walks[wto-1],walks[wto],sep="_")
          #   edges_walk = table(edges_walk0)
          #   order_edges_walk = match(edges_temp,names(edges_walk))
          #   order_edges_walk[is.na(order_edges_walk)] = length(edges_walk) + 1
          #   edges_walk = append(edges_walk,0)
          #   edge_tally[startind:endind] = edges_walk[order_edges_walk]
          #   
          #   to = to+1
          #   from_layers_startind = from_layers_startind[-1]
          # }
          # x11(); plot(edge_tally, edges$weight) # random walk feature vs original feature; much like layer weighted feature
          
          # 
          # random_edges_all_files[[mf[i,id_col]]] = edge_tally
          # random_paths_all_files[[mf[i,id_col]]] = path_tally
          
          #   save(edge_tally, file=lname)
          # } else {
          #   edge_tally = get(load(lname))
          # }
          
          return(edge_tally) #, path_tally=path_tally))
        }
        # }
        
        time_output(start2)
        print("putting together edge feature")
        
        # if (!file.exists(paste0(ename,".Rdata"))) {
        #   
        #   random_edges_all_files = foreach(x = list.files(ename, full.names=T)) %dopar% { get(load(x)) }
        #   
        #   edge_names = Reduce('union',lapply(random_edges_all_files, function(x) names(x)))
        #   random_edges_all_files1 = Reduce('rbind', foreach(x=random_edges_all_files) %dopar% {
        #     x0 = x[match(edge_names, names(x))]
        #     x0[is.na(x0)] = 0
        #     return(x0)
        #   })
        #   rownames(random_edges_all_files1) = mf[,id_col]
        #   colnames(random_edges_all_files1) = edge_names
        #   random_edges_all_files1 = Matrix(random_edges_all_files1, sparse=T)
        #   save(random_edges_all_files1, file=paste0(ename,".Rdata"))
        #   # rm(random_edges_all_files1)
        #   # rm(random_edges_all_files)
        #   
        #   # unlink(ename, recursive=T)
        # }
        
        random_edges_all_files_names = Reduce('union', lapply(tallies, function(x) names(x)))
        edge_tallies = ldply(tallies, function(x) {
          xx = x[match(random_edges_all_files_names, names(x))]
          xx[is.na(xx)] = 0
          return(xx)
        })
        rownames(edge_tallies) = rownames(mf[,id_col])
        colnames(edge_tallies) = random_edges_all_files_names
        
        
        time_output(start2)
        # print("putting together path feature")
        # 
        # if (!file.exists(paste0(dname,".Rdata"))) {
        #   feat_type_path_files_full = list.files(path=dname, full.names=T)
        #   feat_type_path_files = list.files(path=dname, full.names=F)
        #   feat_type_path_filenames = gsub(".Rdata","",feat_type_path_files)
        #   
        #   start1 = Sys.time()
        #   mp0 = foreach(x=feat_type_path_files_full) %dopar% {
        #     # data.frame(as.list(get(load(paste0(feat_dir,"/",feat_type_path,"/",x))))) ) #30sec for flowcap
        #     a = get(load(x))
        #     if (grepl("error",a[1],ignore.case=T)) return("NA")
        #     return(a) 
        #   } #30sec for flowcap
        #   time_output(start1)
        #   good_ind = sapply(mp0,function(x) x[1]!="NA")
        #   mp0 = mp0[good_ind]
        #   feat_type_path_files_full = feat_type_path_files_full[good_ind]
        #   feat_type_path_files = feat_type_path_files[good_ind]
        #   feat_type_path_filenames = feat_type_path_filenames[good_ind]
        #   
        #   # start1 = Sys.time()
        #   # mp0 = foreach(x=feat_type_path_files) %dopar% {
        #   #   return( data.frame(as.list(get(load(paste0(feat_dir,"/",feat_type_path,"/",x))))) ) #8 min 14 cores for flowcap
        #   # }
        #   # time_output(start1)
        #   
        #   path_sum = sapply(mp0, function(x) sum(x))
        #   
        #   path_ratio = path_sum/min(path_sum)
        #   mp = lapply(1:length(mp0), function(xi) mp0[[xi]]*path_ratio[xi])
        #   if (min_walks>0) mp = lapply(mp, function(x) x[x>=min_walks])
        #   
        #   start1 = Sys.time()
        #   path_names = lapply(mp0, function(x) names(x))
        #   path_name = unique(unlist(path_names))
        #   # path_order = lapply(path_names, function(x) match(path_name,path_names))
        #   
        #   # mp2 = Matrix(0,ncol=length(path_names),nrow=length(mp), 
        #   #              dimnames=list(gsub(".Rdata","",feat_type_path_files),path_name), sparse=T)
        #   
        #   loop.inds = loopInd(1:length(mp),no_cores)
        #   
        #   mp2 = list()
        #   mp2ind = 1
        #   # mp2 = foreach (ii = loop.inds, .combine=rbind) %dopar% {
        #   for (ii in loop.inds) {
        #     mtemp = Matrix(0,nrow=length(ii),ncol=length(path_name),sparse=T)
        #     for(i in 1:length(ii)) {
        #       mpi = mp[[ii[i]]]
        #       mpiorder = match(path_name,names(mpi))
        #       mpi2 = mpi[mpiorder]
        #       mpi2[is.na(mpiorder)] = 0
        #       # mtempi = c.sparseVector(mpi2)
        #       mtemp[i,] = c.sparseVector(mpi2)
        #       # return(mtempi)
        #     }
        #     # mtemp = lapply(mtempt, as, "sparseMatrix")
        #     # mtemp = Reduce(cbind, mtemp)
        #     
        #     # return(mtemp)
        #     mp2[[mp2ind]] = mtemp
        #     mp2ind = mp2ind + 1
        #   }
        #   mp2 = Reduce("rbind",mp2)
        #   mp2 = Matrix(mp2, sparse=T)
        #   rownames(mp2) = feat_type_path_filenames
        #   colnames(mp2) = path_name
        #   # mp1 = foreach(x=mp) %dopar% {data.frame(as.list(x))}
        #   # mp2 = Reduce(rbind.fill, mp1)
        #   # mp2[is.na(mp2)] = 0
        #   # mp2 = Matrix(mp2, sparse=T)
        #   
        #   save(mp2, file=paste0(feat_dir,"/",dname,".Rdata"))
        #   
        #   time_output(start1)
        #   
        #   # random_paths_all_files1 = Matrix(random_paths_all_files1,sparse=T)
        #   # rownames(random_paths_all_files1) = feat_type_path_files
        #   # colnames(random_paths_all_files1) = path_names
        #   # mp = random_paths_all_files1
        
      }
      
      # time_output(start2)
      print("putting together node feature")
      
      node_tallies_names = sapply(colnames(edge_tallies), function(x) str_split(x,"_")[[1]][2])
      ntnu = unique(node_tallies_names)
      node_tallies = foreach(x=ntnu, .combine='cbind') %dopar% { rowSums(edge_tally[,node_tallies_names==x]) }
      colnames(node_tallies) = ntnu
      
      rownames(node_tallies) = rownames(edge_tallies) = rownames(mf[,id_col])
      
      
      save(edge_tallies, file=paste0(ename, ".Rdata"))
      if (writecsv) write.csv(edge_tallies, file=paste0(ename, ".csv"))
      save(node_tallies, file=paste0(dname, ".Rdata"))
      if (writecsv) write.csv(node_tallies, file=paste0(dname, ".csv"))
      
      time_output(start2)
      
    }, error = function(err) { cat(paste("ERROR:  ",err)) })#; return(T) })
    # return(F)
    # })
  }
  
  time_output(start)
  
}


