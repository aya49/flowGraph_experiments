## Input: original features --> Output: bicluster & plots
# aya43@sfu.ca 20161220

#Directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)


#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

#Input
sampleMeta_dir = paste(result_dir,  "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir,  "/matrix", sep="")

#Output
rw_dir = paste(result_dir,  "/rw", sep=""); dir.create (rw_dir, showWarnings=F)

#Libraries/Functions
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
libr("igraph")
libr("diffusr")
libr("foreach")
libr("doMC")
libr("stringr")

#Setup Cores
no_cores = 2#detectCores()-3
registerDoMC(no_cores)








#Options for script
# sampletimes = 5
overwrite = F #overwrite biclust?
writecsv = F

plot_size = c(500,500)
plot_size_bar = c(1300,2000)
attributes = c("aml")

cellCountThres = c(200) #(needed if sample x cell pop matrices are used) insignificant if count under
good_samples = c(3)
control = "normal" #control value in target_col column for each centre
date_col = "date"
target_col = "aml"
order_cols = c("specimen")

# good_samples = c(3)

dis = c("rw")
min_steps = 100
tube = 6 #which panel to use (flowCAP-II only); any of tubes 1-7

#data paths
#matrix_type = list.files(path=result_dir,pattern=glob2rx("matrix*.Rdata"))
matrix_type = list(c("Pval_CountAdj","",F,T),c("PvalTRIM_CountAdj","",F,T)) #node, edge, if edge is blank, use node values on edges (of child node); edge weight should be reversed or not; should edge weights be absolute

# matrix_type = c( "Child_entropyTRIM_CountAdjBH"              ,"Child_entropyTRIM_CountAdjbonferroni"      ,"Child_entropyTRIM_CountAdjBY"              ,"Child_entropyTRIM_CountAdjPEERBH"         
#                  ,"Child_entropyTRIM_CountAdjPEERbonferroni"  ,"Child_entropyTRIM_CountAdjPEERBY"          ,"Child_entropyTRIM_CountAdjPEER"            ,"Child_entropyTRIM_CountAdj"               
#                  ,"CountAdjPEER"                              ,"CountAdj"                                  ,"LogFold_CountAdjPEER"                      ,"LogFold_CountAdj"                         
#                  ,"LogFoldTRIM_CountAdjBH"                    ,"LogFoldTRIM_CountAdjbonferroni"            ,"LogFoldTRIM_CountAdjBY"                    ,"LogFoldTRIM_CountAdjPEERBH"               
#                  ,"LogFoldTRIM_CountAdjPEERbonferroni"        ,"LogFoldTRIM_CountAdjPEERBY"                ,"LogFoldTRIM_CountAdjPEER"                  ,"LogFoldTRIM_CountAdj"                     
#                  ,"Parent_entropyTRIM_CountAdjBH"             ,"Parent_entropyTRIM_CountAdjbonferroni"     ,"Parent_entropyTRIM_CountAdjBY"             ,"Parent_entropyTRIM_CountAdjPEERBH"        
#                  ,"Parent_entropyTRIM_CountAdjPEERbonferroni" ,"Parent_entropyTRIM_CountAdjPEERBY"         ,"Parent_entropyTRIM_CountAdjPEER"           ,"Parent_entropyTRIM_CountAdj"              
#                  ,"Pval_CountAdjBH"                           ,"Pval_CountAdjbonferroni"                   ,"Pval_CountAdjBY"                           ,"Pval_CountAdjPEERBH"                      
#                  ,"Pval_CountAdjPEERbonferroni"               ,"Pval_CountAdjPEERBY"                       ,"Pval_CountAdjPEER"                         ,"Pval_CountAdj"                            
#                  ,"PvalTRIM_CountAdjBH"                       ,"PvalTRIM_CountAdjbonferroni"               ,"PvalTRIM_CountAdjBY"                       ,"PvalTRIM_CountAdjPEERBH"                  
#                  ,"PvalTRIM_CountAdjPEERbonferroni"           ,"PvalTRIM_CountAdjPEERBY"                   ,"PvalTRIM_CountAdjPEER"                     ,"PvalTRIM_CountAdj"       
#                  ,"Child_pnratioTRIM_CountAdjBH"              ,"Child_pnratioTRIM_CountAdjbonferroni"      ,"Child_pnratioTRIM_CountAdjBY"              ,"Child_pnratioTRIM_CountAdjPEERBH"         
#                  ,"Child_pnratioTRIM_CountAdjPEERbonferroni"  ,"Child_pnratioTRIM_CountAdjPEERBY"          ,"Child_pnratioTRIM_CountAdjPEER"            ,"Child_pnratioTRIM_CountAdj"               
#                  ,"Child_propTRIM_CountAdjBH"                 ,"Child_propTRIM_CountAdjbonferroni"         ,"Child_propTRIM_CountAdjBY"                 ,"Child_propTRIM_CountAdjPEERBH"            
#                  ,"Child_propTRIM_CountAdjPEERbonferroni"     ,"Child_propTRIM_CountAdjPEERBY"             ,"Child_propTRIM_CountAdjPEER"               ,"Child_propTRIM_CountAdj"                  
#                  ,"Parent_contrib_CountAdjPEER"               ,"Parent_contrib_CountAdj"                   ,"Parent_contribTRIM_CountAdjBH"             ,"Parent_contribTRIM_CountAdjbonferroni"    
#                  ,"Parent_contribTRIM_CountAdjBY"             ,"Parent_contribTRIM_CountAdjPEERBH"         ,"Parent_contribTRIM_CountAdjPEERbonferroni" ,"Parent_contribTRIM_CountAdjPEERBY"        
#                  ,"Parent_contribTRIM_CountAdjPEER"           ,"Parent_contribTRIM_CountAdj"               ,"Parent_effort_CountAdjPEER"                ,"Parent_effort_CountAdj"                   
#                  ,"Parent_effortTRIM_CountAdjBH"              ,"Parent_effortTRIM_CountAdjbonferroni"      ,"Parent_effortTRIM_CountAdjBY"              ,"Parent_effortTRIM_CountAdjPEERBH"         
#                  ,"Parent_effortTRIM_CountAdjPEERbonferroni"  ,"Parent_effortTRIM_CountAdjPEERBY"          ,"Parent_effortTRIM_CountAdjPEER"            ,"Parent_effortTRIM_CountAdj"               
# )
matrix_count = c("CountAdj") #(needed if sample x cell pop matrices are used) count matrix, to get rid of low cell count cells
















start = Sys.time()

start1 = Sys.time()
centre = paste0(panelL," ",centreL)
cat("\n",centre,sep="")

#Prepare data
# matrix_type = list.files(path=paste(result_dir,  sep=""),pattern=glob2rx("matrix*.Rdata"))
# matrix_type = matrix_type[!matrix_type=="matrixCount.Rdata" & !matrix_type=="matrixProp.Rdata"]
# matrix_type = gsub("matrix|.Rdata","",matrix_type)
# matrix_type = matrix_type[!grepl("leaf",matrix_type)]
# matrix_type = matrix_type[!grepl("KO|Max",matrix_type)]
# matrix_type = matrix_type[!grepl("Count_sample",matrix_type)]
# matrix_type = matrix_type[grepl("CountAdj",matrix_type)]

mc = get(load(paste0(matrix_dir, matrix_count,".Rdata")))
sampleMeta = get(load(sampleMeta_dir))

for (mcp in matrix_type) {
  tryCatch({
    cat("\n", mcp, " ",sep="")
    start2 = Sys.time()
    
    ## upload and prep node matrix
    if (file.exists(paste0(matrix_dir, mcp[1],".Rdata"))) { m0 = get(load(paste0(matrix_dir, mcp[1],".Rdata")))
    } else if (file.exists(paste0(matrix_dir, mcp[1],".csv"))) { m0 = read.csv(paste0(matrix_dir, mcp[1],".csv", header=T))
    } else { next } 
    
    ## upload and prep edge matrix
    if (mcp[2]!="") {
      if (file.exists(paste0(matrix_dir, mcp[2],".Rdata"))) { e0 = get(load(paste0(matrix_dir, mcp[2],".Rdata")))
      } else { next }
    }
    
    
    #get matrix colnames/rownames (sample/celltype features)
    m0cn = colnames(m0)
    m0rn = rownames(m0)
    
    #get feature layers if feature names represent cell types
    colsplitlen = cell_type_layers(m0cn)
    # k0 = 0; if (!is.null(colsplitlen)) k0 = c(2,3,4,max(colsplitlen)) #1,4, # how many layers to consider i.e. k=max(phenolevel) only
    
    #get to-delete low count phenotype indices; CountAdj should be first one
    for (countThres in cellCountThres) { cat("\ncountThres: ",countThres," >",sep="")
      
      mm = trimMatrix(m0, TRIM=grepl("TRIM",mcp[1]), mc=mc,
                      sampleMeta=sampleMeta, sampleMeta_to_m1_col="fileName", target_col=target_col, control=str_split(control,"[|]")[[1]], order_cols=NULL, 
                      colsplitlen=NULL, k=max(colsplitlen), countThres=countThres, goodcount=1, good_sample=1)
      m = mm$m
      sm = mm$sm
      
      if (mcp[2]!="") {
        em = trimMatrix(e0, TRIM=grepl("TRIM",mcp[2]), mc=mc,
                        sampleMeta=sampleMeta, sampleMeta_to_m1_col="fileName", target_col=target_col, control=str_split(control,"[|]")[[1]], order_cols=NULL, 
                        colsplitlen=NULL, k=max(colsplitlen), countThres=countThres, goodcount=1, good_sample=1)
        e = em$m
        sme = em$sm
      }
      
      #check matrix; if matrix is a list, merge; if matrix is empty, skip
      if (is.null(m)) next
      if (all(m==0)) next
      
      #trim matrix rows by tube
      m = m[sm$tube==tube,]
      sm = sm[sm$tube==tube,]
      if (mcp[2]!="") {
        tubeind = sme$tube==tube
      e = lapply(e, function(x) {
        x0 = x[tubeind,]
        if (is.null(dim(x0))) {
          x0 = matrix(x0,ncol=1)
          colnames(x0) = colnames(x)
          rownames(x0) = rownames(x)[tubeind]
        }
        return(x0)
      })
      sme = sme[tubeind,]
      }
      
      dname = paste(rw_dir, "/", mcp[1], ifelse(mcp[2]!="","_edge-",""), mcp[2], "_", dis, "_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_goodSample-",good_sample, sep = "")
      
      
      random_walks_all_files = foreach(i = 1:nrow(sm)) %dopar% {
        
        #prepare nodes and edges
        nodes0 = m[i,m[i,]!=0]
        nodes = data.frame(cell=names(nodes0),value=nodes0)
        if (mcp[2]!="") {
          edges0 = lapply(names(e), function(ei) {
            erows = lapply(1:ncol(e[[ei]]), function(jj) {
              j = colnames(e[[ei]])[jj]
              # if (j=="") j = "_"
              return(c(j, ei, e[[ei]][i,jj]))
            })
            return(Reduce("rbind",erows))
          })
          
        } else {
          require(flowType)
          markers = unique(unlist(str_split(nodes$cell,"[+-]")))
          phenocode = sapply(nodes$cell, function(x) encodePhenotype(x,markers[markers!=""]))
          phenolevel = unlist(lapply(phenocode, function(x){return(length(markers) - charOccurences("0", x) - 1) } ))
          pm = data.frame(phenotype=nodes$cell, phenocode=phenocode, phenolevel=phenolevel, stringsAsFactors=F)
          pmchild = getphenoChild(pm)
          edges0 = lapply(1:length(pmchild$phenoChild), function(i) {
            children = lapply(unlist(pmchild$phenoChild[[i]]), function(c) c(names(pmchild$phenoChild)[i], pm$phenotype[c], nodes$value[c]))
            return(Reduce("rbind",children))
          })
        }
        edges = Reduce("rbind",edges0)
        edges = data.frame(from=edges[,1],to=edges[,2],weight=as.numeric(edges[,3]))
        if (mcp[3]!="False") edges$weight = max(edges$weight) - edges$weight
        if (mcp[4]=="100") edges$weight = abs(edges$weight)
        
        unique_cells = intersect(c(edges$from,edges$to),nodes$cell)
        edges = edges[edges$from%in%unique_cells & edges$to%in%unique_cells,]
        from_layers = cell_type_layers(edges$from)
        e_order = order(from_layers)
        edges = edges[e_order,]
        from_layers = from_layers[e_order]
        from_layers_startind = from_layers_startind0 = sapply(0:max(from_layers), function(x) min(which(from_layers==x)))
        
        #add leaf to root node edges + duplicate edges with higher weight as random_walk() doesn't recognize weights
        leaf_nodes = setdiff(edges$to,edges$from)
        leaf_edges = lapply(leaf_nodes, function(x) c(x,"",max(edges$weight))) #connect all leaf nodes to root node
        leaf_edges = Reduce("rbind",leaf_edges)
        
        edges1 = edges
        
        q95 = quantile(edges1$weight, .95)
        edges1$weight = ceiling(100 * edges1$weight / q95)
        edges1$weight[edges1$weight > 100] = 100
        edges1 = edges1[rep(1:nrow(edges1), edges1$weight), 1:2]
        
        edges1 = data.frame(from=append(edges1[,1],leaf_edges[,1]),to=append(edges1[,2],leaf_edges[,2])) #,weight=append(as.numeric(edges[,3]),as.numeric(leaf_edges[,3])))
        nodes = nodes[nodes$cell%in%unique_cells,]
        
        #random walk
        g = graph_from_data_frame(edges1,directed=T,vertices=nodes)
        # walks0 = random_walk(g,start="",steps=max(min_steps,nrow(edges1)*1000),mode="out")
        walks0 = random_walk(g,start="",steps=max(min_steps,nrow(edges1)*500),mode="out")
        walks = as_ids(walks0)
        
        #tally walks on each edge by layer
        edge_tally = rep(0,nrow(edges)) #without leaf to root edges
        names(edge_tally) = apply(edges[,1:2],1,function(x) paste(x,collapse="_"))
        w0ind = which(walks=="")
        to = 1 #to layer 1
        while (length(from_layers_startind)>0) {
          startind = from_layers_startind[1]
          endind = from_layers_startind[2]-1
          if (length(from_layers_startind)==1) endind = nrow(edges)
          
          edges_temp = names(edge_tally)[startind:endind]
          wto = w0ind + to; wtodel = wto>length(walks); wto = wto[!wtodel]
          edges_walk0 = sapply(wto,function(wtoi) paste(walks[c(wtoi-1,wtoi)],collapse="_"))
          edges_walk = table(edges_walk0)
          order_edges_walk = match(edges_temp,names(edges_walk))
          order_edges_walk[is.na(order_edges_walk)] = length(edges_walk) + 1
          edges_walk = append(edges_walk,0)
          edge_tally[startind:endind] = edges_walk[order_edges_walk]
          
          to = to+1
          from_layers_startind = from_layers_startind[-1]
        }
        # x11(); plot(edge_tally, edges$weight) # random walk feature vs original feature; much like layer weighted feature
        return(edge_tally)
      }
      
      edge_names = Reduce('union',lapply(random_walks_all_files, function(x) names(x)))
      random_walks_all_files1 = Reduce('rbind', lapply(random_walks_all_files, function(x) {
        x[match(edge_names, names(x))]
        x[is.na(x)] = 0
        return(x)
      }))
      rownames(random_walks_all_files1) = sm$fileName
      colnames(random_walks_all_files1) = edge_names
      save(random_walks_all_files1, file=paste0(dname,".Rdata"))
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    TimeOutput(start2)
  }, error = function(err) { cat(paste("ERROR:  ",err)); return(T) })
  return(F)
}
TimeOutput(start1)




