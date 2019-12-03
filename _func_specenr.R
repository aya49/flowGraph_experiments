## generic helper functions -----------------------------------

## load libraries
libr = function(pkgs) {
  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) 
    install.packages(setdiff(pkgs, rownames(installed.packages())), verbose=F)
  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
    if (!requireNamespace("BiocManager", quietly=T)) install.packages("BiocManager")
    BiocManager::install(setdiff(pkgs, rownames(installed.packages())), ask=F)
  }
  sapply(pkgs, require, character.only=T)
}


## prepare parallel backend give number of cores no_cores
parallel_backend = function(no_cores=1) {
  parl = F
  if (no_cores>1) {
    parl = T
    require(doMC)
    registerDoMC(no_cores)
  }
  return(parl)
}


## input: Sys.time() value
## output: formatted time as string; used in time_output function
tstr = function(time, tz="GMT") format(.POSIXct(time,tz=tz), "%H:%M:%S")


## input: Sys.time() value
## output: time elapsed since input
time_output = function(start, message="", tz="GMT") {
  start = as.POSIXct(start)
  end = Sys.time()
  time_elapsed = difftime(end, start, units="secs")
  cat(message, ifelse(message=="","",": "), tstr(start,tz=tz), "-", tstr(end,tz=tz), ">", tstr(time_elapsed,tz=tz), "\n")
}


## input: x=loop indices, n=number of vectors to split x into
## output: list of x with n vectors
loop_ind_f = function(x,n) {
  if (n==1) return(list(x))
  return(split(x, ceiling(seq_along(x)/ceiling(length(x)/n))))
}


## flowgraph class definition -----------------

setClass("flowgraph", 
         slots=list(
           feat="list",#list(node="list", edge="list"), 
           feat_summary="list",#list(node="list", edge="list", 
                             #desc=list(node="list", edge="list")),
           meta="data.frame", 
          markers="character", edge_list="list", 
          graph="list", plot_layout="character",
          etc="list"))


## flowgraph class modification ---------------

# replace marker names
gsub_markers = function(object, markers_new, markers_old=NULL) {
  require(plyr)
  
  if (is.null(markers_old)) {
    if (length(markers_new)!=length(object@markers)) {
      cat("incorrect number of markers\n")
      return(object)
    }
    markers_old = object@markers
  } else {
    if (length(markers_old)!=length(markers_new)) {
      cat("incorrect number of markers\n")
      return(object)
    }
  }
  for (mi in 1:length(markers_old)) {
    object@graph$v$phenotype = 
      gsub(markers_old[mi], markers_new[mi], object@graph$v$phenotype)
    object@graph$e$from = 
      gsub(markers_old[mi], markers_new[mi], object@graph$e$from)
    object@graph$e$to = 
      gsub(markers_old[mi], markers_new[mi], object@graph$e$to)
    
    object@edge_list$parent = llply(object@edge_list$parent, function(x)
      gsub(markers_old[mi], markers_new[mi], x))
    names(object@edge_list$parent) = gsub(markers_old[mi], markers_new[mi], names(object@edge_list$parent))
    object@edge_list$child = llply(object@edge_list$child, function(x) 
      llply(x, function(y) {
        gsub(markers_old[mi], markers_new[mi], y)
      }))
    names(object@edge_list$child) = gsub(markers_old[mi], markers_new[mi], names(object@edge_list$child))
    
  }
  object@feat$node = llply(object@feat$node, function(x) {
    colnames(x) = object@graph$v$phenotype
    x
  })
  if (length(object@feat$edge)>0) {
    ecn = paste0(object@graph$e$from, "_", object@graph$e$to)
    object@feat$edge = llply(object@feat$edge, function(x) {
      colnames(x) = ecn
      x
    })
  }
  return(object)
} 


# replace ids
gsub_ids = function(object, ids_new, ids_old=NULL) {
  require(plyr)
  
  if (is.null(ids_old)) {
    if (length(ids_new)!=nrow(object@meta)) {
      cat("incorrect number of ids\n")
      return(object)
    }
    ids_old = object@meta$id
  } else {
    if (length(ids_old)!=length(ids_new)) {
      cat("incorrect number of ids\n")
      return(object)
    }
  }
  
  ids_ind = match(ids_old, object@meta$id)
  object@meta$id[ids_ind] = ids_new
  object@feat$node = llply(object@feat$node, function(x) {
    rownames(x)[ids_ind] = ids_new
    x
  })
  if (length(object@feat$edge)>0) {
    object@feat$edge = llply(object@feat$edge, function(x) {
      rownames(x)[ids_ind] = ids_new
      x
    })
  }
  return(object)
} 


# extract flowgraph for certain sample id's
extract_samples = function(object, sample_ids) {
  require(plyr)
  
  if (!any(sample_ids%in%object@meta$id)) {
    cat("please provide valid sample id's; see @meta$id\n")
    return(object)
  }
  
  id_inds = match(sample_ids,object@meta$id)
  
  object@meta = object@meta[id_inds,]
  object@feat$node = llply(object@feat$node, function(x) x[id_inds,,drop=F])
  if (length(object@feat$edge)>0) {
    object@feat$edge = llply(object@feat$edge, function(x) x[id_inds,,drop=F])
  }
  return(object)
} 


# extract flowgraph for certain cell populations
extract_phenotypes = function(object, phenotypes) {
  require(plyr)
  
  if (!any(phenotypes%in%object@graph$v$phenotype)) {
    cat("please provide valid phenotypes; see @graph$v\n")
    return(object)
  }
  
  id_inds = match(phenotypes,object@graph$v$phenotype)
  id_inds_ = object@graph$e$from%in%phenotypes & object@graph$e$to%in%phenotypes
  
  object@graph$v = object@graph$v[id_inds,,drop=F]
  object@graph$e = object@graph$e[id_inds_,,drop=F]
  object@gr = set_layout_graph(object@gr, as.function(object@plot_layout))
  
  object@feat$node = llply(object@feat$node, function(x) x[,id_inds,drop=F])
  if (length(object@feat$edge)>0) {
    object@feat$edge = llply(object@feat$edge, function(x) x[,id_inds_,drop=F])
  }
  
  object@edge_list$child = 
    object@edge_list$child[phenotypes[phenotypes%in%names(object@edge_list$child)]]
  object@edge_list$child = llply(object@edge_list$child, function(x) 
    llply(x, function(y) {
      a = y[y%in%phenotypes]
      if (length(a)==0) return(NULL)
      a
    }))
  
  object@edge_list$child = plyr::compact(object@edge_list$child)
  object@edge_list$parent = 
    object@edge_list$parent[phenotypes[phenotypes%in%names(object@edge_list$parent)]]
  object@edge_list$parent = llply(object@edge_list$parent, function(x) {
    a = x[x%in%phenotypes]
    if (length(a)==0) return(NULL)
    a
  })
  object@edge_list$parent = plyr::compact(object@edge_list$parent)
  
  return(object)
} 


# merge flowgraph samples
merge_samples = function(obj1, obj2) {
  require(plyr)
  
  object = obj1
  
  meta2 = matrix(NA,nrow=nrow(obj2@meta),ncol=ncol(obj1@meta))
  for (i in 1:ncol(obj1@meta)) {
    j = which(colnames(obj2@meta)==colnames(obj1@meta)[i])
    if (length(j)!=0) 
      meta2[,i] = obj2@meta[,j]
  }
  colnames(meta2) = colnames(obj1@meta)
  object@meta = rbind(obj1@meta, meta2)
  nfs = intersect(names(obj1@feat$node),names(obj2@feat$node))
  object@feat$node = 
    llply(nfs, function(xi) {
      a = obj1@feat$node[[xi]]
      b = obj2@feat$node[[xi]]
      abcol = intersect(colnames(a),colnames(b))
      ab = rbind(a[,match(abcol,colnames(a)),drop=F],
                 b[match(setdiff(rownames(b),rownames(a)),rownames(b)),match(abcol,colnames(b)),drop=F])
    })
  names(object@feat$node) = nfs
  if (length(object@feat$edge)>0) {
    efs = intersect(names(obj1@feat$edge),names(obj2@feat$edge))
    object@feat$edge = llply(efs, function(xi) {
      a = obj1@feat$edge[[xi]]
      b = obj2@feat$edge[[xi]]
      abcol = intersect(colnames(a),colnames(b))
      ab = rbind(a[,match(abcol,colnames(a)),drop=F],
                 b[match(setdiff(rownames(b),rownames(a)),rownames(b)),match(abcol,colnames(b)),drop=F])
    })
    names(object@feat$edge) = efs
  }
  return(object)
} 


# generic flowgraph merge function
merge_flowgraph = function(obj1, obj2, 
                           method_sample=c("union","intersect","setdiff","none"), 
                           method_phenotype=c("intersect","setdiff","none")) {
  method_sample = match.arg(method_sample)
  method_phenotype = match.arg(method_phenotype)
  
  if (method_sample=="union" & method_phenotype!="intersect") {
    print("if method_sample=union, method_phenotype=intersect;\n
            otherwise, feature become difficult to compare")
    return(NULL)
  }
  if (method_sample=="union") return(merge_samples(obj1, obj2))
  if (method_sample=="none" & method_phenotype=="none") return(obj1)
  
  sample_id1 = obj1@meta$id
  sample_id2 = obj2@meta$id
  sample_id1_int = intersect(sample_id1, sample_id2)
  sample_id1_new = setdiff(sample_id1, sample_id2)
  sample_id1_uni = union(sample_id1, sample_id2)
  
  phen1 = obj1@graph$v$phenotype
  phen2 = obj2@graph$v$phenotype
  phen1_int = intersect(phen1, phen2)
  phen1_new = setdiff(phen1, phen2)
  phen1_uni = union(phen1, phen2)
  
  if (method_sample=="none")
    sample1 = obj1
  
  if (method_sample=="intersect") {
    if (length(sample_id1_int)==0) {
      print("no intersecting samples")
      return(NULL)
    }
    sample1 = extract_samples(obj1, sample_id1_int)
  }
  
  if (method_sample=="setdiff") {
    if (length(sample_id1_new)==0) {
      print("no setdiff samples")
      return(NULL)
    }
    sample1 = extract_samples(obj1, sample_id1_new)
  }
  
  if (method_phenotype=="none") return(sample1)
  
  if (method_phenotype=="intersect") {
    if (length(phen1_int)==0) {
      print("no intersecting phenotypes")
      return(NULL)
    }
    if (all(phen1_int%in%phen2)) return(sample1)
    return(extract_phenotypes(sample1, phen1_int))
  }
  
  if (method_phenotype=="setdiff") {
    if (length(phen1_new)==0) {
      print("no setdiff phenotypes")
      return(NULL)
    }
    if (all(phen1_new%in%phen1)) return(sample1)
    return(extract_phenotypes(sample1, phen1_new))
  }
  
}


## flowgraph initialization helper functions -----------------------------

## input: x matrix (phenotypes on columns, will convert in function) -- x0 is optional, plots according to x0 so if you want to plot more cell populations than those in x
## output: f = normalize factors per sample; fidff = difference between f and peak of count ratio per sample
tmm = function(x, x0=NULL, lib.size, refColumn, cutoff=Inf, pngnames=NULL, mains=NULL, no_cores=1, samplesOnCol=F, 
               # tmm parameters
               Acutoff=-10,
               doWeighting = F,
               logratioTrim = .3,
               sumTrim = 0.05,
               minlogR = 1e-6 #min value of log2((obs/obsn)/(ref/refn))
) {
  require(plyr)
  require(pracma)
  require(foreach)
  
  parl = parallel_backend(no_cores)
  if (no_cores>detectCores()) no_cores = detectCores()
  
  plotimg = F
  if(!is.null(pngnames)) plotimg = T
  
  if (is.null(x0)) x0 = x
  if (!samplesOnCol) { x = t(x); x0 = t(x0) } # original function wants samples on the column
  
  ## Taken from TMM 
  ref = x[,refColumn]
  refn = lib.size[refColumn]
  
  ff = foreach(i=1:ncol(x), .combine = list, .maxcombine = ncol(x), .multicombine = T) %dopar% {
    #for(i in ncol(x):1) { cat(i," ",sep="")
    obs = x[,i]
    obsn = lib.size[i]
    #logR = log2((obs/obsn)/(ref/refn))
    logR = log2(obs/ref)			#log ratio of expression, accounting for libr size
    absE = (log2(obs/obsn) + log2(ref/refn))/2	#absolute expression
    v = (obsn-obs)/obsn/obs + (refn-ref)/refn/ref	 #estimated asymptotic variance
    
    #remove infinite values, cutoff based on A
    fin = is.finite(logR) & is.finite(absE) & (absE > Acutoff)
    logR = logR[fin]
    absE = absE[fin]
    v = v[fin]
    
    if(max(abs(logR)) < minlogR)
      return(list(f=1, fdiff=0)) # f[i] = 1
    
    #taken from the original mean() function
    n = length(logR)
    loL = floor(n * logratioTrim) + 1
    hiL = n + 1 - loL
    loS = floor(n * sumTrim) + 1
    hiS = n + 1 - loS
    
    #keep = (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
    #a fix from leonardo ivan almonacid cardenas, since rank() can return
    #non-integer values when there are a lot of ties
    keep = (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)
    
    if(doWeighting) {
      fi = sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
    } else { fi = mean(logR[keep], na.rm=TRUE) } #f[i] = mean(logR[keep], na.rm=TRUE) }
    
    #Results will be missing if the two libraries share no features with positive counts
    #In this case, return unity
    #if(is.na(f[i])) f[i] = 0
    if(is.na(fi)) fi = 0
    
    #check if close to peak; if not, switch to peak
    d = density(log2((obs)/ref), na.rm=T)
    p = as.matrix(findpeaks(d$y)); if(ncol(p)==1) p = t(p)
    p1 = d$x[p[which.max(p[,1]),2]]
    #fdiff[i] = p1-f[i]
    fdiffi = p1-fi
    
    if (plotimg) {
      pngname = pngnames[i]
      png (file=pngname , width=700, height=1800)
      par(mfrow=c(3,1), mar=(c(5, 5, 4, 2) + 0.1))
      
      #plot(d); abline(v=f[i], col="red"); abline(v=p1, col="blue"); 
      plot(d); abline(v=fi, col="red"); abline(v=p1, col="blue"); 
      
      plot((x0[,i]+x0[,refColumn])/2, log(x0[,i]/x0[,refColumn]), cex=.5, main=paste(mains[i],": f=",fi, sep=""))
      #abline(h=f[i], col="red")
      abline(h=fi, col="red")
    }
    
    #if f[i] too far from peak
    #if (abs(f[i]-p1)>cutoff) {
    if (abs(fi-p1)>cutoff) {
      abline(h=p1, col="blue")
      #f[i] = p1
      fi = p1
    }
    
    #f[i] = 1/2^f[i]
    fi = 1/2^fi
    
    #plot((matrixCount[,i]+matrixCount[,refColumn])/2, log2((matrixCount[,i]*f[i])/matrixCount[,refColumn]), cex=.5, main=paste("AFTER CHANGE: mean count vs. log2 fold change: ", sampleMeta$gene[i]," over refColumn ", sampleMeta$gene[refColumn],": f=",f[i], sep=""))
    if (plotimg) {
      plot((x0[,i]+x0[,refColumn])/2, log((x0[,i]*fi)/x0[,refColumn]), cex=.5, main=paste(mains[i],": f=",fi, sep=""))
      abline(h=0, col="red")
      dev.off()
    }
    
    return(list(f=fi, fdiff=fdiffi))
    
    
  }
  #multiple of 1
  rm(x)
  
  f = rep(NA,length(ff))
  fdiff = rep(NA,length(ff)) # diff between density peak and value (note: logged)
  for (i in 1:length(ff)) {
    f[i] = ff[[i]]$f
    try({ fdiff[i] = ff[[i]]$fdiff })
  }
  
  return(list(f=f,fdiff=fdiff))
}


## input: matrix
## output: if features (rownames) have +/- symbols, returns corresponding feature layers; else returns NULL
cell_type_layers = function(phen) sapply(str_split(phen,"[+-]+"), function(x) sum(x!=""))


## input: Phenotype
## output: Phenotype, PhenoCode, Phenolayer (number of markers per phenotype)
getPhen = function(phen, phenocode=NULL) {
  require(stringr)
  pm = data.frame(phenotype=phen)
  markers = unique(unlist(str_split(phen,"[-+]")))
  markers = markers[markers!=""]
  if (is.null(phenocode)) {
    pm$phenocode = laply(phen, function(x){
      if (x=="") return(paste0(rep(0,length(markers)), collapse=""))
      b = str_split(x,"[+-]+")[[1]]
      b = b[-length(b)]
      bo = match(markers,b)
      phec = as.vector(laply(str_extract_all(x,"[+-]+"), 
                              function(y) str_count(y,"[+]"))+1)
      c = phec[bo]
      c[is.na(c)] = 0
      paste0(c,collapse="")
    })
  } else {
    pm$phenocode = phenocode
  }
  pm$phenolayer = cell_type_layers(phen)
  pm$phenotype = as.character(pm$phenotype)
  
  return(pm)
}


## input: cell pop names or getPhen output meta_cell
## output: child (pos/neg) and parent list
getPhenCP = function(meta_cell=NULL, cp=NULL, no_cores=1) {
  require(plyr)
  parl = parallel_backend(no_cores)
  if (parl) if (no_cores>detectCores()) no_cores = detectCores()
  
  if (is.null(cp) & is.null(meta_cell)) { 
    print("give me something!"); return(NULL) 
  } else if (is.null(meta_cell)) {
    meta_cell = getPhen(cp)
  }
  meta_cell_grid = do.call(rbind, str_split(meta_cell$phenocode,""))
  meta_cell_grid = Matrix(apply(meta_cell_grid,2,as.numeric),sparse=T)
  allcolu = llply(1:ncol(meta_cell_grid), function(j) unique(meta_cell_grid[,j]))
  allcol = llply(1:length(allcolu), function(ci) {
    a = llply(allcolu[[ci]], function(ui) 
      meta_cell_grid[,ci]==ui)
    names(a) = allcolu[[ci]]
    a
  })
  
  pchild = list(meta_cell$phenotype[meta_cell$phenolayer==1])
  names(pchild) = ""
  pparen = llply(1:length(pchild[[1]]), function(x) "")
  names(pparen) = pchild[[1]]
  
  jjl = sort(unique(meta_cell$phenolayer))
  if (length(jjl)>2) {
    jjl = jjl[jjl>0]
    jj_inds = llply(jjl, function(l) which(meta_cell$phenolayer==as.numeric(l)))
    cat("\n")
    for (jjli in 2:length(jj_inds)) {
      start1 = Sys.time()
      
      meta_cell_ = meta_cell[jj_inds[[jjli-1]],,drop=F]
      meta_cell__ = meta_cell[jj_inds[[jjli]],,drop=F]
      meta_cell_grid_ = meta_cell_grid[jj_inds[[jjli-1]],,drop=F]
      meta_cell_grid__ = meta_cell_grid[jj_inds[[jjli]],,drop=F]
      
      cat("- ", nrow(meta_cell_), " pops @ layer ", jjli-1)
      allcol__ = llply(allcol, function(x) 
        llply(x, function(y) y[jj_inds[[jjli]]]))
      allcol_ = llply(allcol, function(x) 
        llply(x, function(y) y[jj_inds[[jjli-1]]]))
      
      #child
      loop_ind = loop_ind_f(1:nrow(meta_cell_), no_cores)
      pchildl = llply(loop_ind, function(jj) 
        llply(jj, function(j) {
          colj1 = which(meta_cell_grid_[j,]>0)
          mcgrow = meta_cell_grid_[j,]
          chi = Reduce("&", llply(colj1, function(coli) 
            allcol__[[coli]][[as.character(mcgrow[coli])]]) )
          meta_cell__$phenotype[chi]
        }), .parallel=parl)
      pchildl = unlist(pchildl,recursive=F)
      names(pchildl) = meta_cell_$phenotype
      pchildl = plyr::compact(pchildl)
      pchild = append(pchild, pchildl)
      
      #paren
      loop_ind = loop_ind_f(1:nrow(meta_cell__), no_cores)
      pparenl = llply(loop_ind, function(jj) 
        llply(jj, function(j) {
          colj1 = which(meta_cell_grid__[j,]>0)
          mcgrow = meta_cell_grid__[j,]
          chidf = do.call(cbind, llply(colj1, function(coli) 
            allcol_[[coli]][[as.character(mcgrow[coli])]]) )
          chi = apply(chidf,1,function(x) sum(!x)==1)
          meta_cell_$phenotype[chi]
        }), .parallel=parl)
      pparenl = unlist(pparenl,recursive=F)
      names(pparenl) = meta_cell__$phenotype
      pchildl = plyr::compact(pparenl)
      pparen = append(pparen, pparenl)
      
      time_output(start1)
    }
  }
  edf = do.call(rbind, llply(1:length(pchild), function(x)
    data.frame(from=names(pchild)[x], to=pchild[[x]])))
  
  return(list(pchild=pchild, pparen=pparen, edf=edf))
}


## input: gr$v $e from flowgraph object, and a layout function from igraph
## output: gr$v $e with x and y coordinates
set_layout_graph = function(gr,FUN=layout.reingold.tilford) { # layout.circle
  # FUN is a layout function from the igraph package
  # assume graph is connected, used internally
  require(igraph)
  
  gr_e = gr$e
  gr_v = gr$v
  
  # edit layout
  gr0 = graph_from_data_frame(gr_e)
  gr_vxy_ = FUN(gr0)
  gr_vxy = as.data.frame(gr_vxy_)
  if (as.character(substitute(FUN))=="layout.reingold.tilford") {
    ## edit layout manually
    gys = sort(unique(gr_vxy[,2]))
    gxns = sapply(gys, function(y) sum(gr_vxy[,2]==y))
    maxlayern = max(gxns)
    maxlayer = gys[which.max(gxns)]
    maxlayertf = gr_vxy[,2]==maxlayer
    gr_vxy[maxlayertf,1] = rank(gr_vxy[maxlayertf,1])-1
    for (gy in gys) {
      if (gy==maxlayer) next()
      layertf = gr_vxy[,2]==gy
      gr_vxy[layertf,1] = 
        rank(gr_vxy[layertf,1]) - 1 + floor((maxlayern-sum(layertf))/2)
    }
    # turn plot sideways
    gr_vxy = gr_vxy[,2:1]
    gr_vxy[,1] = max(gr_vxy[,1])-gr_vxy[,1]
  }
  
  # get node
  colnames(gr_vxy) = c("x","y")
  gr_v_xy = gr_vxy[match(gr_v$phenotype,names(V(gr0)[[]])),]
  gr_v$x = gr_v_xy$x
  gr_v$y = gr_v_xy$y
  
  # get edge
  gr_e$from.x = gr_v$x[match(gr_e$from, gr_v$phenotype)]
  gr_e$from.y = gr_v$y[match(gr_e$from, gr_v$phenotype)]
  gr_e$to.x = gr_v$x[match(gr_e$to, gr_v$phenotype)]
  gr_e$to.y = gr_v$y[match(gr_e$to, gr_v$phenotype)]
  
  return(list(e=gr_e,v=gr_v))
}


# change layout
set_layout = function(object, layout_fun=layout.reingold.tilford) {
  object@plot_layout = as.character(substitute(layout_fun))
  object@graph = set_layout_graph(object@graph, layout_fun)
  return(object)
}


## flowgraph initialization -----------------------------------
flowgraph = function(input_, meta=NULL, no_cores=1, 
           layout_fun=layout.reingold.tilford, 
           markers=NULL,
           cumsumpos=F, # whether to make positives + = +/++ (cumsum)
           # prop=T,
           prop=T,
           specenr=T, 
           normalize=T, 
           norm_ind=0, norm_layer=3, norm_path=NULL, ...) { # layout.circle
  start = Sys.time()
  
  require(Matrix)
  require(stringr)
  require(plyr)
  parl = parallel_backend(no_cores)
  if (parl) if (no_cores>detectCores()) no_cores = detectCores()
  
  ## input -----------------------------------------------
  ## meta data frame must have an "id" column for samples names
  ## input can be 
  # - a matrix or vector of the counts from a Phenotypes object but the names/colnames must be the phenotype
  # - a Phenotypes object from the flowType package
  # - a list of Phenotypes objects
  # - a vector of flowtype object paths (full path)
  # - a matrix of raw flowType count values
  # all Phenotype objects should have been made using the same markers and partitionsPerMarker
  ## norm_ind is the sample index used to normalize count
  # - 0 means we define the sample with the median total cell count as the reference sample
  # - NULL means don't normalize
  ## norm_layer specifies which cell populations (on which layer(s)) is used to normalize count; 
  # - norm_layer=NULL means use all layers
  
  start1 = Sys.time()
  cat("preparing input ")
  
  # convert input_ into a list of Phenotype objects input
  msg = "make sure input type is correct"
  
  phenocode = NULL
  if (class(input_)=="integer" | class(input_)=="matrix") {
    ## feature: count (sample x cell population)
    if (is.null(dim(input_))) {
      mc = mc0 = matrix(input_,nrow=1)
      colnames(mc) = names(input_)
      rownames(mc) = "s1"
    } else {
      mc = mc0 = input_
      if (is.null(rownames(input_))) rownames(input_) = paste0("s",1:nrow(input_))
      rownames(mc) = rownames(input_)
    }
    
    ## make meta for samples
    if (!is.null(meta)) rownames(mc) = meta$id
    sample_id = rownames(mc)
    
    ## make meta for cell populations
    if (!any(grepl("[+-]",colnames(mc)))) {
      print(msg)
      return(NULL)
    } 
    phenotype = colnames(mc)
    if (is.null(markers)) {
      markers = unique(unlist(str_split(phenotype,"[-+]")))
      markers = markers[markers!=""]
    }
    
  } else {
    require(flowType)
    badinput = F
    if (class(input_)=="Phenotypes") {
      ftl = list(input=input_)
    } else if (class(input_)=="list") {
      testclass = laply(input_, function(x) class(x)=="Phenotypes")
      if (!all(testclass)) badinput = T
      ftl = input_
    } else if (class(input_)=="character") {
      ftl = llply(input_, function(x) get(load(x)), .parallel=parl)
      names(ftl) = laply(str_split(input_,"/"), function(x) x[length(x)])
    } else {
      badinput = T
    }
    if (badinput) {
      print(msg)
      return(NULL)
    }
    
    # make sample_id
    if (!is.null(meta)) names(ftl) = meta$id
    sample_id = names(ftl)
    
    ## make meta for cell populations
    if (is.null(markers))
      markers = ftl[[1]]@MarkerNames
    phenotype = NULL
    try({
      phenotype = laply(ftl[[1]]@PhenoCodes, function(x) 
        decodePhenotype(x, markers, ftl[[1]]@PartitionsPerMarker))
    }, silent=T)
    if (is.null(phenotype))
    try({
      phenotype = rownames(ftl[[1]]@MFIs)
    }, silent=T)
    if (is.null(phenotype)) {
      print("no phenotype cell populations labels in Phenotype file")
      return(NULL)
    }
    
    try ({
      phenocode = ftl[[1]]@PhenoCodes
    }, silent=T)
    
    ## feature: count (sample x cell population)
    mc0 = as.matrix(ldply(ftl, function(ft) ft@CellFreqs))
    if (class(mc0[1])=="character"){
      mc0 = as.matrix(mc0[,-1])
      class(mc0) = "numeric"
    } 
    mc = mc0
    
    # remove Phenotype list to save memory
    rm(ftl)
    gc()
  }
  
  time_output(start1)
  
  
  ## feature: count (sample x cell population) -----------
  start1 = Sys.time()
  cat("preparing feature(s): raw count; mapping out cell population relationships ")
  
  colnames(mc) = phenotype
  rownames(mc) = sample_id
  keepinds = apply(mc, 2, function(x) any(x>0))
  mc = Matrix(mc[,keepinds,drop=F], sparse=T) # remove all 0 columns
  phenotype = colnames(mc)
  if (!is.null(phenocode)) phenocode = phenocode[keepinds]
  
  ## make meta for cell populations FINAL
  # markern = length(markers)
  meta_cell = getPhen(phenotype,phenocode)
  
  # make parent list (for each cell popultion, list its parents)
  pccell = getPhenCP(meta_cell=meta_cell, no_cores=no_cores)
  edf = pccell$edf
  pchild = pccell$pchild # not used, returned
  pparen = pccell$pparen
  
  # extract cell populations that have parents to calculate expected proportions for; just to be safe
  cells1 = meta_cell$phenotype[meta_cell$phenolayer==1]
  cells1 = append("",cells1[cells1%in%names(pparen)])
  cells = meta_cell$phenotype[meta_cell$phenolayer>1]
  cells = cells[cells%in%names(pparen)]
  cells_ = append(cells1,cells)
  
  meta_cell = meta_cell[match(cells_,meta_cell$phenotype),]
  root = which(meta_cell$phenolayer==0)
  
  mc = mc[,match(cells_, colnames(mc))] # finalize cell populations
  
  # trim/order list of parents and children
  pparen_ = pparen[names(pparen)%in%cells_]
  pparen_[cells] = llply(pparen_[cells], function(x) {
    a = x[x%in%cells_]
    if (length(a)==0) return(NULL)
    return(a)
  })
  pparen_ = plyr::compact(pparen_)
  
  pchild_ = pchild[names(pchild)%in%cells_]
  pchild_ = llply(pchild_, function(x) {
    a = x
    for (i in 1:length(a)) 
      a[[i]] = x[[i]][x[[i]]%in%cells_]
    if (all(sapply(a, length)==0)) return(NULL)
    return(a)
  })
  pchild_ = plyr::compact(pchild_)
  
  edf = edf[edf$to%in%cells_ & edf$from%in%cells_,,drop=F]
  
  time_output(start1)
  
  ## list of outputs

  gr = list(e=edf,
            v=meta_cell)
  gr = set_layout_graph(gr, layout_fun) # layout cell hierarchy
  
  if (is.null(meta)) meta = data.frame(id=sample_id)
  fg = new("flowgraph", 
           feat=list(node=list(count=mc), edge=list()), 
           markers=markers,
           graph=gr, edge_list=list(child=pchild_,parent=pparen_), 
           meta=meta, plot_layout=as.character(substitute(layout_fun)), etc=list(cumsumpos=F))
  
  if (!any(grepl("3",meta_cell$phenocode))) 
    cumsumpos = cumsumpos & T
  if (cumsumpos)
    fg = flowgraph_cumsum(fg, no_cores=no_cores)
  
  if (prop) { # included in specenr
    fg = flowgraph_prop(fg)
    fg = flowgraph_prop_edge(fg, no_cores=no_cores)
  }
  
  if (normalize)
    fg = flowgraph_normalize(fg, norm_ind=norm_ind, norm_layer=norm_layer, norm_path=norm_path, no_cores=no_cores, ...)
  
  if (specenr)
    fg = flowgraph_specenr(fg, no_cores=no_cores)
  
  time_output(start, "total time used")
  return(fg)
}


## feature: proportion
flowgraph_prop = function(fg, overwrite=F) {
  if (!overwrite & "prop"%in%names(fg@feat$node)) {
    return(fg)
  }
  
  start1 = Sys.time()
  cat("preparing feature(s): proportion ")
  
  mc = fg@feat$node$count
  mp = mc/mc[,colnames(mc)==""]
  dimnames(mp) = dimnames(mc)
  fg@feat$node$prop = mp
  
  time_output(start1)
  return(fg)
}


## feature: edge proportions
flowgraph_prop_edge = function(fg, no_cores=1, overwrite=F) {
  if (!overwrite & "prop"%in%names(fg@feat$edge)) {
    return(fg)
  }

  require(plyr)
  parl = parallel_backend(no_cores)
  if (parl) if (no_cores>detectCores()) no_cores = detectCores()
  
  start1 = Sys.time()
  cat("preparing feature(s): proportion on edges ")
  
  edf = fg@graph$e
  mc = as.matrix(fg@feat$node$count)
  rooti = which(colnames(mc)=="")
  loop_ind = loop_ind_f(1:nrow(edf), no_cores)
  childprop_ = do.call(cbind,llply(loop_ind, function(ii) {
    do.call(cbind,llply(ii, function(i) {
      pname = ifelse(edf$from[i]=="", rooti, edf$from[i])
      mc[,edf$to[i]]/mc[,pname]
    }))
  }, .parallel=parl))
  colnames(childprop_) = paste0(edf$from,"_",edf$to)
  childprop_[is.nan(as.matrix(childprop_))] = 0
  fg@feat$edge$prop = childprop_
  
  time_output(start1)
  return(fg)
}


## feature: specenr
flowgraph_specenr = function(fg, no_cores=1, overwrite=F) {
  if (!overwrite & "specenr"%in%names(fg@feat$node)) 
    return(fg)
  
  require(gsubfn)
  require(stringr)
  require(plyr)
  parl = parallel_backend(no_cores)
  if (parl) if (no_cores>detectCores()) no_cores = detectCores()
  
  meta_cell = fg@graph$v
  pparen = fg@edge_list$parent
  if (!"prop"%in%names(fg@feat$node))
    fg = flowgraph_prop(fg)
  mp = fg@feat$node$prop
  if (!"prop"%in%names(fg@feat$edge))
    fg = flowgraph_prop_edge(fg)
  ep = fg@feat$edge$prop
  rooti = which(colnames(mp)=="")
  
  ## feature: expected proportion
  
  cells1 = append("",meta_cell$phenotype[meta_cell$phenolayer==1])
  cells = meta_cell$phenotype[meta_cell$phenolayer>1]
  cellsn = meta_cell$phenolayer
  pparen = fg@edge_list$parent
  pchild = fg@edge_list$child
  ## calculate the expected proportions for cell populatins with positive marker conditions only
  expe1 = matrix(.5,nrow=nrow(mp), ncol=length(cells1), dimnames=list(rownames(mp),cells1))
  expe1[,1] = 1
  
  cells1_p = str_extract(cells1,"[+]+")
  maxp = max(nchar(cells1_p[!is.na(cells1_p)]))
  
  start1 = Sys.time()
  cat("preparing feature(s): expected proportion; ")
  
  cpind = which(!grepl("[-]",cells))
  loop_ind = loop_ind_f(cpind,no_cores)
  expecp = do.call(cbind, llply(loop_ind, function(ii) {
    do.call(cbind, llply(ii, function(ic) {
      cpop = i = cells[ic]
      il = cellsn[ic]
      
      pnames = pparen[[i]]
      parent = mp[,pnames,drop=F]
      parento = laply(1:nrow(parent), function(xi) order(parent[xi,]))
      if (is.null(dim(parento))) parento = matrix(parento,nrow=nrow(parent))
      grprnt = unique(unlist(pparen[pnames]))
      
      pedges = do.call(cbind,llply(1:length(pnames), function(gi) {
        gname = pparen[[pnames[gi]]]
        cns = paste0(gname,"_",pnames[gi])
        ep[,cns,drop=F]
      }))
      
      expect1 = rep(0, nrow(parent))
      for (colj in unique(parento[,1])) {
        jind = parento[,1]==colj
        
        parentjmin = parent[jind,colj]
        
        pnamej = gsub("[+]","[+]",pnames[colj])
        pnamej = gsub("[-]","[-]",pnamej)
        pedgesnj = pedges[jind,!grepl(pnamej,colnames(pedges)),drop=F]
        
        expect1[jind] = apply(pedgesnj,1,max) * parentjmin
      }
      return(expect1)
    }))
  }, .parallel=parl))
  expecp[is.nan(expecp)] = 0
  colnames(expecp) = cells[cpind]
  
  
  ## infer the expected proportion of cell populations with negative marker condition(s)
  expec0 = expec = as.matrix(cbind(expe1,expecp))
  cpopneg = setdiff(cells,colnames(expec))
  cpopnegl = cell_type_layers(cpopneg)
  
  # replaces whatever pattern only once; e.g. gsubfn("_", p, "A_B_C_")
  p = proto(i=1, j=1, fun=function(this, x) if (count>=i && count<=j) "*" else x) 
  
  csp = fg@etc$cumsumpos
  expecn = do.call(cbind,llply(sort(unique(cpopnegl)), function(lev) {
    cpopl = cpopneg[cpopnegl==lev]
    cpopnegno = str_count(cpopl,"[-]") # number of negative marker conditions
    
    # ## version 1; DELETE ONE VERSION
    # expec_temp = NULL
    # for (negno in sort(unique(cpopnegno))) {
    #   cpopi = cpopl[cpopnegno==negno]
    #   sibs = gsubfn("-", p, cpopi)
    #   for (sib in unique(sibs)) {
    #     sib_ = sapply(1:maxp, function(pn) 
    #       gsub("[*]", paste0(rep("+",pn),collapse=""), sib) )
    #     sib_ = sib_[sib_%in%cells]
    #     sib_sum = rowSums(mp[,sib_,drop=F])
    #     sibpars = pparen[[sib_[1]]]
    #     sibis = which(sibs==sib) # shouldn't be length(sibis)>1, but safe
    #     for (sibi in sibis) {
    #       parsin = intersect(sibpars, pparen[[cpopi[sibi]]])
    #       # if (length(sibpars)>1) #there must be only 1 common parent
    #       #   for (sibpar in sibpars)
    #       #     parsin = intersect(sibpar, parsin)
    #       expec_temp = cbind(expec_temp, mp[,parsin] - sib_sum)
    #       colnames(expec_temp)[ncol(expec_temp)] = cpopi[sibi]
    #     }
    #   }
    # }
    
    ## version 2
    expec_temp = do.call(cbind,llply(cpopl, function(cpi) {
      cpig = gsub("[+-]","",cpi)
      parnames = intersect(pparen[[cpi]],colnames(mp))
      a = NULL
      for (parname in parnames) {
        chilnames = unlist(pchild[[parname]])
        chilnames = chilnames[chilnames!=cpi]
        chilnamesg = gsub("[+-]","",chilnames)
        chilgreat = chilnamesg==cpig & chilnames%in%colnames(expec)
        if (any(chilgreat)) {
          children = chilnames[chilgreat]
          if (csp) 
            a = matrix(mp[,parname]-expec[,children[which.min(nchar(children))],drop=F],ncol=1)
          else 
            a = matrix(mp[,parname]-rowSums(expec[,children,drop=F]),ncol=1)
          
          break
        }
      }
      if (is.null(a)) {
        a = b = mp[,parnames,drop=F]
        if (ncol(b)>1) a = matrix(apply(b,1,min), ncol=1)
      }
      colnames(a) = cpi
      return(a)
    }))
    if (is.null(dim(expec_temp))) expec_temp = matrix(expec_temp,ncol=1)
    
    return(expec_temp)
  }))
  expec = cbind(expec,expecn)
  
  exp1 = cbind(expe1,expec[,match(cells,colnames(expec)),drop=F])
  a = mp/exp1
  aa = as.matrix(a)
  a[is.infinite(aa)] = max(a[is.finite(aa)])
  a[is.nan(aa)] = 0
  specenr1 = log(a)
  specenr1[is.nan(as.matrix(specenr1))] = 0
  suppressWarnings({ specenr1[exp1==0] = log(mp[exp1==0]) })
  specenr1[mp==0] = 0
  
  fg@feat$node$expect_prop = exp1
  fg@feat$node$specenr = specenr1
  
  if ("count_norm"%in%names(fg@feat$node)) {
    mca = fg@feat$node$count_norm
    expc1 = exp1*mca[,1]
    fg@feat$node$expect_count_norm = expc1
  }
  
  return(fg)
  
  time_output(start1)
}


## feature(optional): normalized (adjusted) count 
flowgraph_normalize = function(fg, norm_ind=0, norm_layer=3, norm_path=NULL, no_cores=1, overwrite=F, ...) {
  require(plyr)
  parl = parallel_backend(no_cores)
  if (parl) if (no_cores>detectCores()) no_cores = detectCores()
  
  if (!overwrite & "count_norm"%in%names(fg@feat$node)) {
    return(fg)
  }
  
  mca = mc = fg@feat$node$count
  fdiff0 = rep(0,nrow(mc))
  f0 = rep(1,nrow(mc))
  rooti = which(colnames(mc)=="")
  
  if (nrow(mca)>1 & !is.null(norm_ind)) {
    start1 = Sys.time()
    cat("preparing feature(s): normalized count ")
    
    sample_id = fg@meta$id
    meta_cell = fg@graph$v
    maxlayer = max(meta_cell$phenolayer)
    
    # prepare feat_file_cell_counts
    x = x0 = as.matrix(mc)[,-rooti,drop=F] # take out total cell count
    maxx = max(x0[is.finite(x0)])
    rootc = mc[,rooti]
    if (norm_ind==0) {
      norm_ind = which.min(abs( rootc-median(rootc) )) # reference column: median total count out of all control files
      if ("class"%in%colnames(fg@meta))
        if ("control"%in%fg@meta$class) {
          rootcc = which(fg@meta$class=="control")
          norm_ind = rootcc[which.min(abs( rootc[rootcc]-median(rootc[rootcc]) ))]
        }
    } 
    
    # extract cell populations that would define TMM (layer/count)
    if (!is.null(norm_layer)) {
      norm_layer[norm_layer>maxlayer] = maxlayer
      x = x[,colnames(x0)%in%meta_cell$phenotype[meta_cell$phenolayer%in%norm_layer], drop=F]
    }
    
    # prepare paths to plot to
    pngnames = NULL
    dir.create(norm_path, showWarnings=F, recursive=T)
    if (!is.null(norm_path)) {
      pngnames = paste0(norm_path,"/", sample_id,".png")
      mains = paste0("mean count vs. ln fold change over ref sample ", 
                     sample_id[norm_ind], " on layer ", 
                     ifelse(is.null(sample_id[1]),"all",paste(norm_layer,collapse="-")))
    }
    
    # calculate absolute count TMM
    fresult = tmm(x,x0,rootc,norm_ind,cutoff=Inf,pngnames=pngnames,mains=mains,no_cores=no_cores,samplesOnCol=F, ...)
    f0 = fresult$f
    fdiff0 = fresult$fdiff
    
    mca = do.call(rbind,llply(c(1:nrow(mc)), function(x) mc[x,]*f0[x]))
    dimnames(mca) = dimnames(mc)
    
    # plot difference between TMM and peak for all files
    if (!is.null(pngnames[1])) {
      png(paste0(norm_path,"/all.png") , width=700, height=700)
      plot(sort(abs(fdiff0)), cex=.4, ylim=c(0,3), main="cell-count-norm-factor_f_diff_from_peak_abs")
      lines(sort(abs(fdiff0)), col="blue")
      graphics.off()
    }
    
    fg@etc$count_norm_factor = f0
    fg@etc$count_norm_factor_diffpeak = fdiff0
    fg@feat$node$count_norm = mca
    
    if ("specenr"%in%names(fg@feat$node))
      fg@feat$node$expect_count_norm = fg@feat$node$specenr*mca[,1]
    
    time_output(start1)
  } # return f0 and fdiff0
  
  return(fg)
}


## feature(s): get mean of a class for each feature
flowgraph_mean_class = function(fg, class, no_cores=1, node_features=NULL, edge_features=NULL) {
  parl = parallel_backend(no_cores)
  if (parl) if (no_cores>detectCores()) no_cores = detectCores()
  
  if (!class%in%colnames(fg@meta)) {
    cat("invalid class name, choose one from @meta\n")
    return(fg)
  }
  if (length(unique(fg@meta[,class]))<2) return(fg)
  label = paste0("_MEAN_",class)
  
  meandiff = function(m0, classcol) {
    m = m0
    for (pi in unique(classcol)) {
      pii = classcol==pi
      m_ = m[pii,,drop=F]
      m_m = colMeans(m_)
      m[pii,] = do.call(rbind, llply(1:nrow(m_), function(i) m_[i,]-m_m ))
    }
    dimnames(m) = dimnames(m0)
    return(m)
  }
  
  if (is.null(node_features)) node_features = names(fg@feat$node)
  node_features = node_features[node_features%in%names(fg@feat$node)]
  node_features = node_features[!grepl("MEAN", names(node_features))]
  if (length(node_features)>0) {
    fg_nodefs0 = fg@feat$node[node_features]
    fg_nodefs = llply(fg_nodefs0, function(m0) 
      meandiff(m0, fg@meta[,class]), .parallel=parl)
    names(fg_nodefs) = paste0(names(fg_nodefs),label)
    fg@feat$node = append(fg@feat$node, fg_nodefs)
  }
  
  if (is.null(edge_features)) edge_features = names(fg@feat$edge)
  edge_features = edge_features[edge_features%in%names(fg@feat$edge)]
  edge_features = edge_features[!grepl("MEAN", names(edge_features))]
  if (length(edge_features)>0) {
    fg_edgefs0 = fg@feat$edge[edge_features]
    fg_edgefs = llply(fg_edgefs0, function(m0) 
      meandiff(m0, fg@meta[,class]), .parallel=parl)
    names(fg_nodefs) = paste0(names(fg_edgefs),label)
    fg@feat$edge = append(fg@feat$edge, fg_edgefs)
  }
  
  return(fg)
}


## get p values 
flowgraph_p = function(
  fg, no_cores=1, class, control=NULL,
  node_features=NULL, edge_features=NULL, 
  overwrite=F, 
  test_name="t_BY", # name your statistical test / adjustment method; if same name, overwrite=T will overwrite, overwrite=F will return nothing
  diminish=F, # don't test if all parents are insignificant, stricter the lower the layer
  p_thres=.05, p_rate=2, # only used if diminish=T
  test=function(x,y) tryCatch(t.test(x,y)$p.value, error=function(e) 1), 
  adjust=function(x) p.adjust(x, method="BY")
) {
  # if control is NULL, compare all
  # if node/edge_features is NULL, do all
  parl = parallel_backend(no_cores)
  if (parl) if (no_cores>detectCores()) no_cores = detectCores()
  
  if (!class%in%colnames(fg@meta)) {
    cat("invalid class, choose one from @meta\n")
    return(fg)
  }
  
  classes_ = fg@meta[,class]
  classes = unique(classes_)
  if (length(classes)<2) {
    cat("invalid class, choose one from @meta\n")
    return(fg)
  }
  
  # class list lists what classes to compare
  if (!is.null(control)) {
    if (!control%in%fg@meta[,class]) {
      cat("invalid control, choose one from @meta\n")
      return(fg)
    }
    class_list = llply(classes[classes!=control], function(x) c(control, x))
  } else {
    class_list = llply(2:length(classes), function(i) 
      llply(1:(i-1), function(j) c(classes[i], classes[j])))
  }
  
  cat("\n")
  for (ftype in c("node","edge")) {
    if (ftype=="node") features = node_features
    if (ftype=="edge") features = edge_features
    
    if (is.null(features)) features = names(fg@feat[[ftype]])
    features = features[features%in%names(fg@feat[[ftype]])]
    if (length(features)==0) next
    
    for (feature in features) {
      if (!feature%in%names(fg@feat[[ftype]])) next
      for (classl in class_list) {
        desc = fg@feat_summary$desc[[ftype]]
        
        idsname = paste0(classl,collapse="_")
          if (!overwrite & !is.null(desc[[ftype]][[feature]][[test_name]][[class]][[idsname]])) next

        start1 = Sys.time()
        cat("- claculating (", class, ":", paste0(classl), ")",
            test_name,"for", ftype, "feat", feature)
        
        id1 = classes_==classl[1]
        id2 = classes_==classl[2]
        ids = list(fg@meta$id[id1],fg@meta$id[id2]); 
        names(ids) = classl
        fg@feat_summary$desc[[ftype]][[feature]][[test_name]][[class]][[idsname]] = ids

        m = as.matrix(fg@feat[[ftype]][[feature]])
        m1 = m[id1,,drop=F]
        m2 = m[id2,,drop=F]
        
        if (!diminish) {
          loop_ind = loop_ind_f(1:ncol(m), no_cores)
          p = unlist(llply(loop_ind, function(ii) 
            llply(ii,function(i) test(m1[,i], m2[,i])), .parallel=parl ))
          p[is.nan(p)] = 1
          if (!is.null(adjust)) p = adjust(p)
          names(p) = colnames(m)
        } else {
          pparen = fg@edge_list$parent
          if (ftype=="node") layers_ = fg@graph$v$phenolayer
          if (ftype=="edge") {
            phen_to = fg@graph$e$to
            phen_from = fg@graph$e$from
            layers_ = fg@graph$v$phenolayer[match(phen_to, fg@graph$v$phenotype)]
          } 
          layers = sort(unique(layers_))
          root = which(layers==0) # should be 1
          if (ftype=="node") loop_ind_ = llply(layers[-root], function(l)
            return(which(layers_==l)))
          if (ftype=="edge") loop_ind_ = llply(layers, function(l) which(layers_==l))
          
          p = rep(NA,ncol(m))
          if (ftype=="node") {
            root_ = which(colnames(m)=="")
            p[root_] = test(m1[,root_], m2[,root_])
            if (is.nan(p[root_])) p[root_] = 1
          }
          names(p) = colnames(m)
          
          if (ftype=="node") p_thress = rep(p_rate, length(layers)-1)
          if (ftype=="edge") p_thress = rep(p_rate, length(layers))
          p_thress[1] = p_thres
          p_thress = cumprod(p_thress)
          p_thress = sort(p_thress, decreasing=T)
          
          for (lvl in layers[-root]) {
            loop_ind = loop_ind_f(loop_ind_[[lvl]], no_cores)
            pl = unlist(llply(loop_ind, function(ii) {
              if (ftype=="node")
                return(llply(ii, function(i) {
                  pheni = fg@graph$v$phenotype[i]
                  parnis = pparen[[pheni]]
                  ppar = p[parnis]
                  if (lvl==1) ppar = p[root_]
                  if (all(ppar>p_thress[lvl])) return(NA)
                  pt = test(m1[,i], m2[,i])
                  if (is.nan(pt)) return(1)
                  return(pt)
                }))
              if (ftype=="edge")
                return(llply(ii, function(i) {
                  pheni = fg@graph$v$phenotype[i]
                  parnis = phen_to==phen_from[i]
                  if (sum(parnis)>0)
                    if (all(p[parnis]>p_thress[lvl])) return(NA)
                  pt = test(m1[,i], m2[,i])
                  if (is.nan(pt)) return(1)
                  return(pt)
                }))
            }, .parallel=parl
            ))
            pai = !is.na(pl)
            pail = sum(pai)
            if (pail>1)
              if (!is.null(adjust)) 
                pl[pai] = adjust(rep(pl[pai],sum(pai)))[1:pail]
            pl[!pai | is.nan(pl)] = 1
            p[loop_ind_[[lvl]]] = pl
          }
        }
        fg@feat_summary[[ftype]][[feature]][[test_name]][[class]][[idsname]][["p"]] = p
        fg@feat_summary[[ftype]][[feature]][[test_name]][[class]][[idsname]][["m1"]] = colMeans(as.matrix(m1))
        fg@feat_summary[[ftype]][[feature]][[test_name]][[class]][[idsname]][["m2"]] = colMeans(as.matrix(m2))
        
        time_output(start1)
      }
    }
  }
  
  return(fg)
}

flowgraph_clear_p = function(fg) {
  fg@feat_summary = list()
  return(fg)
}


## flowgraph_cumsum ----------------------------
flowgraph_cumsum = function(fg, no_cores) {
  # check if already cumsum
  if (fg@etc$cumsumpos) return(fg)
  
  # check if do-able (there exists multple ++)
  if (!any(grepl("3",fg@graph$v$phenocode))) return(fg)
  
  # phenotype=c("","A+","A++","A-","B-","B+","B++","A-B-","A-B+","A-B++","A+B-","A+B+","A+B++","A++B-","A++B+","A++B++")
  # mc=matrix(c(90,30,30,30,30,30,30,10,10,10,10,10,10,10,10,10),nrow=1)
  # colnames(mc)=phenotype
  # meta_cell = getPhen(phenotype)
  # pparen = getPhenCP(meta_cell=meta_cell)$pparen
  
  mc = fg@feat$node$count
  pparen = fg@edge_list$parent
  meta_cell = fg@graph$v
  meta_cell_grid = do.call(rbind, str_split(meta_cell$phenocode,""))
  meta_cell_grid = apply(meta_cell_grid, 2, as.numeric)
  meta_cell_grid = Matrix(meta_cell_grid,sparse=T)
  
  meta_cell_gridTF = apply(meta_cell_grid,2, function(x) {
    xmax = max(x)
    if (xmax==2) return(rep(F,length(x)))
    x>1 & x<xmax
  })
  meta_cell_gridTF_sum = apply(meta_cell_gridTF,1,sum)
  meta_cell_gridTF_any = meta_cell_gridTF_sum>0
  ccol = which(apply(meta_cell_gridTF, 2, any))
  
  for (marker in ccol) {
    coldo = which(meta_cell_gridTF[,marker])
    coldo = coldo[order(meta_cell_grid[coldo,marker], decreasing=T)]
    loop_ind = loop_ind_f(coldo, no_cores)
    mc[,coldo] = do.call(cbind,llply(loop_ind, function(jj) 
      do.call(cbind,llply(jj, function(j) {
      mcgi = mcgis = mcgip = meta_cell_grid[j,]
      mcgis[marker] = mcgis[marker]+1
      jsib = meta_cell$phenocode==paste0(mcgis,collapse="")
      mc[,j] + mc[,jsib]
    })), .parallel=T))
  }
  fg@feat$node$count_original = fg@feat$node$count
  fg@feat$node$count = mc
  fg@etc$cumsumpos = T
  return(fg)
}


## flowgraph interpretation functions ------------------------------

setMethod(
  "show", "flowgraph",
  function(object) 
    cat("flowgraph object with ", length(object@feat$node)," cell population and ", length(object@feat$edge), " edge feature(s)", 
        "\n- markers: ", paste0(object@markers, collapse=", "),
        "\n- contains: ", nrow(object@graph$v), " cell populations and ", nrow(object@graph$e), " edges\n", sep="")
)


setMethod(
  "summary", "flowgraph",
  function(object) {
    show(object)
    cat("\n")
    
    summary_table = function(m, feat_type) {
      m = as.matrix(m)
      data.frame(feat=feat_type, nrow=nrow(m), ncol=ncol(m), 
                 inf=sum(is.infinite(m)), neginf=sum(m==-Inf), na=sum(is.na(m)), nan= sum(is.nan(m)), 
                 neg=sum(m<0), pos=sum(m>0), zero=sum(m==0), max=max(m[is.finite(m)]), min=min(m[is.finite(m)]))
      
    }
    
    result1 = ldply(names(object@feat$node), function(x) 
      summary_table(object@feat$node[[x]],x))
    result1 = data.frame(type=rep("node",nrow(result1)),result1)
    
    result2 = NULL
    if (length(object@feat$edge)>0) {
      result2 = ldply(names(object@feat$edge), function(x)
        summary_table(object@feat$edge[[x]],x))
      result2 = data.frame(type=rep("edge",nrow(result2)),result2)
    }
    
    tab = rbind(result1,result2)
    return(tab)
  }
)

ggdf = function(gr0) {
  list(e=data.frame(gr0$e, width=1,color="unchanged", e_ind=F),
       v=data.frame(gr0$v, 
                    size=1, color="", sizeb=1, colorb="", fill="", 
                    label=gr0$v$phenotype, 
                    label_ind=F, v_ind=F, vb_ind=F))
}

flowgraph_summary_plot = function(
  fg,
  method, class, 
  nodeft="specenr", edgeft="prop", 
  nodeftlabel="prop", label_max=30,
  p_thres=.01, show_bgedges=T,
  path, width=9, height=9) { # width in inches feat must be list with names node and edge
  # cellhierarchy plots p values
  
  require(ggplot2)
  require(ggrepel)
  
  sum = fg@feat_summary
  desc = sum$desc
  gr0 = ggdf(fg@graph)
  if (!is.null(nodeftlabel)) 
    if (nodeftlabel!=nodeft) 
      lft = T
  
  for (idsname in names(desc$node[[nodeft]][[method]][[class]])) {
    if (is.null(idsname)) next
    dsc = desc$node[[nodeft]][[method]][[class]][[idsname]]
    pms = sum$node[[nodeft]][[method]][[class]][[idsname]]
    
    p = pms$p
    p_ = p<p_thres
    
    m1 = pms$m1; #names(dsc)[1]
    m2 = pms$m2; #names(dsc)[2]
    
    main = paste0("feat: ",nodeft,"; class: ",class,
                  "; pthres: ",p_thres,
                  "\nsize = -ln(p value)",
                  "\nlabel(",names(dsc)[1],"/",names(dsc)[2],") = ",nodeft)
    
    gr = gr0
    gr$v$label_ind = gr$v$v_ind = p_
    if (!is.null(label_max)) {
      if (sum(p_)>label_max) {
        gr$v$label_ind = rep(F, nrow(gr$v))
        gr$v$label_ind[order(p,decreasing=T)[1:30]] = T
      } 
    }
    gr$e$e_ind = gr$e[,1]%in%gr$v$name[p_] & gr$e[,2]%in%gr$v$name[p_]
    gr$v$color = ifelse(m2>m1,"increased","reduced")
    gr$v$label = paste0(gr$v$name,":",round(m1,3),"/",round(m2,3))
    if (lft) {
      pms_ = sum$node[[nodeftlabel]][[method]][[class]][[idsname]]
      m1_ = pms_$m1
      m2_ = pms_$m2
      gr$v$label = paste0(gr$v$label," (",round(m1_,3),"/",round(m2_,3),")")
      main = paste0(main,"(",nodeftlabel,")")
    }
    gr$v$size = -log(p)
    gr$v$size[is.infinite(gr$v$size)] = max(gr$v$size[!is.infinite(gr$v$size)])
    if (!is.null(edgeft)) {
      e1 = sum$edge[[nodeft]][[method]][[class]][[idsname]]$m1
      e2 = sum$edge[[nodeft]][[method]][[class]][[idsname]]$m2
      pe = sum$edge[[nodeft]][[method]][[class]][[idsname]]$p
      pe_ = pe<p_thres
      gr$e$color[pe_] = ifelse(e2>e1,"increased","reduced")
      main = paste0(main,"\nedge size = non(1)/sig(2) p values of ", edgeft)
    }
    
    gp = gggraph(gr, main=main, bgedges=show_bgedges)
    path_ = paste0(path,"/",ifelse(is.null(nodeftlabel),nodeft,paste0(nodeft,"_label-",nodeftlabel)),"/",class)
    dir.create(path_, recursive=T, showWarnings=F)
    ggsave(paste0(path_,"/",idsname,".png"), plot=gp, scale=1, width=width, height=height, units="in", dpi=500, limitsize=T)
  }
}




## graph plot functions
ggblank = function(gr_v=NULL) {
  require(ggplot2)
  if (is.null(gr_v)) {
    gp = ggplot()
  } else {
    gp = ggplot(gr_v,aes(x=x, y=y, colour=color))
  }
  
  gp +
    scale_x_continuous(expand=c(0,1)) +  # expand x limits
    scale_y_continuous(expand=c(0,1)) + # expand y limits
    # theme_bw()+  # use the ggplot black and white theme
    theme(
      axis.text.x = element_blank(),  # rm x-axis text
      axis.text.y = element_blank(), # rm y-axis text
      axis.ticks = element_blank(),  # rm axis ticks
      axis.title.x = element_blank(), # rm x-axis labels
      axis.title.y = element_blank(), # rm y-axis labels
      panel.background = element_blank(), 
      panel.border = element_blank(), 
      panel.grid.major = element_blank(),  #rm grid labels
      panel.grid.minor = element_blank(),  #rm grid labels
      plot.background = element_blank())
}

gggraph = function(gr, main="", bgedges=T) { # indices of whether to apply color size etc
  # gr_v: name x y label size color sizeb colorb
  # gr_e: from to from.x from.y to.x to.y color
  require(ggrepel)
  require(ggplot2)
  
  gr_v = gr$v
  gr_e = gr$e
  label_ind = gr_v$label_ind
  v_ind = gr_v$v_ind
  vb_ind = gr_v$vb_ind
  e_ind = gr_e$e_ind
  
  # base graph
  if (bgedges) { # keep greyed out edges on
    gp = ggblank() + ggtitle(main) +
      scale_fill_brewer(palette="Pastel2") +
      geom_segment(data=gr_e[!e_ind,], color="grey",
                   aes(x=from.x,xend=to.x, y=from.y,yend=to.y)) +
      geom_segment(data=gr_e[e_ind,], 
                   aes(x=from.x,xend=to.x, y=from.y,yend=to.y, color=color)) +
      # geom_point(data=gr_v[vb_ind,],aes(x=x,y=y, color=colorb),size=gr_v[vb_ind,"size"]+1) +
      # geom_point(data=gr_v[!v_ind,],aes(x=x,y=y), size=1, color="grey")+
      geom_point(data=gr_v[v_ind,],aes(x=x,y=y, color=color, size=size)) +
      geom_label_repel(data=gr_v[label_ind,],
                       aes(x=x,y=y,label=label, color=color),
                       nudge_x=-.1, direction="y", hjust=1, segment.size=0.2)
  } else {
    gp = ggblank() + ggtitle(main) +
      scale_fill_brewer(palette="Pastel2") +
      # geom_segment(data=gr_e[!e_ind,], color="grey",
      #              aes(x=from.x,xend=to.x, y=from.y,yend=to.y)) +
      geom_segment(data=gr_e[e_ind,], 
                   aes(x=from.x,xend=to.x, y=from.y,yend=to.y), color="grey50") +
      # geom_point(data=gr_v[vb_ind,],aes(x=x,y=y, color=colorb),size=gr_v[vb_ind,"size"]+1) +
      # geom_point(data=gr_v[!v_ind,],aes(x=x,y=y), size=1, color="grey")+
      geom_point(data=gr_v[v_ind,],aes(x=x,y=y, color=color, size=size)) +
      geom_label_repel(data=gr_v[label_ind,],
                       aes(x=x,y=y,label=label, color=color),
                       nudge_x=-.1, direction="y", hjust=1, segment.size=0.2)
  }

  return(gp)
}
