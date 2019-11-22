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
         slots=list(node_features="list", edge_features="list", 
                    markers="character", graph="list", edge_list="list", 
                    meta="data.frame", plot_layout="character", etc="list",
                    node_p="list", edge_p="list", node_p_class="list", edge_p_class="list"))


## flowgraph class modification ---------------

# replace marker names
gsub_markers = function(object, markers_new, markers_old=NULL) {
  require(plyr)
  
  if (is.null(markers_old)) {
    if (length(markers_new)!=length(object@markers)) {
      cat("incorrect number of markers")
      return(object)
    }
    markers_old = object@markers
  } else {
    if (length(markers_old)!=length(markers_new)) {
      cat("incorrect number of markers")
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
  object@node_features = llply(object@node_features, function(x) {
    colnames(x) = object@graph$v$phenotype
    x
  })
  if (length(object@edge_features)>0) {
    ecn = paste0(object@graph$e$from, "_", object@graph$e$to)
    object@edge_features = llply(object@edge_features, function(x) {
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
      cat("incorrect number of ids")
      return(object)
    }
    ids_old = object@meta$id
  } else {
    if (length(ids_old)!=length(ids_new)) {
      cat("incorrect number of ids")
      return(object)
    }
  }
  
  ids_ind = match(ids_old, object@meta$id)
  object@meta$id[ids_ind] = ids_new
  object@node_features = llply(object@node_features, function(x) {
    rownames(x)[ids_ind] = ids_new
    x
  })
  if (length(object@edge_features)>0) {
    object@edge_features = llply(object@edge_features, function(x) {
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
    cat("please provide valid sample id's; see @meta$id")
    return(object)
  }
  
  id_inds = match(sample_ids,object@meta$id)
  
  object@meta = object@meta[id_inds,]
  object@node_features = llply(object@node_features, function(x) x[id_inds,,drop=F])
  if (length(object@edge_features)>0) {
    object@edge_features = llply(object@edge_features, function(x) x[id_inds,,drop=F])
  }
  return(object)
} 


# extract flowgraph for certain cell populations
extract_phenotypes = function(object, phenotypes) {
  require(plyr)
  
  if (!any(phenotypes%in%object@graph$v$phenotype)) {
    cat("please provide valid phenotypes; see @graph$v")
    return(object)
  }
  
  id_inds = match(phenotypes,object@graph$v$phenotype)
  id_inds_ = object@graph$e$from%in%phenotypes & object@graph$e$to%in%phenotypes
  
  object@graph$v = object@graph$v[id_inds,,drop=F]
  object@graph$e = object@graph$e[id_inds_,,drop=F]
  object@gr = set_layout_graph(object@gr, as.function(object@plot_layout))
  
  object@node_features = llply(object@node_features, function(x) x[,id_inds,drop=F])
  if (length(object@edge_features)>0) {
    object@edge_features = llply(object@edge_features, function(x) x[,id_inds_,drop=F])
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
  nfs = intersect(names(obj1@node_features),names(obj2@node_features))
  object@node_features = 
    llply(nfs, function(xi) {
      a = obj1@node_features[[xi]]
      b = obj2@node_features[[xi]]
      abcol = intersect(colnames(a),colnames(b))
      ab = rbind(a[,match(abcol,colnames(a)),drop=F],
                 b[match(setdiff(rownames(b),rownames(a)),rownames(b)),match(abcol,colnames(b)),drop=F])
    })
  names(object@node_features) = nfs
  if (length(object@edge_features)>0) {
    efs = intersect(names(obj1@edge_features),names(obj2@edge_features))
    object@edge_features = llply(efs, function(xi) {
      a = obj1@edge_features[[xi]]
      b = obj2@edge_features[[xi]]
      abcol = intersect(colnames(a),colnames(b))
      ab = rbind(a[,match(abcol,colnames(a)),drop=F],
                 b[match(setdiff(rownames(b),rownames(a)),rownames(b)),match(abcol,colnames(b)),drop=F])
    })
    names(object@edge_features) = efs
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


# change layout
set_layout = function(object, layout_fun) {
  object@plot_layout = as.character(substitute(layout_fun))
  object@gr = set_layout_graph(object@gr, layout_fun)
  return(object)
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
    
    if(max(abs(logR)) < minlogR) { return(list(f=1, fdiff=0)) # f[i] = 1
    } else {
      
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
## output: Phenotype, PhenoCode, Phenolevel (number of markers per phenotype)
getPhen = function(phen) {
  require(stringr)
  pm = data.frame(phenotype=phen)
  markers = unique(unlist(str_split(phen,"[-+]")))
  markers = markers[markers!=""]
  pm$phenocode = sapply(phen, function(x){
    if (x=="") return(paste0(rep(0,length(markers)), collapse=""))
    b = str_split(x,"[+-]+")[[1]]
    b = b[-length(b)]
    bo = match(markers,b)
    phec = as.vector(sapply(str_extract_all(x,"[+-]+"), 
                            function(y) str_count(y,"[+]"))+1)
    c = phec[bo]
    c[is.na(c)] = 0
    paste0(c,collapse="")
  })
  pm$phenolevel = cell_type_layers(phen)
  pm$phenotype = as.character(pm$phenotype)
  
  return(pm)
}


## input: cell pop names or getPhen output meta_cell
## output: child (pos/neg) and parent list
getPhenCP = function(cp=NULL, meta_cell=NULL, no_cores=1) {
  require(plyr)
  parl = parallel_backend(no_cores)
  if (parl) if (no_cores>detectCores()) no_cores = detectCores()
  
  if (is.null(cp) & is.null(meta_cell)) { 
    print("give me something!"); return(NULL) 
  } else if (is.null(meta_cell)) {
    meta_cell = getPhen(cp)
  }
  meta_cell_grid = ldply(1:nrow(meta_cell), function(i) {
    pc = meta_cell$phenocode[i]
    pc = as.numeric(laply(seq(1, nchar(pc), 1), function(i) substr(pc, i, i)))
  })#, .parallel=parl)
  mcsum = rowSums(meta_cell_grid)
  lvl1ind = which(meta_cell$phenolevel==1 & rowSums(meta_cell_grid)==1)
  lvl1order = lvl1ind[sapply(lvl1ind, function(x) which(meta_cell_grid[x,]==1))]
  colnames(meta_cell_grid) = unique(gsub("[+-]","",meta_cell$phenotype[lvl1order]))
  rownames(meta_cell_grid) = meta_cell$phenotype
  mcpls = llply(unique(meta_cell$phenolevel), function(x) meta_cell$phenolevel==x)
  names(mcpls) = unique(meta_cell$phenolevel)
  maxl = max(meta_cell$phenolevel)
  minl = min(meta_cell$phenolevel)
  
  ilevel = meta_cell$phenolevel==minl
  iparen = NULL; ichild = meta_cell$phenolevel==minl+1
  res = NULL
  
  cat("\n")
  for (pl in sort(unique(meta_cell$phenolevel))) {
    start2 = Sys.time()
    cat("- ",sum(ilevel)," pops @ ")
    
    if(!is.null(ichild)) ccand = meta_cell_grid[ichild,,drop=F]
    if(!is.null(iparen)) pcand = meta_cell_grid[iparen,,drop=F]
    
    loop_ind = loop_ind_f(which(ilevel),no_cores)
    result_ = llply(loop_ind, function(ii) { #for (ii in loop_ind) {
      resulti = NULL
      
      # children
      if (pl==0) {
        pos_ind = apply(ccand,1,function(x) 2%in%x)
        resulti[[""]]$pos = meta_cell$phenotype[which(ichild)[pos_ind]]
        resulti[[""]]$neg = meta_cell$phenotype[which(ichild)[!pos_ind]]
      } else if(!is.null(ichild)) {
        for (i in ii) {
          ci = cn = cp = NULL
          mci = meta_cell_grid[i,,drop=T]
          ci_ind = llply(which(mci>0), function(j) ccand[,j]==mci[j])
          ci = which(ichild)[Reduce("&",ci_ind)] # rows with the same stuff as i
          names(ci) = meta_cell$phenotype[ci]
          if (length(ci)>0) {
            negcind = mcsum[ci]-1==mcsum[i] 
            cn = ci[ negcind]; if(sum( negcind)==0) cn = NULL # A-
            cp = ci[!negcind]; if(sum(!negcind)==0) cp = NULL # A+, A++, A+++ ...
          }
          resulti[[meta_cell$phenotype[i]]]$pos = names(cp)
          resulti[[meta_cell$phenotype[i]]]$neg = names(cn)
          
        }
      } else {
        for (i in ii) resulti[[meta_cell$phenotype[i]]]$pos = NULL
        for (i in ii) resulti[[meta_cell$phenotype[i]]]$neg = NULL
      }
      
      # parent
      if (pl==1) {
        for (i in ii) resulti[[meta_cell$phenotype[i]]]$parent = ""
      } else if(!is.null(iparen)) {
        for (i in ii) {
          mci = meta_cell_grid[i,,drop=T]
          pi = NULL
          pi_ind = sapply(which(mci>0), function(j) pcand[,j]==mci[j])
          pi = which(iparen)[apply(pi_ind, 1, function(x) sum(!x)==1)]
          names(pi) = meta_cell$phenotype[pi]
          if (length(pi)>0)
            resulti[[meta_cell$phenotype[i]]]$parent = names(pi)
        }
      } else {
        for (i in ii) resulti[[meta_cell$phenotype[i]]]$parent = NULL
      }
      
      return(resulti)
    }, .parallel=parl)
    result = Reduce("append",result_)
    res = append(res, result)
    
    iparen = ilevel
    ilevel = ichild
    ichild = NULL; if ((pl+1)!=maxl) ichild = meta_cell$phenolevel==pl+2
    time_output(start2,paste0("layer",pl))
  }
  
  pchild = llply(res, function(x) {
    a = NULL
    if (!is.null(x$neg)) a$neg = x$neg
    if (!is.null(x$pos)) a$pos = x$pos
    return(a)
  })
  pchild = Filter(Negate(is.null), pchild)
  
  pparen = llply(res, function(x) return(x$parent))
  pparen = Filter(Negate(is.null), pparen)
  
  return(list(pchild=pchild, pparen=pparen))
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
    # edit layout manually
    gys = sort(unique(gr_vxy[,2]))
    gxns = sapply(gys,function(y) length(gr_vxy[gr_vxy[,2]==y,1]))
    names(gxns) = gys
    gxnmax = max(gxns)
    gxnmaxl = which(gxns==gxnmax)
    minwidthmid = 2
    maxwidthmid = 4
    gxnmaxwidth = gxnmax*minwidthmid-1
    gxos = unlist(llply(gys, function(gy) {
      gxtf = which(gr_vxy[,2]==gy)
      gx = gr_vxy[gxtf,1]
      gxtf[order(gx)]
    }))
    gr_vxy[gxos,1] = unlist(llply(1:length(gys), function(gyi) {
      if (gyi%in%gxnmaxl) return( seq(0,gxnmaxwidth,by=minwidthmid) )
      if (gxns[gyi]==1) return( gxnmaxwidth/2 )
      by = min(maxwidthmid,(gxnmaxwidth+1)/(gxns[gyi]+1))
      a = seq(0,gxns[gyi]-1)*by
      a + (gxnmaxwidth-(a[length(a)]-1))/2
    }))
    # switch sideways
    gr_vxy = gr_vxy[,2:1]
    gr_vxy[,1] = max(gr_vxy[,1])-gr_vxy[,1]
  }
  
  # get node
  colnames(gr_vxy) = c("x","y")
  gr_v = cbind(gr_v,gr_vxy[match(gr_v$phenotype,names(V(gr0)[[]])),])
  
  # get edge
  gr_e$from.x = gr_v$x[match(gr_e$from, gr_v$phenotype)]
  gr_e$from.y = gr_v$y[match(gr_e$from, gr_v$phenotype)]
  gr_e$to.x = gr_v$x[match(gr_e$to, gr_v$phenotype)]
  gr_e$to.y = gr_v$y[match(gr_e$to, gr_v$phenotype)]
  
  return(list(e=gr_e,v=gr_v))
}


## flowgraph initialization -----------------------------------
flowgraph = function(input_, meta=NULL, layout_fun=layout.reingold.tilford, markers=NULL,
                     specenr=T, normalize=T, norm_ind=0, norm_layer=3, norm_path=NULL, no_cores=1, ...) { # layout.circle
  start = Sys.time()
  
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
      testclass = llply(input_, function(x) class(x)=="Phenotypes")
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
    try({
      phenotype = rownames(ftl[[1]]@MFIs)
    }, silent=T)
    if (is.null(phenotype)) {
      print("no phenotype cell populations labels in Phenotype file")
      return(NULL)
    }
    
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
  mc = mc[,apply(mc, 2, function(x) any(x>0)),drop=F] # remove all 0 columns
  phenotype = colnames(mc)
  
  ## make meta for cell populations FINAL
  # markern = length(markers)
  meta_cell = getPhen(phenotype)
  
  # make parent list (for each cell popultion, list its parents)
  pccell = getPhenCP(meta_cell=meta_cell, no_cores=no_cores)
  pchild = pccell$pchild # not used, returned
  pparen = pccell$pparen
  
  # extract cell populations that have parents to calculate expected proportions for; just to be safe
  cells1 = meta_cell$phenotype[meta_cell$phenolevel==1]
  cells1 = append("",cells1[cells1%in%names(pparen)])
  cells = meta_cell$phenotype[meta_cell$phenolevel>1]
  cells = cells[cells%in%names(pparen)]
  cells_ = append(cells1,cells)
  
  meta_cell = meta_cell[match(cells_,meta_cell$phenotype),]
  root = which(meta_cell$phenolevel==0)
  
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
  
  time_output(start1)
  
  time_output(start, "total time used")
  
  ## list of outputs
  # mc # count matrix
  # mp # proportion matrix
  # childprop_ # edge proportion matrix (colnames = parent_child cell population)
  
  
  
  childprop_names = Reduce("rbind", llply(names(pparen_), function(i) Reduce(rbind,llply(pparen_[[i]], function(j) c(j,i))) ))
  colnames(childprop_names) = c("from","to")
  childprop_names = as.data.frame(childprop_names)
  rownames(childprop_names) = 1:nrow(childprop_names)
  
  gr = list(e=data.frame(childprop_names,width=1,color="", e_ind=F),
            v=data.frame(
              meta_cell, 
              size=1, color="", sizeb=1, colorb="", fill="", 
              label=meta_cell$phenotype, 
              label_ind=F, v_ind=F, vb_ind=F))
  gr = set_layout_graph(gr, layout_fun) # layout cell hierarchy
  
  if (is.null(meta)) meta = data.frame(id=sample_id)
  fg = new("flowgraph", 
           node_features=
             list(count=mc),
           edge_features=list(), markers=markers,
           graph=gr, edge_list=list(child=pchild_,parent=pparen_), 
           meta=meta, plot_layout=as.character(substitute(layout_fun)), etc=list())
  
  if (normalize)
    fg = flowgraph_normalize(fg, norm_ind=norm_ind, norm_layer=norm_layer, norm_path=norm_path, no_cores=no_cores, ...)
  
  if (specenr)
    fg = flowgraph_specenr(fg, no_cores=no_cores)
  
  return(fg)
}


## feature(optional): normalized (adjusted) count 
flowgraph_normalize = function(fg, norm_ind=0, norm_layer=3, norm_path=NULL, no_cores=1, ...) {
  require(plyr)
  parl = parallel_backend(no_cores)
  if (parl) if (no_cores>detectCores()) no_cores = detectCores()
  
  mca = mc = fg@node_features$count
  fdiff0 = rep(0,nrow(mc))
  f0 = rep(1,nrow(mc))
  root = which(colnames(mc)=="")
  
  if (nrow(mca)>1 & !is.null(norm_ind)) {
    start1 = Sys.time()
    cat("preparing feature(s): normalized count ")
    
    sample_id = fg@meta$id
    meta_cell = fg@graph$v
    maxlayer = max(meta_cell$phenolevel)
    
    # prepare feat_file_cell_counts
    x = x0 = as.matrix(mc)[,-root,drop=F] # take out total cell count
    maxx = max(x0[is.finite(x0)])
    rootc = mc[,root]
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
      x = x[,colnames(x0)%in%meta_cell$phenotype[meta_cell$phenolevel%in%norm_layer], drop=F]
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
    
    mca = Reduce(rbind,llply(c(1:nrow(mc)), function(x) mc[x,]*f0[x]))
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
    fg@node_features$count_norm = mca
    
    if ("specenr"%in%names(fg@node_features))
      fg@node_features$expect_count_norm = fg@node_features$specenr*mca[,1]
    
    time_output(start1)
  } # return f0 and fdiff0
  
  return(fg)
}


## feature: proportion
flowgraph_prop = function(fg) {
  start1 = Sys.time()
  cat("preparing feature(s): proportion ")
  
  mc = fg@node_features$count
  mp = fg@node_features$prop = mc/mc[,colnames(mc)==""]
  dimnames(mp) = dimnames(mc)
  fg@node_features$prop = mp
  
  time_output(start1)
  return(fg)
}


## feature: edge proportions
flowgraph_prop_edge = function(fg, no_cores=1) {
  require(plyr)
  parl = parallel_backend(no_cores)
  if (parl) if (no_cores>detectCores()) no_cores = detectCores()
  
  start1 = Sys.time()
  cat("preparing feature(s): proportion on edges ")
  
  mc = fg@node_Features$count
  cells_ = colnames(mc)
  loop_ind = loop_ind_f(which(cells_!=""), no_cores)
  childprop_ = Reduce("cbind",llply(loop_ind, function(ii) {
    Reduce("cbind",llply(ii, function(ic) {
      i = cells_[ic]
      pnames = pparen_[[i]]
      if (pnames[1]=="") {
        parent = mc[,colnames(mc)=="",drop=F]
      } else {
        parent = mc[,pnames,drop=F]
      }
      childprop2 = mc[,i]/parent
      colnames(childprop2) = paste0(pnames, "_", i)
      return(childprop2)
    }))
  }, .parallel=parl))
  
  childprop_[is.nan(childprop_)] = 0
  fg@edge_features$prop = childprop_
  
  time_output(start1)
  return(fg)
}


## feature: specenr
flowgraph_specenr = function(fg, no_cores=1) {
  require(gsubfn)
  require(stringr)
  require(plyr)
  parl = parallel_backend(no_cores)
  if (parl) if (no_cores>detectCores()) no_cores = detectCores()
  
  mc = fg@node_features$count
  meta_cell = fg@graph$v
  pparen = fg@edge_list$parent
  root = which(colnames(mc)=="")
  if (!"prop"%in%names(fg@node_features))
    fg = flowgraph_prop(fg)
  mp = fg@node_features$prop
  if (!"prop"%in%names(fg@edge_features))
    fg = flowgraph_prop_edge(fg)
  
  ## feature: expected proportion
  
  cells1 = meta_cell$phenotype[meta_cell$phenolevel==1]
  cells1 = append("",cells1[cells1%in%names(pparen)])
  cells = meta_cell$phenotype[meta_cell$phenolevel>1]
  cells = cells[cells%in%names(pparen)]
  cellsn = meta_cell$phenolevel
  pparen_ = fg@edge_list$parent
  
  ## calculate the expected proportions for cell populatins with positive marker conditions only
  expe1 = matrix(.5,nrow=nrow(mp), ncol=length(cells1), dimnames=list(rownames(mp),cells1))
  expe1[,1] = 1
  
  cells1_p = str_extract(cells1,"[+]+")
  maxp = max(nchar(cells1_p[!is.na(cells1_p)]))
  
  start1 = Sys.time()
  cat("preparing feature(s): expected proportion; ")
  
  cpind = which(!grepl("[-]",cells))
  loop_ind = loop_ind_f(cpind,no_cores)
  expecp = Reduce("cbind", llply(loop_ind, function(ii) {
    Reduce("cbind", llply(ii, function(ic) {
      cpop = i = cells[ic]
      il = cellsn[ic]
      
      pnames = pparen_[[i]]
      parent = mp[,pnames,drop=F]
      parento = laply(1:nrow(parent), function(xi) order(parent[xi,]))
      if (is.null(dim(parento))) parento = matrix(parento,nrow=nrow(parent))
      grprnt = unique(unlist(pparen_[pnames]))
      
      pedges = Reduce(cbind,llply(1:length(pnames), function(gi) {
        gname = pparen_[[pnames[gi]]]
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
  expec = as.matrix(cbind(expe1,expecp))
  cpopneg = setdiff(cells,colnames(expec))
  cpopnegl = cell_type_layers(cpopneg)
  
  # replaces whatever pattern only once; e.g. gsubfn("_", p, "A_B_C_")
  p = proto(i=1, j=1, fun=function(this, x) if (count>=i && count<=j) "*" else x) 
  
  expecn = Reduce("cbind",llply(sort(unique(cpopnegl)), function(lev) {
    cpopl = cpopneg[cpopnegl==lev]
    cpopnegno = str_count(cpopl,"[-]") # number of negative marker conditions
    
    ## version 1; DELETE ONE VERSION
    expec_temp = NULL
    for (negno in sort(unique(cpopnegno))) {
      cpopi = cpopl[cpopnegno==negno]
      sibs = gsubfn("-", p, cpopi)
      for (sib in unique(sibs)) {
        sib_ = sapply(1:maxp, function(pn) 
          gsub("[*]", paste0(rep("+",pn),collapse=""), sib) )
        sib_ = sib_[sib_%in%cells]
        sib_sum = rowSums(mp[,sib_,drop=F])
        sibpars = pparen_[[sib_[1]]]
        sibis = which(sibs==sib) # shouldn't be length(sibis)>1, but safe
        for (sibi in sibis) {
          parsin = intersect(sibpars, pparen_[[cpopi[sibi]]])
          # if (length(sibpars)>1) #there must be only 1 common parent
          #   for (sibpar in sibpars)
          #     parsin = intersect(sibpar, parsin)
          expec_temp = cbind(expec_temp, mp[,parsin] - sib_sum)
          colnames(expec_temp)[ncol(expec_temp)] = cpopi[sibi]
        }
      }
    }
    
    # ## version 2
    # expec_temp = Reduce(cbind,llply(cpopl, function(cpi) {
    #   cpig = gsub("[+-]","",cpi)
    #   parnames = intersect(pparen_[[cpi]],colnames(mp))
    #   a = NULL
    #   for (parname in parnames) {
    #     chilnames = unlist(pchild[[parname]])
    #     chilnames = chilnames[chilnames!=cpi]
    #     chilnamesg = gsub("[+-]","",chilnames)
    #     chilgreat = chilnamesg==cpig & chilnames%in%colnames(expec)
    #     if (any(chilgreat)) {
    #       a = matrix(mp[,parname]-rowSums(expec[,chilnames[chilgreat],drop=F]),ncol=1)
    #       break
    #     }
    #   }
    #   if (is.null(a)) {
    #     a = b = mp[,parnames,drop=F]
    #     if (ncol(b)>1) a = matrix(apply(b,1,min), ncol=1)
    #   }
    #   colnames(a) = cpi
    #   return(a)
    # }))
    # if (is.null(dim(expec_temp))) expec_temp = matrix(expec_temp,ncol=1)
    
    return(expec_temp)
  }))
  expec = cbind(expec,expecn)
  
  exp1 = cbind(expe1,expec[,match(cells,colnames(expec)),drop=F])
  a = mp/exp1
  a[is.infinite(a)] = max(a[is.finite(a)])
  a[is.nan(a)] = 0
  lnpropexpect1 = log(a)
  lnpropexpect1[is.nan(lnpropexpect1)] = 0
  lnpropexpect1[exp1==0] = log(mp[exp1==0])
  lnpropexpect1[mp==0] = 0
  
  fg@node_features$expect_prop = exp1
  fg@node_features$specenr = lnpropexpect1
  
  if ("count_norm"%in%names(fg@node_features)) {
    mca = fg@node_features$count_norm
    expc1 = exp1*mca[,1]
    fg@node_features$expect_count_norm = expc1
  }
  
  time_output(start1)
}


## feature(s): get mean of a class for each feature
flowgraph_mean_class = function(fg, class, no_cores=1, node_features=NULL, edge_features=NULL) {
  parl = parallel_backend(no_cores)
  if (parl) if (no_cores>detectCores()) no_cores = detectCores()
  
  if (!class%in%colnames(fg@meta)) {
    cat("invalid class name, choose one from @meta")
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
      m[pii,] = Reduce(rbind, llply(1:nrow(m_), function(i) m_[i,]-m_m ))
    }
    dimnames(m) = dimnames(m0)
    return(m)
  }
  
  if (is.null(node_features)) node_features = names(fg@node_features)
  node_features = node_features[node_features%in%names(fg@node_features)]
  node_features = node_features[!grepl("MEAN", names(node_features))]
  if (length(node_features)>0) {
    fg_nodefs0 = fg@node_features[node_features]
    fg_nodefs = llply(fg_nodefs0, function(m0) 
      meandiff(m0, fg@meta[,class]), .parallel=parl)
    names(fg_nodefs) = paste0(names(fg_nodefs),label)
    fg@node_features = append(fg@node_features, fg_nodefs)
  }
  
  if (is.null(edge_features)) edge_features = names(fg@edge_features)
  edge_features = edge_features[edge_features%in%names(fg@edge_features)]
  edge_features = edge_features[!grepl("MEAN", names(edge_features))]
  if (length(edge_features)>0) {
    fg_edgefs0 = fg@edge_features[edge_features]
    fg_edgefs = llply(fg_edgefs0, function(m0) 
      meandiff(m0, fg@meta[,class]), .parallel=parl)
    names(fg_nodefs) = paste0(names(fg_edgefs),label)
    fg@edge_features = append(fg@edge_features, fg_edgefs)
  }
  
  return(fg)
}


## get p values 
flowgraph_p = function(
  fg, class, no_cores=1, 
  node_features=NULL, edge_features=NULL, 
  control=NULL,
  overwrite=F, 
  test_name="t_BY", # name your statistical test / adjustment method; if same name, overwrite=T will overwrite, overwrite=F will return nothing
  diminish=F, # don't test if all parents are insignificant, stricter the lower the layer
  p_thres=.05, # only used if diminish=T
  diminish_level=3, # only used if diminish=T; >3 levels won't be tested if parent not significant
  test=function(x,y) tryCatch(t.test(x,y)$p.value, error=function(e) 1), 
  adjust=function(x) p.adjust(x, method="BY")
) {
  # if control is NULL, compare all
  # if node/edge_features is NULL, do all
  parl = parallel_backend(no_cores)
  if (parl) if (no_cores>detectCores()) no_cores = detectCores()
  
  if (!class%in%colnames(fg@meta)) {
    cat("invalid class, choose one from @meta")
    return(fg)
  }
  
  classes_ = fg@meta[,class]
  classes = unique(classes_)
  if (length(classes)<2) {
    cat("invalid class, choose one from @meta")
    return(fg)
  }
  
  # class list lists what classes to compare
  if (!is.null(control)) {
    if (!control%in%fg@meta[,class]) {
      cat("invalid control, choose one from @meta")
      return(fg)
    }
    class_list = llply(classes[classes!=control], function(x) c(control, x))
  } else {
    class_list = llply(2:length(classes), function(i) 
      llply(1:(i-1), function(j) c(classes[i], classes[j])))
  }
  
  if (is.null(node_features)) node_features = names(fg@node_features)
  node_features = node_features[node_features%in%names(fg@node_features)]
  if (length(node_features)>0) {
    for (classl in class_list) {
      for (node_feature in node_features) {
        # prep index
        nfi = 1
        if (test_name%in%names(fg@node_p_class)) {
          if (class%in%names(fg@node_p_class[[test_name]])) {
            if (node_feature%in%fg@node_p_class[[test_name]][[class]]) {
              nfi = which(sapply(fg@node_p_class[[test_name]][[class]][[node_feature]], function(x)
                setequal(names(x),classl)))
              if (length(nfi)!=0) {
                if (!overwrite) next
              } else {
                nfi = length(fg@node_p_class[[test_name]][[class]][[node_feature]]) + 1
              }
            }
          }
        }
        id1 = fg@node_p_class[[test_name]][[class]][[node_feature]][[nfi]][[1]] = 
          fg@meta$id[classes_==classl[1]]
        id2 = fg@node_p_class[[test_name]][[class]][[node_feature]][[nfi]][[2]] = 
          fg@meta$id[classes_==classl[2]]
        names(fg@node_p_class[[test_name]][[class]][[node_feature]][[nfi]]) = classl
        
        m = fg@node_features[[node_feature]]
        m1 = m[id1,]
        m2 = m[id2,]
        
        if (!diminish) {
          loop_ind = loop_ind_f(1:ncol(m), no_cores)
          p = unlist(llply(loop_ind, function(ii) 
            llply(ii, function(i) test(m1[,i], m2[,i])), .parallel=parl ))
          p = adjust(p)
        } else {
          root = which(levels==0) # should be 1
          pparen = fg@edge_list$parent
          
          levels_ = fg@graph$v$phenolevel
          levels = order(unique(levels_))
          loop_ind_ = llply(levels[-root], function(l) which(levels_==l))
          
          root_ = which(colnames(m)=="")
          p = rep(NA,ncol(m))
          p[root_] = test(m1[,root_], m2[,root_])
          names(p)[root_] = ""
          
          for (lvl in levels[-root]) {
            lvlover = lvl>diminish_level
            loop_ind = loop_ind_f(loop_ind_[[lvl]], no_cores)
            p[loop_ind_[[lvl]]] = unlist(llply(loop_ind, function(ii) 
              llply(ii, function(i) {
                pheni = fg@graph$v$phenotype[i]
                parnis = pparen[[pheni]]
                if (all(p[parnis]>p_thres) & lvlover) return(NA)
                test(m1[,i], m2[,i])
              }), .parallel=parl
            ))
          }
          pai = !is.na(p)
          if (sum(pai)>1) p[pai] = adjust(p[pai])
          p[pai] = 1
        }
        names(p) = colnames(m)
        
        fg@node_p[[test_name]][[class]][[node_feature]][[nfi]] = p
      }
    }
  }
  
  if (is.null(edge_features)) edge_features = names(fg@edge_features)
  edge_features = edge_features[edge_features%in%names(fg@edge_features)]
  if (length(edge_features)>0) {
    for (classl in class_list) {
      for (edge_feature in edge_features) {
        # prep index
        nfi = 1
        if (test_name%in%names(fg@edge_p_class)) {
          if (class%in%names(fg@edge_p_class[[test_name]])) {
            if (edge_feature%in%fg@edge_p_class[[test_name]][[class]]) {
              nfi = which(sapply(fg@edge_p_class[[test_name]][[class]][[edge_feature]], function(x)
                setequal(names(x),classl)))
              if (length(nfi)!=0) {
                if (!overwrite) next
              } else {
                nfi = length(fg@edge_p_class[[test_name]][[class]][[edge_feature]]) + 1
              }
            }
          }
        }
        id1 = fg@edge_p_class[[test_name]][[class]][[edge_feature]][[nfi]][[1]] = 
          fg@meta$id[classes_==classl[1]]
        id2 = fg@edge_p_class[[test_name]][[class]][[edge_feature]][[nfi]][[2]] = 
          fg@meta$id[classes_==classl[2]]
        names(fg@edge_p_class[[test_name]][[class]][[edge_feature]][[nfi]]) = classl
        
        m = fg@edge_features[[edge_feature]]
        m1 = m[id1,]
        m2 = m[id2,]
        
        if (!diminish) {
          loop_ind = loop_ind_f(1:ncol(m), no_cores)
          p = unlist(llply(loop_ind, function(ii) 
            llply(ii, function(i) test(m1[,i], m2[,i])), .parallel=parl ))
          p = adjust(p)
        } else {
          pparen = fg@edge_list$parent
          
          phen_to = fg@graph$e$to
          phen_from = fg@graph$e$from
          levels_ = fg@graph$v$phenolevel[match(phen_to, fg@graph$v$phenotype)]
          levels = order(unique(levels_))
          loop_ind_ = llply(levels, function(l) which(levels_==l))
          
          p = rep(NA,ncol(m))
          
          for (lvl in levels) {
            lvlover = lvl>diminish_level
            loop_ind = loop_ind_f(loop_ind_[[lvl]], no_cores)
            p[loop_ind_[[lvl]]] = unlist(llply(loop_ind, function(ii) 
              llply(ii, function(i) {
                pheni = fg@graph$v$phenotype[i]
                parnis = phen_to==phen_from[i]
                if (sum(parnis)>0)
                  if (all(p[parnis]>p_thres) & lvlover) return(NA)
                return(test(m1[,i], m2[,i]))
              }), .parallel=parl
            ))
          }
          pai = !is.na(p)
          if (sum(pai)>1) p[pai] = adjust(p[pai])
          p[pai] = 1
        }
        names(p) = colnames(m)
        
        fg@edge_p[[test_name]][[class]][[edge_feature]][[nfi]] = p
      }
    }
  }
  
  return(fg)
}


## flowgraph interpretation functions ------------------------------


setMethod(
  "show", "flowgraph",
  function(object) 
    cat("flowgraph object with ", length(object@node_features)," cell population and ", length(object@edge_features), " edge feature(s)", 
        "\n- markers: ", object@markers,
        "\n- contains: ", nrow(object@graph$v), " cell populations and ", nrow(object@graph$e), " edges", sep="")
)


setMethod(
  "summary", "flowgraph",
  function(object) {
    show(object)
    
    summary_table = function(m) 
      data.frame(data=data, feat=feat_type, nrow=nrow(m), ncol=ncol(m), 
                 inf=sum(is.infinite(m)), na=sum(is.na(m)), nan= sum(is.nan(m)), 
                 neg=sum(m<0), pos=sum(m>0), zero=sum(m==0), max=max(m[is.finite(m)]))
    
    result1 = ldply(object@node_features, function(m) summary_table(m))
    result2 = NULL
    if (length(object@edge_features)>0) 
      result2 = ldply(object@edge_features, function(m) summary_table(m))
    tab = rbind(result1,result2)
    colnames(tab)[1] = "feature_type"
    return(tab)
  }
)



