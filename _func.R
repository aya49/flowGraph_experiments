## helper functions
## aya43@sfu.ca
## created 20180517
## last modified 20180517


## default options -----------------------------------------
options(stringsAsFactors=F)


## system functions ------------------------------

## input: list of package names to load
## output: none; load/install package
libr = function(pkgs) {
  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) 
    install.packages(setdiff(pkgs, rownames(installed.packages())), verbose=F)
  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
    if (!requireNamespace("BiocManager", quietly=T)) install.packages("BiocManager")
    BiocManager::install(setdiff(pkgs, rownames(installed.packages())), ask=F)
  }
  sapply(pkgs, require, character.only=T)
}



## unload all pkgs
unloadpkg = function() 
  a = lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""), detach, character.only=T, unload=T)


## helper functions ----------------------------


## input: Sys.time() value
## output: formatted time as string; used in time_output function
tstr = function(time, tz="GMT") return( format(.POSIXct(time,tz=tz), "%H:%M:%S") )


## input: Sys.time() value
## output: time elapsed since input
time_output = function(start, message="", tz="GMT") {
  start = as.POSIXct(start)
  end = Sys.time()
  time_elapsed = difftime(end, start, units="secs")
  cat(message, ifelse(message=="","",": "), tstr(start,tz=tz), "-", tstr(end,tz=tz), ">", tstr(time_elapsed,tz=tz), "\n")
}


## input: x=vector of indices; n=cores or number of parts you want to split x into
## ouput: list of n indices vector
## for parallel processing, one list item per core
loop_ind_f = loopInd = function(x,n) {
  if (n==1) return(list(x))
  return(split(x, ceiling(seq_along(x)/ceiling(length(x)/n))))
}


## input: file path and a file extension
## output: List of all file names in given path with the given extension
fileNames = function(pathList, ext="fcs") {
  temp.str = sapply(strsplit(pathList, "/"), function(x) x[length(x)])
  pathList = sub(paste0(".", ext,"$"), "", temp.str, ignore.case=T)
  return(pathList)
}

# List of all last folder names in given path (Genotype)
folderNames = function(pathList) {
  folders = c()
  for (i in 1:length(pathList)) {
    temp.str = strsplit(pathList[i],split="/")
    folders[i] = temp.str[[1]][length(temp.str[[1]])-1]
  }
  return(folders)
}




#remove NULL elements from list
remove_null = function(x) return(Filter(Negate(is.null), x))


## input: matrix
## output: matrix without all NA col/row
delna = function(m)
  m[apply(m,1,function(x) !all(is.na(x))), apply(m,2,function(x) !all(is.na(x)))]


# deletes columns where all values are the same
delcols = function(x) x[sapply(x, function(y) length(unique(y))>1)]


## input: matrix
## output: returns u1 (col index with only 1 unique element), 
##                 ua (col index where every row is a unique element), prints colnames and its unique elements if there is less than n unique elements
col_probe = function(m,n=15) {
  require(data.table)
  u1 = NULL
  ua = NULL
  nm = nrow(m)
  for (col in 1:ncol(m)) {
    coln = colnames(m)[col]
    if (is.data.table(m)) {
      a = unique(as.matrix(m[,..col]))
    } else {
      a = as.vector(unique(as.matrix(m[,col])))
    }
    la = length(a)
    if (la == 1) u1 = append(u1,col)
    if (la == nm) ua = append(ua,col)
    
    cat("\n", coln, ": ", la)
    if (la<n) cat("; ",a)
  }
  c("\n")
  return(list(u1=u1,ua=ua))
}


## Input: Clt, a matrix
## Output: Outputs column inds of duplicate columns
duplicateindM = function(clt) {
  delcol_dup = duplicated(t(clt))
  col_dup = 1:ncol(clt); col_dup[delcol_dup] = sapply(which(delcol_dup),function(x) {
    for (xi in which(!delcol_dup)) {
      if (identical(clt[,xi],clt[,x]) ) { a = xi; break }
    }
    return(a)
  })
  return(col_dup)
}


## arithmic functions ----------------------------

## input: vector
## outpute: mode, most frequently occuring element in vector
mode_set = Mode = function(y) unique(y)[which.max(tabulate(match(y, unique(y))))]



## input: vector
## outpute: geometric mean
## from package psych
mean_geo = geometric.mean = function(x,na.rm=TRUE) {
  if (is.null(nrow(x))) {
    exp(mean(log(x),na.rm=TRUE))
  } else {
    exp(apply(log(x),2,mean,na.rm=na.rm))
  }
}


# whether or not x is in between two values rng
'%between%' = function(x,rng) (x>rng[1] & x<rng[2]) | (x<rng[1] & x>rng[2])



## Jensen-shannon-divergence (half of kullback divergence both ways)
# http://stackoverflow.com/questions/11226627/jensen-shannon-divergence-in-r

jsd = function(p,q, list=T) { #p,q are two distributions, two lists of distributions if to avg over
  if (list) {
    JS = 0
    for (i in 1:length(p)) {
      m = (p[[i]] + q[[i]])/2
      JS = ((sum(p[[i]]*log(p[[i]] / m)) + sum(q[[i]]*log(q[[i]] / m))) /2) + JS
    }
    JS = JS/length(p)
  } else {
    m = (p + q)/2
    JS = ((sum(p * log(p / m)) + sum(q * log(q / m))) /2) + JS
  }
  return(JS)
}

# plotting functions ------------------------------

## input: flowframe or data frame; 2 columns
## ouput: density intensity plot
# references plotDens from flowDensity
plot_int = function(dat, col, main, pch = ".", ...) {
  if (missing(col)) {
    colPalette = colorRampPalette(c("grey", "black"))
    col = densCols(dat, colramp = colPalette)
  }
  if (nrow(dat) < 2) {
    graphics::plot(1, type = "n", axes = F, ...)
  } else {
    graphics::plot(dat, col = col, pch = pch, ...)
  }
}


## graph plot functions
ggblank = function() {
  require(ggplot2)
  ggplot() +
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

gggraph = function(a,main="", label_ind=NULL, v_ind=NULL, e_ind=NULL) { # indices of whether to apply color size etc
  # gr_v: name x y label size color sizeb colorb
  # gr_e: from to from.x from.y to.x to.y color
require(ggrepel)
  
  gr_v = a$v
  gr_e = a$e
  if (is.null(label_ind)) label_ind = rep(T,nrow(gr_v))
  if (is.null(v_ind)) v_ind = rep(T,nrow(gr_v))
  if (is.null(e_ind)) e_ind = rep(T,nrow(gr_e))
  # base graph
  gp = ggblank() + ggtitle(main) +
    geom_segment(data=gr_e[!e_ind,], color="grey",
                 aes(x=from.x,xend=to.x, y=from.y,yend=to.y)) +
    geom_segment(data=gr_e[e_ind,], 
                 aes(x=from.x,xend=to.x, y=from.y,yend=to.y,
                     color=color)) +
    geom_point(data=gr_v[!v_ind,],aes(x=x,y=y), size=1, color="grey")+
    geom_point(data=gr_v[v_ind,],aes(x=x,y=y, color=colorb,size=sizeb))+
    geom_point(data=gr_v[v_ind,],aes(x=x,y=y, color=color, size=size)) +
    geom_label_repel(
      data=gr_v[v_ind,],
      aes(x=x,y=y,label=label, color=color), nudge_y = .3)
  
  return(gp)
}

gpdf = function(a) {
  gr_e = a$e
  gr_v = a$v
  return(list(e=data.frame(gr_e,width=1,color=""),
              v=data.frame(gr_v,size=1, color="",sizeb=1, colorb="")))
}

# image functions -----------------------------------

## input: gray scale 1d matrix
## output: rgb... repeat 3 times
gray2rgb = function(im) {
  if (length(dim(im)==2)) {
    im = array(im, dim=c(dim(im),1))
  }
  if (length(dim(im)==3)) {
    if (dim(im)[3]==1) {
      img = abind(im,im)
      img = abind(im,im)
    }
  } else {
    return(NULL)
  }
  return(img)
}



## flowtype functions -----------------------------

## Input: Phenotype
## Output: Phenotype, PhenoCode, Phenolevel (number of markers per phenotype)
getPhen = function(phen) {
  require(flowType)
  require(stringr)
  pm = data.frame(phenotype=phen)
  markers = unique(unlist(str_split(phen,"[+-]")))
  markers = markers[!markers==""]
  pm$phenocode = sapply(phen, function(x) encodePhenotype(x, markers))
  pm$phenolevel = cell_type_layers(phen)
  
  return(pm)
}

##Input: matrix
##Output: if features (rownames) have +/- symbols, returns corresponding feature layers; else returns NULL
cell_type_layers <- function(phen) str_count(phen, "[+-]")



## input: cell pop names or getPhen output meta_cell
## output: child (pos/neg) and parent list
getPhenCP = function(cp=NULL, meta_cell=NULL, no_cores=1) {
  require(plyr)
  if (is.null(cp) & is.null(meta_cell)) { 
    print("give me something!"); return(NULL) 
  } else if (is.null(meta_cell)) {
    meta_cell = getPhen(cp)
  }
  meta_cell_grid = ldply(1:nrow(meta_cell), function(i) {
    pc = meta_cell$phenocode[i]
    pc = as.numeric(laply(seq(1, nchar(pc), 1), function(i) substr(pc, i, i)))
  }, .parallel=T)
  colnames(meta_cell_grid) = unique(gsub("[+-]","",meta_cell$phenotype[meta_cell$phenolevel==1]))
  rownames(meta_cell_grid) = meta_cell$phenotype
  mcsum = rowSums(meta_cell_grid)
  mcpls = llply(unique(meta_cell$phenolevel), function(x) meta_cell$phenolevel==x)
  names(mcpls) = unique(meta_cell$phenolevel)
  maxl = max(meta_cell$phenolevel)
  minl = min(meta_cell$phenolevel)
  
  ilevel = meta_cell$phenolevel==minl
  iparen = NULL; ichild = meta_cell$phenolevel==minl+1
  res = NULL
  
  for (pl in unique(meta_cell$phenolevel)) {
    # for (pl in 0:5) {
    start2 = Sys.time()
    cat(sum(ilevel)," pops ")
    
    if(!is.null(ichild)) ccand = meta_cell_grid[ichild,,drop=F]
    if(!is.null(iparen)) pcand = meta_cell_grid[iparen,,drop=F]
    
    loop_ind = loopInd(which(ilevel),no_cores)
    result_ = llply(loop_ind, function(ii) { #for (ii in loop_ind) {
      resulti = NULL
      if (pl==0) {
        pos_ind = apply(ccand,1,function(x) 2%in%x)
        resulti[[""]]$pos = meta_cell$phenotype[which(ichild)[pos_ind]]
        resulti[[""]]$neg = meta_cell$phenotype[which(ichild)[!pos_ind]]
      } else if(!is.null(ichild)) {
        for (i in ii) {
          ci = cn = cp = NULL
          mci = meta_cell_grid[i,,drop=T]
          ci_ind = llply(which(mci>0), function(j) ccand[,j]==mci[j])
          ci = which(ichild)[Reduce("&",ci_ind)]
          names(ci) = meta_cell$phenotype[ci]
          if (length(ci)>0) {
            cn = ci[mcsum[ci]-1==mcsum[i]]; if(length(cn)==0) cn = NULL
            cp = ci[mcsum[ci]-2==mcsum[i]]; if(length(cp)==0) cp = NULL
          }
          resulti[[meta_cell$phenotype[i]]]$pos = names(cp)
          resulti[[meta_cell$phenotype[i]]]$neg = names(cn)
          
        }
      } else {
        for (i in ii) resulti[[meta_cell$phenotype[i]]]$pos = NULL
        for (i in ii) resulti[[meta_cell$phenotype[i]]]$neg = NULL
      }
      
      if (pl==1) {
        for (i in ii) resulti[[meta_cell$phenotype[i]]]$parent = ""
      } else if(!is.null(iparen)) {
        for (i in ii) {
          mci = meta_cell_grid[i,,drop=T]
          pi = NULL
          pi_ind = sapply(which(mci>0), function(j) pcand[,j]==mci[j])
          pi = which(iparen)[rowSums(pi_ind)==(sum(mci>0)-1)]
          names(pi) = meta_cell$phenotype[pi]
          if (length(pi)>0)
            resulti[[meta_cell$phenotype[i]]]$parent = names(pi)
        }
      } else {
        for (i in ii) resulti[[meta_cell$phenotype[i]]]$parent = NULL
      }
      
      return(resulti)
    }, .parallel=T)
    result = Reduce("append",result_)
    # result = unlist(result, recursive=F)
    # names(result) = meta_cell$phenotype[which(ilevel)]
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
  
  
  # graph edge list
  gc = ldply(1:length(pchild), function(x) 
    data.frame(from=names(pchild)[x], to=unlist(pchild[[x]])) )
  gp = ldply(1:length(pparen), function(x) 
    data.frame(from=pparen[[x]], to=names(pparen)[x]) )
  gr_e = unique(rbind(gc,gp))
  gr_v = meta_cell$phenotype

  gr_vp = gr_v[!grepl("[-]",gr_v)]
  gr_ep = gr_e[!grepl("[-]",gr_e[,1]) & !grepl("[-]",gr_e[,2]),]
  
  gr_v = data.frame(name=gr_v)
  gr_vp = data.frame(name=gr_vp)

  return(list(pchild=pchild, pparen=pparen,
              gr=list(e=gr_e,v=gr_v),
              grp=list(e=gr_ep,v=gr_vp)))
}

layout_gr = function(gr_e,gr_v,FUN=layout.reingold.tilford) {
  # FUN is a layout function from the igraph package
  # assume graph is connected, used internally
  require(igraph)
  
  # edit layout
  gr_vxy_ = FUN(gr = graph_from_data_frame(gr_e)) # layout.circle
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
  }
  
  # get node
  colnames(gr_vxy) = c("x","y")
  gr_v = cbind(gr_v,gr_vxy)

  # get edge
  gr_e$from.x <- gr_vxy$x[match(gr_e$from, gr_v$name)]
  gr_e$from.y <- gr_vxy$y[match(gr_e$from, gr_v$name)]
  gr_e$to.x <- gr_vxy$x[match(gr_e$to, gr_v$name)]
  gr_e$to.y <- gr_vxy$y[match(gr_e$to, gr_v$name)]
  
  return(list(e=gr_e,v=gr_v))
}




# other functions ---------------------------------

##Input: matrix with rows in order of date/time etc.
##Output: kalman filter for each column in matrix
kmf = function(m,cols) {
  fkffitall = NULL
  statsfitall = NULL
  for (i in cols) {
    yall = as.numeric(m[,i])
    if (length(unique(yall))==1) {
      fkffitall[[i]] = statsfitall[[i]] = yall
    } else {
      ## Set constant parameters:
      dt = ct = matrix(0) 
      Zt = Tt = matrix(1)
      a0 = yall[1]           # Estimation of the first sample count
      P0 = matrix(100)     # Variance of 'a0'
      ## Estimate parameters 23min if TS
      fit.fkf = optim(c(HHt = var(yall, na.rm = TRUE) * .5,
                         GGt = var(yall, na.rm = TRUE) * .5),
                       fn = function(par, ...) -fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik,
                       yt = rbind(yall), a0 = a0, P0 = P0, dt = dt, ct = ct,
                       Zt = Zt, Tt = Tt, check.input = FALSE)
      ## Filter Nile data with estimated parameters:
      fkf.obj = fkf(a0, P0, dt, ct, Tt, Zt, HHt = matrix(fit.fkf$par[1]), GGt = matrix(fit.fkf$par[2]), yt = rbind(yall))
      fkffitall[[i]] = fkf.obj$att[1,]
      ## Compare with the stats' structural time series implementation: 5min if TS
      statsfitall[[i]] = fitted(StructTS(yall, type = "level"))
    }
  }
  return(list(fkffitall=fkffitall,statsfitall=statsfitall))
}




















# distance scores ------------------------------------

## Input: Distance matrix, class list (numeric)
## Output: NCA score
NCA_score = function(x,y,delta=1, fast=F, doUnderflow=T, no_cores=1) {
  require(Brobdingnag)
  require(foreach)
  require(doMC)
  
  registerDoMC(no_cores)
  
  #preprocess matrix
  x = as.matrix(x)
  if (nrow(x)!=length(y)) {
    cat("x and y dimensions don't match")
    return(NULL) #check if x pairs with y
  }
  power = 2
  if (fast) power = 1
  if (!fast) delta = 1
  xe = -(x^power)/delta ## PREVENT UNDERFLOW
  underflow = F
  if (length(which(xe<(-700)))>0) {
    if (!doUnderflow) {
      cat(length(which(xe<700)), "values too big")
      return(NULL)
    }
    underflow = T
    diag(xe) = -Inf
  } else {
    xe = exp(xe)
    diag(xe) = 0
  }
  
  #pij = prob of xi picking xj as neighbour (fill top half)
  if (underflow) {
    pij = list()
    for (i in 1:length(y)) {
      pij[[i]] = brob(xe[i,])/sum(brob(xe[i,-i])) ## PREVENT UNDERFLOW
    }
  } else {
    pij = matrix(0, nrow=nrow(xe), ncol=ncol(xe))
    for (i in 1:length(y)) {
      pij[i,] = xe[i,]/sum(xe[i,-i])
    }
  }
  
  
  #pi = prob of classifying xi correctly; piy = prob of classifying points of class y correctly
  yt = table(y)
  yf = as.numeric(as.factor(y)) # as factor
  pyt = rep(0,length(yt))
  pi = rep(0,length(y))
  yi = lapply(1:length(yt), function(yn) return(which(y==names(yt)[yn])) )
  
  
  if (underflow) {
    # pi = sum(pij[[1]][yi[[yf[1]]]])
    # if (length(yf)==2) pi = cbrob(pi,sum(pij[[2]][yi[[yf[2]]]]))
    # if (length(yf)>2) {
    #   for (i in 2:length(yf)) {
    #     pi = cbrob(pi,sum(pij[[i]][yi[[yf[i]]]]))
    #   }
    # }
    # for (yn in 1:length(yt)) {
    #   pyt[yn] = as.numeric(sum(pi[yi[[yn]]]))
    # }
    for (yn in 1:length(yt)) { #for each class yn
      pi[yi[[yn]]] = sapply(yi[[yn]], function(i) as.numeric(sum(pij[[i]][yi[[yn]]])) )
      pyt[yn] = sum(pi[yi[[yn]]])
    }
  } else {
    for (yn in 1:length(yt)) { #for each class yn
      pi[yi[[yn]]] = sapply(yi[[yn]], function(i) sum(pij[i,yi[[yn]]]) )
      pyt[yn] = sum(pi[yi[[yn]]])
    }
  }
  names(pyt) = names(yt)
  return(list(pi=pi,yt=yt,pyt=pyt,p=sum(pyt))) #p is total score
}


## Input: 2 vector of labels
## Output: Precision, Recall, Specificity, F, Accuracy
f.measure.comembership = function(la,cl) {
  # require(clusteval)
  # ftpn = comembership_table(la,cl)
  # 
  # tn = ftpn$n_00
  # tp = ftpn$n_11
  # fn = ftpn$n_10
  # fp = ftpn$n_01
  
  tn = tp = fn = fp = 0
  for (lai in 1:(length(la)-1)) {
    for (laj in (lai+1):length(la)) {
      posl = la[lai] == la[laj]
      posc = cl[lai] == cl[laj]
      if (posl) {
        if (posc) {
          tp = tp+1
        } else {
          fn = fn+1
        }
      } else {
        if (posc) {
          fp = fp+1
        } else {
          tn = tn+1
        }
      }
    }
  }
  
  p=tp/(tp+fp)
  r=tp/(tp+fn)
  sp=tn/(tn+fp)
  return(list(p=p, r=r, sp=sp, f_comember=2*p*r/(p+r), a=(tp+tn)/(tn+tp+fn+fp)))
}



# pvalue -------------------------------------

## Significance; take note to put in same size control
## tstat = (a-mean(control))/se
## se = sd(control)/sqrt(length(control))
t.test.single = function(control,a) { #input a vector and a single number
  se = sd(control) #/ sqrt(length(control))
  tstat = (a-mean(control)) / se
  p = 2*pt(-abs(tstat), df=length(control)-1) # 2 tailed (less stringent)
  return(p)
}


# normalize ---------------------------------------

## Input: x matrix (phenotypes on columns, will convert in function) -- x0 is optional, plots according to x0 so if you want to plot more cell populations than those in x
## Output: f = normalize factors per sample; fidff = difference between f and peak of count ratio per sample
tmm = function(x,x0=NULL,lib.size,refColumn,cutoff=Inf,plotimg=T,pngnames=NULL,mains=NULL,no_cores=detectCores()-1,samplesOnCol=F,Acutoff=-10) {
  require(foreach)
  require(doMC)
  registerDoMC(no_cores)
  
  if (is.null(x0)) x0 = x
  if (!samplesOnCol) { x = t(x); x0 = t(x0) }
  
  ## Taken from TMM 
  
  #f = rep(NA,ncol(x))
  #fdiff = rep(NA,ncol(x)) #diff between density peak and value (note: logged)
  ref = x[,refColumn]
  refn = lib.size[refColumn]
  doWeighting = F
  logratioTrim = .3
  sumTrim = 0.05
  minlogR = 1e-6 #min value of log2((obs/obsn)/(ref/refn))
  
  
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
  fdiff = rep(NA,length(ff)) #diff between density peak and value (note: logged)
  for (i in 1:length(ff)) {
    f[i] = ff[[i]]$f
    try({ fdiff[i] = ff[[i]]$fdiff })
  }
  
  return(list(f=f,fdiff=fdiff))
}








# classify/cluster ---------------------------------------

## Input: distance matrix, parameter values
## Output: table of labels/cluster assigned to each item in distance matrix

knntable = function(dd,knn,label) {
  testind = which(is.na(label))
  clt = t(matrix(sapply(testind, function(x) {
    topkind = order(dd[x,])
    topkind = topkind[-(topkind==x)]
    topkfn = colnames(dd)[topkind]
    topkla = label[match(topkfn,rownames(dd))]
    topkla = topkla[!is.na(topkla)]
    return(sapply(knn, function(y) return(Mode(topkla[1:y])) ))
  }),nrow=length(knn)))
  if (!length(knn)>1) clt = as.matrix(clt,ncol=1)
  colnames(clt) = knn
  rownames(clt) = rownames(dd)[testind]
  return(clt)
}


pamtable = function(dd,nclass,pamtries) {
  require(cluster)
  pp0t = sapply(1:length(pamtries), function(j) pam(dd, k=(nclass),diss=T)$clustering)
  if (is.null(dim(pp0t))) pp0t = matrix(pp0t,ncol=1)
  colnames(pp0t) = rep("none",ncol(pp0t))
  rownames(pp0t) = rownames(dd)
  return(pp0t)
}


spectable = function(mm,nclass,label,methods=c("rbf"),tries=1,savedist=NULL,savesim=NULL,replace='kerns') {
  #methods=c("rbf","poly","vanilla","tanh","laplace","bessel","anova","spline")
  require(kernlab)
  require(stringr)
  tryCatch({
    
    parlist = pp0t = NULL
    for (method in methods) {
      parlist0 = sil = NCA = cl = sim = dist = NULL
      for (i in 1:tries) {
        sp = specc(as.matrix(mm),kernel=paste0(method,"dot"),centers=nclass)
        cl[[i]] = sp@.Data
        parlist0[[i]] = unlist(sp@kernelf@kpar)
        sim[[i]] = kernelMatrix(kernel=match.fun(paste0(method,"dot"))(as.numeric(parlist0[[i]])),x=as.matrix(mm))
        dist[[i]] = get_graphd(sim[[i]])
        sil = append(sil,median(silhouette(label,dist[[i]])[,3]))
        NCA = append(NCA,NCA_score(as.matrix(dist[[i]]), label, doUnderflow=doUnderflow)$p)
      }
      maxsilind = which.max(sil)
      pp0t = cbind(pp0t, cl[[maxsilind]])
      parlist = c(parlist, paste0(paste0(names(sp@kernelf@kpar),"-",signif(parlist0[[maxsilind]],4),collapse="_"),"_silmed-",signif(max(sil),5),"_NCA-",signif(NCA[maxsilind],5),collapse=""))
      s = sim[[maxsilind]]
      d = dist[[maxsilind]]
      if (!is.null(savesim)) save(s,file=paste0(gsub(replace,paste0("kern_",method,"_",paste0(names(sp@kernelf@kpar),"-",signif(parlist0[[maxsilind]],4),collapse="_")),savesim),"_simmatrix.Rdata"))
      if (!is.null(savedist)) save(d,file=paste0(gsub(replace,paste0("kern_",method,"_",paste0(names(sp@kernelf@kpar),"-",signif(parlist0[[maxsilind]],4),collapse="_")),savesim),"_dist.Rdata"))
    }
    colnames(pp0t) = parlist
    rownames(pp0t) = rownames(mm)
    return(pp0t)
  }, error = function(e) {
    return(NULL)
  })
  
}

spec1table = function(sim,nclass) {
  require(kernlab)
  sim = sim/max(sim)
  diag(sim) = 1
  pp0t = matrix(specc(as.kernelMatrix(sim),centers=nclass)@.Data, ncol=1)
  colnames(pp0t) = "none"
  rownames(pp0t) = rownames(sim)
  return(pp0t)
}

hctable = function(dd,nclass,links=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")) {
  pp0t = sapply(links,function(x) return(cutree(hclust(as.dist(dd),method=x),k=nclass)))
  colnames(pp0t) = links
  rownames(pp0t) = rownames(dd)
  return(pp0t)
}


lvtable = function(sim,rwThres) {
  require(igraph)
  pp0t = sapply(rwThres, function(rwt) {
    d2 = sim
    tops = quantile(as.vector(sim),rwt)
    d2[d2<tops] = 0
    gr = graph_from_adjacency_matrix(d2, weighted=T,mode='undirected', diag=F)
    return(cluster_louvain(gr)$membership)
  })
  colnames(pp0t) = rwThres
  rownames(pp0t) = rownames(sim)
  return(pp0t)
}

dctable = function(mm,alpha=.85,nu=seq(0.0, 1.0, by=0.05)) {
  require(densitycut)
  clt = DensityCut(X=mm, alpha=alpha, nu=nu, show.plot=F)$cluster
  clt = matrix(clt,ncol=1)
  rownames(clt) = rownames(mm)
  colnames(clt) = "none"
  return(clt)
}

dc1table = function(dd,k=3,alpha=.85,nu=seq(0.0, 1.0, by=0.05)) {
  require(densitycut)
  require(FastKNN)
  
  #create knn.ind, knn.dist
  dd = as.matrix(dd)
  ki = t(sapply(1:nrow(dd), function(x) k.nearest.neighbors(x,dd,k=k)))
  kd = t(sapply(1:nrow(dd), function(x) dd[x,ki[x,]]))
  rownames(ki) = rownames(kd) = rownames(dd)
  colnames(ki) = colnames(kd) = 1:k
  
  clt = DensityCut(knn.dist=kd,knn.index=ki, alpha=alpha, nu=nu, show.plot=F)$cluster
  clt = matrix(clt,ncol=1)
  rownames(clt) = rownames(dd)
  colnames(clt) = "none"
  return(clt)
}


rw1table = function(sim,rwThres) {
  require(igraph)
  gr = graph_from_adjacency_matrix(sim, weighted=T,mode='undirected', diag=F)
  
  pp0t = sapply(rwThres, function(rwt) {
    tops = quantile(as.vector(sim),rwt)
    gr1 = delete.edges(gr, which(E(gr)$weight<tops))
    return(components(gr1)$membership)
  })
  colnames(pp0t) = rwThres
  rownames(pp0t) = rownames(sim)
  return(pp0t)
  
}








##Input: Row x Bicluster (TRUE/FALSE)
##Output: cluster array for each row; put row into largest cluster its in
rowxcluster_to_cluster = function(rowxcluster) {
  clusters_size = apply(rowxcluster, 2, function(x) sum(x) )
  rclusters = apply(rowxcluster, 1, function(x) {
    cl = which(x)
    if (length(cl) == 0) cl = 0
    if (length(cl) > 1) cl = cl[which.max(clusters_size[cl])]
    return(cl)
  })
}


##Input: Bicluster x Col (TRUE/FALSE)
##Output: cluster array for each row; put row into largest cluster its in
clusterxcol_to_cluster = function(clusterxcol) {
  clusters_size = apply(clusterxcol, 1, function(x) sum(x) )
  return(apply(clusterxcol, 2, function(x) {
    cl = which(x)
    if (length(cl) == 0) cl = 0
    if (length(cl) > 1) cl = cl[which.max(clusters_size[cl])]
    return(cl)
  }))
}


