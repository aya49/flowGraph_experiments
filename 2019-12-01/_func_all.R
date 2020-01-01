## helper functions
## aya43@sfu.ca
## created 20180517
## last modified 20180517


## default options -----------------------------------------
options(stringsAsFactors=F)


## generic functions ------------------------------

## input: list of package names to load
## output: none; load/install package
libr <- function(pkgs) {
  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0)
    install.packages(setdiff(pkgs, rownames(installed.packages())), verbose=F)
  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
    if (!requireNamespace("BiocManager", quietly=T)) install.packages("BiocManager")
    BiocManager::install(setdiff(pkgs, rownames(installed.packages())), ask=F)
  }
  sapply(pkgs, require, character.only=T)
}


## input: Sys.time() value
## output: time elapsed since input
time_output = function(start, message="", tz="GMT") {
  start = as.POSIXct(start)
  end = Sys.time()
  time_elapsed = difftime(end, start, units="secs")
  cat(message, ifelse(message=="","",": "), ft(start,tz=tz), "-", ft(end,tz=tz), ">", ft(time_elapsed,tz=tz), "\n")
}

## unload all pkgs
unloadpkg = function()
  a = lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)



## input: vector
## outpute: mode, most frequently occuring element in vector
mode_set <- function(y) unique(y)[which.max(tabulate(match(y, unique(y))))]


## input: vector
## outpute: geometric mean
## from package psych
geometric.mean <- function(x,na.rm=TRUE) {
  if (is.null(nrow(x))) {
    exp(mean(log(x),na.rm=TRUE))
  } else {
    exp(apply(log(x),2,mean,na.rm=na.rm))
  }
}


## input: Sys.time() value
## output: formatted time as string; used in time_output function
ft = function(time, tz="GMT") return( format(.POSIXct(time,tz=tz), "%H:%M:%S") )


## input: x=vector of indices; n=cores or number of parts you want to split x into
## ouput: list of n indices vector
loop_ind_f <- loopInd <- function(x,n) {
  if (n==1) return(list(x))
  return(split(x, ceiling(seq_along(x)/ceiling(length(x)/n))))
}


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

## input: matrix
## output: returns u1 (col index with only 1 unique element), ua (col index where every row is a unique element), prints colnames and its unique elements if there is less than n unique elements
col_probe = function(m,n=15) {
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
  return(list(u1=u1,ua=ua))
}


## input: file path and a file extension
## output: List of all file names in given path with the given extension
fileNames = function(pathList, ext="fcs") {
  temp.str = sapply(strsplit(pathList, "/"), function(x) x[length(x)])
  pathList = sub(paste0(".", ext,"$"), "", temp.str, ignore.case=T)
  return(pathList)
}


## input: matrix
## output: matrix without all NA col/row
delna <- function(m)
  m[apply(m,1,function(x) !all(is.na(x))), apply(m,2,function(x) !all(is.na(x)))]


## input: gray scale 1d matrix
## output: rgb... repeat 3 times
gray2rgb <- function(im) {
  if (length(dim(im)==2)) {
    im = array(im, dim=c(dim(im),1))
  }
  if (length(dim(im)==3)) {
    if (dim(im)[3]==1) {
      img = abind(im,im)
      img = abind(im,im)
    }
  } else {
    rturn(NULL)
  }
  return(img)
}


#Written 20151026 by Alice Yue
#Last modified 20151026

#Note: All parameters are case sensative!


## generic functions ------------------------------

## input: list of package names to load
## output: none; load/install package
libr <- function(pkgs) {
  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0)
    install.packages(setdiff(pkgs, rownames(installed.packages())), verbose=F)
  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
    if (!requireNamespace("BiocManager", quietly=T)) install.packages("BiocManager")
    BiocManager::install(setdiff(pkgs, rownames(installed.packages())), ask=F)
  }
  sapply(pkgs, require, character.only=T)
}


## input: Sys.time() value
## output: time elapsed since input
time_output = function(start, message="", tz="GMT") {
  start = as.POSIXct(start)
  end = Sys.time()
  time_elapsed = difftime(end, start, units="secs")
  cat(message, ifelse(message=="","",": "), ftime(start,tz=tz), "-", ftime(end,tz=tz), ">", ftime(time_elapsed,tz=tz), "\n")
}

## input: Sys.time() value
## output: formatted time as string; used in time_output function
ftime = function(time, tz="GMT") return( format(.POSIXct(time,tz=tz), "%H:%M:%S") )


## unload all pkgs
unloadpkg = function()
  a = lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)


# combines two vectors in every combination
pastevector <- function(x,y,collapse="") apply(expand.grid(x, y), 1, paste, collapse=collapse)

# deletes columns where all values are the same
delcols <- function(x) x[sapply(x, function(y) length(unique(y))>1)]

# https://stackoverflow.com/questions/3443687/formatting-decimal-places-in-r
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))


# generates n colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


##Input: colour scale
##Output: adds a colour scale to plot
##https://menugget.blogspot.com/2011/08/adding-scale-to-image-plot.html#more
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}


#remove NULL elements from list
remove_null <- function(x) return(Filter(Negate(is.null), x))

# INPUT: x=vector of indices; n=cores or number of parts you want to split x into
# Ouput: list of n indices vector
loopInd <- function(x,n) {
  if (n==1) return(list(x))
  return(split(x, ceiling(seq_along(x)/ceiling(length(x)/n))))
}

# INPUT: clustering index vector (numerical)
# OUTPUT: clustering index matrix
cluster_v2m = function(x) {
  ncoll = unique(x)
  xx = Matrix(0, ncol=length(ncoll), nrow=length(x), sparse=T)
  for (i in 1:length(x)) {
    xx[i,ncoll==x[i]] = 1
  }
  colnames(xx) = ncoll
  rownames(xx) = names(x)
  return(xx)
}

# INPUT: clustering index matrix
# OUTPUT: clustering index vector (numerical)
cluster_m2v = function(xx) {
  x = NULL
  for (i in 1:nrow(xx)) {
    x[i] = which(xx[i,]>0)
  }
  return(x)
}

# INPUT: val=vector of values; group=whether these values are in different groups
# OUTPUT: density histogram plot; saves as html
dens_plot <- function(val,group=NULL,filename=NULL,title="Density with Histogram overlay",binwidth=(max(val,na.rm=T)-min(val,na.rm=T))/50,fillh="#FF7C8E",filld="#CAE2D0",multipleoverlay="stack", hist=T) {
  require(plotly)
  if (is.null(group)) {
    if (hist) {
      p <- ggplot(data.frame(x=val,group=rep("A",length(val))), aes(x)) +
        geom_histogram(aes(y = ..density..), alpha = 0.7, fill = fillh, na.rm=T, binwidth=binwidth) +
        geom_density(fill = filld, alpha = 0.5, na.rm=T) +
        theme(panel.background = element_rect(fill = '#ffffff')) +
        ggtitle(title)
    } else {
      p <- ggplot(data.frame(x=val,group=rep("A",length(val))), aes(x)) +
        geom_density(fill = filld, alpha = 0.5, na.rm=T) +
        theme(panel.background = element_rect(fill = '#ffffff')) +
        ggtitle(title)

    }
  } else {
    if (hist) {
      p <- ggplot(data.frame(x=val,group=group), aes(x, fill = group)) +
        geom_histogram(aes(y = ..density..), alpha = 0.7, position=multipleoverlay, na.rm=T, binwidth=binwidth) +
        geom_density(position = multipleoverlay, alpha = 0.3, na.rm=T) +
        theme(panel.background = element_rect(fill = '#ffffff')) +
        ggtitle(title)
    } else {
      p <- ggplot(data.frame(x=val,group=group), aes(x, fill = group)) +
        geom_density(position = multipleoverlay, alpha = 0.3, na.rm=T) +
        theme(panel.background = element_rect(fill = '#ffffff')) +
        ggtitle(title)

    }
  }
  if (is.null(filename)) {
    ggplotly(p)
  } else {
    #devtools::install_github("ropensci/plotly@fix/nse")
    htmlwidgets::saveWidget(ggplotly(p), filename)
  }
}

# INPUT: n= number of cores; x=length of equilateral triangle (e.g. distance matrix, only lower/upper triangle is unique) you want to calculate
# OUTPUT: list of inices to do in parallel
# INCOMPLETE FUNCTION
triLoopInd <- function(x,n) {
  if (n==1) return(list(c(1,x-1)))
  loopInd = 1
  perloop = 0
  if ((x-1)<=n) {
    perloop = 1:min((x-1),n)
  } else if ((x-1)<=(n*2)) {
    perloop = rep(1,n)+c(rep(0,n-((x-1)-n)),rep(1,(x-1)-n))
  } else { loopInd = ceiling(c(1,x*triSplit(n))) }

  if (length(loopInd)==1 & length(perloop)>1) {
    for (i in 2:length(perloop)) {
      loopInd[i] = loopInd[i-1]+perloop[i-1]
    }
  }
  loop.ind = NULL
  for (i in 1:length(loopInd)) {
    if (i==length(loopInd)) {
      loop.ind[[i]] = c(loopInd[i],x-1)
    } else {
      loop.ind[[i]] = c(loopInd[i],loopInd[i+1]-1)
    }
  }
  return(loop.ind)
}

#Output: proportions that split a triangle into n equal pieces
triSplit <- function(n) 1-c(sapply((n-1):1, function(x) sqrt((x)/(n))))

#Output: boolean, whether or not x is in between two values rng
'%between%' <- function(x,rng) (x>rng[1] & x<rng[2]) | (x<rng[1] & x>rng[2])


#Output: columns where all values are lower than or equal to a threshold
colIndBelow <- function(m,thres)
  which(apply(m,2,function(x) all(x<=thres)))


#Output: Random matrix
randomMatrix <- function(nrow,ncol)
  matrix(rexp(ncol*nrow, rate=.1), ncol=ncol)


#Input: package
#Output: TRUE if package is installed, FALSE otherwise
is.installed <- function(pkg){
  is.element(pkg, installed.packages()[,1])
}

charOccurences <- function(char, str) {
  str2 <- gsub(char,"",str)
  return (nchar(str) - nchar(str2))
}

#Input: file paths of matrices with column names you want to compare
#Output: list of unique column names; last vector in list corresponds to which filepath the unique column names belong to
colNamesSame <- function(pathList) {
  count <- 1
  colNameSame <- NULL
  colNameSame[[1]] <- colnames(read.csv(pathList[1]))
  colNameSameIndex <- c(1)
  for (i in 2:length(pathList)) {
    s <- colnames(read.csv(pathList[i]))
    if(!TRUE %in% (colNameSame%in%list(s)) ) {
      count <- count+1
      colNameSameIndex[count] <- i
      colNameSame[[count]] <- s
    }
  }
  colNameSame[[count+1]] <- colNameSameIndex
  return(colNameSame)
}

#Input: file path and a file extension
#Output: List of all file names in given path with the given extension
fileNames <- function(pathList, ext="fcs") {
  for (i in 1:length(pathList)) {
    temp.str <- strsplit(pathList[i], split="[/]")
    pathList[i] <- gsub(paste(".", ext, sep=""), "", temp.str[[1]][ length(temp.str[[1]]) ], ignore.case=TRUE)
  }
  return(pathList)
}

#Input: file path
#Output: List of all last folder names in given path (Genotype)
folderNames <- function(pathList) {
  folders <- c()
  for (i in 1:length(pathList)) {
    temp.str <- strsplit(pathList[i],split="/")
    folders[i] <- temp.str[[1]][length(temp.str[[1]])-1]
  }
  return(folders)
}

loadFrames <- function(pathList) {
  myFrames <- new.env()
  cat("loading", length(pathList), "files:")
  for (i in 1:length(pathList)) {
    myFrames[[as.character(i)]] <- get(load(pathList[i]))
    cat(" ", i, sep="")
  }
  return(myFrames)
}

#input: a phenotype with all markers
#output: vector of markers
getMarkers <- function(phenotype) {
  return (unlist(strsplit(as.character(phenotype), '[-+]')))
}



#given number of markers, output number of nodes, edges
getGraphInfo <- function(m) {
  npl = epl = rep(0,m)
  for (i in 1:m) {
    npl[i] = choose(m,i)*(2^i)
    epl[i] = npl[i]*i
  }
  n = 3^m
  e = n*2*m/3
  return (list(n=n, npl=npl, e=e, epl=epl))
}





## Input: PhenoMeta
## Output: Index of leaves (phenotypes with no children)
getleaves <- function (phenoMeta, no_cores) {
  require(foreach)
  require(doMC)

  registerDoMC(no_cores)

  finalLevel = which(phenoMeta$phenolevel==max(phenoMeta$phenolevel))
  notFinalLevel = setdiff(1:nrow(phenoMeta), finalLevel)
  nflleaves = foreach (i=1:length(notFinalLevel), .combine="c") %dopar% {
    pheni = unlist(strsplit(phenoMeta$phenocode[i],""))
    phenind = which(pheni!="0")
    zeroind = setdiff(1:length(pheni),phenind) #must have zeros because notFinalLevel
    childphenocode = as.vector(sapply(zeroind, function(x) { pi1=pi2=pheni; pi1[x]="1"; pi2[x] = "2"; return(c(paste(pi1,collapse=""),paste(pi2,collapse=""))) } ))
    childrenind = match(childphenocode, phenoMeta$phenocode)
    if (all(is.na(childrenind))) return(i)
    return(NA)
  }
  leaves = notFinalLevel[nflleaves[!is.na(nflleaves)]]
  return(list(leaves=leaves, finalLevel=finalLevel))
}










## Input: Phenotype
## Output: Phenotype, PhenoCode, Phenolevel (number of markers per phenotype)
getPhen <- function(phen) {
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




## Jensen-shannon-divergence (half of kullback divergence both ways)
# http://stackoverflow.com/questions/11226627/jensen-shannon-divergence-in-r

jsd <- function(p,q, list=T) { #p,q are two distributions, two lists of distributions if to avg over
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



#Output: vector of 3iTcell gating strategy phenotypes
GatingStrategyPhenotypes <- function(){
  GatingStrategyPop <- c(NA, NA, "")

  GatingStrategyPop[4] <- "TCRb-CD8a-CD161+TCRd-" #2.1
  GatingStrategyPop[5] <- "CD44+CD62L-TCRb-CD8a-CD161+TCRd-" #2.11
  GatingStrategyPop[6] <- "CD44+CD62L+TCRb-CD8a-CD161+TCRd-" #2.12
  GatingStrategyPop[7] <- "CD44+TCRb-KLRG1+CD8a-CD161+TCRd-" #2.13

  GatingStrategyPop[8] <- "TCRb-TCRd+" #1
  GatingStrategyPop[9] <- "TCRb-KLRG1+GITR-TCRd+" #1.2
  GatingStrategyPop[10] <- "CD44+CD62L-TCRb-TCRd+" #1.3
  GatingStrategyPop[11] <- "CD44+CD62L+TCRb-TCRd+" #1.4

  GatingStrategyPop[12] <- "TCRb+CD161+CD4-TCRd-" #3.1
  GatingStrategyPop[14] <- "CD44+CD62L-TCRb+CD161+CD4-TCRd-" #3.11
  GatingStrategyPop[15] <- "CD44+CD62L+TCRb+CD161+CD4-TCRd-" #3.12
  GatingStrategyPop[16] <- "TCRb+KLRG1+CD161+CD4-TCRd-" #3.13

  GatingStrategyPop[13] <- "TCRb+CD161+CD4+TCRd-" #3.2
  GatingStrategyPop[17] <- "CD44+CD62L-TCRb+CD161+CD4+TCRd-" #3.21
  GatingStrategyPop[18] <- "CD44+CD62L+TCRb+CD161+CD4+TCRd-" #3.22
  GatingStrategyPop[19] <- "TCRb+KLRG1+CD161+CD4+TCRd-" #3.23

  GatingStrategyPop[24] <- "TCRb+CD8a-CD161-CD4+TCRd-" #3.31
  GatingStrategyPop[20] <- "CD25+TCRb+CD8a-CD161-CD4+GITR+TCRd-" #3.311
  GatingStrategyPop[21] <- "CD44+CD62L-CD25+TCRb+CD8a-CD161-CD4+GITR+TCRd-" #3.3112
  GatingStrategyPop[22] <- "CD44+CD62L+CD25+TCRb+CD8a-CD161-CD4+GITR+TCRd-" #3.3111
  GatingStrategyPop[23] <- "CD44+CD25+TCRb+KLRG1+CD8a-CD161-CD4+GITR+TCRd-" #3.3113

  GatingStrategyPop[26] <- "CD25-TCRb+CD8a-CD161-CD4+TCRd-" #3.312
  GatingStrategyPop[27] <- "CD44+CD62L-CD25-TCRb+CD8a-CD161-CD4+TCRd-" #3.3121
  GatingStrategyPop[28] <- "CD62L+CD25-TCRb+CD8a-CD161-CD4+TCRd-" #3.3122
  GatingStrategyPop[29] <- "CD44+CD25-TCRb+KLRG1+CD8a-CD161-CD4+TCRd-" #3.3123

  GatingStrategyPop[25] <- "TCRb+CD8a+CD161-CD4-TCRd-" #3.32
  GatingStrategyPop[30] <- "CD44+CD62L-TCRb+CD8a+CD161-CD4-TCRd-" #3.322
  GatingStrategyPop[31] <- "CD44+CD62L+TCRb+CD8a+CD161-CD4-TCRd-" #3.321
  GatingStrategyPop[32] <- "CD44+TCRb+KLRG1+CD8a+CD161-CD4-TCRd-" #3.324

  GatingStrategyPop[33] <- "TCRb+TCRd-" #3
  GatingStrategyPop[34] <- "TCRb+CD161+TCRd-" #3.1/2 CD4?

  GatingStrategyPop[35] <- "TCRb-KLRG1-GITR+TCRd+" #1.1
  GatingStrategyPop[36] <- "CD44-CD62L+TCRb-TCRd+" #1.5
  GatingStrategyPop[37] <- "TCRb-TCRd-" #2
  GatingStrategyPop[38] <- "CD44-CD62L+TCRb+CD8a+CD161-CD4-TCRd-" #3.323

  GatingStrategyPop[39] <- "TCRb+CD161-TCRd-" #3.3
  return(GatingStrategyPop)
}


#001_Clean: Create Global Frame-------------------------------------------------

#Input: FCS, start & end index of marker columns, path of filePrefixWithDir
clean3iTcell <- function(f, MarkerStartCol, MarkerEndCol, path) {
  fi <- clean(f, vectMarkers=c(MarkerStartCol:MarkerEndCol), filePrefixWithDir=path, ext="fcs", diagnostic=TRUE)
  if (!is.character(fi)) {
    fi@exprs <- fi@exprs[which(fi@exprs[,MarkerEndCol+2] < 10000),]
    f <- fi
  }
  remove(fi)
  return(f)
}

#Note: flowDensity (used to get gating thresholds after clean) requires FSC & SSC therefore don't use this function...
#Input: FCS, start & end index of marker columns
#Output: parameters@data of FCS's marker columns only
extractextractMarkerColPar <- function(f, MarkerStartCol, MarkerEndCol) {
  f@parameters <- f@parameters[MarkerStartCol:MarkerEndCol]
  row.names(f@parameters@data) <- paste("$P", c(1:(1+MarkerEndCol-MarkerStartCol)), sep="")
  parData <- NULL
  parData <- f@parameters@data
  return(parData)
}

#Note: flowDensity (used to get gating thresholds after clean) requires FSC & SSC therefore don't use this function...
#Input: FCS, start & end index of marker columns, parameter data
#Output: FCS with exprs & parameters of only marker columns
extractMarkerCol <- function(f, MarkerStartCol, MarkerEndCol, parData=NULL) {
  f@exprs <- f@exprs[,c(MarkerStartCol:MarkerEndCol)]
  #store parameter data of 1st FCS file only (consistent) i.e. $P1 column name, desc, range, minRange, maxRange
  if(is.null("parData")) {
    parData <- extractextractMarkerColPar(f, MarkerStartCol, MarkerEndCol)
  }
  f@parameters@data <- parData

  return(f)
}

#Input: list of FCS
#Output: global frame FCS with (smpl) 1000 rows from each FCS
globalFrame <- function(fcleanFrames=NULL, fCEFilePath=NULL, smpl=1000) {
  if (!is.null(fcleanFrames)) {
    g <- fcleanFrames[[as.character(1)]]
    g@exprs <- g@exprs[sample(1:nrow(g), smpl),]
    cat("Sampling", smpl, "rows from", length(fcleanFrames), "fcs files: 1")
    for (i in 2:length(fcleanFrames)) {
      f <- fcleanFrames[[as.character(i)]]
      g@exprs <- rbind(g@exprs, f@exprs[sample(1:nrow(f), smpl),])
      cat(", ", i, sep="")
    }
  }
  else if (!is.null(fCEFilePath)) {
    g <- get(load(fCEFilePath[1]))
    g@exprs <- g@exprs[sample(1:nrow(g), smpl),]
    cat("Sampling", smpl, "rows from", length(fCEFilePath), "fcs files: 1")
    for (i in 2:length(fCEFilePath)) {
      f <- get(load(fCEFilePath[i]))
      g@exprs <- rbind(g@exprs, f@exprs[sample(1:nrow(f), smpl),])
      remove(f)
      cat(", ", i, sep="")
    }
  }
  return(g)
}

#002_Gate----------------------------------

#Input: a vector of numbers
#Output: the input vector with outliers replaced by the average of non-outlier values
replaceOutliers <- function(vec) {
  #Outlier detection using arrayMvout
  # 4 rows of random variables numbers between 5 and 5.5
  #  Create fake data that has a a small standard deviation and no outliers to use to trick ArrayOutliers to find outliers from a vector.
  # if (!exists("GoodFakeData")) {
  #   GoodFakeData <- NULL
  #   for(i in 1:4) {
  #     GoodFakeData <- cbind(GoodFakeData, runif(length(vec),5,5.5))
  #   }
  # }
  # isNA <- which(is.na(vec))
  # asDF <- data.frame(cbind(vec, GoodFakeData))
  # if(length(isNA)!=0) {
  #   asDF <- asDF[-isNA,]
  # }
  # outlier <- ArrayOutliers(asDF, alpha=0.5)$outl # find all outliers for each marker
  # outlierNo <- as.numeric(c(isNA, rownames(outlier)))
  # outlierNo <- outlierNo[c(order(outlierNo))]

  #assumes normal distribution; can use either method 1 or 2
  outlier <- getOutliers(vec, method="I", distribution="normal") #Method 1: detects outliers by checking which are below (above) the limit where according to the model distribution less then rho[1] (rho[2]) observations are expected (given length(y) observations)
  # outlier <- getOutliers(vec, method="II") #Method 2:  detects outliers by finding the observations (not used in the fit) who's fit residuals are below (above) the estimated confidence limit alpha[1] (alpha[2]) while all lower (higher) observations are outliers too
  outlierNo <- c(outlier$iRight, outlier$iLeft)
  avg <- mean(vec[-outlierNo])
  if (length(outlierNo)!= 0) {
    cat(paste("\nChanging FCS#", outlierNo, " from ", vec[outlierNo], " to ", round(avg, 4), sep=""))
    vec[outlierNo] <- avg
  }
  return(vec)
}

#Input: list of cleaned FCS
#Output: gthres matrix of thresholds for each FCS and its columns/markers
gthres3iTcell <- function(fcleanFrames=NULL, fcleanFilePath=NULL) {
  if (!is.null(fcleanFrames)) { gthres <- matrix(0, nrow=length(fcleanFrames), ncol=ncol(fcleanFrames[[as.character(1)]]@exprs) ) } #matrix of thresholds
  else if (!is.null(fcleanFilePath)) { gthres <- matrix(0, nrow=length(fcleanFilePath), ncol=16 ) }
  # if ( Genotype[i] != "+_+" && Genotype[i] != "+_Y" ){ # skip wildtypes

  #Col 10, 16---------------------------------
  cat("\nColumn 10 (TCRb), 16 (TCRd): ")
  for(i in 1:max(length(fcleanFrames), length(fcleanFilePath))){
    cat(i, " ", sep="")
    if (!is.null(fcleanFrames)) { f <- fcleanFrames[[as.character(i)]] }
    else if (!is.null(fcleanFilePath)) { f <- get(load(fcleanFilePath[i])) }
    gthres[i,10] <- deGate(f,channel=10)
    gthres[i,16] <- deGate(f,channel=16, upper=T)
  }
  #Find and replace outliers using avg; not channel 16 because the outliers of TCRd are actually fine for the 186XX files
  gthres[,10] <- replaceOutliers(gthres[,10])


  #Col 12, 13, 14---------------------------------
  gthres13 <- NULL
  cat("\nColumn 12 (CD8a), 13 (CD161), 14 (CD4): ")
  for(i in max(length(fcleanFrames), length(fcleanFilePath)):1){
    cat(i, " ", sep="")
    if (!is.null(fcleanFrames)) { f <- fcleanFrames[[as.character(i)]] }
    else if (!is.null(fcleanFilePath)) { f <- get(load(fcleanFilePath[i])) }
    F2 <- getflowFrame(flowDensity(f, channels=c(10,16), position=c(F,F), gates=c(gthres[i,10],gthres[i,16]))) #TCR-
    F3 <- getflowFrame(flowDensity(f, channels=c(10,16), position=c(T,F), gates=c(gthres[i,10],gthres[i,16]))) #ab T-cells
    gthres[i,12] <- deGate(F2,channel=12)
    gthres[i,14] <- deGate(F3,channel=14)
    # gthres[i,13] <- deGate(F3,channel=13)
    ClosestTo2p4 <- c(
      deGate(F3,channel=13, tinypeak.removal=0.001, all.cut=T),
      deGate(F3,channel=13, upper=T, tinypeak.removal=0.9))
    gthres13[i] <- which(rank(abs(ClosestTo2p4 - 2.4)) == 1)
    gthres[i,13] <- ClosestTo2p4[gthres13[i]]
  }
  #Find and replace outliers using avg
  gthres[,12] <- replaceOutliers(gthres[,12])
  gthres[,13] <- replaceOutliers(gthres[,13])
  gthres[,14] <- replaceOutliers(gthres[,14])


  #Col 7, 8, 11---------------------------------
  cat("\nColumn 7 (CD44), 8 (CD62L), 11 (KLRG1): ")
  for(i in 1:max(length(fcleanFrames), length(fcleanFilePath))){
    cat(i, " ", sep="")
    if (!is.null(fcleanFrames)) { f <- fcleanFrames[[as.character(i)]] }
    else if (!is.null(fcleanFilePath)) { f <- get(load(fcleanFilePath[i])) }
    F2 <- getflowFrame(flowDensity(f, channels=c(10,16), position=c(F,F), gates=c(gthres[i,10],gthres[i,16]))) #TCR-
    F3 <- getflowFrame(flowDensity(f, channels=c(10,16), position=c(T,F), gates=c(gthres[i,10],gthres[i,16]))) #ab T-cells
    F2.1 <- getflowFrame(flowDensity(F2, channels=c(13,12), position=c(T,F), gates=c(gthres[i,13],gthres[i,12]))) #TCR-; NK-cells
    F3.3 <- getflowFrame(flowDensity(F3, channels=c(13,14), position=c(F,NA), gates=c(gthres[i,13],gthres[i,14]))) #ab; P4
    gthres[i,7] <- deGate(F2.1,channel=7, use.percentile=T, percentile=0.001)
    # gthres[i,8] <- deGate(F2.1,channel=8)
    gthres[i,8] <- 2
    gthres[i,11] <- deGate(F2.1,channel=11)
  }
  #Find and replace outliers using avg; do channel 2 seperately,
  gthres[,7] <- replaceOutliers(gthres[,7])
  gthres[,11] <- replaceOutliers(gthres[,11])


  #Col 9, 15---------------------------------
  gthres15a <- 0
  gthres15b <- 0
  cat("\nColumn 9 (CD25), 15 (GITR): ")
  for(i in max(length(fcleanFrames), length(fcleanFilePath)):1){
    cat(i, " ", sep="")
    if (!is.null(fcleanFrames)) { f <- fcleanFrames[[as.character(i)]] }
    else if (!is.null(fcleanFilePath)) { f <- get(load(fcleanFilePath[i])) }
    F3 <- getflowFrame(flowDensity(f, channels=c(10,16), position=c(T,F), gates=c(gthres[i,10],gthres[i,16]))) #ab T-cells
    F3.3 <- getflowFrame(flowDensity(F3, channels=c(13,14), position=c(F,NA), gates=c(gthres[i,13],gthres[i,14]))) #ab; P4
    F3.31 <- getflowFrame(flowDensity(F3.3, channels=c(12,14), position=c(F,T), gates=c(gthres[i,12],gthres[i,14]))) #ab; P4; CD4+ T-cells
    gthres[i,9] <- deGate(F3.31,channel=9)
    # gthres[i,15] <- deGate(F3.31,channel=c(15))
    gthres15temp1a <- deGate(F3.31, channel=15, tinypeak.removal=0.01)
    gthres15temp1b <- deGate(F3.31, channel=15, upper=T)
    tempD <- density(F3.31@exprs[,15]) # have
    peakValue <- tempD$x[which.max(tempD$y)]
    if(gthres15temp1a > peakValue && gthres15temp1b > peakValue) {
      gthres15temp1 <- min(gthres15temp1a, gthres15temp1b)
    } else {
      gthres15temp1 <- max(gthres15temp1a, gthres15temp1b)
    }
    if(gthres15temp1==gthres15temp1a) {
      gthres15a <- gthres15a+1
    } else {
      gthres15b <- gthres15b+1
    }
    gthres[i,15] <- gthres15temp1
  }
  #Find and replace outliers using avg
  gthres[,9] <- replaceOutliers(gthres[,9])
  gthres[,15] <- replaceOutliers(gthres[,15])

  colnames(gthres) <- c(1:ncol(gthres))
  return(gthres)
}

#003_ScatterPlots -----------------------------------------------

#By Justin Meskas; modified 20151221 by Alice Yue
#3iTCell project specific (see gating strategy)
#Input: FCS, gthres[i,], path where png scatter plot should be saved
#Output: saves a png scatter lot into path (does so by first extracting required cell populations)
ScatterPlots3iTCell <- function(f, gthresi, path) {
  #Create cell populations, one for each plot; Only cell populations used for plotting are isolated here
  F1 <- getflowFrame(flowDensity(f, channels=c(10,16), position=c(F,T), gates=c(gthresi[10],gthresi[16]))) #gd T-cells
  F2 <- getflowFrame(flowDensity(f, channels=c(10,16), position=c(F,F), gates=c(gthresi[10],gthresi[16]))) #TCR-
  F3 <- getflowFrame(flowDensity(f, channels=c(10,16), position=c(T,F), gates=c(gthresi[10],gthresi[16]))) #ab T-cells

  F2.1 <- getflowFrame(flowDensity(F2, channels=c(13,12), position=c(T,F), gates=c(gthresi[13],gthresi[12]))) #TCR-; NK-cells
  F3.1 <- getflowFrame(flowDensity(F3, channels=c(13,14), position=c(T,F),  gates=c(gthresi[13],gthresi[14]))) #ab; NKT-cells
  F3.2 <- getflowFrame(flowDensity(F3, channels=c(13,14), position=c(T,T),  gates=c(gthresi[13],gthresi[14]))) #ab; iNKT-cells
  F3.3 <- getflowFrame(flowDensity(F3, channels=c(13,14), position=c(F,NA), gates=c(gthresi[13],gthresi[14]))) #ab; P4

  F3.31 <- getflowFrame(flowDensity(F3.3, channels=c(12,14), position=c(F,T), gates=c(gthresi[12],gthresi[14]))) #ab; P4; CD4+ T-cells
  F3.32 <- getflowFrame(flowDensity(F3.3, channels=c(12,14), position=c(T,F), gates=c(gthresi[12],gthresi[14]))) #ab; P4; CD8+ T-cells

  F3.311 <- getflowFrame(flowDensity(F3.31, channels=c(9,15), position=c(T,T),  gates=c(gthresi[9],gthresi[15]))) #ab; P4; CD4+ T-cells; KLRG1+
  F3.312 <- getflowFrame(flowDensity(F3.31, channels=c(9,15), position=c(F,NA), gates=c(gthresi[9],gthresi[15]))) #ab; P4; CD4+ T-cells;

  #Plot! Error :(
  tryCatch({
    png(file=path) #, width=1800, height=1200)
    #par(mar(c(5, 5, 10, 2) + 0.1))
    layout(rbind(c(1:6),c(7:11,0),c(12,13,0,0,0,0)))

    plotDens(f, c(10,16), main=paste("0_TCRb_TCRd"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[10]); abline(h=gthresi[16]);
    plotDens(F1, c(11,15), main=paste("1_KLRG1_GITTR"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[11]); abline(h=gthresi[15]);
    plotDens(F1, c(8,7), main=paste("1_CD62L_CD44"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[8]); abline(h=gthresi[7]);

    plotDens(F2, c(13,12), main=paste("2_CD161_CD8a"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[13]); abline(h=gthresi[12]);
    plotDens(F2.1, c(8,7), main=paste("2.1_CD62L_CD44"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[8]); abline(h=gthresi[7]);
    plotDens(F2.1, c(7,11), main=paste("2.1_CD44_KLRG1"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[7]); abline(h=gthresi[11]);

    plotDens(F3, c(13,14), main=paste("3_CD161_CD4"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[13]); abline(h=gthresi[14]);
    plotDens(F3.1, c(8,7), main=paste("3.1_CD62L_CD44"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[8]); abline(h=gthresi[7]);
    plotDens(F3.1, c(11,15), main=paste("3.1_KLRG1_GITR"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[11]); abline(h=gthresi[15]);
    plotDens(F3.2, c(8,7), main=paste("3.2_CD62L_CD44"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[8]); abline(h=gthresi[7]);
    plotDens(F3.2, c(11,15), main=paste("3.2_KLRG1_GITR"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[11]); abline(h=gthresi[15]);

    plotDens(F3.3, c(12,14), main=paste("3.3_CD8a_CD4"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[12]); abline(h=gthresi[14]);
    plotDens(F3.31, c(9,15), main=paste("3.31_CD25_GITR"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[9]); abline(h=gthresi[15]);
    plotDens(F3.311, c(8,7), main=paste("3.311_CD62L_CD44"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[8]); abline(h=gthresi[7]);
    plotDens(F3.311, c(11,7), main=paste("3.311_KLRG1_CD44"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[11]); abline(h=gthresi[7]);
    plotDens(F3.312, c(8,7), main=paste("3.312_CD62L_CD44"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[8]); abline(h=gthresi[7]);
    plotDens(F3.312, c(11,7), main=paste("3.312_KLRG1_CD44"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[11]); abline(h=gthresi[7]);

    plotDens(F3.32, c(8,7), main=paste("3.32_CD62L_CD44"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[8]); abline(h=gthresi[7]);
    plotDens(F3.32, c(11,7), main=paste("3.32_KLRG1_CD44"), cex.lab = 2, cex.axis = 2, cex.main=2, devn=FALSE); abline(v=gthresi[11]); abline(h=gthresi[7]);

    dev.off()
  }, error = function(err) {
    cat(" Error in plotting\n")
  })
}

flowtypeDF <- function(f, PropMarkers=7:16, Thresholds) {


  #   FP-Growth
  #
  #   Input: A database DB, represented by FP-tree constructed according to Algorithm 1
  #   Output: The complete set of frequent patterns.
  #   Method: call FP-growth(FP-tree, null).
  #   Procedure FP-growth(Tree, a) {
  #     (01) if Tree contains a single prefix path then // Mining single prefix-path FP-tree {
  #       (02) let P be the single prefix-path part of Tree;
  #       (03) let Q be the multipath part with the top branching node replaced by a null root;
  #       (04) for each combination (denoted as ??) of the nodes in the path P do
  #       (05) generate pattern ?? ?? a with support = minimum support of nodes in ??;
  #       (06) let freq pattern set(P) be the set of patterns so generated;
  #     }
  #     (07) else let Q be Tree;
  #     (08) for each item ai in Q do { // Mining multipath FP-tree
  #       (09) generate pattern ?? = ai ?? a with support = ai .support;
  #       (10) construct ????s conditional pattern-base and then ????s conditional FP-tree Tree ??;
  #       (11) if Tree ?? ?? ?? then
  #       (12) call FP-growth(Tree ?? , ??);
  #       (13) let freq pattern set(Q) be the set of patterns so generated;
  #     }
  #     (14) return(freq pattern set(P) ?? freq pattern set(Q) ?? (freq pattern set(P) ?? freq pattern set(Q)))
  #   }

  #   FP-tree construction
  #   Input: A transaction database DB and a minimum support threshold ?.
  #   Output: FP-tree, the frequent-pattern tree of DB.
  #   Method: The FP-tree is constructed as follows.
  #   Scan the transaction database DB once. Collect F, the set of frequent items, and the support of each frequent item. Sort F in support-descending order as FList, the list of frequent items.
  #   Create the root of an FP-tree, T, and label it as ??null??. For each transaction Trans in DB do the following:
  #     Select the frequent items in Trans and sort them according to the order of FList. Let the sorted frequent-item list in Trans be [ p | P], where p is the first element and P is the remaining list. Call insert tree([ p | P], T ).
  #   The function insert tree([ p | P], T ) is performed as follows. If T has a child N such that N.item-name = p.item-name, then increment N ??s count by 1; else create a new node N , with its count initialized to 1, its parent link linked to T , and its node-link linked to the nodes with the same item-name via the node-link structure. If P is nonempty, call insert tree(P, N ) recursively.

  #   FPTree
  #   One root labeled as ??null?? with a set of item-prefix subtrees as children, and a frequent-item-header table (presented in the left side of Figure 1);
  #   Each node in the item-prefix subtree consists of three fields:
  #     Item-name: registers which item is represented by the node;
  #   Count: the number of transactions represented by the portion of the path reaching the node;
  #   Node-link: links to the next node in the FP-tree carrying the same item-name, or null if there is none.
  #   Each entry in the frequent-item-header table consists of two fields:
  #     Item-name: as the same to the node;
  #   Head of node-link: a pointer to the first node in the FP-tree carrying the item-name.
  #   Additionally the frequent-item-header table can have the count support for an item. The Figure 1 below show an example of a FP-tree.
}


#004_CellPopMatrix ?? Phenotypes not working----------------------------------------------------

#By Justin Meskas; modified 20151221 by Alice Yue
#3iTcell specific (MarkerNo & markers 7:16)
#Input: Number of Markers/col, path list of flowtype-ed ft, list of phenotypes, vector of genotypes corresponding to ftFrames/ftFilePath
#Output: Matrix with all cell population proporitions row(phenotypes) col(samples; col names=genotypes)
matrixCellCountF <- function(MarkersNo, ftFrames=NULL, ftFilePath=NULL, markers, Genotype=Genotype) {
  matrixP <- matrix(0, nrow=calcNumPops(rep(2,MarkersNo), MarkersNo), ncol=length(ftFrames) )
  #Fill Matrix by col(FCS file)
  cat("reading & adding ", max(length(ftFilePath), length(ftFrames)), " flowtypes/samples: ", sep="")
  for (i in 1:length(ftFrames)) { cat(" ", i)
    if (!is.null(ftFrames)) { ft <- ftFrames[[as.character(i)]] }
    if (!is.null(ftFilePath)) { load(ftFilePath[i]) }
    matrixP[,i] <- ft@CellFreqs
  }
  Phenotypes <- unlist(lapply(ft@PhenoCodes, function(x){return( decodePhenotype(x, markers, ft@PartitionsPerMarker) )}))
  rownames(matrixP) <- Phenotypes
  colnames(matrixP) <- Genotype

  # -- Merge flowType Res begin; project specific 3iTcells ------------------------------
  # sub("CD44[+][+]","CD44+",(dude[which(regexec("CD44[+][+]", dude) !=-1)])[7:16])
  # dude[which(regexec("CD44-", dude) !=-1)]
  #
  # dude[which(regexec("CD44-", dude) !=-1)][7:16]
  # dude[which(regexec("CD44[+][+]", dude) !=-1)][7:16]
  # setdiff(dude[which(regexec("CD44[+]", dude) !=-1)], dude[which(regexec("CD44[+][+]", dude) !=-1)])[7:16]
  #
  # Proportions[which(regexec("CD44-", dude) !=-1)][7:16]
  # Proportions[which(regexec("CD44[+][+]", dude) !=-1)][7:16]
  # setdiff(Proportions[which(regexec("CD44[+]", dude) !=-1)], Proportions[which(regexec("CD44[+][+]", dude) !=-1)])[7:16]
  # -- merge flowType Res end

  #remove files that look bad ------------------------------
  # reducedIndices <- union(which(regexec("L00001865", FTKOFilePath) !=-1), which(regexec("L00001866", FTKOFilePath) !=-1))

  # FTKOFilePath_reduced <- FTKOFilePath[-reducedIndices]
  # Matrix3iTcell_reduced <- Matrix3iTcell[,-reducedIndices]
  # FTKOFilePath <- FTKOFilePath[-reducedIndices]
  # Matrix3iTcell <- Matrix3iTcell[,-reducedIndices]
  # FTGenotype <- FTGenotype[-reducedIndices]
  # FTGenotypeLong <- FTGenotypeLong[-reducedIndices]

  return(matrixP)
}




#005_Barcode --------------------------------------------------------------------


#Input: Matrix of cell counts (Phenotypes (row) vs. Samples/FCS files (colNames=Genotype)),
matrixSigINPROGRESS <- function(matrixCellCount, colCombos, compareColCombos, method, #comparedCols=list of columns to compare
                                test="wilcox", pAdj="BH", # pValThres=.05, cellCountThres=300, method="pVal"; test=c(), pAdj=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
                                colLabel) {
  if (method=="pVal") {
    #Produce Matrix of p values (Phenotypes (row) vs. Genotype (col))
    matrixP <- 1
  }
  return (matrixP)
}



#Make sure to name your columns and rows if you are going to specify cellCount and pVal Thresholds! Function will cut matrix down
matrixPValBC <- function(matrixCellCount, colCombos, colLabel,
                         cellCountThres=0, pValThres=NULL, #don't set these if you want the full matrix, these trim things down, a bit faster
                         test="wilcox", adjust=NULL) {
  #ensure all parameters exist
  if ( !(test %in% c("wilcox", "ttest")) ) {
    cat("specify whether you would like to do the wilcox test or the ttest")
    return()
  }
  #start!
  matrixPVal <- matrix(1, nrow=nrow(matrixCellCount), ncol=length(colLabel), dimnames=list(rownames(matrixCellCount), colLabel))
  cat(paste("getting pValues of ", nrow(matrixCellCount), " possible cell populations for ", ncol(matrixPVal), " genotypes: ", sep="")) #3iTCell specific
  for (j in ncol(matrixPVal):1){ cat(paste(" ", j, sep="")) #for every KO genotype (# of compares that need to be done)
    for (i in nrow(matrixPVal):1) {
      compare1 <- as.numeric(matrixCellCount[i, colCombos[[j]][[1]] ])
      compare2 <- as.numeric(matrixCellCount[i, colCombos[[j]][[2]] ])
      if (median(compare1)<cellCountThres & median(compare2)<cellCountThres) { #TRIM: if both WT and KO medians (to avoid outler influence) is < cell count threshold then change is not significant (set as 1)
        matrixPVal[i,j] <- 1
      } else {
        if (test=="wilcox") matrixPVal[i,j] <- wilcox.test(compare1, compare2)$p.value
        else if (test=="ttest") try( matrixPVal[i,j] <- t.test(compare1, compare2)$p.value )
      }
    }
  }
  matrixPVal[is.nan(matrixPVal)] <- 1
  if (!is.null(pValThres)) { #TRIM: remove rows and columns with only insignificant p values
    matrixPVal <- matrixPValTrim(matrixPVal, pValThres)
    # if (barcode==TRUE) {
    # matrixBC <- matrixBC[rownames(matrixBC)%in%rownames(matrixPVal), colnames(matrixBC)%in%colnames(matrixPVal)]
    # matrixBC[matrixPVal>=pValThres] <- 0
    # }
  }
  if (!is.null(adjust)) {
    matrixPValAdj <- NULL
    if (adjust=="BH") {
      for (i in 1:ncol(matrixPVal)) {
        matrixPValAdj <- cbind(matrixPValAdj, p.adjust(matrixPVal[,i], method= "BH"))
      }
    } else if (adjust=="bon") {
      for (i in 1:ncol(matrixPVal)) {
        matrixPValAdj <- cbind(matrixPValAdj, p.adjust(matrixPVal[,i], method= "bonferroni"))
      }
    }
    return(list(matrixPVal, matrixPValAdj))
  }

  return(matrixPVal)
}

#Remove rows and columns with only insignificant p values
matrixTrim <- function(m, thres=.05, trimGreater=T) {
  if (trimGreater) {
    trimRowIndex <- which(apply(m[,-1], 1, function(x) all(x>=thres))==T)
    trimColIndex <- which(apply(m[-1,], 2, function(x) all(x>=thres))==T)
  } else {
    trimRowIndex <- which(apply(m[,-1], 1, function(x) all(x<=thres))==T)
    trimColIndex <- which(apply(m[-1,], 2, function(x) all(x<=thres))==T)
  }
  if (length(trimRowIndex)>0) m <- m[-trimRowIndex,]
  if (length(trimColIndex)>0) m <- m[,-trimColIndex]
  return(m)
}

#given matrix, outputs matrix that means/median columns with same column name
matrixM <- function(matrixCellCount, matrixType="mean", densn=1e5) {
  uniqueCol <- unique(colnames(matrixCellCount))
  matrixM <- matrix(0, nrow=nrow(matrixCellCount), ncol=length(uniqueCol), dimnames=list(rownames(matrixCellCount), uniqueCol))
  if (matrixType=="mean") {
    require(Matrix)
    for(j in length(uniqueCol):1) { cat(" ", j)
      #del outliers/row
      matrixM[,j] <- rowMeans(matrixCellCount[,which(colnames(matrixCellCount)==uniqueCol[j])]) }
  } else if(matrixType=="med") {
    for(j in length(uniqueCol):1) { cat(" ", j)
      matrixM[,j] <- apply(matrixCellCount[,which(colnames(matrixCellCount)==uniqueCol[j])], 1, function(x) median(x)) }
  } else if(matrixType=="density") {
    for(j in length(uniqueCol):1) { cat(" ", j)
      matrixM[,j] <- apply(matrixCellCount[,which(colnames(matrixCellCount)==uniqueCol[j])], 1, function(x) which.max(density(x, n=densn)$y)) }
  }
  return(matrixM)
}

matrixMDiff <- function(matrixM) {
  M <- do.call(cbind, lapply(split(matrixM, rep(1:ncol(matrixM), each = nrow(matrixM))), function(x) return(x-matrixM[,1])) ) [,-1]
  colnames(M) <- colnames(matrixM[2,ncol(matrixM)])
  return(M)
}

#nodeInfo variable must be a global variable! k=root phenocode
brokenLatticeAdjList <- function(k) {
  kIndex <- which(nodeInfo$phenoCode==k)
  childIndex <- grep(glob2rx(gsub("0", "?", k)), nodeInfo$phenoCode)
  if (length(childIndex)==0) { edgeList <- c() } #stop condition
  else { #if we're not at leaf, keep calling
    child <- nodeInfo$phenoCode[childIndex] #source/child phenotype index
    edgeList <- data.frame(child)
    edgeList$parent <- rep(k, nrow(edgeList))
    if (nodeInfo$phenoLevel[kIndex]<(max(nodeInfo$phenoLevel)-1)) {
      for(l in intersect(edgeList[,1], phenoMeta$phenoCode[which(nodeInfo$phenoLevel==nodeInfo$phenoLevel[kIndex]+1)]) ) { # call all direct child and get their adjacency lists
        edgeListChild <- brokenLatticeAdjList(l)
        if (length(edgeListChild)>0) {
          commonChild <- intersect(edgeList[,1], edgeListChild[,1])
          if (length(commonChild>0)) { edgeList[-which(edgeList[,1]%in%commonChild==TRUE),] } #if direct child has same source, delete from current nodes' source
          edgeList <- rbind(edgeList, edgeListChild) #return list
        }
      }
    }
  }
  return(edgeList) #return list
}



#Input: Reference phenocode & its markers; phenocode and markers that need to be changed to reference
#Output: phenoMatch order (p2[phenoMatch]) & markerMatch (m2[markerMatch])
phenoCC <- function(p1,m1,p2,m2) {

  markerMatch <- match(m1, m2)
  m2new <- m2[markerMatch] #CD44&CD45 different (10), will put in same position here # match(markeri, markers) #index of marker[[i]] in markers e.g. c(1,2,5,7)
  m2NA <- which((m2%in%m1)==F)
  m2newAll <- append(m2new, m2[m2NA])

  # in p2, make NA phenocodes not in p1 & Reorder those in p1
  if (length(m2NA)>0) {
    p2excludedi <- NULL
    for (i in 1:length(m2NA)) { p2excludedi <- append(p2excludedi, which( substr( p2, start=m2NA[i], stop=m2NA[i] ) !="0")) }
    p2excludedi <- unique(p2excludedi)
    p2temp <- p2
    p2temp[p2excludedi] <- NA
    p2includedi <- c(1:length(p2))[-p2excludedi]
  } else {
    p2temp <- p2
    p2includedi <- c(1:length(p2))
  }
  for (i in p2includedi) {
    pc2new <- strsplit(p2[i],"")[[1]][markerMatch]
    pc2new[which(is.na(pc2new))] <- "0"
    p2temp[i] <- paste(pc2new, collapse="")
  }

  phenoMatch <- match(p1, p2) #p1 <- p2[phenoMatch]
  return(list(phenoMatch, markerMatch))
}




#006_RchyOptymix ??Change for more tests ---------------------------------------

#Input: Matrix with actual cell count (row(phenotypes/cell pops) vs col(samples/FCS)), Genotypes vector, unique Genotypes vector that should be tested (optional)
#Output: Matrix with p values of row(phenotypes), col(genotypes - combo samples)
matrixPValFull <- function(matrixPR, Genotype, compareCol, test="wilcox") {
  #ensure all parameters exist
  if ((test!="wilcox")&(test!="ttest")) {
    cat("please specify test = as either wilcox or ttest")
    return()
  }
  WTIndex <- which(is.na(match(Genotype,uniqueWTGT))==FALSE)

  #start
  Phenotypes <- rownames(matrixPR)
  matrixPVal <- matrix(NA, nrow=nrow(matrixPR), ncol=length(compareCol), dimnames=list(Phenotypes, uniqueKOGT)) #get rid of uniqueKOGT
  matrixPVal[1,] <- rep(1, length(compareCol))
  cat(paste("getting pValues of ", nrow(matrixPR), " possible cell populations for ", length(uniqueKOGT), " genotypes which contain 3+ samples: ", sep="")) #3iTCell specific
  for (i in 1:length(compareCol)){ #for every KO genotype (# of compares that need to be done)
    cat(paste(" ", i, sep=""))
    if (test=="wilcox") {
      for (j in 2:nrow(matrixPR)){ #2 because first row is just total cell count; for every cell population #High p values = likely no difference
        matrixPVal[j,i] <- wilcox.test(as.numeric(matrixPR[j, compareCol[[i]][[1]] ]), as.numeric(matrixPR[j, compareCol[[i]][[2]] ])) $p.value
      }
    } else if (test=="ttest") {
      for (j in 2:nrow(matrixPR)){
        try( matrixPVal[j,i] <- t.test(matrixPR[j, WTIndex], matrixPR[j, iKOIndex])$p.value )
      }
    }
  }
  matrixPVal[is.nan(matrixPVal)] <- 1
  return(matrixPVal)
}

#Output: Unique values of vec that appears >=threshold times
frequentValues <- function(vec, threshold, display=FALSE) {
  for (i in 1:(threshold-1)) {
    vec <- vec[duplicated(vec)]
    if (display==TRUE) { cat(length(unique(vec)), "phenotypes significant in ", i+1, " KO genotypes\n", sep="") }
  }
  return(vec)
}

#Plots RchyOptimyx; 3iTcell specific Marker Names
plotRchyPDF <- function(path, markersName=NULL, ylab='-log10(Pvalue)', maxScore=max(-log10(pVal)), pVal=NULL, pheno.codes, phenotypeScores=-log10(pVal), phenotypeScoresPlot=phenotypeScores, startPhenotype, pathCount=1, trimPaths=FALSE, trim.level=length(markersName)) {
  if (is.null(markerNames)) { markersName= c("CD44", "CD62L", "CD25", "TCRb", "KLRG1", "CD8a", "CD161", "CD4", "GITR", "TCRd") }
  res <- RchyOptimyx(pheno.codes, phenotypeScores, startPhenotype[1], pathCount, trimPaths, trim.level)
  if (length(startPhenotype)>1) {
    for (j in 2:length(startPhenotype)){
      res <- merge(res, RchyOptimyx(pheno.codes, phenotypeScores, startPhenotype[j], pathCount, trimPaths, trim.level))
    }
  }
  if (maxScore=="repeat") {
    maxScoreIndices <- NULL
    for (j in 1:length(res@nodes[1,])){
      if (length(which(pheno.codes==res@nodes[1,j])) == 0){
        maxScoreIndices[j] <- maxScoreIndices[j-1] #same as prev
      } else { maxScoreIndices[j] <- which(pheno.codes==res@nodes[1,j]) }
    }
    maxScore <- c(0,max(-log10(pVals[maxScoreIndices])))
  }
  pdf(file=path, height=18, width=28)
  plot(res, phenotypeScoresPlot, phenotypeCodes=pheno.codes, markersName, ylab, maxScore)
  dev.off();
}





#007_FPM ----------------------------------------------------

#Input: 2 rows to compare, distance measure
#Output: Ranking of farthes to nearest features
featDist <- function(m, dis) {
  require(vegan)
  m[1,which(m[1,]==0 & m[2,]!=0)] <- 1
  m[2,which(m[1,]!=0 & m[2,]==0)] <- 1

  return(unlist( lapply(c(1:ncol(m)), function(x) { return(vegdist(as.matrix(m[,x]), method=dis))}) ))
}

featDistList <- function(m1, m2=NULL, dis) {
  d <- NULL
  if (is.null(m2)) {
    for (i in 1:nrow(m1)) {
      for (j in (i+1):nrow(m1)) {
        d[[i]][[j]] <- featDist(m1[c(i,j),],dis)
      }
    }
  } else {
    for (i in 1:nrow(m1)) {
      for (j in 1:nrow(m2)) {
        d[[i]][[j]] <- featDist(rbind(m1[i,],m2[j,]),dis)
      }
    }
  }
  return(d)
}

createTransactionMatrix <- function(phenoCode, markers) {
  m1 <- matrix(0,nrow=length(phenoCode),ncol=length(markers))
  m2 <- matrix(0,nrow=length(phenoCode),ncol=length(markers))
  for (i in 1:length(phenoCode)) {
    pc <- as.numeric(strsplit(as.character(phenoCode[i]),"")[[1]])
    pc1 <- pc2 <- pc
    pc1[which(pc1==2)] <- 0
    pc2[which(pc2==1)] <- 0
    pc2[which(pc2==2)] <- 1
    m1[i,] <- pc1
    m2[i,] <- pc2
  }
  colnames(m1) <- paste0(markers,"-",sep="")
  colnames(m2) <- paste0(markers,"+",sep="")
  return(cbind(m1,m2))
}





getGTindex <- function(attlist,control,sampleCountThres,explist=NULL, all=F) {
  if (is.null(explist)) explist <- attlist
  if (length(control)==1) {
    controlIndex <- grep(control, attlist)
  } else {
    controlIndex <- which(attlist%in%control)
  }
  if (all | length(controlIndex)==0) {
    exp <- unique(explist)
  } else {
    exp <- unique(explist[-controlIndex]) #Unique KO Genotypes
  }
  expIndex <- list() #Index of unique KO Genotypes (one genotype per vector in the list)
  for (i in length(exp):1) { expIndex[[i]] <- which(explist==exp[i]) }
  goodexpIndex <- which( unlist(lapply(expIndex, length)) >=sampleCountThres ) #Index of unique KO genotypes that has 3+ samples available expIndex(expIndexIndex)
  return(list(attlist=attlist, control=control, controlIndex=controlIndex, exp=exp, expIndex=expIndex, goodexpIndex=goodexpIndex))
}


# from https://aurelienmadouasse.wordpress.com/2012/01/13/legend-for-a-continuous-color-scale-in-r/
legend.col <- function(col, lev){
  opar <- par
  n <- length(col)
  bx <- par("usr")
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000, bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n

  xx <- rep(box.cx, each = 2)

  par(xpd = TRUE)
  for(i in 1:n){
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25)
  par <- opar
}












# http://stackoverflow.com/questions/4787332/how-to-remove-outliers-from-a-dataset
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}





































# Gating Strategy populations --------------------------------------------

#SangerTCellSPLEEN
getGSsanger <- function() {
  GatingStrategyPop <- NULL

  GatingStrategyPop[3] <- (("") )

  GatingStrategyPop[4]  <- (("TCRb-CD8a-CD161+TCRd-") )
  GatingStrategyPop[5]  <- (("CD44+CD62L-TCRb-CD8a-CD161+TCRd-") )
  GatingStrategyPop[6]  <- (("CD44+CD62L+TCRb-CD8a-CD161+TCRd-") )
  GatingStrategyPop[7]  <- (("CD44+TCRb-KLRG1+CD8a-CD161+TCRd-") )

  GatingStrategyPop[9]  <- (("TCRb-KLRG1+GITR-TCRd+") )
  GatingStrategyPop[10] <- (("CD44+CD62L-TCRb-TCRd+") )
  GatingStrategyPop[11] <- (("CD44+CD62L+TCRb-TCRd+") )

  GatingStrategyPop[12] <- (("TCRb+CD161+CD4-TCRd-") )
  GatingStrategyPop[14] <- (("CD44+CD62L-TCRb+CD161+CD4-TCRd-") )
  GatingStrategyPop[15] <- (("CD44+CD62L+TCRb+CD161+CD4-TCRd-") )
  GatingStrategyPop[16] <- (("TCRb+KLRG1+CD161+CD4-TCRd-") )

  GatingStrategyPop[13] <- (("TCRb+CD161+CD4+TCRd-") )
  GatingStrategyPop[17] <- (("CD44+CD62L-TCRb+CD161+CD4+TCRd-") )
  GatingStrategyPop[18] <- (("CD44+CD62L+TCRb+CD161+CD4+TCRd-") )
  GatingStrategyPop[19] <- (("TCRb+KLRG1+CD161+CD4+TCRd-") )

  GatingStrategyPop[24] <- (("TCRb+CD8a-CD161-CD4+TCRd-") )
  GatingStrategyPop[20] <- (("CD25+TCRb+CD8a-CD161-CD4+GITR+TCRd-") )
  GatingStrategyPop[21] <- (("CD44+CD62L-CD25+TCRb+CD8a-CD161-CD4+GITR+TCRd-") )
  GatingStrategyPop[22] <- (("CD44+CD62L+CD25+TCRb+CD8a-CD161-CD4+GITR+TCRd-") )
  GatingStrategyPop[23] <- (("CD44+CD25+TCRb+KLRG1+CD8a-CD161-CD4+GITR+TCRd-") )

  GatingStrategyPop[26] <- (("CD25-TCRb+CD8a-CD161-CD4+TCRd-") )
  GatingStrategyPop[27] <- (("CD44+CD62L-CD25-TCRb+CD8a-CD161-CD4+TCRd-") )
  GatingStrategyPop[28] <- (("CD62L+CD25-TCRb+CD8a-CD161-CD4+TCRd-") )
  GatingStrategyPop[29] <- (("CD44+CD25-TCRb+KLRG1+CD8a-CD161-CD4+TCRd-") )

  GatingStrategyPop[25] <- (("TCRb+CD8a+CD161-CD4-TCRd-") )
  GatingStrategyPop[30] <- (("CD44+CD62L-TCRb+CD8a+CD161-CD4-TCRd-") )
  GatingStrategyPop[31] <- (("CD44+CD62L+TCRb+CD8a+CD161-CD4-TCRd-") )
  GatingStrategyPop[32] <- (("CD44+TCRb+KLRG1+CD8a+CD161-CD4-TCRd-") )

  GatingStrategyPop[34] <- (("TCRb+CD161+TCRd-") )

  GatingStrategyPop[35] <- (("TCRb-KLRG1-GITR+TCRd+") )
  GatingStrategyPop[36] <- (("CD44-CD62L+TCRb-TCRd+") )
  GatingStrategyPop[38] <- (("CD44-CD62L+TCRb+CD8a+CD161-CD4-TCRd-") )

  GatingStrategyPop[39] <- (("TCRb+CD161-TCRd-") )

  GatingStrategyPop <- GatingStrategyPop[-which(is.na(GatingStrategyPop))]
  GatingStrategyPop <- gsub("TCRb[-+]","",GatingStrategyPop)
  GatingStrategyPop <- gsub("CD8a","CD8",GatingStrategyPop)
  GatingStrategyPop <- gsub("KLRG1","KLRG",GatingStrategyPop)


  return (GatingStrategyPop)
}
































getPhenIndex <- function(phenName, phenoMeta, markers) {
  phenIndex <- which(phenoMeta$phenotype==phenName)
  if (length(phenIndex)==0) {
    phenNameAlph <- strsplit(phenName,"[-+]")[[1]]
    phenNameEtc <- strsplit(phenName,"[0-9A-z]")[[1]]
    phenNameEtc <- phenNameEtc[which(phenNameEtc!="")]
    orders <- NULL
    for (i in 1:length(phenNameAlph)) {
      orders[i] <- which(markers==phenNameAlph[i])
    }
    orders <- order(orders)
    phenName2 <- NULL
    for (i in 1:length(orders)) {
      phenName2 <- paste(phenName2, phenNameAlph[orders[i]], phenNameEtc[orders[i]], sep="")
    }
    phenIndex <- which(phenoMeta$phenotype==phenName2)
  }

  return(phenIndex)
}


sortByDate <- function(phenData) {
  require(lubridate)
  phenData$date <- ymd(phenData$date)
  barcodes <- str_extract(phenData$fileName,"L[0-9]+")
  barcodes <- as.numeric(str_extract( barcodes, "[0-9]+" ))
  timeOrder <- order(phenData$date,barcodes)
  all(order(barcodes)==c(1:length(barcodes)))
  phenData <- phenData[timeOrder,]
}



#008_Dist_plot --------------------------------------------------


plotTsne = function(x, continuous=F, main, colBarWidth=.08, colBarTitle="", leg=T, text=F) {
  require(lubridate)
  x1 = rownames(x)
  if (is.Date(x1) | is.timepoint(x1)) {
    x1 = as.numeric(x1-min(x1))
  }
  x = x[order(x1),]
  x1 = sort(x1)
  if (continuous) {
    ts = as.numeric(factor(x1))
    colour = heat.colors(length(unique(ts)))
    plot(x, col=colour[ts], pch=20, main=main)
    legend.col(col = colour, lev = ts)
  } else {
    c1 = 7
    if (length(unique(x1))/7>25) {
      c1 = ceiling(length(unique(x1))/25)
    }
    colour = rainbow(c1)
    cp = expand.grid(c(1:c1),c(20,1:19,21:25))
    plot(x, t='n', main=main)
    if (text) {
      text(x,rownames(x))
    } else {
      for (g in 1:length(unique(x1))) {
        points(x[which(x1%in%unique(x1)[g]),1],x[which(x1%in%unique(x1)[g]),2], col=colour[cp[g%%nrow(cp),1]], pch=cp[g%%nrow(cp),2])
      }
    }
    if (leg) legend("topleft", legend=unique(x1), fill=colour[rep(c(1:c1),g)], pch=rep(c(20,1:19,21:25),each=c1)[c(1:g)])
  }
}










## check matrix for anything weird, add onto name
checkm <- function(m,fname) {
  d = as.matrix(m)
  if (sum(is.na(d))>0) fname = paste0(fname,"_NA")
  if (sum(is.nan(d))>0) fname = paste0(fname,"_NaN")
  isinf = sum(d==Inf)
  if (!is.na(isinf)) { if (isinf>0) fname = paste0(fname,"_Inf")}
  return(fname)
}

## Input: Matrix
## Output: coordinates of na in a matrix http://stackoverflow.com/questions/4833300/get-positions-for-nas-only-in-the-middle-of-a-matrix-column
# x1 is an indexing sequence which increments at non-NA entries, counting top-to-bottom. x2 is the same, counting backward (bottom-to-top). They are only both nonzero at internal entries enclosed by non-NAs on both top and bottom, hence their non-NA indices counting in both directions are >0. Finally gate that with a & to filter out just the internal NAs. x1 = nonNA_idx.tb, x2 = nonNA_idx.bt
checkmna <- function(m) {
  return(which(is.na(m),arr.ind=T))
  # z = as.matrix(m)
  # isNA <- is.na(z)
  # # Vertical index which increments at non-NA entries, counting top-to-bottom:
  # nonNA_idx.tb <- apply(!isNA, 2, cumsum)
  # # Vertical index which increments at non-NA entries, counting bottom-to-top:
  # nonNA_idx.bt <- apply(!isNA, 2, function(x) { rev(cumsum(rev(x))) })
  # which(isNA & nonNA_idx.tb>0 & nonNA_idx.bt>0, arr.ind=TRUE)
}


## trim a list of matrices according to reference matrix, colnames and rownames must match.
## OLD
trimlistmatrix = function(refmatrix,refdel=0, targetlist,targetdel=0, ncore=1) {
  require(foreach)
  require(doMC)
  registerDoMC(ncore)

  pt0 = union(names(targetlist), Reduce('union',lapply(targetlist,function(y) return(colnames(y)))) )
  lpi = intersect(pt0, colnames(refmatrix))
  targetlistTRIM = targetlist[names(targetlist)%in%colnames(refmatrix)]
  nrows = nrow(targetlist[[1]])#original
  targetlistTRIM = foreach (i = 1:length(targetlistTRIM)) %dopar% {
    a = as.matrix(targetlistTRIM[[i]],nrow=nrows)
    dimnames(a) = dimnames(targetlistTRIM[[i]])
    row = !is.na(match(rownames(a),rownames(refmatrix)))
    col = !is.na(match(colnames(a),colnames(refmatrix)))
    if (sum(col)>0) {
      a0 = a
      a = as.matrix(a0[row,col],ncol=sum(col))
      colnames(a) = colnames(a0)[col]
      rownames(a) = rownames(a0)[row]
      a[ matrix(refmatrix[,!is.na(match(colnames(refmatrix),colnames(a)))],ncol=sum(col)) ==refdel ] = targetdel
      # if (sum(col)>0) a[ refmatrix[!is.na(match(rownames(refmatrix),rownames(a)))m, !is.na(match(colnames(refmatrix),colnames(a)))]==0 ] = targetdel
    } else { return(NULL) }
    return(Matrix(a,sparse=T))
  }
  names(targetlistTRIM) = names(targetlist)[names(targetlist)%in%colnames(refmatrix)]
  targetlistTRIM = targetlistTRIM[!sapply(targetlistTRIM, is.null)]
  return(targetlistTRIM)
}

## Input : list of matrix/matrixList with phenotype on colummn, samples on rows (filenames)
## Output: Trimmed version of matrix
## OLD
Loadintermatrices <- function(mfilenames,verbose=T) {
  require(Matrix)
  mind = 1
  mml = list()
  mmlname = c()
  if (verbose) cat("loading feature matrices, ",sep="")
  for (mind in 1:length(mfilenames)) {
    mcp = mfilenames[mind]
    if (!file.exists(mcp)) {cat("doesn't exist"); break}
    mml[[mind]] = get(load(mcp))
    if (data.class(mml[[mind]])=="dist") mml[[mind]] = as.matrix(mml[[mind]])
    mmlname[mind] = mcp
    if (mind==1) {
      if (is.null(dim(mml[[1]]))) {
        pt = union(names(mml[[1]]), Reduce('union',lapply(mml[[1]],function(y) return(colnames(y)))) )
        gt = rownames(mml[[1]][[1]])
      } else {
        pt = colnames(mml[[1]])
        gt = rownames(mml[[1]])
      }
    } else {
      if (is.null(dim(mml[[mind]]))) {
        pt0 = union(names(mml[[mind]]), Reduce('union',lapply(mml[[mind]],function(y) return(colnames(y)))) )
        gt0 = rownames(mml[[mind]][[1]])
      } else {
        pt0 = colnames(mml[[mind]])
        gt0 = rownames(mml[[mind]])
      }
      pt = intersect(pt, pt0)
      gt = intersect(gt, gt0)
    }
  }

  #trim matrices (sample) based on samples/phenotypes in common
  if (length(mml)>1) {
    if (verbose) cat("trimming, ",sep="")
    for (i in 1:length(mml)) {
      b = b0 = mml[[i]]
      if (is.null(dim(b))) {
        gt0 = rownames(b[[1]])
        pt0 = union(names(b), Reduce('union',lapply(b,function(y) return(colnames(y)))) )
        ind = match(gt,gt0)
        b = foreach (j = 1:length(b)) %dopar% {
          ptind = match(pt,colnames(b[[j]]))
          if(sum(!is.na(ptind))>0){
            ptind = ptind[!is.na(ptind)]
            # if (is.null(dim(b[[j]]))) {
            #   a = as.matrix(b[[j]],ncol=1)[ind,ptind]
            #   dimnames(a) = dimnames(b[[j]][ind,ptind])
            # } else {
            a = b[[j]][ind,ptind]
            if (is.null(dim(a))) a = matrix(a,nrow=length(ind))
            rownames(a) = rownames(b[[j]])[ind]
            colnames(a) = colnames(b[[j]])[ptind]
            # }
          } else { a = NULL }
          return(a)
        }
        names(b) = names(b0)
        map = match(pt,names(b0))
        b1 = b[map[!is.na(map)]]
        b1ind = !sapply(b1, is.null)
        if (sum(b1ind)>0) {
          mml[[i]] = b1[which(b1ind)]
        } else {
          mml[[i]] = NULL
        }
      } else {
        gt0 = rownames(b)
        pt0 = colnames(b)
        mml[[i]] = b[match(gt,gt0),match(pt,pt0)]
      }
    }
    if (length(mml)<length(mmlname)) { names(mml) = mmlname[1:length(mml)]
    } else { names(mml) = mmlname }
    mml = mml[!sapply(mml, is.null)]
  } else { names(mml) = mfilenames }
  if (verbose) cat("done")
  return(list(mml=mml,pt=pt,gt=gt))
}


## OLD
trimMatrices <- function(mml0,m0,pt,gt,phenoMeta,leavePhenotype,doneAll, countThres,k) {
  lowCountpt = colnames(m0)[apply(m0,2,function(x) all(x<=countThres))]
  if (!length(lowCountpt)>0) lowCountpt = c()

  if (!is.null(dim(phenoMeta))) phenolevel = phenoMeta$phenolevel

  highLevelpt = colnames(m0)[phenolevel>k]
  if (!length(highLevelpt)>0) highLevelpt = c()

  # cat("length of countthres: ", length(lowCountInd))
  # cat("\nlength of highlevelind: ", length(highLevelInd))
  ## Load & fix cell count/countAdj/proportion matrix -----------------------------------------------------
  dpi0 = union(lowCountpt,highLevelpt)
  # cat("\nlength of dpi: ", length(dpi0), "\n")

  lpi = setdiff(pt,dpi0)

  #check if indices already calculated for on this matrix
  if (Position(function(x) identical(x, lpi), leavePhenotype, nomatch = 0) > 0) {cat("-skipped ", sep=""); return(0)}
  leavePhenotype[[paste0(k="k.",k,"_countThres.",countThres)]] = lpi

  if (!doneAll & length(lpi)==length(pt)) doneAll = T


  #trim matrices (phenotype)
  mmlname = names(mml0)
  mml = mml0
  for (i in 1:length(mml)) {
    b0 = b = mml[[i]]
    if (is.null(dim(b))) {
      pt0 = union(names(b), Reduce('union',lapply(b,function(y) return(colnames(y)))) )
      b = foreach (j = 1:length(b)) %dopar% {
        ptind = match(lpi,colnames(b[[j]]))
        if(sum(!is.na(ptind))>0){
          ptind = ptind[!is.na(ptind)]
          # if (is.null(dim(b[[j]]))) {
          #   a = as.matrix(b[[j]],ncol=1)[ind,ptind]
          #   dimnames(a) = dimnames(b[[j]][ind,ptind])
          # } else {
          a = b[[j]][,ptind]
          if (is.null(dim(a))) {
            a = matrix(a,nrow=nrow(b[[j]]))
            rownames(a) = rownames(b[[j]])
            colnames(a) = colnames(b[[j]])[ptind]
          }
        } else { a = NULL }
        return(a)
      }
      names(b) = names(b0)
      map = match(lpi,names(b0))
      b1 = b[map[!is.na(map)]]
      b1ind = !sapply(b1, is.null)
      if (sum(b1ind)>0) {
        mml[[i]] = b1[which(b1ind)]
      } else {
        mml[[i]] = NULL
      }
    } else {
      pt0 = colnames(b)
      mml[[i]] = b[,match(lpi,pt0)]
    }
  }
  if (!length(mml)>0) return(NULL)
  if (length(mml)<length(mmlname)) { names(mml) = mmlname[1:length(mml)]
  } else { names(mml) = mmlname }
  mml = mml[!sapply(mml, is.null)]

  if (!is.null(dim(phenoMeta))) {
    pm = phenoMeta[match(lpi,phenoMeta$phenotype),]
    return(list(mml=mml,pm=pm,leavePhenotype=leavePhenotype,doneAll=doneAll))
  } else {
    return(list(mml=mml,leavePhenotype=leavePhenotype,doneAll=doneAll))
  }
}





































Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}



## Input: dissimilarity Matrix
## Output: similarity matrix
get_graph = function(d0) {
  d0 = as.matrix(d0)
  d1 = d0
  rownames(d1) = colnames(d0) = rownames(d0)
  d1 = 1/d0
  d1[d0==0] = 1/min(d0[d0>0])
  diag(d1) = 0
  return(d1)
}

get_graphd = function(s0) {
  df = 1/s0
  df[s0==0] = 1/min(s0[s0>0])
  diag(df) = 0
  df = as.dist(df)
  return(df)
}




## Input: matrices with items labelled on row + variables corresponding to those rows to split by
## Output: list of original matrix
splitmatrix = function(m,s=NULL,split=T,defaultlabel="all",missym=T) { #only choose one of avgall=F,includeall=F, else will do includeall=T only
  mm = list()
  if (length(unique(s))>1 & split) {
    for (tubei in 1:length(unique(s))) {
      index = which(s==unique(s)[tubei])
      if (missym) { mm[[tubei]] = m[index,index]
      } else { mm[[tubei]] = m[index,] }
    }
    names(mm) = unique(s)
    # if (avgall) {
    #   mm1 = Reduce("+",mm)/length(mm)
    #   mm[[defaultlabel]] = mm1
    # }
    # if (includeall) {
    #   mm[[defaultlabel]] = m
    # }
  } else {
    mm[[1]] = m
    names(mm) = defaultlabel
  }
  return(mm)
}





## Input: distance matrix filenames
## Output: distMeta
distMetafun = function(distmfile, dis,features=NULL) {
  distmfilenames = fileNames(distmfile)
  dmfns = strsplit(distmfilenames,"_")
  dmfns = lapply(1:length(dmfns),function(y) {
    x = dmfns[[y]]
    disind = which(x%in%dis)
    if (suppressWarnings(is.na(as.numeric(x[disind-1]))) ) { typeind = 1:(disind-1); rand = 0
    } else { typeind = 1:(disind-2); rand = x[disind-1] }
    return(c(type=paste(x[typeind],collapse="_"), rand=rand, x[(disind):length(x)]))
  })
  path = distmfile
  distMeta = as.data.frame(path)
  distMeta$type = sapply(dmfns,function(x) x[1])
  distMeta$layer = as.numeric(gsub("layer","",str_extract(distmfilenames, "layer[0-9]+")))
  distMeta$norm = gsub("normalize-","",str_extract(distmfilenames, "normalize-[a-z]+"))
  distMeta$count = as.numeric(gsub("countThres-","",str_extract(distmfilenames, "countThres-[0-9]+")))
  distMeta$weighted = sapply(distmfilenames,function(x) return(grepl("weighted",x)))
  distMeta$weightedorig = sapply(distmfilenames,function(x) return(grepl("orig-weighted",x)))
  distMeta$dist = sapply(dmfns,function(x) return(x[x%in%dis]))
  distMeta$feature = sapply(dmfns,function(x) x[1])
  distMeta$rand = sapply(dmfns,function(x) as.numeric(x[2]))
  distMeta$rw = sapply (dmfns, function(x) any(grepl("_rw",x)))
  distMeta$sim = grepl("simmatrix",distmfilenames)
  if(!is.null(features)) {
    feature = paste0("^",features)
    distMeta$feature = sapply(distMeta$type, function(x) features[sapply(feature,function(y) grepl(y,x,ignore.case=T))] )
    distMeta$featend = sapply(1:nrow(distMeta), function(i) gsub(distMeta$feature[i],"",distMeta$type[i]) )
  }

  return(distMeta)
}








#last n characters
substrRight = function(x, n) substr(x, nchar(x)-n+1, nchar(x))














### Matrix (sample x phenotype) Preprocessing ########################################

##Input: matrix / list of matrices with columnname = cell population;;; layer to prune at, each cell population must have atleast goodcount number of samples with more than countThres cells, each class must have more thatn good_sample amount of samples or else prune
##Output: trimmed matrix / list of matrices
trimMatrix <- function(m0,TRIM=T,mc=NULL,sampleMeta,sampleMeta_to_m1_col="id", target_col="class" ,control="control", order_cols=NULL, colsplitlen=NULL, k=NULL, konly=F,countThres=1000, goodcount=3, good_sample=3) {

  TRIM=T;mc=mc;sampleMeta=meta_file;sampleMeta_to_m1_col="id"; target_col="class" ;control="control"; order_cols=NULL; colsplitlen=NULL; k=4; konly=F;countThres=1000; goodcount=3; good_sample=3

  require(stringr)
  #get matrix colnames/rownames (sample/celltype features)
  m0cn = colnames(m0)
  m0cn = sapply(str_split(m0cn,"_"), function(x) ifelse(length(x)>1, x[2], x[1]) )
  m0rn = rownames(m0)

  #get feature layers if feature names represent cell types
  colsplitlen = cell_type_layers(m0cn)
  if (is.null(colsplitlen)) print("missing cell type features or column names; columns untrimmed.")

  #get to-delete low count phenotype indices; CountAdj should be first one
  colIndexC = rep(T,length(m0cn))
  rowIndexC = rep(T,length(m0rn))
  if(countThres>0 & !is.null(mc)) {
    mcordercol = match(m0cn,colnames(mc))
    mcordercol = mcordercol[!is.na(mcordercol)]
    mcorderrow = match(m0rn,rownames(mc))
    mcorderrow = mcorderrow[!is.na(mcorderrow)]
    mc0 = mc[mcorderrow,mcordercol]
    colIndexC = apply(as.matrix(mc0), 2, function(x) sum(x>=countThres)>=goodcount)
    rowIndexC = apply(as.matrix(mc0), 1, function(x) sum(x>=countThres)>=goodcount)
  }

  #get to-delete high no of marker phenotypes
  colIndexL = rep(T,length(m0cn))
  if (!is.null(k) & k != 0) {
    if (k<max(colsplitlen)) {
      if (konly) {
        colIndexL = colsplitlen == k
      } else {
        colIndexL = colsplitlen <= k
      }
    }
  }



  #trim matrix based on cell count and layer
  # if (m0list) {
  #   m = m0[colIndexC & colIndexL]
  #   if (TRIM & F) {
  #     m = lapply(m, function(y0) {
  #       colIndex0 = apply(y0, 2, function(x) sum(x!=0)>0 )
  #       # rowIndex0 = apply(y, 1, function(x) sum(x!=0)>0 )
  #       # y = y[rowIndex0,colIndex0]
  #       y = y0[,colIndex0]
  #       if (length(colIndex0)==1) {
  #         y = matrix(y,ncol=1)
  #         rownames(y) = rownames(y0)
  #         colnames(y) = colnames(y0)[colIndex0]
  #       }
  #       if (sum(dim(y))>0) return(y)
  #       return(NULL)
  #     })
  #     m = Filter(Negate(is.null), m)
  #   }
  #
  # } else {
  m = m0[rowIndexC, colIndexC & colIndexL]
  if (TRIM) {
    colIndex0 = apply(m, 2, function(x) any(x>0) )
    rowIndex0 = apply(m, 1, function(x) any(x>0) )
    m = m[rowIndex0,colIndex0]
  }
  # }

  #check matrix; if matrix is a list, merge; if matrix is empty, skip
  if (is.null(m)) return (NULL)
  if (all(m==0)) return (NULL)


  sm = sampleMeta
  if (!is.null(sampleMeta)) {
    # if (m0list) {
    #   smorder = match(rownames(m[[1]]),sampleMeta[,sampleMeta_to_m1_col])
    # } else {
    smorder = match(rownames(m),as.matrix(sampleMeta)[,sampleMeta_to_m1_col])
    # }
    smorder = smorder[!is.na(smorder)]
    sm = sampleMeta[smorder,]

    #order samples by date etc.
    if (!is.null(order_cols)) {
      for (order_col in order_cols) { # different centres have different filename colnames
        if (length(grep(order_col, colnames(sm)))>0) sm = sm[order(sm[,order_col]),]
      }
      # if (m0list) {
      #   mrowind = match(sm[,sampleMeta_to_m1_col],rownames(m[[1]]))
      #   morder = mrowind[!is.na(mrowind)]
      #   m = lapply(m, function(x) {
      #     xo = x[morder,]
      #     if (is.null(dim(xo))) xo = matrix(xo, ncol=1, dimnames=list(names(xo),colnames(x)))
      #     return(xo)
      #   })
      # } else {
      mrowind = match(sm[,sampleMeta_to_m1_col],rownames(m))
      m = m[mrowind[!is.na(mrowind)],]
      # }
    }

    #exclude genotypes with less than 3 samples
    if (good_sample>1 & !is.null(target_col)) {
      g = getGTindex(sm[,target_col], control, good_sample)
      goodind = as.vector(c(unlist(g$expIndex[g$goodexpIndex]), g$controlIndex))
      goodind = sort(goodind)
      sm = sm[goodind,]
      # if(m0list) {
      #   m = lapply(m, function(x) {
      #     xo =  x[goodind,]
      #     if (is.null(dim(xo))) xo = matrix(xo, ncol=1, dimnames=list(names(xo),colnames(x)))
      #     return(xo)
      #   })
      # } else {
      m = m[goodind,]
      # }
    }
  }

  return(list(m=m,sm=sm))
}







##Input: matrix
##Output: if features (rownames) have +/- symbols, returns corresponding feature layers; else returns NULL
cell_type_layers <- function(cell_types) {
  colsplit = str_split(cell_types,"[+-]")
  colsplitlen = sapply(colsplit, function(x) length(x)) - 1
  return(colsplitlen)
}






##Input: matrix with rows in order of date/time etc.
##Output: kalman filter for each column in matrix
kmf <- function(m,cols) {
  fkffitall = NULL
  statsfitall = NULL
  for (i in cols) {
    yall = as.numeric(m[,i])
    if (length(unique(yall))==1) {
      fkffitall[[i]] <- statsfitall[[i]] <- yall
    } else {
      ## Set constant parameters:
      dt <- ct <- matrix(0)
      Zt <- Tt <- matrix(1)
      a0 <- yall[1]           # Estimation of the first sample count
      P0 <- matrix(100)     # Variance of 'a0'
      ## Estimate parameters 23min if TS
      fit.fkf <- optim(c(HHt = var(yall, na.rm = TRUE) * .5,
                         GGt = var(yall, na.rm = TRUE) * .5),
                       fn = function(par, ...) -fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik,
                       yt = rbind(yall), a0 = a0, P0 = P0, dt = dt, ct = ct,
                       Zt = Zt, Tt = Tt, check.input = FALSE)
      ## Filter Nile data with estimated parameters:
      fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, HHt = matrix(fit.fkf$par[1]), GGt = matrix(fit.fkf$par[2]), yt = rbind(yall))
      fkffitall[[i]] <- fkf.obj$att[1,]
      ## Compare with the stats' structural time series implementation: 5min if TS
      statsfitall[[i]] <- fitted(StructTS(yall, type = "level"))
    }
  }
  return(list(fkffitall=fkffitall,statsfitall=statsfitall))
}









## Input: clustering filepath (currently only biclust)
## Output: original feature matrix
get_feat_matrix <- function(clust_fileName, feat_dir, mc, meta_file, id_col, target_col, control, order_cols=NULL, good_count=NULL, good_sample=NULL) {
  # get original feature matrix and meta file
  x = str_split(clust_fileName,"_")[[1]]
  # bcmethod = x[1]
  feature = x[2]
  layer = as.numeric(gsub("layer","",x[4]))
  countThres = as.numeric(str_split(x[5],"-")[[1]][2])
  m0 = get(load(paste0(feat_dir,"/",feature,".Rdata")))

  mm = trimMatrix(m0,TRIM=T, mc=mc, sampleMeta=meta_file, sampleMeta_to_m1_col=id_col, target_col=target_col, control=control, order_cols=order_cols, colsplitlen=NULL, k=layer, countThres=countThres, goodcount=good_count, good_sample=good_sample)
  m_ordered = mm$m
  meta_file_ordered = mm$sm

  #split up analysis by tube etc.
  split_col = str_split(str_split(x[3],"-")[[1]][2],"[.]")[[1]][1]
  if (split_col=="") {
    split_ind = list(all = 1:nrow(meta_file_ordered))
  } else {
    split_ids = unique(meta_file_ordered[,split_col])
    split_ids = split_ids[!is.na(split_ids)]
    split_ind = lapply(split_ids, function(split_id) which(meta_file_ordered[,split_col]==split_id) )
    names(split_ind) = split_ids
  }
  tube = str_split(x[3],"[.]")[[1]][2]
  m = m_ordered[split_ind[[tube]],]
  sm = meta_file_ordered_split = meta_file_ordered[split_ind[[tube]],]
  return(list(m=m,sm=sm))
}

get_feat_matrix2 <- function(clust_fileName=NULL, feat_dir=NULL, meta_file=NULL, id_col=NULL, row_names, col_names=NULL, getcsv=F) {
  # get original feature matrix and meta file
  getm = F
  if (!is.null(clust_fileName) & !is.null(feat_dir)) {
    getm = T

    x = str_split(clust_fileName,"_")[[1]]
    # bcmethod = x[1]
    feature = x[2]
    if (getcsv) {
      m0 = read.csv(paste0(feat_dir,"/",feature,".csv"),row.names=1, check.names=F)
    } else {
      m0 = get(load(paste0(feat_dir,"/",feature,".Rdata")))
    }
    morder = match(row_names,rownames(m0))
    morder = morder[!is.na(morder)]
    m = m0[morder,]
  }

  if (!is.null(col_names) & getm) {
    mcorder = match(col_names,colnames(m0))
    mcorder = mcorder[!is.na(mcorder)]
    m = m[,mcorder]
  }

  getsm = F
  if (!is.null(meta_file) & !is.null(id_col)) {
    getsm = T

    if (getm) {
      smorder = match(rownames(m),meta_file[,id_col])
    } else {
      smorder = match(row_names,meta_file[,id_col])
    }
    smorder = smorder[!is.na(smorder)]
    sm = meta_file[smorder,]
  }
  if (getm & getsm) return(list(m=m,sm=sm))
  if (!getm & getsm) return(list(sm=sm))
  if (getm & !getsm) return(list(m=m))
}












### Dist Score ########################################################################

## Input: Distance matrix, class list (numeric)
## Output: NCA score
NCA_score <- function(x,y,delta=1,fast=F,doUnderflow=T) {
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









### Feature Generation ########################################################################



## Input: meta_cell ($phenolevel, $phenotype, $phenocode)
## Output: Index of list of children
getphenoChild <- function (meta_cell, no_cores=1) {
  require(foreach)
  require(doMC)

  registerDoMC(no_cores)

  notFinalLevel = which(meta_cell$phenolevel!=max(meta_cell$phenolevel))
  ppcc = foreach (i=1:nrow(meta_cell)) %dopar% {
    if (!i%in%notFinalLevel) return(list(pc=NULL,pcpn=NULL))
    pheni = unlist(strsplit(meta_cell$phenocode[i],""))
    phenind = which(pheni!="0")
    zeroind = setdiff(1:length(pheni),phenind)
    childphenocode = as.vector(sapply(zeroind, function(x) { pi1=pi2=pheni; pi1[x]="1"; pi2[x] = "2"; return(c(paste(pi1,collapse=""),paste(pi2,collapse=""))) } ))
    childrenind = match(childphenocode, meta_cell$phenocode)
    childrenind = childrenind[!is.na(childrenind)]
    if (length(childrenind)>0) { #if non-leaf node; this condition is only used when full hierarchy not available; otherwise, can take this out
      childsum = sapply(meta_cell$phenocode[childrenind], function(x) sum(as.numeric(unlist(strsplit(x,"")))) ) #sum children phenocodes up, larger ones have +, smaller ones have - in new marker, split them up
      childplus = which(childsum==max(childsum))
      #split between positive and negative
      # phenoChild[[i]] = list()
      # phenoChild[[i]][[1]] = childrenind[-childplus]
      # phenoChild[[i]][[2]] = childrenind[childplus]
      pc = list()
      pc[[1]] = childrenind[-childplus]
      pc[[2]] = childrenind[childplus]

      ## if there are missing counterparts
      nopm1 = gsub("[+-]","",meta_cell$phenotype[pc[[1]]])
      nopm2 = gsub("[+-]","",meta_cell$phenotype[pc[[2]]])
      # unique1 = which(nopm1%in%setdiff(nopm1,nopm2))
      # unique2 = which(nopm2%in%setdiff(nopm2,nopm1))
      # if (length(unique1)>0) {
      # }
      # if (length(unique2)>0) {
      # }
      negposintersect1 = which(nopm1%in%intersect(nopm1,nopm2))
      negposintersect2 = which(nopm2%in%intersect(nopm2,nopm1))
      if (length(negposintersect1)>0) {
        # phenoChildpn[[i]] = list()
        # phenoChildpn[[i]][[1]] = phenoChild[[i]][[1]][negposintersect1]
        # phenoChildpn[[i]][[2]] = phenoChild[[i]][[2]][negposintersect2]
        pcnp = list()
        pcnp[[1]] = pc[[1]][negposintersect1]
        pcnp[[2]] = pc[[2]][negposintersect2]

        return(list(pc=pc, pcnp=pcnp))
      } else {
        return(list(pc=pc, pcnp=NULL))
      }
    } else {
      return(list(pc=NULL, pcnp=NULL))
    }
  }

  phenoChild = phenoChild_names = phenoChild_ind = phenoChildpn = phenoChildpn_names = phenoChildpn_ind = list()
  for (i in 1:length(ppcc)) {
    if (!is.null(ppcc[[i]]$pc)) {
      phenoChild[[i]] = ppcc[[i]]$pc
      phenoChild_names[[i]] = list()
      phenoChild_names[[i]][[1]] = meta_cell$phenotype[ppcc[[i]]$pc[[1]]]
      phenoChild_names[[i]][[2]] = meta_cell$phenotype[ppcc[[i]]$pc[[2]]]
    }
    if (!is.null(ppcc[[i]]$pcnp)) {
      phenoChildpn[[i]] = ppcc[[i]]$pcnp
      phenoChildpn_names[[i]] = list()
      phenoChildpn_names[[i]][[1]] = meta_cell$phenotype[ppcc[[i]]$pcnp[[1]]]
      phenoChildpn_names[[i]][[2]] = meta_cell$phenotype[ppcc[[i]]$pcnp[[2]]]
    }
  }

  phenoChild_ind <- which(vapply(phenoChild, Negate(is.null), NA))
  phenoChildpn_ind <- which(vapply(phenoChildpn, Negate(is.null), NA))

  phenoChild = meta_cell$phenotype[phenoChild_ind]
  phenoChildpn = meta_cell$phenotype[phenoChildpn_ind]

  phenoChild = phenoChild[phenoChild_ind]
  phenoChild_names = phenoChild_names[phenoChild_ind]
  phenoChildpn = phenoChildpn[phenoChildpn_ind]
  phenoChildpn_names = phenoChildpn_names[phenoChildpn_ind]

  names(phenoChild) = names(phenoChild_names) = meta_cell$phenotype[phenoChild_ind]
  names(phenoChildpn) = names(phenoChildpn_names) = meta_cell$phenotype[phenoChildpn_ind]

  return(list(phenoChild=phenoChild, phenoChild_ind=phenoChild_ind, phenoChild_names=phenoChild_names, phenoChildpn=phenoChildpn, phenoChildpn_ind=phenoChildpn_ind, phenoChildpn_names=phenoChildpn_names))
}

getphenoParent <- function(meta_cell, no_cores=1) {
  require(foreach)
  require(doMC)

  registerDoMC(no_cores)

  phenoParent = foreach (i=c(1:nrow(meta_cell))[meta_cell$phenolevel!=0]) %dopar% {
    pheni = unlist(strsplit(meta_cell$phenocode[i],""))
    phenind = which(pheni!="0")
    zeroind = which(pheni=="0")
    parentphenocode = as.vector(sapply(phenind, function(x) {
      pi0=pheni;
      pi0[x]="0"
      return(paste(pi0,collapse=""))
    } ))
    parentsind = match(parentphenocode, meta_cell$phenocode)
    parentsind = parentsind[!is.na(parentsind)]
    if (length(parentsind)>0) {
      return(parentsind)
    } else {
      return(NULL)
    }
  }
  phenoParent_ind <- which(vapply(phenoParent, Negate(is.null), NA))

  phenoParent = phenoParent[phenoParent_ind]
  phenoParent_names = lapply(phenoParent, function(x) meta_cell$phenotype[x])

  names(phenoParent) = names(phenoParent_names) = meta_cell$phenotype[phenoParent_ind]

  return(list(phenoParent=phenoParent, phenoParent_names=phenoParent_names, phenoParent_ind=phenoParent_ind))
}

# getphenoParent <- function(meta_cell, phenoChildpn, phenoChildpn_ind, no_cores) {
#   require(foreach)
#   require(doMC)
#
#   registerDoMC(no_cores)
#
#   phenoParent = foreach (i=(1:nrow(meta_cell))) %dopar% {
#     pheni = unlist(strsplit(meta_cell$phenocode[i],""))
#     phenind = which(pheni!="0")
#     #enumerate all possible parents
#     parentphenocode = sapply(phenind, function(x) { pi=pheni; pi[x]="0"; return(paste(pi,collapse="")) } )
#     parentind = match(parentphenocode, meta_cell$phenocode)
#     parentind = parentind[!is.na(parentind)]
#     if (length(parentind)>0) return(parentind)
#     return(NULL)
#   }
#
#   phenoParent_ind = which(vapply(phenoParent, Negate(is.null), NA))
#   phenoParent = phenoParent[phenoParent_ind]
#   names(phenoParent) = meta_cell$phenotype[phenoParent_ind]
#
#   #delete parents from child whom doesn't have a twin under that parent.
#   phenoParentpn = phenoParent
#   phenoParentpn_ind = phenoParent_ind
#   delind = c()
#   for (i in 1:length(phenoParent)) {
#     child = phenoParent_ind[i]
#     parents = phenoParent[[i]]
#     delparents = c()
#     for (j in 1:length(phenoParent[[i]])) {
#       parent = phenoParent[[i]][j]
#       ind = (phenoChildpn_ind==parent)
#       if (sum(ind)>0) {
#         if (!sum(unlist(phenoChildpn[[which(ind)]])==child)>0) delparents = append(delparents, j)
#       }
#     }
#     if (length(delparents)>0) parents = parents[-delparents]
#     if (!length(parents)>0) {
#       delind = append(delind, i)
#     } else {
#       phenoParentpn[[i]] = parents
#     }
#   }
#
#   if (length(delind)>0) {
#     phenoParentpn[delind] = NULL
#     phenoParentpn_ind = phenoParentpn_ind[-delind]
#   }
#
#   return(list(phenoParent=phenoParent, phenoParent_ind=phenoParent_ind, phenoParentpn=phenoParentpn, phenoParentpn_ind=phenoParentpn_ind))
#
# }




## Input: feature matrix (row = sample, col = phenotypes), phenocode or phenotype (prefer phenocode)
## Output: feature matrix with proportions of parent
prop_parent <- function(m, phen, list=F, par=T, no_cores=detectCores()-3) {
  if (ncol(m)!=length(phen)) {
    cat("ncol(matrix) doesn't match length(phen)")
    return(NULL)
  }

  #list of child phenotypes for each non-leaf phenotype
  children = phen_children(phen)
  parentind = which(!is.null(children))

  if (par) {
    mplist = foreach(i = 1:length(children), .combine='list') %dopar% { #for each phenotype
      if (is.null(children[[i]])) return(NULL)
      return(m[,children[[i]]]/m[,i])
    }
  } else {
    mplist = list()
    for (i in parentind) {
      mplist[[i]] = m[,children[[i]]]/m[,i]
    }
  }
  if (list) return(list(mp=mplist,parentind=parentind,children=children))
  mp = do.call(cbind, mplist)
  return(list(mp=mp,parentind=parentind,children=children))
}



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
  gc = ldply(names(pchild), function(x)
    data.frame(from=x, to=unlist(ldply[[x]])) )
  gp = ldply(names(pparen), function(x)
    data.frame(from=pparen[[x]], to=x) )
  gr_e = unique(rbind(gc,gp))
  gr_v = meta_cell$phenotype
  gr_e = rbind(ldply(meta_cell$phenotype[meta_cell$phenolevel==1], function(x) data.frame(from="",to=x)), gr_e)

  gr_vp = gr_v[!grepl("[-]",gr_v)]
  gr_ep = gr_e[!grepl("[-]",gr_e[,1]) & !grepl("[-]",gr_e[,2]),]

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
s


## IN PROGRESS
## Input: normalized count matrix m, phenoChild list of children, phenoChild_ind indices of parents in phenoChild, whether this list is negative only (with positive counterpart), weighting function
child_entrop <- function(m, phenoChild, phenoChild_ind, pc_overlapneg=T) {
  #get entropy for each cell population
  require(entropy)
  m[m<1] = 1 #prevent Inf, NaN's
  for (i in 1:(nrow(m)-1)) {
    for (j in (i+1):nrow(m)) { #all combination of rows
      ratio = m[,i]/m[,j] #10/5=2, 5/10=.5;; log(2)=.69, log(.5)=-.69
      growth = m[,i]-m[,j]
      weight = matrix(0,nrow=2,ncol=length(ratio)) #row 1 = negative, row 2 = positive (more than parent)
      for (s in 1:length(phenoChild)) { # all features
        l = log(ratio[phenoChild[[s]]]/ratio[phenoChild_ind[s]])
      }
    }
  }
}



## Input: normalized count matrix m, phenoChild list of children, phenoChild_ind indices of parents in phenoChild, whether this list is negative only (with positive counterpart), weighting function
child_entropeasy <- function(m, phenoChild, phenoChild_ind, pc_overlapneg=T) { # use only overlapping indices, easier to normalize
  #get entropy for each cell population
  require(entropy)
  m[m<1] = 1 #prevent Inf, NaN's
  for (i in 1:(nrow(m)-1)) {
    for (j in (i+1):nrow(m)) { #all combination of rows
      ratio = m[,i]/m[,j] #10/5=2, 5/10=.5;; log(2)=.69, log(.5)=-.69
      growth = m[,i]-m[,j]
      weight = matrix(0,nrow=2,ncol=length(ratio)) #row 1 = negative, row 2 = positive (more than parent)
      for (s in 1:length(phenoChild)) { # all features
        l = log(ratio[phenoChild[[s]]]/ratio[phenoChild_ind[s]])
      }
    }
  }
}


# JS Divergence; 0<JSD<log(2)

# input: normalized count matrix m, phenoChild list of children, phenoChild_ind indices of parents in phenoChild, whether this list is negative only (with positive counterpart), weighting function
child_JSD <- function(m, phenoChildpn, phenoChildpn_ind) { # use only overlapping indices, easier to normalize
  #get entropy for each cell population
  m[m<1] = 1 #prevent Inf, NaN's
  for (i in 1:(length(phenoChildpn_ind)-1)) {
    for (j in (i+1):length(phenoChildpn_ind)) { #all combination of rows
      coli = phenoChildpn_ind[i]
      colj = phenoChildpn_ind[j]
      dist = rep(0,nrow(m)) #jsd for all features
      dist[m[,coli]==0 & m[,colj]==0] = 0
      dist[(m[,coli]!=0 & m[,colj]==0) | m[,coli]==0 & m[,colj]!=0] = log(2)
      non0parents = (m[,coli]!=0 & m[,colj]!=0) #parents that are not 000

      parenti = m[,coli]
      childreni = cbind(m[,unlist(phenoChildpn[[i]][[1]])])
      parentj = m[,colj]
      childrenj = cbind(m[,unlist(phenoChildpn[[j]][[1]])])

      dist[non0parents] = sapply(non0parents, function(x) return(jsd(
        mapply(c,as.list(childreni[x,]),as.list(parenti[x]-childreni[x,]),SIMPLIFY=F),
        mapply(c,as.list(childrenj[x,]),as.list(parentj[x]-childrenj[x,]),SIMPLIFY=F) )))

    }
  }
}







## Input: matrix (row = objects)
## Output: distance matrix










### Pvalue #########################################################################

## Significance; take note to put in same size control
## tstat = (a-mean(control))/se
## se = sd(control)/sqrt(length(control))
t.test.single <- function(control,a) { #input a vector and a single number
  se = sd(control) #/ sqrt(length(control))
  tstat = (a-mean(control)) / se
  p = 2*pt(-abs(tstat), df=length(control)-1) # 2 tailed (less stringent)
  return(p)
}


## find lv such that b= exp(log(a,lv))
getlv = function(a,b=100) return(exp(log(a)/(log(b))))









### Normalize ######################################################################

## Input: x matrix (phenotypes on columns, will convert in function) -- x0 is optional, plots according to x0 so if you want to plot more cell populations than those in x
## Output: f = normalize factors per sample; fidff = difference between f and peak of count ratio per sample
tmm <- function(x,x0=NULL,lib.size,refColumn,cutoff=Inf,plotimg=T,pngnames=NULL,mains=NULL,no_cores=detectCores()-1,samplesOnCol=F,Acutoff=-10) {
  require(foreach)
  require(doMC)
  registerDoMC(no_cores)

  if (is.null(x0)) x0 = x
  if (!samplesOnCol) { x = t(x); x0 = t(x0) }

  ## Taken from TMM

  #f <- rep(NA,ncol(x))
  #fdiff <- rep(NA,ncol(x)) #diff between density peak and value (note: logged)
  ref <- x[,refColumn]
  refn <- lib.size[refColumn]
  doWeighting <- F
  logratioTrim <- .3
  sumTrim <- 0.05
  minlogR <- 1e-6 #min value of log2((obs/obsn)/(ref/refn))


  ff <- foreach(i=1:ncol(x), .combine = list, .maxcombine = ncol(x), .multicombine = T) %dopar% {
    #for(i in ncol(x):1) { cat(i," ",sep="")
    obs <- x[,i]
    obsn <- lib.size[i]
    #logR <- log2((obs/obsn)/(ref/refn))
    logR <- log2(obs/ref)			#log ratio of expression, accounting for libr size
    absE <- (log2(obs/obsn) + log2(ref/refn))/2	#absolute expression
    v <- (obsn-obs)/obsn/obs + (refn-ref)/refn/ref	 #estimated asymptotic variance

    #remove infinite values, cutoff based on A
    fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
    logR <- logR[fin]
    absE <- absE[fin]
    v <- v[fin]

    if(max(abs(logR)) < minlogR) { return(list(f=1, fdiff=0)) # f[i] <- 1
    } else {

      #taken from the original mean() function
      n <- length(logR)
      loL <- floor(n * logratioTrim) + 1
      hiL <- n + 1 - loL
      loS <- floor(n * sumTrim) + 1
      hiS <- n + 1 - loS

      #keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
      #a fix from leonardo ivan almonacid cardenas, since rank() can return
      #non-integer values when there are a lot of ties
      keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)

      if(doWeighting) {
        fi <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
      } else { fi <- mean(logR[keep], na.rm=TRUE) } #f[i] <- mean(logR[keep], na.rm=TRUE) }

      #Results will be missing if the two libraries share no features with positive counts
      #In this case, return unity
      #if(is.na(f[i])) f[i] <- 0
      if(is.na(fi)) fi <- 0

      #check if close to peak; if not, switch to peak
      d <- density(log2((obs)/ref), na.rm=T)
      p <- as.matrix(findpeaks(d$y)); if(ncol(p)==1) p <- t(p)
      p1 <- d$x[p[which.max(p[,1]),2]]
      #fdiff[i] <- p1-f[i]
      fdiffi <- p1-fi

      if (plotimg) {
        pngname <- pngnames[i]
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
        #f[i] <- p1
        fi = p1
      }

      #f[i] <- 1/2^f[i]
      fi <- 1/2^fi

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

  f <- rep(NA,length(ff))
  fdiff <- rep(NA,length(ff)) #diff between density peak and value (note: logged)
  for (i in 1:length(ff)) {
    f[i] <- ff[[i]]$f
    try({ fdiff[i] <- ff[[i]]$fdiff })
  }

  return(list(f=f,fdiff=fdiff))
}








## Classify/Cluster ############################################

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
rowxcluster_to_cluster <- function(rowxcluster) {
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
clusterxcol_to_cluster <- function(clusterxcol) {
  clusters_size = apply(clusterxcol, 1, function(x) sum(x) )
  rclusters = apply(clusterxcol, 2, function(x) {
    cl = which(x)
    if (length(cl) == 0) cl = 0
    if (length(cl) > 1) cl = cl[which.max(clusters_size[cl])]
    return(cl)
  })
}


