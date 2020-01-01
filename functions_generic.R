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


## time ----------------------------

#' @title Formats time into string.
#' @description Formats time into a string HH:MM:SS given time zone.
#' @param time A time variable of class \code{POSIXct}, \code{POSIXt}.
#' @param tz A string indicating time zone.
#' @return Time formatted as a string; used in \code{time_output} function.
#' @examples
#' \dontrun{
#'  tstr(Sys.time())
#' }
#' @rdname tstr
tstr <- function(time, tz="GMT")
    format(.POSIXct(time, tz=tz), "%H:%M:%S")


#' @title Outputs elapsed time.
#' @description Given a time, prints the time elapsed from that time until now.
#' @param start A time variable of class \code{POSIXct}, \code{POSIXt}.
#' @param msg A string with a message to print out after the elapsed time.
#' @param tz A string indicating time zone.
#' @return Prints to console, the time from which process
#'  started \code{start} - ended, and > time elapsed from \code{start} until now.
#' @examples
#' \dontrun{
#'  start <- Sys.time()
#'  time_output(start,'start - now > time elapsed')
#' }
#' @rdname time_output
#' @export
time_output <- function(start, msg="", tz="GMT") {
    start <- as.POSIXct(start)
    end <- Sys.time()
    time_elapsed <- difftime(end, start, units="secs")
    message(msg, ifelse(msg == "", "", ": "),
        tstr(start, tz=tz), "-",
        tstr(end, tz=tz), " > ",
        tstr(time_elapsed, tz=tz))
}


## parallel -------------------------------------------

#' @title Prepares number of cores needed for parallel backend.
#' @description Prepares number of cores needed for parallel backend.
#' @param no_cores An integer indicating how many cores to parallelize on.
#' @return An integer indicating how many cores to parallelize on.
#' @details Given the number of cores a user wishes to parallelize processes on
#'  \code{no_cores}, \code{parallel_backend} ensures this value does not
#'  exceed the actual number of cores the user's computer contains.
#' @examples
#' \dontrun{
#'  no_cores <- 100
#'  parallel_backend(no_cores)
#' }
#' @seealso
#'  \code{\link[parallel]{detectCores}}
#' @rdname parallel_backend
#' @importFrom parallel detectCores
parallel_backend <- function(no_cores=1)
    max(1, min(no_cores, parallel::detectCores()))


#' @title Prepares parallel loop indices.
#' @description \code{loop_ind_f} is a helper function that splits
#'  a vector of loop indices into a list of multiple loop indices
#'  for use in parallel processes within the flowGraph package.
#' @param x A vector of loop indices.
#' @param n An integer, or the number of vectors to split \code{x} into.
#' @return list of \code{n} vectors with elements from \code{x}.
#' @examples
#' \dontrun{
#'  old_loop_inds <- 1:10
#'  no_cores <- 5
#'  new_loop_inds <- loop_ind_f(old_loop_inds, no_cores)
#'  example_indices <- plyr::llply(new_loop_inds, function(ii) {
#'      plyr::llply(ii, function(i) {
#'          i
#'      })
#'  }, .parallel=T)
#' }
#' @rdname loop_ind_f
#' @export
loop_ind_f <- function(x, n) {
    if (n == 1)
        return(base::list(x))
    return(base::split(x, ceiling(seq_along(x)/ceiling(base::length(x)/n))))
}

#Output: proportions that split a triangle into n equal pieces
triSplit <- function(n) 1-c(sapply((n-1):1, function(x) sqrt((x)/(n))))


## matrix analysis -------------------------------------------------

#' @title Summarizes a numeric matrix.
#' @description Summarizes a numeric matrix.
#' @param m A numeric matrix.
#' @param feat_type Name of the matrix \code{m}.
#' @return A data frame containing one row summarizing \code{m}; see \code{\link[flowGraph]{fg_get_feat_desc}}.
#' @examples
#' \dontrun{
#'  summary_table(matrix(rnorm(12),nrow=3), feat_type='random')
#' }
#' @rdname summary_table
#' @export
summary_table <- function(m, feat_type="") {
    m <- as.matrix(m)
    base::data.frame(feat=feat_type, nrow=base::nrow(m),
        ncol=base::ncol(m), inf=sum(is.infinite(m)),
        neginf=sum(m == -Inf), na=sum(base::is.na(m)),
        nan=sum(is.nan(m)), neg=sum(m < 0),
        pos=sum(m > 0), zero=sum(m == 0), max=max(m[is.finite(m)]),
        min=min(m[is.finite(m)]))
}

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
                                   
#Output: columns where all values are lower than or equal to a threshold
colIndBelow <- function(m,thres) 
  which(apply(m,2,function(x) all(x<=thres)))







## matrix modifiers --------------------------------------------

              
#Output: Random matrix
randomMatrix <- function(nrow,ncol) 
  matrix(rexp(ncol*nrow, rate=.1), ncol=ncol)
              
              
## input: matrix
## output: matrix without all NA col/row
delna = function(m)
  m[apply(m,1,function(x) !all(is.na(x))), apply(m,2,function(x) !all(is.na(x)))]


# deletes columns where all values are the same
delcols = function(x) x[sapply(x, function(y) length(unique(y))>1)]

                               
#' @title Normalizes matrix values by class.
#' @description Used only in the \code{\link[flowGraph]{flowGraph_mean_class}} function;
#'  for each class in the \code{classes} vector, \code{meandiff} takes the column mean
#'  of the rows in the given matrix associated with that class;
#'  it then takes the difference point by point between these means and
#'  the original rows for that class.
#' @param m0 A numeric matrix.
#' @param classes A vector whose length is equal to the number of rows in the given matrix.
#' @return A numeric matrix whose dimensions equate to that of the input
#'  and whose values are normalized per class.
#' @examples
#' \dontrun{
#'  classes <- append(rep('apples',4), rep('oranges',3))
#'  m0 <- matrix(rnorm(35), nrow=7)
#'  m <- meandiff(m0, classes)
#' }
#' @seealso
#'  \code{\link[flowGraph]{flowGraph_mean_class}}
#' @rdname meandiff
#' @export
#' @importFrom plyr llply
mean_diff <- function(m0, classes) {
    m <- m0
    for (pi in base::unique(classes)) {
        pii <- classes == pi
        m_ <- m[pii, , drop=FALSE]
        m_m <- base::colMeans(m_)
        m[pii, ] <- base::do.call(rbind,
            plyr::llply(base::seq_len(base::nrow(m_)), function(i) m_[i, ] - m_m))
    }
    base::dimnames(m) <- base::dimnames(m0)
    return(m)
}


## vector functions ------------------------------------
                        
# combines two vectors in every combination
pastevector <- function(x,y,collapse="") apply(expand.grid(x, y), 1, paste, collapse=collapse)

                        
## numeric functions -----------------------------------------------

# https://stackoverflow.com/questions/3443687/formatting-decimal-places-in-r
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

                        
## arithmetic functions ------------------------------------------------

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


## distance scoring functions ------------------------------------

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


## pvalue -------------------------------------

## Significance; take note to put in same size control
## tstat = (a-mean(control))/se
## se = sd(control)/sqrt(length(control))
t.test.single = function(control,a) { #input a vector and a single number
  se = sd(control) #/ sqrt(length(control))
  tstat = (a-mean(control)) / se
  p = 2*pt(-abs(tstat), df=length(control)-1) # 2 tailed (less stringent)
  return(p)
}


## string functions ------------------------------
                            
charOccurences <- function(char, str) {
  str2 <- gsub(char,"",str)
  return (nchar(str) - nchar(str2))
}
                        

## cluster output functions -------------------------

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

