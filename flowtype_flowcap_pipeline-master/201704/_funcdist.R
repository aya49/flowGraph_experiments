



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












## Significance; take note to put in same size control
## tstat = (a-mean(control))/se
## se = sd(control)/sqrt(length(control))
t.test.single <- function(control,a) { #input a vector and a single number
  se = sd(control) #/ sqrt(length(control))
  tstat = (a-mean(control)) / se
  p = 2*pt(-abs(tstat), df=length(control)-1) # 2 tailed (less stringent)
  return(p)
}















## Distance
## Input: m (row=object, col=attributes), (optional) m0=list(m1,m2,m3,...), lcombw=linear combo weight, w (weight vector for each object), wlog=T (log weights first), wpercentage=T (weights add up to 1)
dist_a <- function(m,m0=NULL,lcombow=NULL,w=NULL,wlog=T,wpercentage=T,method="canberra") {
  d = matrix(Inf,nrow=nrow(m),ncol=nrow(m))
  weight = rep(1,nrow(m))
  if (!is.null(m0) & is.null(lcombow)) lcombow = rep(1,length(m0))
  if (is.null(m0)) 
  for (i in 1:(nrow(m)-1)) {
    for (j in (i+1):nrow(m)) {
      if (!is.null(w)) {
        weight = pmax(w[i,],w[j,])
        if (wlog) weight = log(weight)
        if (wpercentage) weight = weight/sum(weight)
      }
      
      if (grep("canberra",method)) {
      }
      d[i,j] = d[j,i] = weight
    }
  }
  return(as.dist(d))
}

