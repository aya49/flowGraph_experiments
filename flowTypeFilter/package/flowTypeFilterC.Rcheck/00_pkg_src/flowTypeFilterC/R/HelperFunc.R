#Author: Mehrnoush Malek  
#Date Modified: Oct 19,2017
#R source code for calculating counts for all permutation, when thresholds are filters


find.channels <- function(Thresholds,Frame)
{
  inds <-which(names(Thresholds)=="gate")
  chs <- unlist(lapply (inds, function(x) {
    ch <-colnames(Thresholds[[x]])
    if (is.null(ch))
      stop("Filters should have colnames that matches FCS channel names.")
    ch <-match(ch,colnames(exprs(Frame)))
    if (length(ch)!=2)
      stop("Cannot match names of filters with channel names")
    return(ch)}))
  return(chs)
}
###
###Things to do:
#You need to think of PropMarkers for cases when filter is involved, The current code is
#based on PropMArkers only include the channel that have thresholds
#You need to update flowType.R on Tuesday

####
###Generating all permutations

PopsGenerator<- function(l,PartitionsPerMarker,MaxMarkersPerPop )
{
  if (length(PartitionsPerMarker)==1)
  {
    vect <-lapply(1:l, function(x)  rep(as.character(c(0:PartitionsPerMarker)),1))
  }else if (length(PartitionsPerMarker)==l){
    vect <- lapply(1:length(PartitionsPerMarker), function( x) as.character(c(0:PartitionsPerMarker[x])))
  }else{
    stop("PartitionsPerMarker is either an integer or a vector of same length as MarkerNames")
  }
  pops <- expand.grid(vect)
  if (MaxMarkersPerPop<l)
  {
    pops <-pops[apply(pops,1, function(x) return(length(which(x!=0))<=MaxMarkersPerPop))==TRUE,]
  }
  return(pops)
}


######
is.in.range <- function(vect,threshold, segment)
{
  segment<-as.numeric(segment)
  if (segment == 1)
  {
    
    return (vect <= threshold[1])
  }
  else
  {
    return(vect> threshold[(segment-1)])
  }
  
}


######
is.in.filter <- function(mat,threshold, segment)
{
  segment<-as.numeric(segment)
  pts <- sp:: point.in.polygon(point.x = mat[,1],point.y = mat[,2],pol.x = threshold[,1],pol.y = threshold[,2] )
  if(segment==1)
    return(pts==0)
  if (segment==2)
    return(pts!=0)
  if (segment > 2)
    stop("Cannot split data to more than 2, when polygon filter is provided.")
}

#####
#terminalphenotype<-function(phenotype, data,Thresholds)
  
#{
  
  
 # ind <- rep(T, nrow(data))
  #phenotype2<-unlist(strsplit(phenotype,split=""))
  ###Indices problem
  #p1<-which(phenotype2!="0")
  #if (length(p1)<1)
   # stop("Something is wrong with coding partitions.")
    #  if (names(Thresholds)[p1]=="gate")
     # {
      #  ind <- is.in.filter(mat=data[,colnames(Thresholds[[p1]])],threshold = Thresholds[[p1]],segment =  phenotype2[p1])
      #}else{
       # ind <- is.in.range(vect= data[,names(Thresholds)[p1]],threshold = as.vector(Thresholds[[p1]]), segment = phenotype2[p1])
      #}
  #return(ind)
#}

terminalphenotype<-function(phenotype, data,Thresholds)
  
{
  
  
  ind <- rep(T, nrow(data))
  phenotype2<-unlist(strsplit(phenotype,split=""))
  ###Indices problem
  p1<-which(phenotype2!="0")
  if (length(p1)<1)
    stop("Something is wrong with coding partitions.")
  if (names(Thresholds)[p1]=="gate")
  {
    ind <- is.in.filter(mat=data[,colnames(Thresholds[[p1]])],threshold = Thresholds[[p1]],segment =  phenotype2[p1])
  }else{
    ind <- is.in.range(vect= data[,names(Thresholds)[p1]],threshold = as.vector(Thresholds[[p1]]), segment = phenotype2[p1])
  }
  return(which(ind==T))
}


count.cells <-function(Thresholds,PartitionsPerMarker,MaxMarkersPerPop,EXPRS,cores,verbose,...)
{
  
  l<- length(Thresholds)
  pops <- as.matrix(PopsGenerator(l,PartitionsPerMarker,MaxMarkersPerPop))
  #First level
  start.time <- Sys.time()
  d<-apply(pops,1,function(x) return (length(which(as.numeric(x)==0))==(l-1)))
  first.level.pheno <- apply(pops[d==T,],1, paste0,sep="",collapse = '')
  first.level.ind<- vector("list",length = length(first.level.pheno))
  first.level.ind <- parallel::mclapply(first.level.pheno,function(ph) 
                              return(terminalphenotype (phenotype = ph,  data =  EXPRS,Thresholds = Thresholds)),mc.cores = cores,...)
  counts <-nrow(EXPRS)
  codes<-paste0(pops[1,],collapse="")
  names(first.level.ind)<- first.level.pheno
  if(verbose)
{
    print(paste0("Level 1 done. Identified phenotypes: ", length(first.level.pheno)))
  each.time<-Sys.time()

    print(paste0("Running time for level",i1,": " ,each.time-start.time))
}
 counts<-c(counts,unlist(parallel::mclapply(first.level.ind,length)))

  #counts<-c(counts,unlist(parallel::mclapply(first.level.ind,sum,mc.cores = cores,...)))
  codes<-c(codes,first.level.pheno)
  prev.level.ind<-first.level.ind
  if (MaxMarkersPerPop==1)
    return(counts)
for(i1 in 2:MaxMarkersPerPop)
#  lapply(2:MaxMarkersPerPop, function(i1)
  {

    d<-apply(pops,1,function(x) return (length(which(as.numeric(x)==0))==(l-i1)))
    curr.level.pheno <-apply(pops[which(d==T),],1, paste0,sep="",collapse = '')
    curr.level.ind<- vector("list",length = length(curr.level.pheno))
    curr.level.ind <- parallel::mclapply(curr.level.pheno,function(ph) {
      ph2<-unlist(strsplit(ph,split=""))
      ind <-which(ph2!="0")
      temp <- rep("0",l)
      temp[-ind[1]]<-ph2[-ind[1]]
      #temp[ind[-length(ind)]]<-ph2[ind[-length(ind)]]
      prev.pheno <- paste0(temp, sep="",collapse = "")
      temp <- rep("0",l)
      #temp[ind[length(ind)]]<-ph2[ind[length(ind)]]
      temp[ind[1]]<-ph2[ind[1]]
      
      first.pheno <- paste0(temp, sep="",collapse = "")
      #colnames(prev.level.ind)==prev.pheno changed to just prev.pheno
 inds <- intersect( prev.level.ind[[prev.pheno]] ,
      #which(first.level.pheno[,tail(ind,1)]==ph[tail(ind,1)]) changed to first.pheno
      first.level.ind[[first.pheno]])
     # inds <- prev.level.ind[[prev.pheno]] &
       # #which(first.level.pheno[,tail(ind,1)]==ph[tail(ind,1)]) changed to first.pheno
      #  first.level.ind[[first.pheno]]
      return(inds)
    },mc.cores=cores)
    if(verbose){
      print(paste0("Level ", i1, " done. Identified phenotypes: ", length(curr.level.pheno)))
      each.time<-Sys.time()

    print(paste0("Running time for level",i1,": " ,each.time-start.time))
    }
system.time({
    names(curr.level.ind)<- curr.level.pheno
    codes<-c(codes,curr.level.pheno)
 counts<-c(counts,unlist(parallel::mclapply(curr.level.ind,length)))
   # counts<-c(counts,unlist(parallel::mclapply(curr.level.ind,sum)))
    prev.level.ind <- curr.level.ind
})
gc(verbose = getOption("verbose"), reset = T)
  }
  end.time<-Sys.time()
  if (verbose)
    print(paste0("Running time:" ,end.time-start.time))
  return(list(counts=counts,codes=codes,mfi=NULL))
}
