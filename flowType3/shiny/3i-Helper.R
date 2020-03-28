find.farthest.pheno <- function(pheno,max.limit=300, min.limit=5)
{
  pheno.mat <-sapply(pheno,function(x) as.numeric(unlist(strsplit(x, ""))))
  rownames(pheno.mat)<-colnames(pheno.mat) <- NULL
  distance <- fields::rdist(t(pheno.mat))
  dmax <-apply(distance,2,which.max)
  dmax <- unique(dmax)
  mat <- pheno.mat[,dmax]
  prev.ind <- length(dmax)
  while(ncol(mat)>min.limit)
  {
    distance <- rdist(t(mat))
    dmax <-unique(apply(distance,2,which.max))
    ind <-sort(dmax,decreasing = T,index.return=T)$ix[1:min(max.limit,length(dmax))]
    ind2<-unique(apply(distance[ind,],1, which.max))
    mat <-mat[,ind2]
    if (length(ind2)==1)
      return(paste(mat, sep="",collapse=""))
    if ( ncol(mat)==prev.ind){
      break
    }else{
      prev.ind <- ncol(mat)
    }
  }
  return(apply(mat,2,paste,sep="",collapse=""))
}
##########################
##########################
########################## 

##########################
##########################
########################## 
flowtype.settingup<-function(panel,threshold,which.file=1)
{
  th1<-threshold[[which.file]]
 s.names <-unlist(lapply(threshold,function(x) return(x[[1]])))
 class.type <- unlist(lapply(threshold[[which.file]],class))
  if (panel=="Panel_BM-cell")
  {
    for ( i in 1:length(threshold))
    {
      print(i)
      if (!all(unlist(lapply(threshold[[i]],class))==class.type))
      {
        threshold[[i]]<-1
      }else{
      }
      threshold[[i]] <- threshold[[i]][c("gr1.gate","CD3Tcells.filter","plasma.filter","b220.gate","cd43.gate","cd24.gate","bp1.gate","igm.gate","igd.gate")]
      colnames(threshold[[i]][[2]]) <- c( "PE-Cy7-A",  "BV786-A")
      colnames(threshold[[i]][[3]]) <- c( "BV650-A","PE-Cy7-A")
    }
    names(threshold)<-s.names
    PropMarkers <- c(10,18,8,9,17,12,7)
    PartitionsPerMarker <- 2
    MaxMarkersPerPop<-7
    partition.vec <- rep(2,9)
    filter.vec <- c(F,T,T,F,F,F,F,F,F)
    marker.names<-c("GR1","Tcells","Plasma","B220","CD43","CD24","BP1","IgM","IgD")
    return(list(prp=PropMarkers,part=PartitionsPerMarker,max=MaxMarkersPerPop,
                part.v=partition.vec,filt.v=filter.vec,names=marker.names,threshold=threshold))
  }else if (panel=="Panel_M-cell")
  {
    for ( i in 1:length(threshold))
    {
      print(i)
      if (!all(unlist(lapply(threshold[[i]],class))==class.type) | any(is.na(names(threshold[[i]]))) | length(threshold[[i]])<11)
      {
        threshold[[i]]<-1
      }else{
        threshold[[i]] <- threshold[[i]][c("lin.neg.filter","lin.mac.filter","eosinophils.filter", "Ly6g.neg.filter",
                                           "Ly6c.gate","cd11b.gate","pdc.flow.filter", "cd11c.gate", "mhcII.gate","cd8a.filter","cd103.gate" )]
        colnames(threshold[[i]][[1]]) <- c( "BV605-A","BV421-A")
        colnames(threshold[[i]][[2]]) <- c( "PerCP-Cy5-5-A","BV510-A")
        colnames(threshold[[i]][[3]]) <- c( "APC-A","SSC-A")
        colnames(threshold[[i]][[4]]) <- c( "APC-A","SSC-A")
        colnames(threshold[[i]][[7]]) <- c( "Alexa Fluor 700-A","BV510-A")
        colnames(threshold[[i]][[10]]) <- c( "BV510-A","PE-Cy7-A")
      }
    }
    PropMarkers <- c(10,13,16,7,17)
    PartitionsPerMarker <- 2
    names(threshold)<-s.names
    
    l<- length(th1)
    MaxMarkersPerPop<-11
    partition.vec <- rep(2,11)
    filter.vec <- c(T,T,T,T,F,F,T,F,F,T, F)
    marker.names<-c("Lin-","Lin-Mac-","Eosinophils","Ly6G-","Ly6C","CD11b","PDC","CD11c","MHCII","CD8A","CD103") 
    return(list(prp=PropMarkers,part=PartitionsPerMarker,max=MaxMarkersPerPop,
                part.v=partition.vec,filt.v=filter.vec,names=marker.names,threshold=threshold))
  }else if (panel=="Panel_B-cell")
  {
    for ( i in 1:length(threshold))
    {
      
      if (!all(unlist(lapply(threshold[[i]],class))==class.type))
      {
        threshold[[i]]<-1
      }else{
        threshold[[i]] <- threshold[[i]][ c("b220.gate","plasma.cells.filter","b1.bcells.filter", "cd95.gate","gl7.gate","igg1.gate","mzp.filter","transitional.filter",
                                            "cd.23.gate","igm.gate","igd.gate" ,"iggp.mem.bcells.filter")]
        
        colnames(threshold[[i]][[2]]) <- c( "Alexa Fluor 700-A",  "BV650-A")
        colnames(threshold[[i]][[3]]) <- c( "PE-A",  "BV510-A")
        colnames(threshold[[i]][[7]]) <- c( "BV786-A","FITC-A")
        colnames(threshold[[i]][[8]]) <- c( "BV786-A","FITC-A")
        colnames(threshold[[i]][[12]]) <- c( "BV786-A","PE-A")
      }
    }
    names(threshold)<-s.names
    
    PropMarkers <- c(10,18,9,17,12,16,8)
    PartitionsPerMarker <- 2
    l<- length(th1)
    MaxMarkersPerPop<-6
    partition.vec <- rep(2,12)
    filter.vec <- c(F,T,T,F,F,F,T,T,F,F,F,T)
    marker.names<-c("b220","Plasma","B1Bcells", "CD95","GL7","IgG1","MZP","Transitional cells",
                    "CD23","IgM","IgD" ,"IgG MemoryBcells")
    return(list(prp=PropMarkers,part=PartitionsPerMarker,max=MaxMarkersPerPop,names=marker.names,
                part.v=partition.vec,filt.v=filter.vec,threshold=threshold))
  }
}

##########################
##########################
########################## 

#Naive way of generating WT rank based on what I understood from the call
rank.generator <- function(ranks)
{
  dup<-duplicated(ranks)
  wt.ranks<-setdiff(1:140, unique(ranks))
  for ( i in 1:length(ranks))
  {
    if(dup[i])
      wt.ranks <-setdiff(wt.ranks,wt.ranks[wt.ranks-ranks[i]>0][1])
    
  }
  return(wt.ranks)
}


##########################
##########################
########################## 

pval.calc <- function(props, label,geneders=NULL,varaince, WT=F,which.sex=NULL,wilcox=F, rank.wilcox=F,method="BH",rank.given=F)
{
  
  if(rank.given)
  {
    
    
    
  }
  if (WT)
  {
    print("Comparing wildtype to KO")
    label.1  <- which(label==which.sex)
    label.2 <- which(label=="WT") 
  }else{
    label.1  <- which(label=="Female")
    label.2 <- which(label=="Male")
  }
  pvals <- apply(props,2, function(prop.x) {
    if (all(prop.x==1))
      return(1)
    if (wilcox)
    { 
      p<-tryCatch(wilcox.test(prop.x[label.1],prop.x[label.2],var.equal = varaince)$p.value, error=function(x) {return("wrong")})
      
    }else{
      p<-tryCatch(t.test(prop.x[label.1],prop.x[label.2],var.equal = varaince)$p.value, error=function(x) {return("wrong")})
    }
    if (mode(p)=="character" | is.nan(p))
    {
      #print("t.test failed, as data was almost constant.")
      p<-1
    }
    if (rank.wilcox)
    { 
      m_rank<- rank(c(prop.x[which(genders$ko=="Male")],prop.x[which(genders$wt=="Male")]))
      f_rank <- rank(c(prop.x[which(genders$ko=="Female")],prop.x[which(genders$wt=="Female")]))
      
      #And this is significant: 
      p2<-tryCatch(wilcox.test(x = c( m_rank[1:length(which(genders$ko=="Male"))], f_rank[1:length(which(genders$ko=="Female"))] ),
                               y= c( m_rank[(length(which(genders$ko=="Male"))+1): length(m_rank)], f_rank[(length(which(genders$ko=="Female"))+1): length(f_rank)] ),
                               var.equal = varaince)$p.value, error=function(x) {return("wrong")})
    }
    
    if (mode(p2)=="character"| is.nan(p2))
    {
      #print("t.test failed, as data was almost constant.")
      p2<-1
    }
    return(c(p,p2))
  })
  p.w<- unlist(lapply(pvals, function(x) x[1]))
  p.r<- unlist(lapply(pvals, function(x) x[2]))
  mean.1 <- apply(props[label.1,],2,median)
  mean.2 <- apply(props[label.2,],2,median)
  names(pvals)<- colnames(props)
  pvalsw.adjust <- p.adjust(p.w, method=method)
  pvalsr.adjust <- p.adjust(p.r, method=method)
  return(cbind(p.w,pvalsw.adjust, mean.1,mean.2,p.r,pvalsr.adjust))
}

##########################
##########################
##########################
Find.WTsamples<-function(g1, allProportions,ko.ind,wt.ind, sex)
{
  ind1 <- which(allProportions[ko.ind,"Gender"]==sex)
  group <- c("Female","Male")
  if (length(ind1)==0)
  {
    warning(paste("There were not any", sex, "KO for this gene. Taking the range from the opposite sex instead.", sep=" "))
    ind1 <- which(allProportions[ko.ind,"Gender"]==group[which(group!=sex)])
  }  
  ind2 <- which(allProportions[wt.ind,"Gender"]==sex)
  assay.date<-as.Date(allProportions[ko.ind[ind1],"Assay Date"],format = "%d-%b-%Y")
  wt.date <- as.Date(allProportions[wt.ind[ind2],"Assay Date"],format = "%d-%b-%Y")
  t1<-assay.date[order(assay.date)][1]
  t2 <- tail(assay.date[order(assay.date)],1)
  # on.range1 <-which( wt.date ==t1 )
  within.range <-  which( wt.date >=t1 & wt.date  <=t2)
  #on.range2 <-  which(wt.date  ==t2)
  out.range.left <- which(wt.date  <t1)
  out.range.right <- which(wt.date  >t2)
  
  if (length(within.range)>=70)
  {
    set.seed(240)
    wt.range <- sample(size =70,within.range)
  }else{
    needed <- 70-length(within.range)
    # if (needed%%2==0)
    # {
    if (length(out.range.left)>=needed/2)
    {
      wt.l.range <-tail(out.range.left,floor(needed/2))
      if (length(out.range.right)<needed/2)
      {
        wt.r.range <-out.range.right
        wt.l.range <- tail(out.range.left,needed-length(out.range.right))
      }else{
        wt.r.range <-tail(out.range.right,ceiling(needed/2))
      }
    }else
    {
      wt.l.range <-out.range.left
      if (length(out.range.right)<(needed/2-length(out.range.left)))
      {
        print("Not enough WT sample")
        next;
      }
      wt.r.range <- tail(out.range.right,needed-length(out.range.left))
    }
    wt.range <- c( wt.l.range,within.range,wt.r.range)
  }
  
  return(wt.range)
}


