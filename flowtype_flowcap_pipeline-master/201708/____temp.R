libr(Matrix)
libr(foreach)
libr(doMC)
source("~/projects/IMPC/code/_funcAlice.R")

options(stringsAsFactors=FALSE)

no_cores = detectCores()-1
registerDoMC(no_cores)

layerlimit = 2 #any node 3 layers apart or more gets score of 0


m=get(load("/home/maruko/projects/IMPC/result/P1/Sanger_SPLEEN/matrixLogFoldTRIM_CountAdj.Rdata"))
load("/home/maruko/projects/IMPC/result/P1/Sanger_SPLEEN/sampleMeta.Rdata")
load("/home/maruko/projects/IMPC/result/P1/Sanger_SPLEEN/phenoMeta.Rdata")
markers = unique(gsub("[+-]","",phenoMeta$phenotype[phenoMeta$phenolevel==1])) #used to normalize edge matches
sampleMeta=sampleMeta[match(rownames(m),sampleMeta$fileName),]

#sampleMeta$gene[!is.na(match(sampleMeta$fileName,rownames(m)))] #view all genes
aind = which(sampleMeta$gene[!is.na(match(sampleMeta$fileName,rownames(m)) )]%in%c("1700007K13Rik(b)_1700007K13Rik(b)","1700067K01Rik_+"))
# aind = which(sampleMeta$gene[!is.na(match(sampleMeta$fileName,rownames(m)) )]=="Tmem127_Tmem127")
# aind = which(sampleMeta$gene[!is.na(match(sampleMeta$fileName,rownames(m)) )]=="Cog6_Cog6")

maind0 = m[aind,]
colswsig = which(apply(as.matrix(maind0)[-1,], 2, function(x) any(x>0))==T)
# colswsig = which(apply(1:ncol(maind0), function(x) any(as.matrix(maind0)[,x]>0))==T)
maind = maind0[,colswsig]
maindlist = list()
for (i in 1:nrow(maind)) {
  maindlist[[i]] = colnames(maind)[which(maind[i,]>0)] #-log(.025))]
}
maindl = unlist(maindlist)
# maind[,names(sort(table(maindl),decreasing=T))]
# sampleMeta[!is.na(match(sampleMeta$fileName,rownames(maind))),]



## isolate 2 samples to compare pairwise
compareind = c(5,6)
maind2 = maind[compareind,]
colsanysig = which(sapply(1:ncol(maind2), function(x) any(as.matrix(maind2)[,x]>0))==T)
maind2 = maind2[,colsanysig]
pm2 = phenoMeta[!is.na(match(phenoMeta$phenotype,colnames(maind2))),] #same order as columns

#seperate out trimmed values for the 2 samples
colsallsig = which(sapply(1:ncol(maind2), function(x) all(as.matrix(maind2)[,x]>0))==T)
if (length(colsallsig)>0) maind20 = maind2[,colsallsig]
maind21 = maind2[1,maind2[1,]>0]
maind22 = maind2[2,maind2[2,]>0]

#list of marker vectors of each node
if (length(colsallsig)>0) pt0 = lapply(colnames(maind20), function(x) return(unlist(strsplit(x, "(?<=[+-])", perl=TRUE))))
pt1 = lapply(names(maind21), function(x) return(unlist(strsplit(x, "(?<=[+-])", perl=TRUE))))
pt2 = lapply(names(maind22), function(x) return(unlist(strsplit(x, "(?<=[+-])", perl=TRUE))))

#phenoMeta
if (length(colsallsig)>0) pm20 = phenoMeta[which(!is.na(match(pm2$phenotype,colnames(maind20)))),]
pm21 = pm2[which(!is.na(match(pm2$phenotype,names(maind21)))),]
pm22 = pm2[which(!is.na(match(pm2$phenotype,names(maind22)))),]

#vector of layers
# if (length(colsallsig)>0) pm20layer = sapply(pt0,function(x) return(length(x)))
# pm21layer = sapply(pt1,function(x) return(length(x)))
# pm22layer = sapply(pt2,function(x) return(length(x)))

#edge list
from = c()
to = c()
addedm = c()
ei = 1
for (i in sort(unique(pm21$phenolevel))) {
  if (sum(grepl(i+1,pm21$phenolevel))>0){
    for (j in which(pm21$phenolevel==i)) {
      for (k in which(pm21$phenolevel==i+1)) {
        if (all(!is.na(match(pt1[[j]],pt1[[k]])))) { to[ei]=pm21$phenotype[k]; from[ei]=pm21$phenotype[j]; addedm[ei]=setdiff(pt1[[k]],pt1[[j]]); ei=ei+1}
      }
    }
  }
}
evalue = rep(0,length(addedm))
e1 = data.frame(from,to,addedm,evalue)

from = c()
to = c()
addedm = c()
ei = 1
for (i in sort(unique(pm22$phenolevel))) {
  if (sum(grepl(i+1,pm22$phenolevel))>0){
    for (j in which(pm22$phenolevel==i)) {
      for (k in which(pm22$phenolevel==i+1)) {
        if (all(!is.na(match(pt2[[j]],pt2[[k]])))) { to[ei]=pm22$phenotype[k]; from[ei]=pm22$phenotype[j]; addedm[ei]=setdiff(pt2[[k]],pt2[[j]]); ei=ei+1}
      }
    }
  }
}
evalue = rep(0,length(addedm)) #difference of 1 is max difference (0 value)
e2 = data.frame(from,to,addedm,evalue)

if (length(colsallsig)>0) e0 = e1[intersect( which(e1$from%in%e2$from), which(e1$to%in%e2$to) ),]
















#matrix of node matches based on their markers (1=match, 0 = no markers in common):: match/all; or #indelswaps (not proportional to original)
#nodematch = Matrix( foreach(y=pt1, .combine='rbind') %dopar% {return( sapply(pt2, function(x) return(length(intersect(x,y)) / max(length(y),length(x))) ) )}, sparse=T)
nodematch0 = Matrix( foreach(y=pt1, .combine='rbind') %dopar% {return( sapply(pt2, function(x) {
  if (length(intersect(x,y))>0) { return(max(length(y),length(x))-length(intersect(x,y)))
  } else { return(Inf) }
}) )}, sparse=T)
rownames(nodematch0) = names(maind21)
colnames(nodematch0) = names(maind22)

#+1 with log
nodematch1 = log(nodematch0+1) * (1/log(layerlimit+2)) #layerlimit+1 = 1 (score of 0)
nodematch1[nodematch1>1] = Inf
nodematch2 = matrix(1,nrow=nrow(nodematch1),ncol=ncol(nodematch1)) #scores for noninf values in nodematch0
nodematch = nodematch2-(nodematch2*nodematch1) #2+1 node score
emdefault = 0
edgescore = 1 #base score for each edge below node, then minus its percentage (absolute edgematch difference)
# edgematch = matrix(emdefault,nrow=nrow(nodematch0),ncol=ncol(nodematch0)) #edge score




alpha = .2 #percentage of topology
beta = .8 #percentage of node similarity
theta = 0 # [0,1] penalty for multiple matches on one node; not used yet


# exact matches only









start = Sys.time()
# for (l in sort(unique(pm21$phenolevel),decreasing=T)) {
#   pind = which(pm21$phenolevel==l & pm21$phenotype%in%e1$from)
#   if (!length(pind)>0) next
#   imatched = which(nodematch[i,]>=0 & pm22$phenolevel>=l-layerlimit & pm22$phenolevel<=l+layerlimit & pm22$phenotype%in%e2$from)

#need to made between [0,1];; unused
edgematch = foreach (i=1:nrow(nodematch),.combine='rbind') %dopar% { #find all edges up from i
  iedges = e1[e1$to==pm21$phenotype[i],]
  edgematchi = rep(emdefault,nrow(pm22))
  if (!nrow(iedges)>0) return(edgematchi)

  imatched = which(nodematch[i,]>=0 & pm22$phenolevel>=pm21$phenolevel[i]-layerlimit & pm22$phenolevel<=pm21$phenolevel[i]+layerlimit & pm22$phenotype%in%e2$to)
  # edgematch[i,imatched] = sapply(imatched, function(j) {
  edgematchi[imatched] = sapply(imatched, function(j) {
    jedges = e2[e2$to==pm22$phenotype[j],]
    ije = match(iedges$addedm,jedges$addedm)
    ie = which(!is.na(ije))
    if (length(ie)>0) {
      je = ije[ie]
      return(sum( rep(edgescore,length(ie)) * (1-abs(iedges$evalue[ie]-jedges$evalue[je])) ))
      # edgematch[i,j] = sum( rep(edgescore,length(ie)) * (1-abs(iedges$evalue[ie]-jedges$evalue[je])) )
      # nonmatch = nrow(jedges)+nrow(iedges)-2*length(ie)+edgematch[i,j]
      # edgematch[i,j] = sum( rep(edgescore,length(ie)) * (1-abs(iedges$evalue[ie]-jedges$evalue[je])) ) / nonmatch #over # of non-matching edges
    } else { return(emdefault) }
  })

  return(edgematchi)
}
# }
TimeOutput(start)

dimnames(nodematch0) = dimnames(nodematch1) = dimnames(nodematch) = dimnames(edgematch) = dimnames(nodematch)


# sort(unique(as.vector(nodematch)))
# sort(unique(exp(as.vector(nodematch))))
# sort(unique(log(as.vector(nodematch))))









rind = which(!pm21$phenotype%in%e1$to) #root nodes
lind = which(!pm21$phenotype%in%e1$from) #leaf nodes
match = alpha*nodematch #+ beta*edgematch #default for all nodes, do no alpha for leaf nodes?
parentno = matrix(0,nrow=nrow(nodematch),ncol=ncol(nodematch))

if (length(unique(pm21$phenolevel))>1) {
  matchitb = NULL
  for (l in sort(unique(pm21$phenolevel),decreasing=T)[-1]) {
    pind = which(pm21$phenolevel==l & pm21$phenotype%in%e1$from) #for all nonleaf nodes
    if (!length(pind)>0) next
    enorm = 1-log(l+1,length(markers)+1) #used to normalize edge scores

    # result = foreach (i=pind,.combine='rbind') %dopar% { #find all edges up from i
    for (i in pind) {
      iedges = e1[e1$from==pm21$phenotype[i],]
      # matchi = nodematch[i,]
      imatched = which(matchi>=0 & pm22$phenolevel>=pm21$phenolevel[i]-layerlimit & pm22$phenolevel<=pm21$phenolevel[i]+layerlimit)# & pm22$phenotype%in%e2$from)
      # edgematch[i,imatched] = sapply(imatched, function(j) {
      
      matchitb[[as.character(i)]] = NULL #traceback, lists children chosen
      for (j in imatched) {
        jedges = e2[e2$from==pm22$phenotype[j],]
        ije = match(iedges$addedm,jedges$addedm)
        if (!all(is.na(ije))) {
          ie = iedges[!is.na(ije),]
          je = jedges[ije[!is.na(ije)],]
          for (ij in 1:nrow(ie)) {
            k = which(rownames(nodematch)==ie$to[ij])
            escore = edgescore * (1-abs(ie$evalue[ij]-je$evalue[ij]))
            escore = escore * enorm #normalized over level
            # escore = escore / (escore - 2*nrow(ie) + nrow(jedges) + nrow(iedges)) #normalized over all edges
            
            kupedges = e1[e1$to==ie$to[ij]]
            nescore = (match[i,j]+beta*escore)/nrow(kupedges)
            if (nescore > max(match[i,]) | j==which.max(match[i,])) {
              matchitb[[as.character(i)]][[as.character(j)]] = append(matchitb[[as.character(i)]],k)
              match[i,j] = match[i,j] + nescore/nrow(kupedges)
            } else { match[i,j] = match[i,j] + max(match[i,])/nrow(kupedges) }
          }
        }
      }
      # return(list(matchi=matchi,matchitb=matchitb))
    }
    
    match[pind,]
    # for i, match with all on pm22, record max score for each; that way next level can do each too
    for (i in pind) {
      iedges = e1[e1$from==pm21$phenotype[i],] #nrow=0
      imatched = which(match[i,]>=0 & pm22$phenolevel>=pm21$phenolevel[i]-layerlimit & pm22$phenolevel<=pm21$phenolevel[i]+layerlimit) #& pm22$phenotype%in%e2$from)
      for (j in imatched) {
        jedges = e2[e2$from==pm22$phenotype[j],] #nrow=0
        if (nrow(jedges)>0) {
          ije = match(iedges$addedm,jedges$addedm)
          ie = which(!is.na(ije))
          if (length(ie)>0) {
            je = ije[ie]
  
                      #iedges$evalue[ie]-jedges$evalue[je]
            # add up edge and node scores according to alpha beta: T = 
            match[i,j] = match[i,j] + sum( sapply(1:length(ie), function(x) return(match[ rownames(match)==iedges$to[ie[x]], colnames(match)==jedges$to[je[x]]])) )
          }
        }
      }
    }
  }
}
#add up scores for all parentless nodes










