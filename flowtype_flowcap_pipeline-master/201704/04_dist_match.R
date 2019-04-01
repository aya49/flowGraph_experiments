sampleMeta[match(rownames(m),sampleMeta$fileName),]

compareind = c(6,7)
maind2 = m[compareind,]
colsanysig = which(sapply(1:ncol(maind2), function(x) any(as.matrix(maind2)[,x]>0))==T)
maind2 = maind2[,colsanysig]
pm2 = phenoMeta[!is.na(match(phenoMeta$phenotype,colnames(maind2))),] #same order as columns


## seperate out trimmed values for the 2 samples
colsallsig = which(sapply(1:ncol(maind2), function(x) all(as.matrix(maind2)[,x]>0))==T)

maind20 = maind2[,colsallsig]
maind21 = maind2[1,maind2[1,]>0]
maind22 = maind2[2,maind2[2,]>0]

pt0 = lapply(colnames(maind20), function(x) return(unlist(strsplit(x, "(?<=[+-])", perl=TRUE))))
pt1 = lapply(names(maind21), function(x) return(unlist(strsplit(x, "(?<=[+-])", perl=TRUE))))
pt2 = lapply(names(maind22), function(x) return(unlist(strsplit(x, "(?<=[+-])", perl=TRUE))))

pm20 = pm2[which(!is.na(match(pm2$phenotype,colnames(maind20)))),]
pm21 = pm2[which(!is.na(match(pm2$phenotype,names(maind21)))),]
pm22 = pm2[which(!is.na(match(pm2$phenotype,names(maind22)))),]


## candidate matches for each phenotype
matchfor1 = list()
for (i in 1:length(maind21)) {
  matchfor1[[i]] = sapply(1:length(pt2), function(x) {
    inter = intersect(pt2[[x]],pt1[[i]])
    a = -Inf
    if (length(inter)>0) a= -(length(pt2[[x]])+length(pt1[[i]])-(2*length(inter))) / length(pt1[[i]])
    return(a) #if no overlap
  })
}
names(matchfor1) = names(maind21)

matchfor2 = list()
for (i in 1:length(maind22)) {
  matchfor2[[i]] = sapply(1:length(pt1), function(x) {
    inter = intersect(pt1[[x]],pt2[[i]])
    a = -Inf
    if (length(inter)>0) a= -(length(pt1[[x]])+length(pt2[[i]])-(2*length(inter))) / length(pt2[[i]])
    return(a) #if no overlap
  })
}
names(matchfor2) = names(maind22)

## find best match layer by layer
for (l in sort(unique(pm21$phenolevel),decreasing=T)) {
  lind = which(pm21$phenolevel==l)
  matchfor1
}

