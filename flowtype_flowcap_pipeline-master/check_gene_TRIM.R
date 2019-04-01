sampleMeta$gene[!is.na(match(sampleMeta$fileName,rownames(m))]
aind = which(sampleMeta$gene[!is.na(match(sampleMeta$fileName,rownmes(m)))])=="")
maind = m[aind,]
maindlist = list()
for (i in 1:nrow(maind)) {
maindlist[[i]] = colnames(maind)[which(maind[i,]> -log(.025))]
}
maindl = unlist(maindlist)
sort(table(maindl))
