#extract gene names from sampleMeta
#INCORPORATED INTO network/data.r; DELETE!


root = "~/projects/IMPC"
result_dir = "result"
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN") #"Sanger_SPLEEN","Sanger_MLN","CIPHE","TCP",
ci=1

gene_dir = paste(result_dir, "/", panelL, "/", centreL, "/gene.txt", sep="")
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
sampleMeta = get(load(sampleMeta_dir[ci]))

libr(stringr)


gene0 = unique(sampleMeta$gene)
gene = gene0[!grepl("failed", gene0)]
gene = gene0[!grepl(" KO$", gene0)]
gene = str_split(gene,"_")
gene = sapply(gene, function(x) x[1])
gene = str_split(gene,"[(]")
gene = sapply(gene, function(x) x[1])
gene = gene[!gene%in%"+"]
write.table(gene,file=gene_dir,sep="\t",row.names=F,col.names=F,quote=F)







