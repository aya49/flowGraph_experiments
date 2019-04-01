# aya43@sfu.ca 20161220
# Distance metric learning

root = "~/projects/flowCAP-II"
result_dir = "result"
setwd(root)

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
matrixCount_dir = paste(result_dir, "/matrixCount.Rdata", sep="")
matrixCountAdj_dir = paste(result_dir, "/matrixCountAdj.Rdata", sep="")
matrixProp_dir = paste(result_dir, "/matrixProp.Rdata", sep="")

libr(MASS)
libr(scatterplot3d)
libr(dml)
libr(ggplot2)
libr(readr)
libr(Rtsne)
libr(plot3D)
libr(rgl)
set.seed(1)

load(sampleMeta_dir)
load(phenoMeta_dir)
m = get(load(matrixCountAdj_dir))

##Order by Class
target = as.numeric(sampleMeta$aml)
amlorder = order(sampleMeta$aml)
sampleMeta = sampleMeta[amlorder,]
m = m[amlorder,]

# generate simulated Gaussian data
k = length(which(sampleMeta$aml=="aml")) #number of aml patients
n = nrow(m)

# define similar constrains
simi = rbind(t(combn(1:k, 2)), t(combn((k+1):n, 2)))

temp =  as.data.frame(t(simi))
tol = as.data.frame(combn(1:n, 2))

# define disimilar constrains
dism = t(as.matrix(tol[!tol %in% simi]))

# trim matrix

# transform data using GdmDiag
print(Sys.time())
result = GdmDiag(m, simi, dism)
print(Sys.time())
m2 = result$newData
# plot original data
color = gl(2, k, labels = c("red", "blue"))

tsne = Rtsne(m, check_duplicates = FALSE, pca = TRUE, perplexity=30, theta=0.5, dims=3)
tsne2 = Rtsne(m2, check_duplicates = FALSE, pca = TRUE, perplexity=30, theta=0.5, dims=3)
print(Sys.time())

tsne$Y$class = target
tsne2$Y$class = target

plot = plot3d(x=tsne$Y$V1, y=tsne$Y$V2, z=tsne$Y$V3, col=tsne$Y$class, pch = ".")
plot2 = plot3d(x=tsne2$Y$V1, y=tsne2$Y$V2, z=tsne2$Y$V3, col=tsne2$Y$class, pch = ".")