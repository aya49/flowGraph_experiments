# aya43@sfu.ca 20161220
# Tries to clean out surrogate variable effects from normalized cell count matrix

root <- "~/projects/flowCAP-II"
result_dir <- "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

#Input
ft_dir <- "data/FT"
csv_dir <- "data/AML.csv"
phenoMeta_dir <- paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir <- paste(result_dir, "/sampleMeta.Rdata", sep="")
matrixCount_dir <- paste(result_dir, "/matrixCount.Rdata", sep="")
matrixCountAdj_dir <- paste(result_dir, "/matrixCountAdj.Rdata",sep="")
matrixProp_dir <- paste(result_dir, "/matrixProp.Rdata", sep="")

#Output


libr(stringr)
libr(RUVSeq)
libr(edgeR)
libr(EDASeq)
libr(scater) #pData
source("code/_funcAlice.R")



start <- Sys.time()

sampleMeta <- get(load(sampleMeta_dir))
m <- get(load(matrixCountAdj_dir))
#phenotype on rows
if (ncol(m)!=nrow(sampleMeta)) { m <- t(m) }















libr(zebrafishRNASeq)
libr(RColorBrewer)
data(zfGenes)

## Exploratory Analysis

#get rid of non-expressing genes, at least 5 reads in at least 2 samples
filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
filtered <- zfGenes[filter,]
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
spikes <- rownames(filtered)[grep("^ERCC*", rownames(filtered))]


#store object as class seqExpressionSet
x <- as.factor(rep(c("Ctl", "Trt"), each=3))
set <- newSeqExpressionSet( as.matrix(filtered), phenoData=data.frame(x, row.names=colnames(filtered)))
set

#boxplots
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

#normalize the data using upper-quartile (UQ) normalization
set0 <- betweenLaneNormalization(set, which="upper")
plotRLE(set0, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set0, col=colors[x], cex=1.2)



## Unwanted Var using Control Genes: genes not influenced by vocariates of interest (purely unwanted)
set1 <- RUVg(set0, spikes, k=1)
pData(set1)

plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set1, col=colors[x], cex=1.2)



## Unwanted Var without using Control Genes: find them empirically

design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[ which(!(rownames(set) %in% rownames(top)[1:5000]))]

#use top 5000 empiricals ranked by edgeR pvalues
set2 <- RUVg(set, empirical, k=1)
pData(set2)

plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set2, col=colors[x], cex=1.2)


## differential analysis
libr(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=counts(set2), colData=pData(set2), design=~W_1+x)
dds <- DESeq(dds)
res <- results(dds)
res

design <- model.matrix(~x + W_1, data=pData(set1))
y <- DGEList(counts=counts(set1), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)




















## Analysis of simulated data

libr(zebrafishRNASeq)
libr(RSkittleBrewer)
libr(genefilter)
libr(polyester)
libr(RUVSeq)
libr(edgeR)
libr(sva)
libr(ffpe)
libr(RColorBrewer)
libr(corrplot)
libr(limma)
trop = RSkittleBrewer('tropical')
set.seed(3532333)


data(zfGenes)
filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
counts <- zfGenes[filter,]

plot(rowMeans(log(counts+1)),rowVars(log(counts+1)),pch=19,col=trop[1])

## Estimate the zero inflated negative binomial parameters
params = get_params(counts)

plot(rowMeans(log(counts+1)),rowVars(log(counts+1)),pch=19,col=trop[1],main="zebrafish data w/fit")
lines(params$fit,col=trop[2])





group = rep(c(-1,1),each=10)
batch = 1 - 2*rbinom(20,size=1,prob=0.5)
gcoeffs = c(rnorm(400),rep(0,600))
nullindex=401:1000
bcoeffs = c(rep(0,100),rnorm(400,sd=1),rep(0,500))
coeffs = cbind(bcoeffs,gcoeffs)
controls = (bcoeffs != 0) & (gcoeffs==0)
mod = model.matrix(~-1 + batch + group)

dat0 = create_read_numbers(params$mu,params$fit,
                           params$p0,m=dim(counts)[1],n=dim(counts)[2],
                           beta=coeffs,mod=mod,seed=4353)
rm(mod)


## Set null and alternative models (ignore batch)
mod1 = model.matrix(~group)
mod0 = cbind(mod1[,1])

## Estimate batch with svaseq (unsupervised)
batch_unsup_sva = svaseq(dat0,mod1,mod0)$sv

## Estimate batch with svaseq (supervised)
batch_sup_sva = svaseq(dat0,mod1,mod0,controls=controls)$sv

## Estimate batch with pca
ldat0 = log(dat0 + 1)
batch_pca = svd(ldat0 - rowMeans(ldat0))$v[,1]

## Estimate batch with ruv (controls known)
batch_ruv_cp <- RUVg(dat0, cIdx= controls, k=1)$W

## Estimate batch with ruv (residuals)
## this procedure follows the RUVSeq vignette
## http://www.bioconductor.org/packages/devel/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf

x <- as.factor(group)
design <- model.matrix(~x)
y <- DGEList(counts=dat0, group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")
seqUQ <- betweenLaneNormalization(dat0, which="upper")
controls = rep(TRUE,dim(dat0)[1])
batch_ruv_res = RUVr(seqUQ,controls,k=1,res)$W

## Estimate batch with ruv empirical controls
## this procedure follows the RUVSeq vignette
## http://www.bioconductor.org/packages/devel/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf

y <- DGEList(counts=dat0, group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

controls =rank(lrt$table$LR) <= 400
batch_ruv_emp <- RUVg(dat0, controls, k=1)$W


## Plot the results
plot(batch,col=trop[1],pch=19,main="batch")

