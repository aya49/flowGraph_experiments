# Refine distances using pagerank (multiply edges)
# aya43@sfu.ca 20170419

root = "~/projects/flowCAP-II"
result_dir = "result"
setwd(root)

options(stringsAsFactors=FALSE)
#options(device="cairo")
options(na.rm=T)

#Input
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
dist_dir = paste(result_dir, "/dist", sep="")

KOonly = F #delete all WT samples
adjustD = T #divide by a number if the numbers are too big e.g. manhattan
overwrite = T #redo and overwrite all past scores
rwThres = c(3,10,50,100,1000,2000)
steps = 10000

ignoredist = ".csv|_rw|_pr"


#Output

libr(stringr)
libr(foreach)
libr(doMC)
libr(igraph)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")



start = Sys.time()

no_cores = 7#detectCores()-1
registerDoMC(no_cores)

# ftWTGT = "normal"

sampleMeta0 = get(load(sampleMeta_dir))
distmfile = list.files(dist_dir, recursive=F, full.names=T, pattern=".Rdata")
delpvalind = grepl(ignoredist,distmfile,ignore.case=T)
if (sum(delpvalind)>0) distmfile = distmfile[!delpvalind]
distmfilenames = fileNames(distmfile)
distmfilelayer = str_extract(distmfile, "layer[0-9]+")

#calculate score_nac for each distance matrix & interested classes
#result = list()
#for(i in 1:length(distmfile)) {
foreach(i=order(distmfilelayer)) %dopar% {
  tryCatch({
    #loop.ind = 1:length(distmfile)
    #for(i in 1:length(distmfile)) {
    #result = foreach(i=loop.ind) %dopar% { 
    start1 = Sys.time()
    cat("\n", distmfile[i], ", loading dist matrix", sep="")
    #result[[i]] = list()
    d0 = as.matrix(get(load(distmfile[i])))
    # if (length(d)==1) { d=d[[1]]; save(d,file=distmfile[i])}
    if (sum(is.na(as.matrix(d0)))>0 | sum(as.matrix(d0)==Inf)>0) return(NULL)
    
    d = d0
    rownames(d) = colnames(d) = 1:nrow(d)
    d = 1/d
    d[d0==0] = 1/min(d0[d0>0])
    gr = graph_from_adjacency_matrix(d, weighted=T,mode='undirected', diag=F)
    
    
    #get line graph by multiplying weights :: http://stackoverflow.com/questions/28890790/how-to-computer-new-edge-weight-when-converting-vertices-to-edges-using-line-gra
    g.edges <- get.edges(gr, E(gr))
    enames <- paste(g.edges[,1], g.edges[,2], sep=",")
    ewts <- E(gr)$weight
    new.edges <- do.call(rbind, sapply(1:vcount(gr), function(x) {
      incident <- which(g.edges[,1] == x | g.edges[,2] == x)
      if (length(incident) <= 1) {
        return(NULL)
      } else {
        all.comb <- combn(incident, 2)
        return(data.frame(x=enames[all.comb[1,]], y=enames[all.comb[2,]], weight=apply(all.comb, 2, function(x) prod(ewts[x]))))
      }
    }))
    gre = graph.data.frame(new.edges)
    
    
    rw = page_rank(gre)$vector
    fromto = strsplit(names(rw),",")
    from = as.integer(sapply(fromto,function(y) y[1]))
    to = as.integer(sapply(fromto,function(y) y[2]))
    d1 = matrix(0,nrow=nrow(d),ncol=ncol(d))
    for (fromtoi in 1:length(rw)) {
      d1[to[fromtoi],from[fromtoi]] = d1[from[fromtoi],to[fromtoi]] = d1[to[fromtoi],from[fromtoi]]+rw[fromtoi]
    }
    d10 = d1
    d1 = 1/d10
    d1[d10==0] = max(d10[d10>0])
    diag(d1) = 0
    
    save(d1,file=paste0(gsub(".Rdata","",distmfile[i]),"_pr.Rdata"))
    
  }, error = function(e) {
    cat(paste("ERROR:  ",e)); graphics.off()
  })
  
}


TimeOutput(start)








