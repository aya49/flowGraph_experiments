# Refine distances using random walk
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
rwThres = c(0,.05,.1,.2,.5,.75,1)
steps = 5 # * number of edges

ignoredist = ".csv|_rw|_pr"



#Output

libr(stringr)
libr(foreach)
libr(doMC)
libr(igraph)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")



start = Sys.time()

no_cores = detectCores()-1
registerDoMC(no_cores)

# ftWTGT = "normal"

sampleMeta0 = get(load(sampleMeta_dir))
distmfile = list.files(dist_dir, recursive=F, full.names=T, pattern=".Rdata")
distmfile = distmfile[grep("layer07",distmfile)]
distmfile = distmfile[grep("manhattan",distmfile)]
distmfile = distmfile[grep("cellpop",distmfile)]
delpvalind = grepl(ignoredist,distmfile,ignore.case=T)
if (sum(delpvalind)>0) distmfile = distmfile[!delpvalind]
distmfilenames = fileNames(distmfile)
distmfilelayer = str_extract(distmfile, "layer[0-9]+")

#calculate score_nac for each distance matrix & interested classes
#result = list()
for(i in order(distmfilelayer,decreasing=T)) {
#a = foreach(i=order(distmfilelayer,decreasing=T)) %dopar% {
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
    
    d1 = get_graph(d0)
    # vector on row; going out, sums to 1
    d = d1/rowSums(d1)
    gr = graph_from_adjacency_matrix(d, weighted=T,mode='directed', diag=F)
    
    cat(", random walk ")
    rw = as.vector(random_walk(gr,start=sample(as.vector(V(gr)),1), steps=steps*length(E(gr))))
    # df = matrix(0,nrow=nrow(d),ncol=ncol(d))
    df = foreach (v = sort(as.vector(V(gr))), .combine='rbind') %dopar% {
    # for (v in sort(as.vector(V(gr)))) {
      u1 = rw[c(rw[-1]==v,F)]
      u2 = rw[c(F,rw[-length(rw)]==v)]
      u = append(u1,u2)
      ut = table(u)
      
      # df[as.integer(names(ut))] = ut
      a = rep(0,length(as.vector(V(gr))))
      a[as.integer(names(ut))] = ut
      return(a)
    }
    
    rownames(df) = colnames(df) = rownames(d0)
    df0 = df

    save(df0,file=paste0(gsub(".Rdata","",distmfile[i]),"_rw_simmatrix.Rdata"))
    for (rwt in rwThres) {
      df1 = df0
      tops = quantile(as.vector(df0),rwt)
      df1[df0<tops] = 0
      df = 1/df1
      df[df1==0] = 2/min(df0[df0>0])
      diag(df) = 0
      dfd1 = as.dist(df)
      save(dfd1,file=paste0(gsub(".Rdata","",distmfile[i]),"_rw_",rwt,".Rdata"))
    }
    TimeOutput(start1)
  }, error = function(e) {
    cat(paste("ERROR:  ",e)); graphics.off()
  })
  
}


TimeOutput(start)








