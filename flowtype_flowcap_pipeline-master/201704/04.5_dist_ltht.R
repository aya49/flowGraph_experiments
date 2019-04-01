# Refine distances using laplacian transform hitting time
# aya43@sfu.ca 20170424

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
rwThres = c(.05,.1,.2,.5,.75,1)
steps = 500 # * number of edges

ignoredist = ".csv|_rw|_pr"

beta = c(.01,.1,1,2,5)



#Output

libr(stringr)
libr(foreach)
libr(doMC)
libr(igraph)
libr(DTMCPack)
libr(e1071)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")



start = Sys.time()

no_cores = 8#detectCores()-1
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
a = foreach(i=order(distmfilelayer,decreasing=T)) %dopar% {
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
    
    for (b in beta) {
      d1 = get_graph(d0)
      # vector on row; going out, sums to 1
      dd = d1/rowSums(d1)
      
      d = 2
      
      #stationary distribution
      pi = statdistr(dd)
      degree = apply(dd,1,function(x) sum(x>0))
      dest = pi^(d/(d+2.0)) * degree^(2/(d+2.0)); dest = dest/sum(dest)
      epsest = pi^(-1.0/(d+2.0)) * degree^(1.0/(d+2.0))
      distmat = diag(epsest)*dd; diag(distmat) = 0
      sps = allShortestPaths(distmat)$length; sps = sps + t(sps)
      spd = (.5*(sps+t(sps)))^2
      nb = lapply(1:nrow(dd), function(i) which(dd[i,]!=0))
      fasthittime = function(orig,nsamples,runtime,nb) {
        klens = append(0,sapply(nb,function(x) length(x)))
        kind = cumsum(klens)
        kflat = unlist(nb)
        j = orig
        t = runtime
        hittimes = matrix(0,nrow=length(nb),ncol=runtime)
        hittemp = rep(0,length(nb))
        for(i in 1:nsamples) {
          for (k in 1:nrow(hittimes)) {
            kl = kind[j+1]-kind[j]
            jn = sample(1:kl,1)
            j = kflat[kind[j]+jn]
            if(hittemp[j]<0) {
              hittemp[j] = t
              
            }
          }
        }
      }
        
      
      save(df,file=paste0(gsub(".Rdata","",distmfile[i]),"_ltht_",b,".Rdata"))
    }

  }, error = function(e) {
    cat(paste("ERROR:  ",e)); graphics.off()
  })
  
}


TimeOutput(start)








