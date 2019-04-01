# Count Stats: significant cell population stats
# aya43@sfu.ca 20161220

#Directory
root = "~/projects/IMPC"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#,"Sanger_MLN","CIPHE","TCP","H")

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
matrixPvalTRIM_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixPvalTRIM_CountAdj.Rdata",sep="")
rchy_dir = paste(result_dir, "/", panelL, "/", centreL, "/rchy", sep="")

#Output
stats_dir = paste(result_dir, "/", panelL, "/", centreL, "/stats", sep=""); for (ri in 1:length(stats_dir)) { suppressWarnings(dir.create(stats_dir[ri])) }
sigGene_dir = paste(stats_dir, "/sigGene", sep="") #sig nodes in common per gene
perFile_dir = paste(stats_dir, "/perFile", sep="") #sig nodes & layers per file
sigFile_dir = paste(stats_dir, "/sigFile", sep="") #sig nodes in common per file
rchyGene_dir = paste(stats_dir, "/rchyGene", sep="") #should be max of all included
rchyFile_dir = paste(stats_dir, "/rchyFile", sep="") #should be max of all included
rchyEdgesFile_dir = paste(stats_dir, "/rchyEdgesFile", sep="") #should be max of all included
rchyNodesFile_dir = paste(stats_dir, "/rchyNodesFile", sep="") #should be max of all included
rchyEdgesFileDiff_dir = paste(stats_dir, "/rchyEdgesFileDiff", sep="") #should be max of all included
rchyNodesFileDiff_dir = paste(stats_dir, "/rchyNodesFileDiff", sep="") #should be max of all included

source("~/projects/IMPC/code/_funcAlice.R")
libr("foreach")
libr("doMC")
libr("RchyOptimyx")

#Setup Cores
no_cores = 1#detectCores()-1
registerDoMC(no_cores)







#Options for script
interestedCols = c("fileName","gene")











start = Sys.time()

for (ci in 1:length(paste0(panelL,centreL))) {
  start1 = Sys.time()
  cat("\n",paste0(panelL," ",centreL)[ci],"; loading matrix ", sep="")
  
  m0 = get(load(matrixPvalTRIM_dir[ci]))
  sampleMeta0 = get(load(sampleMeta_dir[ci]))
  sampleMeta = sampleMeta0[match(rownames(m0),sampleMeta0$fileName),]
  phenoMeta0 = get(load(phenoMeta_dir[ci]))
  phenoMeta = phenoMeta0[match(colnames(m0),phenoMeta0$phenotype),]
  
  
  uniqueVals = sort(unique(sampleMeta[,interestedCols[2]]))
  
  
  #rchy stats
  rchyfolders = list.dirs(rchy_dir[ci],full.names=T,recursive=F)
  rchyfolders = rchyfolders[sapply(rchyfolders, function(x) length(list.files(x))>0)]
  fileNames0 = lapply(rchyfolders, function(rf) list.files(rf,full.names=T,recursive=F,pattern="_merged"))
  fileNames0 = lapply(fileNames0, function(x) gsub("_merged.Rdata",".fcs",fileNames(x)))
  fn = unique(Reduce("union",fileNames0))
  
  rchyNodes0 = rchyEdges0 = matrix(0,nrow=length(fn),ncol=length(rchyfolders),dimnames=list(fn,NULL))

  #break on task 1
  result = foreach (rfi = 1:length(rchyfolders)) %dopar% {
    nodesAll = edgesAll = rchyEdges0 = rchyEdges0 = rchyNodes = rchyEdges = rchyNodesDiff = rchyEdgesDiff = NULL #all nodes/edges, gene #of , file #of
    for (gi in 1:length(uniqueVals)) {
      fnInG = fileNames0[[rfi]][fileNames0[[rfi]] %in% sampleMeta$fileName[sampleMeta$gene==uniqueVals[gi]]]
      if (length(fnInG)>0) {
        nodes = edges = NA
        for (fnii in 1:length(fnInG)) {
          fni = fnInG[fnii]
          
          rchy = get(load(paste0(rchyfolders[rfi],"/",fni,"_merged.Rdata")))
          nodesAll[[fni]] = nodes0 = Reduce("union",sapply(rchy,function(x) x@nodes[1,]))
          edgesAll[[fni]] = edges0 = Reduce("union",sapply(rchy,function(x) x@edges[1,]))
          
          if(!length(nodes)==0) {
            if (is.na(nodes)) { 
              nodes = nodes0
            } else { 
              nodes = intersect(nodes,nodes0)
            }
          }
          rchyNodes0[fni] = length(nodes0)
          
          if(!length(edges)==0) {
            if (is.na(edges)) { edges = edges0
            } else {
              edges = intersect(edges,edges0)
            }
          }
          rchyEdges0[fni] = length(edges0)
        }
        rchyEdges[gi] = length(edges)
        rchyNodes[gi] = length(nodes)
        rchyNodesDiff[gi] = length(setdiff(Reduce('union',nodesAll[fnInG]), nodes))
        rchyEdgesDiff[gi] = length(setdiff(Reduce('union',edgesAll[fnInG]), edges))
      } else {
        rchyEdges[gi] = rchyNodes[gi] = rchyNodesDiff[gi] = rchyEdgesDiff[gi] = 0
      }
    }
    rchyNodes0 = rchyNodes0[match(fn,names(rchyNodes0))]
    rchyEdges0 = rchyEdges0[match(fn,names(rchyEdges0))]
    rchyFile = cbind(rchyNodes0,rchyEdges0)
    colnames(rchyFile) = c(paste0("nodes_Rchy_",fileNames(rchyfolders[rfi])),paste0("edges_Rchy)",fileNames(rchyfolders[rfi])))
    rchyGene = cbind(rchyEdges,rchyNodes,rchyNodesDiff,rchyEdgesDiff)
    colnames(rchyGene) = c(paste0("nodesInCommon_Rchy_",fileNames(rchyfolders[rfi])),paste0("edgesInCommon_Rchy",fileNames(rchyfolders[rfi])),
                           paste0("edgesNotInCommon_Rchy",fileNames(rchyfolders[rfi])),paste0("edgesNotInCommon_Rchy",fileNames(rchyfolders[rfi])))
    
    #inCommon between files
    rchyNodesCommon = rchyEdgesCommon = rchyNodesDiff = rchyEdgesDiff = matrix(0,nrow=length(fileNames0[[rfi]]),ncol=length(fileNames0[[rfi]]))
    for (fni in 2:length(fileNames0[[rfi]])) {
      for (fnj in 1:(fni-1)) {
        rchyNodesCommon[fni,fnj] = rchyNodesCommon[fnj,fni] = length(intersect(nodesAll[[fni]],nodesAll[[fnj]]))
        rchyEdgesCommon[fni,fnj] = rchyEdgesCommon[fnj,fni] = length(intersect(edgesAll[[fni]],edgesAll[[fnj]]))
        rchyNodesDiff[fni,fnj] = rchyNodesDiff[fnj,fni] = length(union(setdiff(nodesAll[[fni]],nodesAll[[fnj]]),setdiff(nodesAll[[fnj]],nodesAll[[fni]])))
        rchyEdgesDiff[fni,fnj] = rchyEdgesDiff[fnj,fni] = length(union(setdiff(nodesAll[[fni]],nodesAll[[fnj]]),setdiff(edgesAll[[fnj]],edgesAll[[fni]])))
      }
    }
    rchyNodesCommon = cbind(sampleMeta$gene[match(fileNames0[[rfi]],sampleMeta$fileName)],fileNames0[[rfi]],rchyNodesCommon)
    rchyEdgesCommon = cbind(sampleMeta$gene[match(fileNames0[[rfi]],sampleMeta$fileName)],fileNames0[[rfi]],rchyEdgesCommon)
    rchyNodesDiff = cbind(sampleMeta$gene[match(fileNames0[[rfi]],sampleMeta$fileName)],fileNames0[[rfi]],rchyNodesDiff)
    rchyEdgesDiff = cbind(sampleMeta$gene[match(fileNames0[[rfi]],sampleMeta$fileName)],fileNames0[[rfi]],rchyEdgesDiff)
    colnames(rchyNodesCommon) = colnames(rchyEdgesCommon) = colnames(rchyNodesDiff) = colnames(rchyEdgesDiff) = c("gene","fileName",fileNames0[[rfi]])
    
    save(rchyNodesCommon,file=paste0(rchyNodesFile_dir,"_",fileNames(rchyfolders[rfi]),".Rdata"))
    write.csv(rchyNodesCommon,file=paste0(rchyNodesFile_dir,"_",fileNames(rchyfolders[rfi]),".csv"))
    save(rchyEdgesCommon,file=paste0(rchyEdgesFile_dir,"_",fileNames(rchyfolders[rfi]),".Rdata"))
    write.csv(rchyEdgesCommon,file=paste0(rchyEdgesFile_dir,"_",fileNames(rchyfolders[rfi]),".csv"))
    save(rchyNodesDiff,file=paste0(rchyNodesFileDiff_dir,"_",fileNames(rchyfolders[rfi]),".Rdata"))
    write.csv(rchyNodesDiff,file=paste0(rchyNodesFileDiff_dir,"_",fileNames(rchyfolders[rfi]),".csv"))
    save(rchyEdgesDiff,file=paste0(rchyEdgesFileDiff_dir,"_",fileNames(rchyfolders[rfi]),".Rdata"))
    write.csv(rchyEdgesDiff,file=paste0(rchyEdgesFileDiff_dir,"_",fileNames(rchyfolders[rfi]),".csv"))
    
    return(list(rchyGene=rchyGene,rchyFile=rchyFile))
  }
  rchyGene0 = cbind(uniqueVals,Reduce('cbind',lapply(result,function(x) x$rchyGene)))
  rchyFile0 = cbind(fn,Reduce('cbind',lapply(result,function(x) x$rchyFile)))
  colnames(rchyGene0)[1] = "gene"
  save(rchyGene0,file=paste0(rchyGene_dir,".Rdata"))
  write.csv(rchyGene0,file=paste0(rchyGene_dir,".csv"))
  colnames(rchyFile0)[1] = "file"
  save(rchyFile0,file=paste0(rchyFile_dir,".Rdata"))
  write.csv(rchyFile0,file=paste0(rchyFile_dir,".csv"))
  
  

  TimeOutput(start1)
  
  
  #file stats
  nLayers = nSig = NULL
  for (fi in 1:nrow(m0)) {
    sigs = m0[fi,]!=0
    nLayers[fi] = max(phenoMeta$phenolevel[sigs])
    nSig[fi] = sum(sigs)
  }
  inFile = cbind(sampleMeta[,interestedCols],nLayers,nSig)
  colnames(inFile) = c(interestedCols,"nlayers","nSigNodes")
  rownames(inFile) = NULL
  save(inFile,file=paste0(perFile_dir,".Rdata"))
  write.csv(inFile,file=paste0(perFile_dir,".csv"))
  
  #sig nodes in common between files
  # sigFile0 = matrix(0,nrow=nrow(m0),ncol=nrow(m0))
  # sigFileNot0 = matrix(0,nrow=nrow(m0),ncol=nrow(m0))
  # 
  # break on task 2
  loop.ind = triLoopInd(nrow(m0),no_cores)
  sigF = foreach (i = 1:length(loop.ind)) %dopar% {
    sigFile0 = sigFileNot0 = list()
    indexx = 1
    for (fi in loop.ind[[i]][1]:loop.ind[[i]][2]) {
      sigFile = rep(0,nrow(m0))
      sigFileNot = rep(0,nrow(m0))
      for (fj in (fi+1):nrow(m0)) {
        sigFile[fj] = sum(apply(m0[c(fi,fj),], 2, function(x) all(x!=0)))
        sigFileNot[fj] = sum(apply(m0[c(fi,fj),], 2, function(x) !all(x!=0) & any(x!=0) ))
      }
      sigFile0[[indexx]] = sigFile
      sigFileNot0[[indexx]] = sigFileNot
      indexx = indexx+1
    }
    return(list(sigFile0=Reduce('cbind',sigFile0),sigFileNot0=Reduce('cbind',sigFileNot0)))
  }
  sigFileNot = cbind(rep(0,nrow(m0)), Reduce('cbind',lapply(sigF,function(x)x$sigFileNot0)) )
  for (fi in 2:nrow(m0)) { for (fj in 1:(fi-1)) { sigFileNot[fj,fi] = sigFileNot[fi,fj] } }
  sigFileNot = cbind(sampleMeta[,interestedCols],sigFileNot)
  
  sigFile = cbind(rep(0,nrow(m0)), Reduce('cbind',lapply(sigF,function(x)x$sigFile0)))
  for (fi in 2:nrow(m0)) { for (fj in 1:(fi-1)) { sigFile[fj,fi] = sigFile[fi,fj] } }
  sigFile = cbind(sampleMeta[,interestedCols],sigFile)
  
  colnames(sigFile) = colnames(sigFileNot) = c(interestedCols,rownames(m0))
  rownames(sigFile) = rownames(sigFileNot) = NULL
  save(sigFile,file=paste0(sigFile_dir,".Rdata"))
  write.csv(sigFile,file=paste0(sigFile_dir,".csv"))
  save(sigFileNot,file=paste0(sigFile_dir,"Not.Rdata"))
  write.csv(sigFileNot,file=paste0(sigFile_dir,"Not.csv"))
  
  TimeOutput(start1)
  
  
  #sig nodes in common per gene
  inCommonDir = inCommon = NULL
  # for (uvi in 1:length(uniqueVals)) {
  result = foreach (uvi = 1:length(uniqueVals)) %dopar% {
    inCommonDir = inCommon = NULL
    uv = uniqueVals[uvi]
    muv = m0[sampleMeta$gene==uv,]
    if (is.null(dim(muv))) { muv = matrix(muv,nrow=1); colnames(muv) = colnames(m0) }
    
    incommonnNOT = apply(muv, 2, function(x) !all(x!=0) & any(x!=0))
    
    incommonn = apply(muv, 2, function(x) all(x!=0))
    incommonDirr = apply(muv, 2, function(x) all(x>0)) | apply(muv, 2, function(x) all(x<0))==T
    
    return(list(inCommon=sum(incommonn),inCommonDir=sum(incommonDirr),inCommonNOT=sum(incommonnNOT)))
    # inCommon[uvi] = sum(incommonn)
    # inCommonDir[uvi] = sum(incommonDirr)
  }
  inCommon = sapply(result, function(x) x$inCommon)
  inCommonDir = sapply(result, function(x) x$inCommonDir)
  incommonnNOT = sapply(result, function(x) x$incommonNOT)
  inComm = cbind(uniqueVals,table(sampleMeta$gene),inCommonDir,inCommon,incommonnNOT)
  rownames(inComm) = NULL
  colnames(inComm) = c("gene","filesPerGene","nodesInCommon_UpDown","nodesInCommon_noUpDown","nodesNotInCommon")
  save(inComm,file=paste0(sigGene_dir,".Rdata"))
  write.csv(inComm,file=paste0(sigGene_dir,".csv"))
  
  
  TimeOutput(start1)
  
}


TimeOutput(start)
