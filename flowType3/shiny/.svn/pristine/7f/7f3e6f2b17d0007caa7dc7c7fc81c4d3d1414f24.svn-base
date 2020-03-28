#Helperfunctions for flowType shiny
#Written by Mehrnoush Malek (MM)
#Revised: March 2020, MM

make.sliding.graph <- function(arguments,ft, counts, method1, method2, p.level,count.lim,parent="CD45+")
  #This is a function that grabs the required arguments from rchy results to use for visNetwork
 # For this used case  labels are created randomly for the samples.
{
  print("start")
  colnames(counts)[1] <- parent
  pop.length <- unlist(lapply(ft@PhenoCodes, function(x) length(which(as.numeric(unlist(strsplit(x,split = "")))!=0))))
  ind <- which (pop.length < (p.level+1))
  counts <- counts[,ind]
  ft.props <- counts/counts[,1]
  colnames(ft.props) <- colnames(counts)
 
  lbls <-sample (x = 0:1,size =nrow(counts) ,replace = T)
  pvals <- apply(ft.props,2, function(prop.x) {
    if (all(prop.x==1))
      return(1)
    if (method1=="Wilcoxon")
    { 
      p<-tryCatch(wilcox.test(prop.x[which(lbls==0)],prop.x[which(lbls==1)])$p.value, error=function(x) {return("wrong")})
      
    }else{
      p<-tryCatch(t.test(prop.x[which(lbls==0)],prop.x[which(lbls==1)])$p.value, error=function(x) {return("wrong")})
    }
    if (mode(p)=="character" | is.nan(p))
    {
      #print("t.test failed, as data was almost constant.")
      p<-1
    }
    return(p)
  })
  names(pvals)<- colnames(ft.props)
  if (method2!="NULL"){
    pvals.adjust <- p.adjust(pvals, method=method2)
  }else{
    pvals.adjust <-pvals
  }
  selected.pheno <- ft@PhenoCodes[which(pvals.adjust<0.08)]
  if (length(selected.pheno)>10)
    selected.pheno <- find.farthest.pheno(pheno = selected.pheno,min.limit = 8)
  Scores <- -log10(pvals.adjust)
  if(length(selected.pheno)<1)
  {
    print("oops")
    rchy.plot <- gplots::textplot(object = "Nothing significant")
  }else{
    rchy.res <- list()
    
    for (i in 1:length(selected.pheno))
    {
      rchy.res[[i]]<-RchyOptimyx(pheno.codes=ft@PhenoCodes[ind], phenotypeScores=Scores,
                                 startPhenotype=selected.pheno[i], pathCount=1,trimPaths=FALSE) 
    }
    merged.res <- rchy.res[[1]]
    i1<-2
    while(i1 < length(rchy.res))
    {
      merged.res <- merge(merged.res,rchy.res[[i1]])
      i1 <- i1+1
    }
    cell.props <-colMeans(ft.props)
    cell.counts <-colMeans(counts)
    pop.names <- sapply(ft@PhenoCodes, flowTypeFilter::decodePhenotype,marker.names =arguments$names , partitions.per.marker=arguments$part.v,arguments$filt.v)
    nodes.lbl <-pop.names[match(merged.res@nodes[1,],
                                names(pop.names))]
    nodes.lbl[1] <- parent
    colors=c('blue','cyan','yellow','red')
    colorFunc <- colorRampPalette(colors)

    z=pretty(c(min(Scores[nodes.lbl]), max(Scores[nodes.lbl])))
    colinds=unlist(lapply(1:length(nodes.lbl), function(i){max(1,ceiling((Scores[nodes.lbl]-min(z))*1000/(max(z)-min(z)))[i])}))
    cols=colorFunc(1000)[colinds]

    names(cols) <-names(Scores[nodes.lbl])
    white.inds<- which(cell.counts[nodes.lbl] < count.lim)
    if (length(white.inds)>0)
    {
      cols[white.inds]<-"#FFFFFF"
      
      nodes.lbl[white.inds]<-sapply(white.inds, function(i1) sprintf(paste0("%-",nchar(as.vector( nodes.lbl[i1])),"s"),""))
    }
    node.levels <- as.vector(sapply(merged.res@nodes[1,],function(x) length(which(unlist(strsplit(x,""))!=0))))+1
    
    nodes <- data.frame(id=merged.res@nodes[1,],label= nodes.lbl,color=cols,shape="ellipse",level=node.levels)
    
    for (i in 1:dim(merged.res@edges)[2])
      merged.res@edges[2,i] <- max(0, as.numeric(merged.res@edges[2,i]))
    edgevalues <- as.numeric(merged.res@edges[2,])
    edge.weights <- as.numeric(unlist(lapply(1:length(edgevalues), function(x) {return(0.75+ 15*(edgevalues[x] - min(edgevalues))/(max(edgevalues) - min(edgevalues)))})))
    #You can add font.color, label, and font.size to edges
    edges <- data.frame(from =unlist(lapply(merged.res@edges[1,],function(x) unlist(strsplit(x,split = "~"))[1])),
                        to =  unlist(lapply(merged.res@edges[1,],function(x) unlist(strsplit(x,split = "~"))[2])),
                        weight=edge.weights,value=edge.weights)
    
    print("end")
    return(list(nodes=nodes,edges=edges))
  }
  
}

making.graph <- function(ft,start.path,nodes.lbl, edges.from,edges.to,frame.list,gates.list,nodes.id,parent="CD45+")
{
  nxt<- tail(start.path,1)
  pop.path <- start.path
  while(nxt!=as.vector(edges.from[1]))
  {
    prnt <- which(edges.to==nxt)[1]
    nxt<- as.vector(edges.from[prnt])
    pop.path <-c( pop.path ,as.vector(nodes.id[which(nodes.id==nxt)]))
  }
  parent.names <- rev(nodes.lbl[match(pop.path,nodes.id)])
  parent.names <- parent.names[-length(parent.names)]
  pop.path <- rev(pop.path)[-1]
  
  pop.path <- lapply(pop.path,function(p) unlist(strsplit(p,"")))
 
  mat <-matrix(c(1:(ceiling(sqrt(length(pop.path)))*2)),nrow =ceiling(sqrt(length(pop.path))) ,byrow = T)
  if (length(frame)==2)
      mat <- cbind(mat, mat+mat[nrow(mat),ncol(mat)])
  layout(mat)

  for (frame in frame.list)
  {
    gates <- gates.list[[identifier(frame)]]
  for (i1 in 1:length(pop.path))
  {
    print(i1)
    pop.vec <- pop.path[[i1]]
    ind <-which(pop.vec!="0")
    if (names(ft@Thresholds)[ind]!="gate")
    {
      plotDens(frame, c(names(ft@Thresholds)[ind],"SSC-A"),main=parent.names[i1])
      abline(v=gates[[ind]],lwd=3)
      f.temp <- frame
      if (pop.vec[ind]=="1")
        exprs(f.temp) <- exprs(f.temp)[which(exprs(f.temp)[,names(ft@Thresholds)[ind]]<=gates[[ind]]),]
      if (pop.vec[ind]=="2")
        exprs(f.temp) <- exprs(f.temp)[which(exprs(f.temp)[,names(ft@Thresholds)[ind]]>gates[[ind]]),]
      frame <- f.temp
      if (i1<length(pop.path))
        for (k1 in (i1+1):length(pop.path))
          pop.path[[k1]][ind]<-"0"

    }else{
      plotDens(frame, colnames(gates[[ind]]),main=parent.names[i1],
               ylim=c(min(c(exprs(frame)[,colnames(gates[[ind]])[2]],
                            gates[[ind]][,2]  ),na.rm = T),
                      max(c(exprs(frame)[,colnames(gates[[ind]])[2]],
                            gates[[ind]][,2]),na.rm = T)),
               xlim=c(min(c(exprs(frame)[,colnames(gates[[ind]])[1]],
                            gates[[ind]][,1]  ),na.rm = T),
                      max(c(exprs(frame)[,colnames(gates[[ind]])[1]],
                            gates[[ind]][,1]),na.rm = T)))
      lines(gates[[ind]], lwd=3,col=1)
      if (pop.vec[ind]=="1")
        frame <- getflowFrame(notSubFrame(frame, colnames(gates[[ind]]),position = c(T,T),filter=gates[[ind]]))
      if (pop.vec[ind]=="2")
        frame <- getflowFrame(flowDensity(frame, colnames(gates[[ind]]),position = c(T,T),filter=gates[[ind]]))
      if (i1<length(pop.path))
        for (k1 in (i1+1):length(pop.path))
          pop.path[[k1]][ind]<-"0"
    }
  }
  }
}
######################################################
plot_results <- function(cell.proportion,outliers,showsample=FALSE,shiny=T,samplestoplotdata=NULL){
  
  
  if(!is.null(samplestoplotdata) & showsample){  #...and user-specified samples to plot, add multiple sets of geom_point_interactive.
    y.scale <- round(range(cell.proportion$props,na.rm = T)*c(.9,1.1),digits = 2)
    r <- ggplot() +
      geom_boxplot_interactive(data=cell.proportion,aes(x = groups, y = props, colour=groups,
                                                tooltip=groups),outlier.shape=NA,na.rm=T)+
      
      geom_point_interactive(position=position_jitter(width=0.05,height=.01),size=.7,data=outliers,
                             aes(x = groups, y = props,colour=groups,
                                 tooltip=paste0(round(props,digits=2),"%, ",files),data_id=files))+
      geom_point_interactive(position=position_jitter(width=0.1,height=.05),data=samplestoplotdata,size=0.9,
                             aes(x = gropus, y = props,color=groups,fill="black",
                                 tooltip=paste0(round(props,digits=2),"%, ",files),data_id=files))+
      scale_y_continuous(name="Proportion % ",limits = y.scale)+
      xlab(label = "Cell Population")+
      theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10),
            axis.text.y = element_text(hjust = 1, size=10,color=1),
            legend.position="none",plot.margin=unit(c(0.8,0,0,0.85),"cm"))
    if (!shiny){
      r <- girafe(ggobj=r)
    
    }
  }else{ 
    y.scale <- round(range(cell.proportion$props,na.rm = T)*c(.9,1.1),digits = 2)
    r <- ggplot() +
      geom_boxplot_interactive(data=cell.proportion,aes(x = groups, y = props, colour=groups,
                                                tooltip=groups),outlier.shape=NA,na.rm=T)+
      
      geom_point_interactive(position=position_jitter(width=0.05,height=.01),size=.7,data=outliers,
                             aes(x = groups, y = props,colour=groups,
                                 tooltip=paste0(round(props,digits=2),"%, ",files),data_id=files))+
      scale_y_continuous(name="Proportion % ",limits = y.scale)+
      xlab(label = "Cell Population")+
      theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10),
            axis.text.y = element_text(hjust = 1, size=10,color=1),
            legend.position="none",plot.margin=unit(c(0.8,0,0,0.85),"cm"))
    if (!shiny){
      r <- girafe(ggobj=r)
                  
    }
  }
  
  return(r)
}
###########################################
make.gs<- function(ft,ids, arguments,pop.path,parent="CD45+")
{
  pop.path <- lapply(pop.path,function(p) unlist(strsplit(p,"")))
  print(pop.path)
  fs <- lapply(ids, function(id){
    temp <- load(paste0(files.path,"/",id,".Rdata"));
    return(get(temp))})
  gs <- GatingSet(as(fs,"flowSet"))
  sampleNames(gs) <- ids
  fs <- gs_pop_get_data(gs)
  gates <-  arguments$threshold[ids]
  scat.chans <- unlist(sapply (c("SSC-A","FSC-A"),function(x) grep(x = colnames(fs),pattern = x,value = T)))[1]
  for (i1 in 1:length(pop.path))
  {
    
    pop.vec <- pop.path[[i1]]
    print(paste(i1,pop.vec))
    ind <-which(pop.vec!="0")
    print(ind)
    if (names(ft@Thresholds)[ind]!="gate")
    {
      
      if (pop.vec[ind]=="1")
      {
        fdens <- fsApply(fs, function(frame) return(flowDensity(frame, c(names(ft@Thresholds)[ind],scat.chans), 
                                                         position = c(F,NA),gates= c(gates[[identifier(frame)]][[ind]],NA))))
        filter.id <- paste0(arguments$names[ind],"-")
      }
      if (pop.vec[ind]=="2")
      {
        fdens <-  fsApply(fs, function(frame) flowDensity(frame, c(names(ft@Thresholds)[ind],scat.chans),
                                                          position = c(T,NA),gates= c(gates[[identifier(frame)]][[ind]],NA)))
        filter.id <- paste0(arguments$names[ind],"+")
        
      }
    }else{
      
      if (pop.vec[ind]=="1"){
        fdens <-  fsApply(fs, function(frame) notSubFrame(frame, colnames(gates[[identifier(frame)]][[ind]]),position = c(T,T),
                                                          filter=as.matrix(gates[[identifier(frame)]][[ind]])))
        filter.id <- paste0("Not ",arguments$names[ind])
        
      }
      if (pop.vec[ind]=="2"){
        fdens <-  fsApply(fs, function(frame) flowDensity(frame, colnames(gates[[identifier(frame)]][[ind]]),position = c(T,T),
                                                          filter=as.matrix(gates[[identifier(frame)]][[ind]])))
        filter.id <- arguments$names[ind]
        
      }
    }
    if (i1<length(pop.path))
      for (k1 in (i1+1):length(pop.path))
        pop.path[[k1]][ind]<-"0"
    poly <- lapply(fdens, function(x) polygonGate(filterId = filter.id,.gate  = x@filter))
    names(poly) <- sampleNames(gs)
    
    nodeID<-add(gs, poly,parent=tail(gs_get_pop_paths(gs),1))
    recompute(gs)
    
  }
  return(gs)
}
###########################################
plot.gates <- function(gs,popstoplot=NULL,upto=NULL,plot.layout=NULL,show.stats=TRUE,contour.pops=NULL,density.overlay.pops=NULL,
                      flowCut.channels=c("Time","BV510-A"),legend.cex=1.5,cell.size=c(1,2.5,3),show.palette=TRUE,offset=.5,
                      ...)
{
  #contour.pops is a vector of population to draw contour lines, it's the children in the plot
  #density.overlay.pops is a vector of population to draw, it's always the parent in the plot
  #cell.size is a vector of length 3, giving the 
  #Popstoplot is either a vector of integer defining which nodes to plot, or the node names itself in "auto" path.
  nodes <- gh_get_pop_paths(gs[[1]],path="auto")
  print(nodes)
  print(paste("Length of gs is: ",length(gs)))
  print(sampleNames(gs))
  if (!is.null(popstoplot))
  { 
    if (is.integer(popstoplot))
      nodes <-nodes[popstoplot]
    else 
      nodes<-popstoplot
  }
  nodes <- setdiff(nodes,"root")
 
    par(mar = c(6,8,3,1)+0.1)
    par(oma=c(0,0,4,0))
    m<-layout(plot.layout)
    for ( j1 in 1:length(gs))
    {
      plotList <- .mergeGates(gh = gs[[j1]],i = nodes,bool = T,merge = T,projections = list())
      to.remove <- which(names(plotList)==names(unlist(lapply(plotList, grep,pattern="_1"))))
      if (length(to.remove)>0)
        plotList <- plotList[-to.remove]
      if (!is.null(upto))
      {
        ind.stop <- unlist(lapply(1:length(plotList), function(k1){
          if (length(plotList[[k1]])==1)
          {
            if (plotList[[k1]]==upto)
              return(k1)
            
          }else{
            if (upto %in% plotList[[k1]]$popIds)
              return(k1)
          }
        }))
        plotList <- plotList[1:ind.stop]
      }
      print(paste("j1 is:", j1))
  tmp <- lapply(plotList,function(p.l) {setup.plot(y=p.l,gh=gs[[j1]],show.stats=show.stats,contour.pops=contour.pops,
                                                   density.overlay.pops=density.overlay.pops,flowCut.channels=flowCut.channels,
                                                   legend.cex=legend.cex,offset=offset,show.palette = show.palette,cell.size=cell.size,...)})

  }
}
