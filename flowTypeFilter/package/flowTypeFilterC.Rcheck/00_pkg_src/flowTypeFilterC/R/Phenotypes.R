setClass("Phenotypes",
         slots = c(CellFreqs="numeric", 
                   MFIs="matrix",
                   PhenoCodes="vector",
                   PropMarkers="vector", 
                   MFIMarkers="vector", 
                   MarkerNames="vector", 
                   Partitions="matrix",
                   MaxPopSize="numeric",
                   PartitionsPerMarker="numeric",
                   Thresholds="list"))

setMethod("summary", signature(object="Phenotypes"), function(object){
  cat(sprintf("Phenotypes object identified by flowTypeFilter with %d cell-frequency-based on %d MFI-based features.\n", length(object@CellFreqs), dim(object@MFIs)[1]*dim(object@MFIs)[2]));
})

setMethod("plot", signature(x="Phenotypes", y="flowFrame"), function(x, y, ...){
  RawData=y
  object=x
  X=exprs(RawData)[,object@PropMarkers]
  if (all (unlist(lapply(object@Thresholds,class))=="matrix"))
     stop("Cannot plot densities as all thresholds are polygone filtes.")
  par(mfrow=c(ceiling(sqrt(length(object@PropMarkers))),ceiling(sqrt(length(object@PropMarkers)))))
  for (i in 1:ncol(X)){
    plot(density(X[,i]), main='', xlab='', ylab='',axes=FALSE, lwd=2, ...);
    axis(1);
    axis(2);
    title(main=object@MarkerNames[object@PropMarkers[i]], ylab='Density', cex=2,cex.lab=1.5);
        abline(v=object@Thresholds[[colnames(X)[i]]],col=2, lwd=2);
  }
})

setMethod("plot", signature(x="Phenotypes", y="numeric"), function(x, y, Frame, ...){
  if (any (unlist(lapply(object@Thresholds,class))!="matrix"))
     plot(data.frame(exprs(Frame)[,x@PropMarkers]), col=getLabels(x, y), ...)
})

setMethod("plot", signature(x="Phenotypes", y="character"), function(x, y, Frame, ...){
       if (any (unlist(lapply(object@Thresholds,class))!="matrix"))
        {
	  phenoCode <- encodePhenotype(pheno.string=y, marker.names=x@MarkerNames[x@PropMarkers])
	  plot(data.frame(exprs(Frame)[,x@PropMarkers]), col=getLabels(x, phenoCode), ...)
         }
})

setMethod("plot", signature(x="Phenotypes", y="flowFrame"), function(x, y,type="density",channel="SSC-A", ...){
  RawData=y
  object=x
  if(type=="density")
  {
  	X=exprs(RawData)[,object@PropMarkers]
 	 if (all (unlist(lapply(object@Thresholds,class))=="matrix"))
     	   stop("Cannot plot densities as all thresholds are polygone filtes.")
	 par(mfrow=c(ceiling(sqrt(length(object@PropMarkers))),ceiling(sqrt(length(object@PropMarkers)))))
  	for (i in 1:ncol(X)){
    		plot(density(X[,i]), main='', xlab='', ylab='',axes=FALSE, lwd=2, ...);
    		axis(1);
    		axis(2);
    		title(main=object@MarkerNames[object@PropMarkers[i]], ylab='Density', cex=2,cex.lab=1.5);
        	abline(v=object@Thresholds[[colnames(X)[i]]],col=2, lwd=2);
       }
  }else if (type=="dotplot")
  {  
      par(mfrow=c(ceiling(sqrt(length(object@Thresholds))),ceiling(sqrt(length(object@Thresholds)))))
     for (i in 1:length(object@Thresholds)){  
          if(class(object@Thresholds[[i]])=="matrix")
              chs <- colnames(object@Thresholds[[i]])
         else
              chs <- c(names(object@Thresholds)[i], channel)
	  colPalette <- colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))
          col <- densCols(exprs(y)[,chs], colramp = colPalette)
        plot(exprs(y[,chs]),pch=".", col = col, ...)
        if(class(object@Thresholds[[i]])=="matrix")
           lines(object@Thresholds[[i]],lty=3, lwd=2)
        else
          abline(v= as.vector(object@Thresholds[[i]]),lwd=2)

  
     }
  }

})
