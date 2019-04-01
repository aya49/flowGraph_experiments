############################################################################################
# removes all cells that exist on the edges of FSC-A, FSC-H, SSC-A, and SSC-H. 
removeMargins <-function(f.temp){
    

for ( q1 in 1:length(colnames(f.temp))){ # 6 because there may be 3 FSC and 3 SSC columns 
        if (length(which(f.temp@exprs[,q1]==262143))>0){f.temp <- f.temp[-which(f.temp@exprs[,q1]==262143)]}
        if (length(which(f.temp@exprs[,q1]==0))>0)     {f.temp <- f.temp[-which(f.temp@exprs[,q1]==0)]}
}
return(f.temp)

}

############################################################################################
# removes all cells that exist on the edges of FSC-A, FSC-H, SSC-A, and SSC-H. 
removeMarginsAndNegatives <-function(f.temp){
    

for ( q1 in 1:length(colnames(f.temp))){ # 6 because there may be 3 FSC and 3 SSC columns 
        if (length(which(f.temp@exprs[,q1]==262143))>0){f.temp <- f.temp[-which(f.temp@exprs[,q1]==262143)]}
        if (length(which(f.temp@exprs[,q1]<=0))>0)     {f.temp <- f.temp[-which(f.temp@exprs[,q1]<=0)]}
}
return(f.temp)

}


############################################################################################
#Remove all cell types on the outside of an elipse from the FSC-A and FSC-H plot.
removeDoublets <-function(f.temp, temp.nameU, q){
if (length(which(f.temp@parameters@data$name=="FSC-H"))>=1) {   
  png ( file = paste("Plots/Figures/",temp.nameU, "_", q, "_B4DoubMarg.png", sep=""))
    plotDens(f.temp, c(1,2))  
  dev.off()
  
  singlet <- flowDensity(f.temp,c(which(f.temp@parameters@data$name=="FSC-A"),which(f.temp@parameters@data$name=="FSC-H")),
                                position=c(T,F),percentile=c(.01,.99),use.percentile=c(T,T),ellip.gate=T)
  f.temp <- getflowFrame(singlet)
  
  png ( file = paste("Plots/Figures/",temp.nameU, "_", q, "_AftDoubMarg.png", sep=""))
    plotDens(f.temp, c(1,2))  
  dev.off()
}
return(f.temp)
}



############################################################################################
# Remove all cell types on the outside of an elipse from the FSC-A and FSC-H plot.
removeDoubletsSSC <-function(f.temp, name1, name2, directory, Plt){
   
if (length(which(f.temp@parameters@data$name=="SSC-H"))>=1) {   
  singlet <- flowDensity(f.temp,c(which(f.temp@parameters@data$name=="SSC-A"),which(f.temp@parameters@data$name=="SSC-H")),
                                position=c(T,F),percentile=c(0.01,.99),use.percentile=c(T,T),ellip.gate=T,scale = 0.999999)

  f.temp <- getflowFrame(singlet)
}
return(f.temp)
}

#################################################################################
# rearranges the order of the flowFrame
reorderfSet <- function(f,commonMarker,columnPlacement){

NumOfChannels <- length(f@parameters@data$desc)    
temp.number <- which(f@parameters@data$desc==commonMarker)
if ( length(temp.number)>1){ temp.number <- temp.number[length(temp.number)]}
row.order <- NULL
for(q in 1:NumOfChannels) {row.order <- c(row.order,q)}
row.order[temp.number]    <- columnPlacement
row.order[columnPlacement]<- temp.number
f <- f[,row.order]
return(f)
}


#################################################################################
# Prints out the time since start_time. Used for optimizing code.
TimeOutput <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}
TimeOutput(Sys.Date())


#################################################################################
#SVD reduction for flowType to be used for grouping patients based on reduced matrix
#Kmeans used for grouping patients, any other method can be used

#Author: Mehrnoush Malek
#Date: April 2012

svd.reduction <- function(cell.prop,kmean.centres=3,kmeans.start=1000)
  #cell.prop is a matrix of size 3^k by m, where k is number of markers used in flowType and m is number of samples
{
  
svdr.Result<-svd(t(cell.prop))
#Find a threshold where samples get far from others
x11();plot(sort(svdr.Result$d))

inds<-which(svdr.Result$d > locator()$y)
acf<-t(cell.prop) %*% (svdr.Result$v[,inds])
kmeans.cluster<-kmeans(acf, centers=kmean.centres, nstart=kmeans.start)
x11();plot(data.frame(acf), pch=16,col=kmeans.cluster$cluster)
return(kmeans.cluster)
}

# length(flowType.res@CellFreqs)
# which(names(flowType.res@CellFreqs)=="")
# flowType.res@CellFreqs[1093]
# flowType.res@CellFreqs[1:2]
# ?heatmap
# d<-plot(flowType.res,flowFrame_fB)
# d
