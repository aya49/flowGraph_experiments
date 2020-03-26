calcNumPops <- function(PartitionsPerMarker, MaxMarkersPerPop)
{
  TotalMarkers <- length(PartitionsPerMarker)
  #PartitionsPerMarker <- PartitionsPerMarker -1
  
  #If the base case, don't do DP
  if(MaxMarkersPerPop ==1)
  	return(TotalMarkers)
  
  #Matrix for storing results, to aid dynamic programming:
  CountsMatrix <- matrix(rep(0, TotalMarkers*MaxMarkersPerPop), nrow=TotalMarkers, ncol=MaxMarkersPerPop)
  
  #Calculate first column of matrix:
  calcForOneMarker <- function (TotalMarkersUsed)
  {
    sum(PartitionsPerMarker[1:TotalMarkersUsed])
  }
  
  CountsMatrix[,1] <- sapply(1:TotalMarkers, calcForOneMarker)
  
  #Calculate the rest of the matrix (excluding cells where m<k)
  calcForKMarkers <- function(TotalMarkersUsed, MaxMarkersUsed)
  {
    lhs <- CountsMatrix[TotalMarkersUsed-1, MaxMarkersUsed]
    #print(lhs)
    rhs <- CountsMatrix[TotalMarkersUsed-1, MaxMarkersUsed-1]*PartitionsPerMarker[TotalMarkersUsed]
    #print(rhs)
    lhs + rhs
  }
  
  #Now calculate counts for all populations 
  for (k in 2:MaxMarkersPerPop)
  {
    for (m in 1:TotalMarkers)
    {
      if(m>k)
        CountsMatrix[m,k] <- calcForKMarkers(m,k)   
      if(m==k)
        CountsMatrix[m,k] <- Reduce(function(x,y){x*y}, PartitionsPerMarker[1:m])
    }
  }
  
  #print(CountsMatrix)
  
  #CountsMatrix is now populated, so just sum the last row to get total populations:
  return (sum(CountsMatrix[TotalMarkers,]) + 1)
  
}
