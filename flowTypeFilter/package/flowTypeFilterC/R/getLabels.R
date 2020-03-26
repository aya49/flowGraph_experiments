##Returns cluster membership labels of the given phenotype number.
getLabels <- function(Phenotypes, PhenotypeNumber){
	
  if (!is.numeric(PhenotypeNumber)){   
    PhenotypeNumber=which(Phenotypes@PhenoCodes==PhenotypeNumber);
  }
  
  Partitions=Phenotypes@Partitions;
  NumMarkers=length(Phenotypes@Thresholds);
  NumCells=length(Partitions[,1]);
  PartitionNumbers <- as.numeric(strsplit(Phenotypes@PhenoCodes[PhenotypeNumber],split='')[[1]])
  index <- as.vector(matrix(1,1,NumCells))
  
  
  calcMarkerMembership <- function(MarkerNum)
  {
  	#If this marker is not included in this phenotype definition, set its inclusion vector to all True:
  	if(PartitionNumbers[MarkerNum]==0)
  		return(rep(T, NumCells))
  	
  	#Otherwise return which cells belong to this partition for this marker:
  	return(Partitions[,MarkerNum] == PartitionNumbers[MarkerNum])
  }
  
  Partition.Membership <- lapply(1:NumMarkers, calcMarkerMembership)
  
  Pheno.Membership <- Partition.Membership[[1]]
  for(Membership in Partition.Membership[2:NumMarkers])
  {
  	Pheno.Membership <- Pheno.Membership& Membership
  }
  
  return(as.numeric(Pheno.Membership));
}
