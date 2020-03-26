#Calculate how many populations you can have from total.markers, only using up to max.markers at a time:

calcMemUse <- function(NumPops, NumPropMarkers, NumMFIMarkers, NumCells, MaxMarkersPerPop, PartitionsPerChannel=2)
{
	double.size <- 16
	int.size <- 8
	#Then calculate memory usage in bytes:
	partitions.size <- NumPropMarkers * NumCells * int.size #size to store single channel per-cell partition information
	mfi.size <- NumPops * NumMFIMarkers * double.size #size to store MFI results
	count.size <- NumPops * int.size
	code.size <- NumPops * ceiling(NumPropMarkers/8)*8 + 88 #size for storing the 1-byte per channel code for each pop
	
	#Also calculate size of intermediate storage used while executing:
	calcLevel <- function(level){choose(NumPropMarkers, level) * PartitionsPerChannel^level}
	levels <- sapply(1:MaxMarkersPerPop, calcLevel)
	level.pops <- max(levels)
	if(which.max(levels)==1)
	{
		level.minus.pops <- 0
	}else
	{
		level.minus.pops <- levels[which.max(levels)-1]	
	}
	
	mfi.index.size <- NumCells * NumPropMarkers * 8 #Size for creating sorted mfi data
	temp.size <- NumCells/8 * (level.pops + level.minus.pops)
	total.size <- mfi.size + count.size + code.size + temp.size + partitions.size + mfi.index.size
	total.size
}

