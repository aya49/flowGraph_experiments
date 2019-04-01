##################################################################
# removes unwanted channels and reorders them to all be the same #
##################################################################
removeAndOrganizeChannels <- function(store.allFCS, ref=15, markers, new_dir, renameChannels=F){
  # Reference FCS
  f <- read.FCS(paste0(store.allFCS[ref,"Path"], "/", store.allFCS[ref,"FCS files"]))
  
  # Finding markers-------------------------------------------------------------------------------------------------
  channels.ind <- Find.markers(frame=f,marker.list=markers)
  names(channels.ind)[grep(names(channels.ind),pattern = "Live*")] <- "Live"

  edited_ind <- c()
  for (h1 in 1:nrow(store.allFCS)){
    f <- read.FCS(paste0(store.allFCS[h1,"Path"], "/", store.allFCS[h1,"FCS files"]))
    
    channels.ind.b4 <- suppressWarnings(Find.markers(frame=f,marker.list=markers)) # It is okay to redefine channel.ind. It will be redefined again below.
    names(channels.ind.b4)[grep(names(channels.ind),pattern = "Live*")] <- "Live"

    FSC.SSC.indices <- c( which(grepl("fsc", tolower(f@parameters@data$desc))), which(grepl("ssc", tolower(f@parameters@data$desc))) )
    channels.ind.b4 <- suppressWarnings(Find.markers(frame=f,marker.list=markers)) # It is okay to redefine channel.ind. It will be redefined again below.
    good.indices <- unique(sort(c(FSC.SSC.indices, channels.ind.b4)))
    remove.indices <- f@parameters@data$desc[setdiff(1:length(f@parameters@data$desc), good.indices)]
    
    if (length(remove.indices)>0) {
      edited_ind <- append(edited_ind, h1)
      f2 <- removeChannels(f=f, remove.indices, 1)
      channels.ind.af <- Find.markers(frame=f2,marker.list=markers)
      suppressWarnings(write.FCS(f2, file=paste0(new_dir, "/", store.allFCS[h1,"FCS files"])))
    }
    
    #if (renameChannels) f@parameters@data[,2] <- markers
    
    #for (h2 in channels.ind.af[-length(channels.ind.af)]){fz[[h1]]@description[paste0("P", h2, "DISPLAY")] <- "LOG"}
    #for (h2 in channels.ind.af[ length(channels.ind.af)]){fz[[h1]]@description[paste0("P", h2, "DISPLAY")] <- NULL}
  }
  return(edited_ind)
}
##################################################################
removeChannels <- function(f, extra.markers, index.FCS.originalMarkers){
  
  ## Getting rid of the extra channels (exprs, description, parameter) and re-organizing the channels according to the original ones
  if(any(f@parameters@data$desc %in% extra.markers == TRUE)){
    index.extraMarkers <- which(f@parameters@data$desc %in% extra.markers)
    index.SPILL.extraMarker <- which(colnames(f@description$SPILL) %in% extra.markers)
    
    # exprs matrix
    ex <- exprs(f)[,-index.extraMarkers]
    # Spill over matrices
    description(f)$SPILL <- description(f)$SPILL[-index.SPILL.extraMarker,-index.SPILL.extraMarker]
    
    # keywords
    for ( x in index.extraMarkers ) {
      if ( x == length(f@parameters@data$desc)) {break} # >=????
      for ( y in x:(length(f@parameters@data$desc) - which(x==index.extraMarkers) + 1) ) {
        index.Desc <- grep(paste0("P",y), names(f@description))
        for(z in 1:length(index.Desc)){
          names(f@description)[index.Desc[z]] <- gsub(paste0("P",y), paste0("P",y-which(x==index.extraMarkers)), names(f@description)[index.Desc[z]])
        }	}	}
    
    # table of data for the FCS file
    o <- parameters(f)@data[-index.extraMarkers,]
    rownames(o) <- paste0("$P",1:length(rownames(o)))
    metaData <- varMetadata(parameters(f))
    # number of markers ( need to change for flowClean )
    f@description$'$PAR' <- as.character(as.numeric(f@description$'$PAR') - length(index.extraMarkers))
    # create the new flowFrame
    f <- new("flowFrame", exprs=ex, parameters=new("AnnotatedDataFrame", o, varMetadata=metaData), description=description(f)[order(names(f@description))])
  }
  return(f)
}

