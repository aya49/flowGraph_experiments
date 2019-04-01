# Written by Justin Meskas Dec 2016
# This code will clean all the files listed in the store.allFCS.Rdata file in the results.dir directory.
# The marker names are required as input.
# If overwrite.channel.names is defined as a location of a file, then the channel names will be overwritten.
# Cleaned flowFrames and indices that are removed will be saved independently.
# A matrix of the results of TvsF is also saved, and is the only thing returned to the user.


SangerCleanFunc <- function(results.dir, marker.names, overwrite.channel.names = NULL){

    libr("flowCore")
    libr("flowDensity")
    libr("foreach")
    libr("parody")
    libr("e1071")
    libr("Cairo")
  libr("doMC")

    load(paste0(results.dir, "/Genotype.Rdata"))
    load(paste0(results.dir, "/uniqueGT.Rdata"))
    load(paste0(results.dir, "/store.allFCS.Rdata"))
    # if (!exists("globalFrame") ) {load(paste0(results.dir, "/globalFrame.Rdata"))} # avoid loading it all the time suring testing (takes 15 seconds)
    load(paste0(results.dir, "/globalFrame.Rdata"))

    if(!is.null(overwrite.channel.names))
        globalFrame@parameters@data$desc <- overwrite.channel.names

    # Removing (max) margin and negative events in scatter channels
    scat.chans <- c(grep (colnames(globalFrame),pattern = "FSC*"),grep (colnames(globalFrame),pattern = "SSC*"))
    res.remove <- removeMargNegs(globalFrame, chans=scat.chans, neg=0, verbose= F)
    globalFrame <- res.remove$f; ind.marg.neg.gf <- res.remove$ind.marg.neg; remove(res.remove)
    globalFrame@exprs <- globalFrame@exprs[-ind.marg.neg.gf,]

    # Channels
    channels.ind <- Find.markers(frame=globalFrame,marker.list=marker.names)
    #names(channels.ind)[grep(names(channels.ind),pattern = "TCRd*")] <- "TCRd"
    names(channels.ind)[grep(names(channels.ind),pattern = "Live*")] <- "Live"
    channels.ind <- sort(channels.ind)

    ## Tranformation
    lgl<-estimateLogicle(globalFrame, channels= colnames(globalFrame)[channels.ind])
    save(channels.ind, file = paste(results.dir, "/channels.ind.Rdata", sep = ""))
    save(lgl,          file = paste(results.dir, "/lgl.Rdata",sep=""))

    # Create folders
    suppressWarnings( dir.create(paste0(results.dir, "/Figures/Clean"),   recursive = T) )

    
    no_cores <- detectCores() - 1
    #no_cores <- detectCores() - 3
    registerDoMC(no_cores)
    start <- Sys.time()

    loop.ind <- 1:nrow(store.allFCS)
    #loop.ind <- c(1:5)
    results <- foreach(i = loop.ind, .combine = list, .maxcombine = nrow(store.allFCS), .multicombine = T) %dopar% {
        res <- tryCatch({
            cat(paste0(i, ": Starting ", store.allFCS[i,"Genotype"],"/",store.allFCS[i,"Barcodes"], ".\n"))
            f <- read.FCS(paste0(store.allFCS[i,"Path"], "/", store.allFCS[i,"FCS files"]))

            if(!is.null(overwrite.channel.names))
                f@parameters@data$desc <- overwrite.channel.names

            # Removing (max) margin and negative events in scatter channels
            scat.chans <- c(grep (colnames(f),pattern = "FSC*"),grep (colnames(f),pattern = "SSC*"))
            res.remove <- removeMargNegs(f, chans=scat.chans, neg=0, verbose= F)
            f <- res.remove$f; ind.marg.neg <- res.remove$ind.marg.neg; remove(res.remove)

            ## Compensation
            cat("Starting Compensation.\n")
            if( det(f@description$SPILL)==1 ){
                message("Check the spillover matrix, it's probably an identity matrix!")
            } else {
                f <- compensate(f, f@description$SPILL)
            }

            f <- transform(f, lgl)

            # directory where TvsF will save its plots
            save.dir_TvsF <- paste0(results.dir, "/Figures/Clean/", store.allFCS[i,"Genotype"])

            ## Cleaning
            res.tvsf <- timeVsFluorescence(f, segment=500, directory=save.dir_TvsF, FileID = store.allFCS[i,"Barcodes"], Plot="Failed Only",
                                           MaxContin=0.2, MeanOfMeans=0.33, MaxOfMeans=0.5, verbose=T)
            f <- res.tvsf$frame
            ind.clean <- res.tvsf$ind;

            return(list(data=res.tvsf$data, ind=list(ind.marg.neg=ind.marg.neg, ind.clean=ind.clean) ))

        },error = function(e) {
            print(e)
            error.data <- rbind(e, store.allFCS[i,"Barcodes"])
            rownames(error.data) <- c('error message', "FileID")

            return(list(data=error.data))
        }) # end of tryCatch

    }

    # figure out which files failed the cleaning step, remove from store.allFCS and results
    # Then save barcodes of FCS files which failed cleaning to .Rdata file
    idx.failed.cleaning <- NULL
    for(i in loop.ind){
        if(row.names(results[[i]]$data)[1] == "error message"){
            #print(results[[i]]$data)
            idx.failed.cleaning <- c(idx.failed.cleaning, i)
        }
    }
    if(length(idx.failed.cleaning) > 0){
        results <- results[-idx.failed.cleaning]
        save(store.allFCS, file =  paste0(results.dir, "store.allFCS.preCleaning.Rdata")) # in case you want to revert it
        FCS.files.failed.cleaning <- store.allFCS[idx.failed.cleaning, 'FCS files']
        save(FCS.files.failed.cleaning, file =  paste0(results.dir, "failedCleaning.FCS.Rdata"))
        store.allFCS <- store.allFCS[-idx.failed.cleaning,]
        save(store.allFCS, file =  paste0(results.dir, "store.allFCS_errors_removed.Rdata"))
    }

    # save TvsF results
    if (length(loop.ind) != 1){
        res.clean <- cbind(sapply(1:length(results), function(x) { results[[x]]$data }))
        row.names(res.clean) <- row.names(results[[1]]$data)
    } else {
        res.clean <- results$data
        row.names(res.clean) <- row.names(results$data)
    }
    colnames(res.clean) <- store.allFCS[loop.ind,"Barcodes"]
    save(res.clean, file = paste0(results.dir, "/res.clean.Rdata"))

    # save indices removed from margins, negatives and TvsF
    if (length(loop.ind) != 1){
        ind.marg.neg.clean.all <- lapply(1:length(results), function(x) { results[[x]]$ind })
    } else {
        ind.marg.neg.clean.all <- results$ind
    }
    save(ind.marg.neg.clean.all, file = paste0(results.dir, "/ind.marg.neg.clean.all.Rdata"))

    cat("Time is: ", TimeOutput(start),"\n",sep="")
    return(res.clean)
}
