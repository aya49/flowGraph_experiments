## Developed by Albina Rahim
## Date: November 08, 2018
## This is the main script for gating the CyTOF datasets for the Immnune Clock of Human Pregnancy

remove(list=ls())

setwd("/code/Projects/Immune_Clock_Human_Pregnancy/Codes/")


###############################################################################################
## NOTE: Run this section functions either separately before you run the rest of the script or you can run the entire script in just one go
## The helperFunc.R prompts user for input: Training data or Validation data


source("helperFunc.R")
if (interactive()){
  dataType <- readDataFunc()
  dataInfo <- NULL
  
  if(dataType == "1"){
    dataInfo <- "Training"
  }else if(dataType == "2"){
    dataInfo <- "Validation"
  }
  
}



#source("flowCut.R")


#############################################################################################
##############################################################################################

library("plyr")
library("doMC")
library('colorRamps')
# library('e1071') # for flowCut
# library('Cairo') # for flowCut
library('flowCore')
library('flowDensity')
library('pracma') # for findpeaks
library('tools')
library('MASS')
library('KernSmooth')
library('stringr')
library('flowViz')

no_cores <- detectCores() - 2
registerDoMC(no_cores)



start <- Sys.time()

## Paths to the FCS files, metadata spreadsheets, and output folder in the Bioinformatics drive
if(dataInfo == "Training"){
  outputPath <- "/mnt/f/FCS data/Immune_Clock_Human_Pregnancy/Results/FR-FCM-ZY3Q/"
}else if(dataInfo == "Validation"){
   outputPath <- "/mnt/f/FCS data/Immune_Clock_Human_Pregnancy/Results/FR-FCM-ZY3R/"
}


load(paste0(outputPath,"/store.allFCS.Rdata"))


# Create directories
suppressWarnings(dir.create(paste0(outputPath,"/Figures")))
suppressWarnings(dir.create(paste0(outputPath,"/Figures/Scatter-Plots/")))
suppressWarnings(dir.create(paste0(outputPath,"/Figures/flowCut-Plots/")))
suppressWarnings(dir.create(paste0(outputPath,"/LeukocytesRdata")))

## We will be focusing on Unstim files only for now (Refer back to email from Nima dated January 03, 2019)
index <- which(store.allFCS[,c('Cell Stimulation')] == "Unstim")
store.allFCS.Unstim <- store.allFCS[index,]
save(store.allFCS.Unstim, file = paste0(outputPath,"/store.allFCS.Unstim.Rdata"))

file.names <- data.frame(store.allFCS.Unstim, stringsAsFactors = F)

#file.names <- data.frame(store.allFCS, stringsAsFactors = F)

print("Starting Gating & Plotting")

#props.events.gates <- llply(15:60, function(i){ 
props.events.gates <- llply(1:nrow(file.names), function(i){ 
  
  x<- file.names[i,]
  
  all.events <- matrix(nrow = 1, ncol = 30, data = NA)# matrix for saving the event counts
  all.props <- matrix(nrow = 1, ncol = 30, data = NA)# matrix for saving the proportions
  all.gthres <- matrix(nrow = 1, ncol = 37, data = NA) # matrix for saving the gating thresholds
  all.freqFeatures <- matrix(nrow = 1, ncol = 25, data = NA) ## matrix for saving proportions based on Mononuclear cells
  
  # all.fClean <-  matrix(nrow = 1, ncol = 1, data = NA) # matrix for saving the outputs from flowCut
  
  ## Saving the flowDensity filter for flowType
  Filters.list <- list("singlets.filter" = 0, "mononuclear.filter" = 0, "granulocytes.filter" = 0,
                       "naiveTregs.filter" = 0, "memoryTregs.filter" = 0)
  
  tryCatch({
    
      f <- try(read.FCS(filename = paste0(x$Path, "/", x$FCS.files)), silent = TRUE)
     
    
      # pregating_channels <- c("Bead", "DNA1", "DNA2", "Dead", "Event_length")
      gating_channels <- c("CD57", "CD19", "CD4", "CD8", "IgD", "CD11c", "CD16", "CD3", "CD38", "CD27", "CD14", "CXCR5", "CCR7", "CD45RA", "CD20", "CD127", "CD33", "CD28", "CD161", "TCRgd", "CD123", "CD56", "HLADR", "CD25", 
                           "CD235ab_CD61", "CD66", "CD45", "Tbet", "CD7", "FoxP3", "CD11b")
      # instrument_channels <- c("Time", "Event_length", "Center", "Offset", "Width", "Residual", "sample")
      
      channels.ind <- sort(Find.markers(f, gating_channels))
   
      if(is.na(channels.ind["Time"])){  
        channels.ind["Time"] <- grep('Time', pData(parameters(f))$name)
      }
    
      pregating.channels <- sort(c(grep(colnames(f),pattern = "Bead|bead*"), grep(f@parameters$desc, pattern = "DNA*"), grep(colnames(f), pattern = "Dead"), grep(colnames(f), pattern = "Event_length")))
      names(pregating.channels) <- colnames(f)[pregating.channels]
    
      ## NOTE: Not removing any margins or negative values since these are CyTOF files and we want to match the results of the manual analysis
      # # Remove margins of pregating channels and compensate ---------------------------------------------
      # f <- removeMargins(f, chans = pregating.channels, verbose = F)
      # #Removing negative values in pregating channels 
      # f <- removeMargins(f, chans = pregating.channels, debris = T, neg = T, verbose = F)
      # 
      # if(!is.null(f@description$SPILL)){
      #   f <- compensate(f, f@description$SPILL)
      # }
      
      ## Running flowCut just to generate plots and see how the markers vs time plots look. If there are any clogs.
      ## fT will be used for running flowCut
      # fT <- transform(f, lgl)
      # 
      # channels.to.clean <- c(pregating.channels[c("Ir191Di", "Ir193Di")], channels.ind)
      # channels.to.clean <- setdiff(channels.to.clean, channels.to.clean['Time'])
      # 
      # 
      # f.Clean <- flowCut(fT, Channels = channels.to.clean, Directory = paste0(outputPath,"Figures/flowCut-Plots/"), Plot = "All")
   
      channels.to.transform <- c(pregating.channels[c("Ir191Di", "Ir193Di")], channels.ind)
      channels.to.transform <- setdiff(channels.to.transform, channels.to.transform['Time'])
      
      ## Arcsinhtransformation on CyTOF data
      asinhTrans <- arcsinhTransform(transformationId="ln-transformation", a=1, b=1, c=1)
      translist <- transformList(colnames(f)[channels.to.transform], asinhTrans) 
     
      f <- transform(f, translist)
    
      #f <- transform(f, lgl)
      
      ## Cell counts/proportions after cleaning and margin event removal
      all.events[1] <- nrow(f)
      all.props[1] <- (nrow(f)/nrow(f))*100
    
    
    ###############################################################################################################
    ###############################################################################################################
    ## PRE-GATING
    
    ## Singlets from All Events---------------------------------------
    
    #dna.gate.low.allEvents <-  deGate(f, channel = pregating.channels["Ir191Di"], use.upper = T, upper = F, tinypeak.removal = 0.9)
    dna.gate.low.allEvents <- deGate(f, channel = pregating.channels["Ir191Di"], use.percentile = T, percentile = 0.015)
    dna.gate.high <- deGate(f, channel = pregating.channels["Ir191Di"], use.percentile = T, percentile = 0.999)
    #eventLength.gate <- deGate(f, channel = pregating.channels["Event_length"], use.percentile = T, percentile = 0.98)
    
    singlets.flowD.temp <- flowDensity(f, channels = c('Event_length', 'Ir191Di'), position = c(NA,T),  gates = c(NA, dna.gate.low.allEvents))
    singlets.flowD <- flowDensity(singlets.flowD.temp, channels = c('Event_length', 'Ir191Di'), position = c(NA,F),  gates = c(NA, dna.gate.high))

    singlets <- getflowFrame(singlets.flowD)
    
  
    all.events[2] <- nrow(singlets)
    all.props[2] <- (nrow(singlets)/nrow(f))*100

    ## Saving the gating thresholds
    all.gthres[1] <- dna.gate.low.allEvents
    all.gthres[2] <- dna.gate.high
    
   
    #plot(f, singlets.flowD)
    
    
    plotDens(f, c(pregating.channels["Event_length"], pregating.channels["Ir191Di"]), main = "Singlets", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(singlets.flowD@filter, lwd = 2)
    
    
    ######################################################################################################

    ## Leukocytes from Singlets population-------------------------------
    
    dna.gate.low.Singlets <- deGate(singlets, channel = pregating.channels["Ir191Di"], use.percentile = T, percentile = 0.1)-0.75
    if(round(dna.gate.low.Singlets,2) <= 4.6){
      dna.gate.low.Singlets <- deGate(singlets, channel = pregating.channels["Ir191Di"])
      if(round(dna.gate.low.Singlets) >= 6){
        dna.gate.low.Singlets <- deGate(singlets, channel = pregating.channels["Ir191Di"], use.percentile = T, percentile = 0.1)-0.75
      }
    }
    
    numPeaks.cd235ab_cd61 <- getPeaks(singlets, channel = channels.ind["CD235ab_CD61"])
    if(length(numPeaks.cd235ab_cd61$Peaks) > 1){
      if(length(deGate(singlets, channel =  channels.ind["CD235ab_CD61"], all.cuts = T)) > 2){
        cd235ab_cd61.gate <- max(deGate(singlets, channel =  channels.ind["CD235ab_CD61"], all.cuts = T))
      }else{
        cd235ab_cd61.gate <- deGate(singlets, channel = channels.ind["CD235ab_CD61"], use.percentile = T, percentile = 0.95)
        if(cd235ab_cd61.gate > 3.75) {
          cd235ab_cd61.gate <- deGate(singlets, channel = channels.ind["CD235ab_CD61"], use.percentile = T, percentile = 0.9)
          if(cd235ab_cd61.gate > 3.5){
            cd235ab_cd61.gate <- deGate(singlets, channel = channels.ind["CD235ab_CD61"], use.upper = T, upper = T)
          }
        }
      }
    }else{
      cd235ab_cd61.gate <- deGate(singlets, channel = channels.ind["CD235ab_CD61"], use.percentile = T, percentile = 0.95)
      if(cd235ab_cd61.gate > 3.75) {
        cd235ab_cd61.gate <- deGate(singlets, channel = channels.ind["CD235ab_CD61"], use.percentile = T, percentile = 0.9)
        if(cd235ab_cd61.gate > 3.5){
          cd235ab_cd61.gate <- deGate(singlets, channel = channels.ind["CD235ab_CD61"], use.upper = T, upper = T)
        }
      }
    }
    
    #cd235ab_cd61.gate <- deGate(singlets, channel = channels.ind["CD235ab_CD61"], use.upper = T, upper = T, use.percentile = T, percentile = 0.95)
      
    leukocytes.flowD.temp <- flowDensity(singlets, channels = c('In113Di', 'Ir191Di'), position = c(NA, T), gates = c(NA, dna.gate.low.Singlets))
    leukocytes.flowD <- flowDensity(leukocytes.flowD.temp, channels = c('In113Di', 'Ir191Di'), position = c(F, F), gates = c(cd235ab_cd61.gate, dna.gate.high))

    leukocytes <- getflowFrame(leukocytes.flowD)

    all.events[3] <- nrow(leukocytes)
    all.props[3] <- (nrow(leukocytes)/nrow(singlets))*100

    ## Saving the gating thresholds
    all.gthres[3] <- dna.gate.low.Singlets
    all.gthres[4] <- cd235ab_cd61.gate
    # ## Saving Leukocytes flowFrame for running flowType later
    # save(leukocytes, file = paste0(outputPath, "/LeukocytesRdata/", x$FCS.files, ".Rdata"))


    min.y <-  min(exprs(singlets)[, c('Ir191Di')])-2
    max.y <- max(exprs(singlets)[, c('Ir191Di')])+5
    min.x <- min(exprs(singlets)[, c('In113Di')])
    max.x <- max(exprs(singlets)[, c('In113Di')])+2
    plotDens(singlets, c(channels.ind["CD235ab_CD61"],pregating.channels["Ir191Di"]), main = "Leukocytes", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(leukocytes.flowD@filter, lwd = 2)

    ####################################################################################################################
    ####################################################################################################################
    ## MAIN GATING

    ## Mononuclear cells & Granulocytes from Leukocytes-----------------------
    if(length(deGate(leukocytes, channel = channels.ind["CD66"], all.cuts = T)) > 1){
      cd66.gate.low <- min(deGate(leukocytes, channel = channels.ind["CD66"], all.cuts = T))
    }else{
      cd66.gate.low <- deGate(leukocytes, channel = channels.ind["CD66"], use.upper = T, upper = F)
    }
    # cd66.gate <- deGate(leukocytes, channel = channels.ind["CD66"])
    # if(cd66.gate.low == cd66.gate){
    #   cd66.gate <- max(deGate(leukocytes, channel = channels.ind["CD66"], all.cuts = T))
    # }
    # 
    # cd45.gate.low <- deGate(leukocytes, channel = channels.ind["CD45"], use.upper = T, upper = F)
    # if(cd45.gate.low < 2){
    #   cd45.gate.low <- deGate(leukocytes, channel = channels.ind["CD45"], use.percentile = T, percentile = 0.01)
    # }
    # 
    # cd45.gate <- deGate(leukocytes, channel = channels.ind["CD45"])
    # if(cd45.gate < 3) {
    #   cd45.gate <- max(deGate(leukocytes, channel = channels.ind["CD45"], all.cuts = T))
    # }
    # 
    # mononuclear.flowD.temp <- flowDensity(leukocytes, channels = c('La139Di','In115Di'), position = c(F,T), gates = c(cd66.gate, cd45.gate-0.5))
    # mononuclear.flowD <- flowDensity(mononuclear.flowD.temp, channels = c('La139Di','In115Di'), position = c(T,NA), gates = c(cd66.gate.low, NA), ellip.gate = T)
    # 
    # #mononuclear.flowD <- flowDensity(leukocytes, channels = c('La139Di','In115Di'), position = c(F,T), gates = c(cd66.gate-0.05, cd45.gate-0.5))
    # 
    # granulocytes.flowD.temp <- flowDensity(leukocytes, channels = c('La139Di','In115Di'), position = c(T,F), gates = c(cd66.gate-0.4, cd45.gate+0.4), ellip.gate = T)
    # granulocytes.flowD <- flowDensity(granulocytes.flowD.temp, channels = c('La139Di','In115Di'), position = c(NA,T), gates = c(NA, cd45.gate.low-0.75), ellip.gate = F)
    # 
    # mononuclear <- getflowFrame(mononuclear.flowD)
    # granulocytes <- getflowFrame(granulocytes.flowD)
    # 
    # all.events[4] <- nrow(mononuclear)
    # all.props[4] <- (nrow(mononuclear)/nrow(leukocytes))*100
    # 
    # all.events[5] <- nrow(granulocytes)
    # all.props[5] <- (nrow(granulocytes)/nrow(leukocytes))*100
    # 
    # ## Saving the gating thresholds
    # all.gthres[5] <- cd66.gate.low
    # all.gthres[6] <- cd66.gate
    # all.gthres[7] <- cd45.gate.low
    # all.gthres[8] <- cd45.gate
    # 
    # plotDens(leukocytes, c(channels.ind["CD66"],channels.ind["CD45"]), main = "Mononuclear cells & Granulocytes", cex.lab = 2, cex.axis = 2, cex.main=2)
    # lines(mononuclear.flowD@filter, lwd = 2)
    # lines(granulocytes.flowD@filter, lwd = 2)

    
    
    theta = atan(atan(tan(pi/4)))
    R <- matrix( c(cos(theta), sin(theta), -sin(theta), cos(theta)),2 ,2)
    
    leukocytes.temp <-rotate.data(leukocytes,c(channels.ind["CD66"],channels.ind["CD45"]),theta = pi/4)$data
    
    cd45.gate <- deGate(leukocytes.temp, channel = channels.ind["CD45"])
    if(cd45.gate > 0.80){
      cd45.gate <- mean(c(deGate(leukocytes.temp, channel = channels.ind["CD45"], tinypeak.removal = 0.9), deGate(leukocytes.temp, channel = channels.ind["CD45"])))
    }
    if(cd45.gate < 0){
      cd45.gate <- median(deGate(leukocytes.temp, channel = channels.ind["CD45"], all.cuts = T))
      if(cd45.gate < 0){
        cd45.gate <- max(deGate(leukocytes.temp, channel = channels.ind["CD45"], all.cuts = T))
      }
    }
    cd66.gate <- deGate(leukocytes.temp, channel = channels.ind["CD66"], use.upper = T, upper = T)
    
    #cd66.gate.low <- deGate()
    mononuclear.flowD.temp <- flowDensity(leukocytes.temp, channels = c('La139Di','In115Di'), position = c(F,T), gates = c(cd66.gate-0.75, cd45.gate+0.5), ellip.gate = T)
    mononuclear.temp <- getflowFrame(mononuclear.flowD.temp)
    mononuclear.temp@exprs[,c(channels.ind["CD66"],channels.ind["CD45"])] <- t(t(R) %*% t(mononuclear.temp@exprs[,c(channels.ind["CD66"],channels.ind["CD45"])]))
    mononuclear.flowD.temp@filter <- rotate.data(mononuclear.flowD.temp@filter,c(channels.ind["CD66"],channels.ind["CD45"]),theta = -pi/4)$data
    mononuclear.flowD.temp@flow.frame <- rotate.data(getflowFrame(mononuclear.flowD.temp),c(channels.ind["CD66"],channels.ind["CD45"]),theta = -pi/4)$data
    mononuclear.flowD.temp@proportion <- (mononuclear.flowD.temp@cell.count/nrow(leukocytes))*100
    
    mononuclear.flowD <- flowDensity(mononuclear.flowD.temp, channels = c('La139Di','In115Di'), position = c(T,NA), gates = c(cd66.gate.low, NA))
    #mononuclear.temp <- getflowFrame(mononuclear.flowD)
    
    mononuclear <- getflowFrame(mononuclear.flowD)
  
    
    granulocytes.flowD <- flowDensity(leukocytes.temp, channels = c('La139Di','In115Di'), position = c(F,F), gates = c(cd66.gate, cd45.gate-0.5), ellip.gate = T)
    granulocytes <- getflowFrame(granulocytes.flowD)
    granulocytes@exprs[,c(channels.ind["CD66"],channels.ind["CD45"])] <- t(t(R) %*% t(granulocytes@exprs[,c(channels.ind["CD66"],channels.ind["CD45"])]))
    granulocytes.flowD@filter <- rotate.data(granulocytes.flowD@filter,c(channels.ind["CD66"],channels.ind["CD45"]),theta = -pi/4)$data
    granulocytes.flowD@flow.frame <- rotate.data(getflowFrame(granulocytes.flowD),c(channels.ind["CD66"],channels.ind["CD45"]),theta = -pi/4)$data
    granulocytes.flowD@proportion <- (granulocytes.flowD@cell.count/nrow(leukocytes))*100
    granulocytes <- getflowFrame(granulocytes.flowD)
    
    all.events[4] <- nrow(mononuclear)
    all.props[4] <- (nrow(mononuclear)/nrow(leukocytes))*100
    
    all.events[5] <- nrow(granulocytes)
    all.props[5] <- (nrow(granulocytes)/nrow(leukocytes))*100
    
    ## Saving the gating thresholds
    all.gthres[5] <- cd66.gate.low
    all.gthres[6] <- cd66.gate
    #all.gthres[7] <- cd45.gate.low
    all.gthres[7] <- cd45.gate
    
    plotDens(leukocytes, c(channels.ind["CD66"],channels.ind["CD45"]), main = "Mononuclear cells & Granulocytes", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(mononuclear.flowD@filter, lwd = 2)
    lines(granulocytes.flowD@filter, lwd = 2)
    
    
    plotDens(leukocytes.temp, c(channels.ind["CD66"],channels.ind["CD45"]), main = "", cex.lab = 2, cex.axis = 2, cex.main=2); abline(h=cd45.gate); abline(v=cd66.gate)
    lines(mononuclear.flowD.temp@filter, lwd=2)
    lines(granulocytes.flowD@filter, lwd=2)
    
    
    #############################################################################
    ## B cells & T cells from Mononuclear cells---------------------

    cd3.gate <- deGate(mononuclear, channel = channels.ind["CD3"])
  
    mononuclear.flowD.temp <- flowDensity(mononuclear, channels = c('Er170Di', 'Nd142Di'), position = c(T, NA), gates = c(cd3.gate, NA))
    cd19.gate <- deGate(mononuclear.flowD.temp, channel = channels.ind["CD19"], tinypeak.removal = 0.001)
    
    if(cd19.gate > 5){
      mononuclear.flowD.temp <- flowDensity(mononuclear, channels = c('Er170Di', 'Nd142Di'), position = c(F, NA), gates = c(cd3.gate, NA))
      cd19.gate <- deGate(mononuclear.flowD.temp, channel = channels.ind["CD19"])
      if(cd19.gate < 4){
        mononuclear.flowD.temp <- flowDensity(mononuclear, channels = c('Er170Di', 'Nd142Di'), position = c(T, NA), gates = c(cd3.gate, NA))
        cd19.gate <- deGate(mononuclear.flowD.temp, channel = channels.ind["CD19"], use.upper = T, upper = T)
    
      }
    }
    Bcells.flowD <- flowDensity(mononuclear, channels = c('Er170Di', 'Nd142Di'), position = c(F, T), gates = c(cd3.gate, cd19.gate))
    Tcells.flowD <- flowDensity(mononuclear, channels = c('Er170Di', 'Nd142Di'), position = c(T, F), gates = c(cd3.gate, cd19.gate))
    NK.LinNeg.flowD <- flowDensity(mononuclear, channels = c('Er170Di', 'Nd142Di'), position = c(F,F), gates = c(cd3.gate, cd19.gate))

    Bcells <- getflowFrame(Bcells.flowD)
    Tcells <- getflowFrame(Tcells.flowD)
    NK.LinNeg <- getflowFrame(NK.LinNeg.flowD) ## Population used in the gating next step to obtain NK cells & lin-

    all.events[6] <- nrow(Bcells)
    all.props[6] <- (nrow(Bcells)/nrow(mononuclear))*100

    all.events[7] <- nrow(Tcells)
    all.props[7] <- (nrow(Tcells)/nrow(mononuclear))*100

    ## Saving proportions as percentage of Mononuclear cells
    all.freqFeatures[1] <- nrow(Bcells)/nrow(mononuclear)
    all.freqFeatures[2] <- nrow(Tcells)/nrow(mononuclear)
    
    ## Saving the gating thresholds
    all.gthres[8] <- cd3.gate
    all.gthres[9] <- cd19.gate
    
    plotDens(mononuclear, c(channels.ind["CD3"],channels.ind["CD19"]), main = "B cells & T cells", cex.lab = 2, cex.axis = 2, cex.main=2); abline(h=cd19.gate, lwd=2); abline(v=cd3.gate, lwd=2)
    lines(Bcells.flowD@filter, lwd = 2)
    lines(Tcells.flowD@filter, lwd = 2)

    #############################################################################
    ## NK cells & lin- from NK.LinNeg------------------

    #cd14.gate.high <- deGate(NK.LinNeg, channel = channels.ind["CD4"], use.percentile = T, percentile = 0.9999)
    cd14.gate.high <- deGate(NK.LinNeg, channel = channels.ind["CD4"], use.percentile = T, percentile = 0.0001)
    cd7.gate <- deGate(NK.LinNeg, channel = channels.ind["CD7"])
    
    NKcells.flowD <- flowDensity(NK.LinNeg, channels = c('Lu175Di', 'Pr141Di'), position = c(T, T), gates = c(cd14.gate.high, cd7.gate))
    LinNeg.flowD <- flowDensity(NK.LinNeg, channels = c('Lu175Di', 'Pr141Di'), position = c(T, F), gates = c(cd14.gate.high, cd7.gate))
   
    NKcells <- getflowFrame(NKcells.flowD)
    LinNeg <- getflowFrame(LinNeg.flowD)

    all.events[8] <- nrow(NKcells)
    all.props[8] <- (nrow(NKcells)/nrow(NK.LinNeg))*100

    all.events[9] <- nrow(LinNeg)
    all.props[9] <- (nrow(LinNeg)/nrow(NK.LinNeg))*100

    ## Saving proportions as percentage of Mononuclear cells
    all.freqFeatures[3] <- nrow(NKcells)/nrow(mononuclear)
    all.freqFeatures[4] <- nrow(LinNeg)/nrow(mononuclear)
    
    ## Saving the gating thresholds
    all.gthres[10] <- cd14.gate.high
    all.gthres[11] <- cd7.gate


   
    min.y <-  min(exprs(NK.LinNeg)[, c('Pr141Di')])
    max.y <- max(exprs(NK.LinNeg)[, c('Pr141Di')])+1
    min.x <- min(exprs(NK.LinNeg)[, c('Lu175Di')])
    max.x <-  max(exprs(NK.LinNeg)[, c('Lu175Di')])+2
    plotDens(NK.LinNeg, c(channels.ind["CD14"],channels.ind["CD7"]), main = "NK cells & lin-", cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x, max.x), ylim = c(min.y, max.y))
    lines(NKcells.flowD@filter, lwd = 2)
    lines(LinNeg.flowD@filter, lwd = 2)


    #############################################################################
    ## CD56low CD16+ & CD56+ CD16- cells from NK cells---------------

    cd56.gate.low <- deGate(NKcells, channel = channels.ind["CD56"], use.percentile = T, percentile = 0.1)
    cd56.gate.mid <- deGate(NKcells, channel = channels.ind["CD56"], use.upper = T, upper = T)
    cd56.gate.high <- deGate(NKcells, channel = channels.ind["CD56"], use.percentile = T, percentile = 1)

    #cd16.gate.low <- deGate(NKcells, channel = channels.ind["CD16"], tinypeak.removal = 0.9999)
    cd16.gate.low <- deGate(NKcells, channel = channels.ind["CD16"], use.upper = T, upper = F)
    cd16.gate.mid <- deGate(NKcells, channel = channels.ind["CD16"], use.percentile = T, percentile = 0.15)
    cd16.gate.high <- deGate(NKcells, channel = channels.ind["CD16"], use.percentile = T, percentile = 0.99)

    ## CD56low CD16+ cells
    cd56low.cd16pos.temp <- flowDensity(NKcells, channels = c('Yb176Di', 'Ho165Di'), position = c(F,T), gates = c(cd56.gate.mid, cd16.gate.mid))
    cd56low.cd16pos.flowD <- flowDensity(cd56low.cd16pos.temp, channels = c('Yb176Di', 'Ho165Di'), position = c(T,F), gates = c(cd56.gate.low-0.15, cd16.gate.high))



    ## CD56+ CD16- cells
    cd56pos.cd16neg.flowD.temp <- flowDensity(NKcells, channels = c('Yb176Di', 'Ho165Di'), position = c(F,F), gates = c(cd56.gate.high+0.25, cd16.gate.mid-0.05))
    cd56pos.cd16neg.flowD <- flowDensity(cd56pos.cd16neg.flowD.temp, channels = c('Yb176Di', 'Ho165Di'), position = c(T,T), gates = c(cd56.gate.low+0.5, cd16.gate.low-1))

    cd56low.cd16pos <- getflowFrame(cd56low.cd16pos.flowD)
    cd56pos.cd16neg <- getflowFrame(cd56pos.cd16neg.flowD)

    all.events[10] <- nrow(cd56low.cd16pos)
    all.props[10] <- (nrow(cd56low.cd16pos)/nrow(NKcells))*100

    all.events[11] <- nrow(cd56pos.cd16neg)
    all.props[11] <- (nrow(cd56pos.cd16neg)/nrow(NKcells))*100

    ## Saving proportions as percentage of Mononuclear cells
    all.freqFeatures[5] <- nrow(cd56low.cd16pos)/nrow(mononuclear)
    all.freqFeatures[6] <- nrow(cd56pos.cd16neg)/nrow(mononuclear)
    
    ## Saving the gating thresholds
    all.gthres[12] <- cd56.gate.low
    all.gthres[13] <- cd56.gate.mid
    all.gthres[14] <- cd56.gate.high
    all.gthres[15] <- cd16.gate.low
    all.gthres[16] <- cd16.gate.mid
    all.gthres[17] <- cd16.gate.high

    
    min.y <-  min(exprs(NKcells)[, c('Ho165Di')])-1
    max.y <- max(exprs(NKcells)[, c('Ho165Di')])+1
    min.x <- min(exprs(NKcells)[, c('Yb176Di')])
    max.x <-  max(exprs(NKcells)[, c('Yb176Di')])+2
    plotDens(NKcells, c(channels.ind["CD56"],channels.ind["CD16"]), main = "CD56low CD16+/CD56+CD16- NKcells", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(cd56low.cd16pos.flowD@filter)
    lines(cd56pos.cd16neg.flowD@filter)

    #############################################################################
    ## cMC, ncMC, & intMC from lin- cells--------------------------

    cd14.gate.low <- deGate(LinNeg, channel = channels.ind["CD14"], use.percentile = T, percentile = 0.3)
    cd16.gate.LinNeg <- deGate(LinNeg, channel = channels.ind["CD16"]) 
    
    
    cMC.flowD <- flowDensity(LinNeg, channels = c('Lu175Di', 'Ho165Di'), position = c(T, F), gates = c(cd14.gate.low, cd16.gate.LinNeg))
    ncMC.flowD <- flowDensity(LinNeg, channels = c('Lu175Di', 'Ho165Di'), position = c(F, T), gates = c(cd14.gate.low, cd16.gate.LinNeg))
    intMC.flowD <- flowDensity(LinNeg, channels = c('Lu175Di', 'Ho165Di'), position = c(T, T), gates = c(cd14.gate.low, cd16.gate.LinNeg))

    # ## CD14+, which is a combined population of cMC and int MC is needed for obtaining M-MDSCs population
    # cd14pos.flowD <- flowDensity(LinNeg, channels = c('Lu175Di', 'Ho165Di'), position = c(T, NA), gates = c(cd14.gate.low, NA))

    ## NOT MC is needed for obtaining pDCs and mDCs
    NOT.MC.flowD <- flowDensity(LinNeg, channels = c('Lu175Di', 'Ho165Di'), position = c(F, F), gates = c(cd14.gate.low, cd16.gate.LinNeg))


    cMC <- getflowFrame(cMC.flowD)
    ncMC <- getflowFrame(ncMC.flowD)
    intMC <- getflowFrame(intMC.flowD)
    NOT.MC <- getflowFrame(NOT.MC.flowD)
    # cd14pos <- getflowFrame(cd14pos.flowD)


    all.events[12] <- nrow(cMC)
    all.props[12] <- (nrow(cMC)/nrow(LinNeg))*100

    all.events[13] <- nrow(ncMC)
    all.props[13] <- (nrow(ncMC)/nrow(LinNeg))*100

    all.events[14] <- nrow(intMC)
    all.props[14] <- (nrow(intMC)/nrow(LinNeg))*100

    ## Saving proportions as percentage of Mononuclear cells
    all.freqFeatures[7] <- nrow(cMC)/nrow(mononuclear)
    all.freqFeatures[8] <- nrow(ncMC)/nrow(mononuclear)
    all.freqFeatures[9] <- nrow(intMC)/nrow(mononuclear)
    
    ## Saving the gating thresholds
    all.gthres[18] <- cd14.gate.low
    all.gthres[19] <- cd16.gate.LinNeg
    
    min.y <-  min(exprs(LinNeg)[, c('Ho165Di')])-1
    max.y <- max(exprs(LinNeg)[, c('Ho165Di')])+1
    min.x <- min(exprs(LinNeg)[, c('Lu175Di')])
    max.x <-  max(exprs(LinNeg)[, c('Lu175Di')])+1
    plotDens(LinNeg, c(channels.ind["CD14"],channels.ind["CD16"]), main = "cMC, ncMC, & intMC", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(cMC.flowD@filter, lwd = 2)
    lines(ncMC.flowD@filter, lwd = 2)
    lines(intMC.flowD@filter, lwd = 2)
    lines(NOT.MC.flowD@filter, lwd = 2)
    # lines(cd14pos.flowD@filter, lwd =2)


    ##############################################################################################################
    ## pDCs from the NOT MC cells

    cd123.gate <- deGate(NOT.MC, channel = channels.ind["CD123"])
    hladr.gate.DCs <- deGate(NOT.MC, channel = channels.ind["HLADR"], use.upper = T, upper = F, tinypeak.removal = 0.9)

    pDC.flowD <- flowDensity(NOT.MC, channels = c('Nd148Di', 'Yb174Di'), position = c(T, T), gates = c(cd123.gate, hladr.gate.DCs))

    pDC <- getflowFrame(pDC.flowD)

    all.events[15] <- nrow(pDC)
    all.props[15] <- (nrow(pDC)/nrow(NOT.MC))*100

    ## Saving proportions as percentage of Mononuclear cells
    all.freqFeatures[10] <- nrow(pDC)/nrow(mononuclear)
    
    
    ## Saving the gating thresholds
    all.gthres[20] <- cd123.gate
    all.gthres[21] <- hladr.gate.DCs

    min.y <-  min(exprs(NOT.MC)[, c('Yb174Di')])
    max.y <- max(exprs(NOT.MC)[, c('Yb174Di')])+1
    min.x <- min(exprs(NOT.MC)[, c('Nd148Di')])
    max.x <-  max(exprs(NOT.MC)[, c('Nd148Di')])+2
    plotDens(NOT.MC, c(channels.ind["CD123"],channels.ind["HLADR"]), main = "pDCs", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(pDC.flowD@filter, lwd=2)


    ##############################################################################################################
    # mDCs from the NOT MC cells

    # cd11c.gate <- deGate(NOT.MC, channel = channels.ind["CD11c"], all.cuts = T)
    # if(length(cd11c.gate) > 1){
    #   cd11c.gate <- cd11c.gate[2]
    # }
    
    cd11c.gate <- deGate(NOT.MC, channel = channels.ind["CD11c"], use.percentile = T, percentile = 0.1)

    mDCs.flowD <- flowDensity(NOT.MC, channels = c('Sm147Di', 'Yb174Di'), position = c(T,T), gates = c(cd11c.gate, hladr.gate.DCs))

    mDC <- getflowFrame(mDCs.flowD)

    all.events[16] <- nrow(mDC)
    all.props[16] <- (nrow(mDC)/nrow(NOT.MC))*100

    ## Saving proportions as percentage of Mononuclear cells
    all.freqFeatures[11] <- nrow(mDC)/nrow(mononuclear)
    
    ## Saving the gating thresholds
    all.gthres[22] <- cd11c.gate
    
    min.y <-  min(exprs(NOT.MC)[, c('Yb174Di')])
    max.y <- max(exprs(NOT.MC)[, c('Yb174Di')])+1
    min.x <- min(exprs(NOT.MC)[, c('Sm147Di')])
    max.x <-  max(exprs(NOT.MC)[, c('Sm147Di')])+2
    plotDens(NOT.MC, c(channels.ind["CD11c"],channels.ind["HLADR"]), main = "mDCs", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(mDCs.flowD@filter, lwd = 2)

    #############################################################################################################
    ## M-MDSCs from cMC population.

    #cd11b.gate <- deGate(cd14pos, channel = channels.ind["CD11b"], use.upper = T, upper = T)
    numPeaks <- getPeaks(cMC, channel = channels.ind["CD11b"])
    if(length(numPeaks) == 1){
      cd11b.gate <- numPeaks$Peaks
    }else{
      cd11b.gate <- numPeaks$Peaks[2]
    }
   

    hladr.gate.MMDSCs <- deGate(cMC, channel = channels.ind["HLADR"], use.upper = T, upper = F, use.percentile = T, percentile = 0.05)
    
    MMDSC.flowD <- flowDensity(cMC, channels = c('Nd144Di', 'Yb174Di'), position = c(T, F), gates = c(cd11b.gate, hladr.gate.MMDSCs))

    MMDSC <- getflowFrame(MMDSC.flowD)

    all.events[17] <- nrow(MMDSC)
    all.props[17] <- (nrow(MMDSC)/nrow(cMC))*100

    ## Saving proportions as percentage of Mononuclear cells
    all.freqFeatures[12] <- nrow(MMDSC)/nrow(mononuclear)
    
    ## Saving the gating thresholds
    all.gthres[23] <- cd11b.gate
    all.gthres[24] <- hladr.gate.MMDSCs
    
    min.y <-  min(exprs(cMC)[, c('Yb174Di')])
    max.y <- max(exprs(cMC)[, c('Yb174Di')])+1
    min.x <- min(exprs(cMC)[, c('Nd144Di')])
    max.x <-  max(exprs(cMC)[, c('Nd144Di')])+2
    plotDens(cMC, c(channels.ind["CD11b"],channels.ind["HLADR"]), main = "M-MDSCs", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(MMDSC.flowD@filter, lwd =2)
    lines(x=c(cd11b.gate, max.x+1), y=c(hladr.gate.MMDSCs, hladr.gate.MMDSCs), lwd = 2)
    lines(x=c(cd11b.gate, cd11b.gate), y=c(min.y-1, hladr.gate.MMDSCs), lwd = 2)


   #############################################################################################################
    ## CD8+ & CD4+ T cells from T cells

    cd4.gate <- deGate(Tcells, channel = channels.ind["CD4"], all.cuts = T)
    if(length(cd4.gate > 1)){
     cd4.gate <- cd4.gate[2]
    }

    cd8a.gate <- deGate(Tcells, channel = channels.ind["CD8"], all.cuts = T)
    if(length(cd8a.gate > 1)){
      cd8a.gate <- cd8a.gate[2]
    }

    cd4.Tcells.flowD <- flowDensity(Tcells, channels = c('Nd145Di','Nd146Di'), position = c(T,F), gates = c(cd4.gate, cd8a.gate))
    cd8.Tcells.flowD <- flowDensity(Tcells, channels = c('Nd145Di','Nd146Di'), position = c(F,T), gates = c(cd4.gate, cd8a.gate+0.5))
    NOT.cd4.cd8.Tcells.flowD <- flowDensity(Tcells, channels = c('Nd145Di','Nd146Di'), position = c(F,F), gates = c(cd4.gate-0.25, cd8a.gate))


    cd4.Tcells <- getflowFrame(cd4.Tcells.flowD)
    cd8.Tcells <- getflowFrame(cd8.Tcells.flowD)
    NOT.cd4.cd8.Tcells <- getflowFrame(NOT.cd4.cd8.Tcells.flowD)

    all.events[18] <- nrow(cd4.Tcells)
    all.props[18] <- (nrow(cd4.Tcells)/nrow(Tcells))*100
    
    all.events[19] <- nrow(cd8.Tcells)
    all.props[19] <- (nrow(cd8.Tcells)/nrow(Tcells))*100

    ## Saving proportions as percentage of Mononuclear cells
    all.freqFeatures[13] <- nrow(cd4.Tcells)/nrow(mononuclear)
    all.freqFeatures[14] <- nrow(cd8.Tcells)/nrow(mononuclear)
    

    ## Saving the gating thresholds
    all.gthres[25] <- cd4.gate
    all.gthres[26] <- cd8a.gate
    
 
    plotDens(Tcells, c(channels.ind["CD4"],channels.ind["CD8"]), main = "CD8+ Tcells / CD4+ T cells", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(cd4.Tcells.flowD@filter, lwd = 2)
    lines(cd8.Tcells.flowD@filter, lwd = 2)
    lines(NOT.cd4.cd8.Tcells.flowD@filter, lwd = 2)


    ##################################################################################################
    ## CD4+ T Naive & CD4+ T memory cells from CD4+ T cells

    #cd45RA.gate.low.cd4 <- min(deGate(cd4.Tcells, channel = channels.ind["CD45RA"], all.cuts = T))
    cd45RA.gate.high.cd4 <- max(deGate(cd4.Tcells, channel = channels.ind["CD45RA"], all.cuts = T))
    # if(length(cd45RA.gate > 1)){
    #   cd45RA.gate <- cd45RA.gate[2]
    # }

    cd4.T.Naive.flowD <- flowDensity(cd4.Tcells, channels = c('Nd145Di', 'Nd143Di'), position = c(NA, T), gates = c(NA, cd45RA.gate.high.cd4))
    cd4.T.Memory.flowD <- flowDensity(cd4.Tcells, channels = c('Nd145Di', 'Nd143Di'), position = c(NA, F), gates = c(NA, cd45RA.gate.high.cd4))
    #cd4.T.Memory.flowD <- flowDensity(cd4.T.Memory.flowD.temp, channels = c('Nd145Di', 'Nd143Di'), position = c(NA, T), gates = c(NA, cd45RA.gate.low.cd4))

    cd4.T.Naive <- getflowFrame(cd4.T.Naive.flowD)
    cd4.T.Memory <- getflowFrame(cd4.T.Memory.flowD)

    all.events[20] <- nrow(cd4.T.Naive)
    all.props[20] <- (nrow(cd4.T.Naive)/nrow(cd4.Tcells))*100
    
    all.events[21] <- nrow(cd4.T.Memory)
    all.props[21] <- (nrow(cd4.T.Memory)/nrow(cd4.Tcells))*100

    ## Saving proportions as percentage of Mononuclear cells
    all.freqFeatures[15] <- nrow(cd4.T.Naive)/nrow(mononuclear)
    all.freqFeatures[16] <- nrow(cd4.T.Memory)/nrow(mononuclear)
    
    ## Saving the gating thresholds
    all.gthres[27] <- cd45RA.gate.high.cd4

    plotDens(cd4.Tcells, c(channels.ind["CD4"],channels.ind["CD45RA"]), main = "CD4+ T naive & CD4+ T memory cells", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(cd4.T.Naive.flowD@filter, lwd = 2)
    lines(cd4.T.Memory.flowD@filter, lwd = 2)

    ##################################################################################################

    ## CD4+ Th1 naive & CD4+ Th1 memory cells from CD4+ T cells

    cd4.Th1.Naive.flowD.temp <- flowDensity(cd4.Tcells, channels = c('Gd160Di', 'Nd143Di'), position = c(NA, T), gates = c(NA, cd45RA.gate.high.cd4))
    
    #plot(cd4.Tcells, cd4.Th1.Naive.flowD.temp)
    
    Tbet.gate.CD4T <- deGate(cd4.Th1.Naive.flowD.temp, channel = channels.ind["Tbet"], use.upper = T, upper = T)

    cd4.Th1.Naive.flowD <- flowDensity(cd4.Tcells, channels = c('Gd160Di', 'Nd143Di'), position = c(T, T), gates = c(Tbet.gate.CD4T, cd45RA.gate.high.cd4))
    cd4.Th1.Memory.flowD <- flowDensity(cd4.Tcells, channels = c('Gd160Di', 'Nd143Di'), position = c(T, F), gates = c(Tbet.gate.CD4T, cd45RA.gate.high.cd4))
    #cd4.Th1.Memory.flowD <- flowDensity(cd4.Th1.Memory.flowD.temp, channels = c('Gd160Di', 'Nd143Di'), position = c(NA, T), gates = c(NA, cd45RA.gate.low.cd4))

    cd4.Th1.Naive <- getflowFrame(cd4.Th1.Naive.flowD)
    cd4.Th1.Memory <- getflowFrame(cd4.Th1.Memory.flowD)

    all.events[22] <- nrow(cd4.Th1.Naive)
    all.props[22] <- (nrow(cd4.Th1.Naive)/nrow(cd4.Tcells))*100
    
    all.events[23] <- nrow(cd4.Th1.Memory)
    all.props[23] <- (nrow(cd4.Th1.Memory)/nrow(cd4.Tcells))*100

    ## Saving proportions as percentage of Mononuclear cells
    all.freqFeatures[17] <- nrow(cd4.Th1.Naive)/nrow(mononuclear)
    all.freqFeatures[18] <- nrow(cd4.Th1.Memory)/nrow(mononuclear)
    
    ## Saving the gating thresholds
    all.gthres[28] <- Tbet.gate.CD4T


    min.y <-  min(exprs(cd4.Tcells)[, c('Nd143Di')])
    max.y <- max(exprs(cd4.Tcells)[, c('Nd143Di')])+1
    min.x <- min(exprs(cd4.Tcells)[, c('Gd160Di')])
    max.x <- max(exprs(cd4.Tcells)[, c('Gd160Di')])+2
    plotDens(cd4.Tcells, c(channels.ind["Tbet"],channels.ind["CD45RA"]), main = "CD4+ Th1 naive & CD4+ Th1 memory cells", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(cd4.Th1.Naive.flowD@filter, lwd=2)
    lines(cd4.Th1.Memory.flowD@filter, lwd=2)


    ##################################################################################################
    ## Tregs naive from CD4+ T naive cells
    ## Implementing rotation to gate Tregs Naive

    theta = atan(atan(tan(pi/3)))
    R <- matrix( c(cos(theta), sin(theta), -sin(theta), cos(theta)),2 ,2)

    cd4.T.Naive.temp <-rotate.data(cd4.T.Naive,c(channels.ind["FoxP3"],channels.ind["CD25"]),theta = -pi/3)$data
    
    FoxP3.gate.Naive.high <-  deGate(cd4.T.Naive.temp, channel = channels.ind["FoxP3"], use.percentile = T, percentile = 0.1)
    FoxP3.gate.Naive.low <- deGate(cd4.T.Naive.temp, channel = channels.ind["FoxP3"], use.percentile = T, percentile = 0.001)
    cd25.gate.Naive <- deGate(cd4.T.Naive.temp, channel = channels.ind["CD25"], use.percentile = T, percentile = 0.95)

    Tregs.Naive.flowD.temp <- flowDensity(cd4.T.Naive.temp, channels = c('Dy162Di', 'Tm169Di'), position = c(F,T), gates = c(FoxP3.gate.Naive.high, cd25.gate.Naive))
    Tregs.Naive.flowD <- flowDensity(Tregs.Naive.flowD.temp, channels = c('Dy162Di', 'Tm169Di'), position = c(T,NA), gates = c(FoxP3.gate.Naive.low, NA))
    
    Tregs.Naive <- getflowFrame(Tregs.Naive.flowD)

    Tregs.Naive@exprs[,c(channels.ind["FoxP3"],channels.ind["CD25"])] <- t(t(R) %*% t(Tregs.Naive@exprs[,c(channels.ind["FoxP3"],channels.ind["CD25"])]))
    # TregsMempoints <- t(t(R) %*% t(Tregs.Naive.flowD@filter))


    Tregs.Naive.flowD@filter <- rotate.data(Tregs.Naive.flowD@filter,c(channels.ind["FoxP3"],channels.ind["CD25"]),theta = pi/3)$data
    Tregs.Naive.flowD@flow.frame <- rotate.data(getflowFrame(Tregs.Naive.flowD),c(channels.ind["FoxP3"],channels.ind["CD25"]),theta = pi/3)$data
    Tregs.Naive.flowD@proportion <- (Tregs.Naive.flowD@cell.count/nrow(cd4.T.Naive))*100

    Tregs.Naive <- getflowFrame(Tregs.Naive.flowD)

    all.events[24] <- nrow(Tregs.Naive)
    all.props[24] <- (nrow(Tregs.Naive)/nrow(cd4.T.Naive))*100

    ## Saving proportions as percentage of Mononuclear cells
    all.freqFeatures[19] <- nrow(Tregs.Naive)/nrow(mononuclear)
    

    ## Saving the gating thresholds (rotated gates)
    all.gthres[29] <- FoxP3.gate.Naive.low
    all.gthres[30] <- FoxP3.gate.Naive.high
    all.gthres[31] <- cd25.gate.Naive

    min.y <-  min(exprs(cd4.T.Naive)[, c('Tm169Di')])
    max.y <- max(exprs(cd4.T.Naive)[, c('Tm169Di')])+2
    min.x <- min(exprs(cd4.T.Naive)[, c('Dy162Di')])
    max.x <-  max(exprs(cd4.T.Naive)[, c('Dy162Di')])+2
    plotDens(cd4.T.Naive, c(channels.ind["FoxP3"],channels.ind["CD25"]), main = "Tregs Naive", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    # abline(h=cd25.gate.Naive); abline(v=FoxP3.gate.Naive)
    lines(Tregs.Naive.flowD@filter, lwd = 2)


    ##################################################################################################
    ## Tregs memory from CD4+ T naive cells
    ## Implementing rotation to gate Tregs Memory

    theta = atan(atan(tan(pi/4)))
    R <- matrix( c(cos(theta), sin(theta), -sin(theta), cos(theta)),2 ,2)

    cd4.T.Memory.temp <-rotate.data(cd4.T.Memory,c(channels.ind["FoxP3"],channels.ind["CD25"]),theta = -pi/4)$data

    # plotDens(cd4.T.Memory.temp, c(channels.ind["FoxP3"],channels.ind["CD25"]), main = "Tregs Memory", cex.lab = 2, cex.axis = 2, cex.main=2)
    # lines(Tregs.Memory.flowD@filter)
    
    FoxP3.gate.Memory.low <-  deGate(cd4.T.Memory.temp, channel = channels.ind["FoxP3"], use.upper = T, upper = F)
    FoxP3.gate.Memory.high <-  deGate(cd4.T.Memory.temp, channel = channels.ind["FoxP3"])
    cd25.gate.Memory <- deGate(cd4.T.Memory.temp, channel = channels.ind["CD25"], use.percentile = T, percentile = 0.95)

    Tregs.Memory.flowD.temp <- flowDensity(cd4.T.Memory.temp, channels = c('Dy162Di', 'Tm169Di'), position = c(T,T), gates = c(FoxP3.gate.Memory.low, cd25.gate.Memory))
    Tregs.Memory.flowD <- flowDensity(Tregs.Memory.flowD.temp, channels = c('Dy162Di', 'Tm169Di'), position = c(F,NA), gates = c(FoxP3.gate.Memory.high, NA))
    
    Tregs.Memory <- getflowFrame(Tregs.Memory.flowD)

    Tregs.Memory@exprs[,c(channels.ind["FoxP3"],channels.ind["CD25"])] <- t(t(R) %*% t(Tregs.Memory@exprs[,c(channels.ind["FoxP3"],channels.ind["CD25"])]))

    Tregs.Memory.flowD@filter <- rotate.data(Tregs.Memory.flowD@filter,c(channels.ind["FoxP3"],channels.ind["CD25"]),theta = pi/4)$data
    Tregs.Memory.flowD@flow.frame <- rotate.data(getflowFrame(Tregs.Memory.flowD),c(channels.ind["FoxP3"],channels.ind["CD25"]),theta = pi/4)$data
    Tregs.Memory.flowD@proportion <- (Tregs.Memory.flowD@cell.count/nrow(cd4.T.Memory))*100

    Tregs.Memory <- getflowFrame(Tregs.Memory.flowD)

    all.events[25] <- nrow(Tregs.Memory)
    all.props[25] <- (nrow(Tregs.Memory)/nrow(cd4.T.Memory))*100
    
    ## Saving proportions as percentage of Mononuclear cells
    all.freqFeatures[20] <- nrow(Tregs.Memory)/nrow(mononuclear)
    
    ## Saving the gating thresholds (rotated gates)
    all.gthres[32] <- FoxP3.gate.Memory.low
    all.gthres[33] <- FoxP3.gate.Memory.high
    all.gthres[34] <- cd25.gate.Memory
    
  
    min.y <-  min(exprs(cd4.T.Memory)[, c('Tm169Di')])
    max.y <- max(exprs(cd4.T.Memory)[, c('Tm169Di')])+2
    min.x <- min(exprs(cd4.T.Memory)[, c('Dy162Di')])
    max.x <-  max(exprs(cd4.T.Memory)[, c('Dy162Di')])+2
    plotDens(cd4.T.Memory, c(channels.ind["FoxP3"],channels.ind["CD25"]), main = "Tregs Memory", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    # abline(h=cd25.gate); abline(v=FoxP3.gate)
    lines(Tregs.Memory.flowD@filter, lwd = 2)
    
    #############################################################################################
    ## Gamma-Delta T-cells from NOT.cd4.cd8.Tcells.flowD
    
    TCRD.gate <- deGate(NOT.cd4.cd8.Tcells, channel = channels.ind["TCRgd"])
    
    gammaDelta.Tcells.flowD <- flowDensity(NOT.cd4.cd8.Tcells, channels = c('Sm152Di', 'Er170Di'), position = c(T,NA), gates = c(TCRD.gate+0.5,NA), ellip.gate=T)
    
    gammaDelta.Tcells <- getflowFrame(gammaDelta.Tcells.flowD)
    
    all.events[26] <- nrow(gammaDelta.Tcells)
    all.props[26] <- (nrow(gammaDelta.Tcells)/nrow(NOT.cd4.cd8.Tcells))*100
    
    ## Saving proportions as percentage of Mononuclear cells
    all.freqFeatures[21] <- nrow(gammaDelta.Tcells)/nrow(mononuclear)
    
    ## Saving the gating thresholds 
    all.gthres[35] <- TCRD.gate
    
  
    min.y <-  min(exprs(NOT.cd4.cd8.Tcells)[, c('Er170Di')])
    max.y <- max(exprs(NOT.cd4.cd8.Tcells)[, c('Er170Di')])+2
    min.x <- min(exprs(NOT.cd4.cd8.Tcells)[, c('Sm152Di')])
    max.x <-  max(exprs(NOT.cd4.cd8.Tcells)[, c('Sm152Di')])+2
    plotDens(NOT.cd4.cd8.Tcells, c(channels.ind["TCRgd"],channels.ind["CD3"]), main = "gamma-delta T-cells", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(gammaDelta.Tcells.flowD@filter, lwd = 2)
    
    
    #################################################################################################################
    ## CD8+ T naive & CD8+ T memory cells from CD8+ T cells
    
    cd45RA.gate.cd8 <- deGate(cd8.Tcells, channel = channels.ind["CD45RA"], all.cuts = T)
    if(length(cd45RA.gate.cd8) > 1){
      cd45RA.gate.cd8 <- cd45RA.gate.cd8[2]
    }
    
    cd8.T.Naive.flowD <- flowDensity(cd8.Tcells, channels = c('Nd146Di', 'Nd143Di'), position = c(NA, T), gates = c(NA, cd45RA.gate.cd8))
    cd8.T.Memory.flowD <- flowDensity(cd8.Tcells, channels = c('Nd146Di', 'Nd143Di'), position = c(NA, F), gates = c(NA, cd45RA.gate.cd8))
    
    cd8.T.Naive <- getflowFrame(cd8.T.Naive.flowD)
    cd8.T.Memory <- getflowFrame(cd8.T.Memory.flowD)
    
    
    all.events[27] <- nrow(cd8.T.Naive)
    all.props[27] <- (nrow(cd8.T.Naive)/nrow(cd8.Tcells))*100
    
    all.events[28] <- nrow(cd8.T.Memory)
    all.props[28] <- (nrow(cd8.T.Memory)/nrow(cd8.Tcells))*100
    
    
    ## Saving proportions as percentage of Mononuclear cells
    all.freqFeatures[22] <- nrow(cd8.T.Naive)/nrow(mononuclear)
    all.freqFeatures[23] <- nrow(cd8.T.Memory)/nrow(mononuclear)
    
    ## Saving the gating thresholds 
    all.gthres[36] <- cd45RA.gate.cd8
    
   
    min.y <-  min(exprs(cd8.Tcells)[, c('Nd143Di')])
    max.y <- max(exprs(cd8.Tcells)[, c('Nd143Di')])+2
    min.x <- min(exprs(cd8.Tcells)[, c('Nd146Di')])-2
    max.x <-  max(exprs(cd8.Tcells)[, c('Nd146Di')])
    plotDens(cd8.Tcells, c(channels.ind["CD8"],channels.ind["CD45RA"]), main = "CD8+ T Naive & CD8+ T Memory", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(cd8.T.Naive.flowD@filter, lwd = 2)
    lines(cd8.T.Memory.flowD@filter, lwd = 2)
    
    #################################################################################################################
    ## CD25+CD8+ T naive & CD25+CD8+ T memory cells from CD8+ T cells
    
    Tbet.gate.CD8T <- deGate(cd8.Tcells, channel = channels.ind["Tbet"])
    
    cd25.cd8.T.Naive.flowD <- flowDensity(cd8.Tcells, channels = c('Gd160Di', 'Nd143Di'), position = c(T, T), gates = c(Tbet.gate.CD8T, cd45RA.gate.cd8))
    cd25.cd8.T.Memory.flowD <- flowDensity(cd8.Tcells, channels = c('Gd160Di', 'Nd143Di'), position = c(T, F), gates = c(Tbet.gate.CD8T, cd45RA.gate.cd8))
    
    cd25.cd8.T.Naive <- getflowFrame(cd25.cd8.T.Naive.flowD)
    cd25.cd8.T.Memory <- getflowFrame(cd25.cd8.T.Memory.flowD)
    
    all.events[29] <- nrow(cd25.cd8.T.Naive)
    all.props[29] <- (nrow(cd25.cd8.T.Naive)/nrow(cd8.Tcells))*100
    
    all.events[30] <- nrow(cd25.cd8.T.Memory)
    all.props[30] <- (nrow(cd25.cd8.T.Memory)/nrow(cd8.Tcells))*100
    
    ## Saving proportions as percentage of Mononuclear cells
    all.freqFeatures[24] <- nrow(cd25.cd8.T.Naive)/nrow(mononuclear)
    all.freqFeatures[25] <- nrow(cd25.cd8.T.Memory)/nrow(mononuclear)
    
    ## Saving the gating thresholds 
    all.gthres[37] <- Tbet.gate.CD8T
    
    
    min.y <-  min(exprs(cd8.Tcells)[, c('Nd143Di')])
    max.y <- max(exprs(cd8.Tcells)[, c('Nd143Di')])+2
    min.x <- min(exprs(cd8.Tcells)[, c('Gd160Di')])
    max.x <-  max(exprs(cd8.Tcells)[, c('Gd160Di')])+2
    plotDens(cd8.Tcells, c(channels.ind["Tbet"], channels.ind["CD45RA"]), main = "CD25+CD8+ T Naive & CD25+CD8+ T Memory", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(cd25.cd8.T.Naive.flowD@filter, lwd = 2)
    lines(cd25.cd8.T.Memory.flowD@filter, lwd = 2)
    
    
    
    #################################################################################################################
    #################################################################################################################

       
    ## Start Plots-----
    png ( file = paste0(outputPath,"/Figures/Scatter-Plots/", x$FCS.files, ".png"), width=2100, height=2100*4/4)
    
    par(mfrow=c(4,5),mar=(c(5, 5, 4, 2) + 0.1))
    #par(mfrow=c(2,2),mar=(c(5, 5, 4, 2) + 0.1))
    
    ## Singlets 
    plotDens(f, c(pregating.channels["Event_length"], pregating.channels["Ir191Di"]), main = "Singlets", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(singlets.flowD@filter, lwd = 2)
    
    ## Leukocytes
    min.y <-  min(exprs(singlets)[, c('Ir191Di')])-2
    max.y <- max(exprs(singlets)[, c('Ir191Di')])+5
    min.x <- min(exprs(singlets)[, c('In113Di')])
    max.x <- max(exprs(singlets)[, c('In113Di')])+2
    plotDens(singlets, c(channels.ind["CD235ab_CD61"],pregating.channels["Ir191Di"]), main = "Leukocytes", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(leukocytes.flowD@filter, lwd = 2)
    
    ## Mononuclear cells & Granulocytes
    plotDens(leukocytes, c(channels.ind["CD66"],channels.ind["CD45"]), main = "Mononuclear cells & Granulocytes", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(mononuclear.flowD@filter, lwd = 2)
    lines(granulocytes.flowD@filter, lwd = 2)
  
    
    ## B cells & T cells 
    plotDens(mononuclear, c(channels.ind["CD3"],channels.ind["CD19"]), main = "B cells & T cells", cex.lab = 2, cex.axis = 2, cex.main=2); abline(h=cd19.gate, lwd=2); abline(v=cd3.gate, lwd=2)
    lines(Bcells.flowD@filter, lwd = 2)
    lines(Tcells.flowD@filter, lwd = 2)
    
    ## NK cells & lin-
    min.y <-  min(exprs(NK.LinNeg)[, c('Pr141Di')])
    max.y <- max(exprs(NK.LinNeg)[, c('Pr141Di')])+1
    min.x <- min(exprs(NK.LinNeg)[, c('Lu175Di')])
    max.x <-  max(exprs(NK.LinNeg)[, c('Lu175Di')])+2
    plotDens(NK.LinNeg, c(channels.ind["CD14"],channels.ind["CD7"]), main = "NK cells & lin-", cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x, max.x), ylim = c(min.y, max.y))
    lines(NKcells.flowD@filter, lwd = 2)
    lines(LinNeg.flowD@filter, lwd = 2)
    
   
    
    ## CD56low CD16+/CD56+CD16- NKcells
    min.y <-  min(exprs(NKcells)[, c('Ho165Di')])-1
    max.y <- max(exprs(NKcells)[, c('Ho165Di')])+1
    min.x <- min(exprs(NKcells)[, c('Yb176Di')])
    max.x <-  max(exprs(NKcells)[, c('Yb176Di')])+2
    plotDens(NKcells, c(channels.ind["CD56"],channels.ind["CD16"]), main = "CD56low CD16+/CD56+CD16- NKcells", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(cd56low.cd16pos.flowD@filter)
    lines(cd56pos.cd16neg.flowD@filter)
    
    ## cMC, ncMC, intMC, & NOT MC (DC gate)
    min.y <-  min(exprs(LinNeg)[, c('Ho165Di')])-1
    max.y <- max(exprs(LinNeg)[, c('Ho165Di')])+1
    min.x <- min(exprs(LinNeg)[, c('Lu175Di')])
    max.x <-  max(exprs(LinNeg)[, c('Lu175Di')])+1
    plotDens(LinNeg, c(channels.ind["CD14"],channels.ind["CD16"]), main = "cMC, ncMC, & intMC", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(cMC.flowD@filter, lwd = 2)
    lines(ncMC.flowD@filter, lwd = 2)
    lines(intMC.flowD@filter, lwd = 2)
    lines(NOT.MC.flowD@filter, lwd = 2)
    
    
    ## pDCs
    min.y <-  min(exprs(NOT.MC)[, c('Yb174Di')])
    max.y <- max(exprs(NOT.MC)[, c('Yb174Di')])+1
    min.x <- min(exprs(NOT.MC)[, c('Nd148Di')])
    max.x <-  max(exprs(NOT.MC)[, c('Nd148Di')])+2
    plotDens(NOT.MC, c(channels.ind["CD123"],channels.ind["HLADR"]), main = "pDCs", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(pDC.flowD@filter, lwd=2)
    
    
    ## mDCs
    min.y <-  min(exprs(NOT.MC)[, c('Yb174Di')])
    max.y <- max(exprs(NOT.MC)[, c('Yb174Di')])+1
    min.x <- min(exprs(NOT.MC)[, c('Sm147Di')])
    max.x <-  max(exprs(NOT.MC)[, c('Sm147Di')])+2
    plotDens(NOT.MC, c(channels.ind["CD11c"],channels.ind["HLADR"]), main = "mDCs", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(mDCs.flowD@filter, lwd = 2)
    
    ## M-MDSCs
    min.y <-  min(exprs(cMC)[, c('Yb174Di')])
    max.y <- max(exprs(cMC)[, c('Yb174Di')])+1
    min.x <- min(exprs(cMC)[, c('Nd144Di')])
    max.x <-  max(exprs(cMC)[, c('Nd144Di')])+2
    plotDens(cMC, c(channels.ind["CD11b"],channels.ind["HLADR"]), main = "M-MDSCs", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(MMDSC.flowD@filter, lwd =2)
    lines(x=c(cd11b.gate, max.x+1), y=c(hladr.gate.MMDSCs, hladr.gate.MMDSCs), lwd = 2)
    lines(x=c(cd11b.gate, cd11b.gate), y=c(min.y-1, hladr.gate.MMDSCs), lwd = 2)
    
    
    ## CD8+ Tcells / CD4+ T cells
    plotDens(Tcells, c(channels.ind["CD4"],channels.ind["CD8"]), main = "CD8+ Tcells / CD4+ T cells", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(cd4.Tcells.flowD@filter, lwd = 2)
    lines(cd8.Tcells.flowD@filter, lwd = 2)
    lines(NOT.cd4.cd8.Tcells.flowD@filter, lwd = 2)
    
    
    ## CD4+ T naive & CD4+ T memory cells
    plotDens(cd4.Tcells, c(channels.ind["CD4"],channels.ind["CD45RA"]), main = "CD4+ T naive & CD4+ T memory cells", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(cd4.T.Naive.flowD@filter, lwd = 2)
    lines(cd4.T.Memory.flowD@filter, lwd = 2)
    
    
    ## CD4+ Th1 naive & CD4+ Th1 memory cells
    min.y <-  min(exprs(cd4.Tcells)[, c('Nd143Di')])
    max.y <- max(exprs(cd4.Tcells)[, c('Nd143Di')])+1
    min.x <- min(exprs(cd4.Tcells)[, c('Gd160Di')])
    max.x <-  max(exprs(cd4.Tcells)[, c('Gd160Di')])+2
    plotDens(cd4.Tcells, c(channels.ind["Tbet"],channels.ind["CD45RA"]), main = "CD4+ Th1 naive & CD4+ Th1 memory cells", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(cd4.Th1.Naive.flowD@filter, lwd=2)
    lines(cd4.Th1.Memory.flowD@filter, lwd=2)
    
    ## Tregs Naive
    min.y <-  min(exprs(cd4.T.Naive)[, c('Tm169Di')])
    max.y <- max(exprs(cd4.T.Naive)[, c('Tm169Di')])+2
    min.x <- min(exprs(cd4.T.Naive)[, c('Dy162Di')])
    max.x <-  max(exprs(cd4.T.Naive)[, c('Dy162Di')])+2
    plotDens(cd4.T.Naive, c(channels.ind["FoxP3"],channels.ind["CD25"]), main = "Tregs Naive", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(Tregs.Naive.flowD@filter, lwd = 2)
    
    ## Tregs Memory
    min.y <-  min(exprs(cd4.T.Memory)[, c('Tm169Di')])
    max.y <- max(exprs(cd4.T.Memory)[, c('Tm169Di')])+2
    min.x <- min(exprs(cd4.T.Memory)[, c('Dy162Di')])
    max.x <-  max(exprs(cd4.T.Memory)[, c('Dy162Di')])+2
    plotDens(cd4.T.Memory, c(channels.ind["FoxP3"],channels.ind["CD25"]), main = "Tregs Memory", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(Tregs.Memory.flowD@filter, lwd = 2)
    
    
    ## gamma-delta T-cells
    min.y <-  min(exprs(NOT.cd4.cd8.Tcells)[, c('Er170Di')])
    max.y <- max(exprs(NOT.cd4.cd8.Tcells)[, c('Er170Di')])+2
    min.x <- min(exprs(NOT.cd4.cd8.Tcells)[, c('Sm152Di')])
    max.x <-  max(exprs(NOT.cd4.cd8.Tcells)[, c('Sm152Di')])+2
    plotDens(NOT.cd4.cd8.Tcells, c(channels.ind["TCRgd"],channels.ind["CD3"]), main = "gamma-delta T-cells", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(gammaDelta.Tcells.flowD@filter, lwd = 2)
    
    
    ## CD8+ T Naive & CD8+ T Memory
    min.y <-  min(exprs(cd8.Tcells)[, c('Nd143Di')])
    max.y <- max(exprs(cd8.Tcells)[, c('Nd143Di')])+2
    min.x <- min(exprs(cd8.Tcells)[, c('Nd146Di')])-2
    max.x <-  max(exprs(cd8.Tcells)[, c('Nd146Di')])
    plotDens(cd8.Tcells, c(channels.ind["CD8"],channels.ind["CD45RA"]), main = "CD8+ T Naive & CD8+ T Memory", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(cd8.T.Naive.flowD@filter, lwd = 2)
    lines(cd8.T.Memory.flowD@filter, lwd = 2)
    
    ## CD25+CD8+ T Naive & CD25+CD8+ T Memory
    min.y <-  min(exprs(cd8.Tcells)[, c('Nd143Di')])
    max.y <- max(exprs(cd8.Tcells)[, c('Nd143Di')])+2
    min.x <- min(exprs(cd8.Tcells)[, c('Gd160Di')])
    max.x <-  max(exprs(cd8.Tcells)[, c('Gd160Di')])+2
    plotDens(cd8.Tcells, c(channels.ind["Tbet"], channels.ind["CD45RA"]), main = "CD25+CD8+ T Naive & CD25+CD8+ T Memory", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    lines(cd25.cd8.T.Naive.flowD@filter, lwd = 2)
    lines(cd25.cd8.T.Memory.flowD@filter, lwd = 2)
    
    
    
    
    dev.off()
    par(mfrow=c(1,1))
    
    ## End Plots-----
    
  },error = function(err) { return(NULL) } ) # end of tryCatch

  #list(x$FCS.files,all.props, all.events, all.gthres, Filters.list)
  list(x$FCS.files,all.props, all.events, all.gthres, all.freqFeatures)
}, .parallel = TRUE) # end llply

print("Finished Gating & Plotting")


# Saving the big dataframe of Proportions, Events, and Gating thresholds 
save(props.events.gates, file = paste0(results.dir,"/Props.Events.Gates.Rdata"))  

# Extracting the FCS Files names from the large list 
all.FCS.store <- as.matrix(sapply(1:length(props.events.gates), function(x){unlist(props.events.gates[[x]][1])}))

# Extracting the Proportions from the large list & finding the files which failed the gating
all.props.store <- as.matrix(t(sapply(1:length(props.events.gates), function(x){unlist(props.events.gates[[x]][2])})))


failedGating.files.index.temp <- sapply(1:nrow(all.props.store), function(x){which(is.na(all.props.store[x,2:ncol(all.props.store)]))})
failedGating.files.index <- which(sapply(1:length(failedGating.files.index.temp), function(x){length(failedGating.files.index.temp[[x]])})!=0)
failedGating.files <- store.allFCS.Unstim[failedGating.files.index,]
save(failedGating.files.index, file = paste0(results.dir,"/failedGating.files.index.Rdata"))


colnames(all.props.store) <- c("All Events", "Singlets", "Leukocytes", "Mononuclear cells", "Granulocytes", "B cells", "T cells",
                               "NK cells", "lin-", "CD56lowCD16+NK", "CD56+CD16-NK", "cMC", "ncMC", "intMC", "pDC", "mDC", "M-MDSC",
                               "CD4+ T cells", "CD8+ T cells", "CD4+ T Naive", "CD4+ T Memory", "CD4+ Th1 Naive", "CD4+ Th1 Memory", 
                               "Tregs Naive", "Tregs Memory", "gamma-delta T cells", "CD8+ T Naive", "CD8+ T Memory", "CD25+CD8+ T Naive", "CD25+CD8+ T Memory")
all.props.store <- cbind(store.allFCS.Unstim, all.props.store)


# Extracting the Event counts from the large list 
all.events.store <-  as.matrix(t(sapply(1:length(props.events.gates), function(x){unlist(props.events.gates[[x]][3])})))
colnames(all.events.store) <- c("All Events", "Singlets", "Leukocytes", "Mononuclear cells", "Granulocytes", "B cells", "T cells",
                                "NK cells", "lin-", "CD56lowCD16+NK", "CD56+CD16-NK", "cMC", "ncMC", "intMC", "pDC", "mDC", "M-MDSC",
                                "CD4+ T cells", "CD8+ T cells", "CD4+ T Naive", "CD4+ T Memory", "CD4+ Th1 Naive", "CD4+ Th1 Memory", 
                                "Tregs Naive", "Tregs Memory", "gamma-delta T cells", "CD8+ T Naive", "CD8+ T Memory", "CD25+CD8+ T Naive", "CD25+CD8+ T Memory")
all.events.store <- cbind(store.allFCS.Unstim, all.events.store)

# Extracting the Gating Thresholds from the large list 
all.gthres.store <-  as.matrix(t(sapply(1:length(props.events.gates), function(x){unlist(props.events.gates[[x]][4])})))
colnames(all.gthres.store) <- c("fsca.live.gate", "live.gate", "ss.high.gate", "cd5.gate", "cd161.gate", "cd62.gate", "cd44.gate", 
                          "cd8a.gate", "cd161.gate.high", "cd161.gate.low", "cd4.gate", "cd25.gate.low", "cd62.gate.NKT", "cd4.gate.high", "cd25.gate.high",
                          "cd62.gate.Tregs", "cd44.gate.Thelper")
all.gthres.store <- cbind(store.allFCS.Unstim, all.gthres.store)


# Extracting the Proportions based on the Mononuclear cells from the large list 
all.freqFeatures.store <-  as.matrix(t(sapply(1:length(props.events.gates), function(x){unlist(props.events.gates[[x]][5])})))
colnames(all.freqFeatures.store) <- c("B cells", "T cells", "NK cells", "lin-", "CD56lowCD16+NK", "CD56+CD16-NK", "cMC", "ncMC", "intMC", "pDC", "mDC", "M-MDSC",
                                "CD4+ T cells", "CD8+ T cells", "CD4+ T Naive", "CD4+ T Memory", "CD4+ Th1 Naive", "CD4+ Th1 Memory", 
                                "Tregs Naive", "Tregs Memory", "gamma-delta T cells", "CD8+ T Naive", "CD8+ T Memory", "CD25+CD8+ T Naive", "CD25+CD8+ T Memory")
all.freqFeatures.store <- cbind(store.allFCS.Unstim, all.freqFeatures.store)

# Extracting the Filters from the large list 
all.filters.store <-  sapply(1:length(props.events.gates), function(x){props.events.gates[[x]][6]})
names(all.filters.store) <- all.FCS.store



## Saving thresholds and filters in a list (will need them later for flowType)
Gates.Filter.list <- list()
for (i in 1:nrow(all.gthres.store)) {
  Gates.Filter.list[[i]] <- list(all.filters.store[i], all.gthres.store[i,c(12:ncol(all.gthres.store))])  
}


save(all.props.store, file = paste0(results.dir,"/all.props.store.Rdata"))
save(all.events.store, file = paste0(results.dir,"/all.events.store.Rdata"))
save(all.gthres.store, file = paste0(results.dir,"/all.gthres.store.Rdata"))
save(all.filters.store, file = paste0(results.dir,"/all.filters.store.Rdata"))
save(Gates.Filter.list, file = paste0(results.dir,"/Gates.Filter.list.Rdata"))

date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
write.csv(all.props.store, file =  paste0(results.dir, "/DCCResults_Proportions_",toupper(centre),"_Panel1", date.time), row.names = FALSE)

date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
write.csv(all.events.store, file =  paste0(results.dir, "/DCCResults_EventCounts_",toupper(centre),"_Panel1", date.time), row.names = FALSE)

cat("Total time is: ",TimeOutput(start),sep="")

