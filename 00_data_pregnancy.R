## extract albina's results and flowtype
#aya43@sfu.ca 20161220

## root directory
root = "~/projects/flowtype_metrics"
setwd(root)

## input
result_dir = paste0(root, "/result/pregnancy"); dir.create(result_dir, showWarnings=F, recursive=T)
data_dir0 = "/mnt/f/FCS data/Immune_Clock_Human_Pregnancy"
data_dir = paste0(data_dir0, "/Results/FR-FCM-ZY3Q") #training
# data_dir = paste0(data_dir0, "/Results/FR-FCM-ZY3R") #validation

## ouput
temp_fcs_dir = "/mnt/f/Brinkman group/current/Alice/temp_pregnancy/leucocytes"; dir.create(temp_fcs_dir, showWarnings=F, recursive=T) # temporary leucocyte fcs files to make flowworkspace work
temp_gates_dir = "/mnt/f/Brinkman group/current/Alice/temp_pregnancy/gates"; dir.create(temp_gates_dir, showWarnings=F, recursive=T) # temporary leucocyte fcs files to make flowworkspace work



## libraries
source("source/_funcAlice.R")
source("source/pregnancy/helperFunc.R")
libr(c("flowCore", "flowType", "flowDensity", "flowViz",
       "CytoML", "flowWorkspace",
       "pracma", "tools", "MASS", "KernSmooth",
       "colorRamps",
       "foreach", "doMC", "plyr", "stringr"))

## cores
no_cores = 8#detectCores()-1
registerDoMC(no_cores)

## options
widthheight = 200 # plot width/height
overwrite = F


## make meta_file ---------------------------------

load(paste0(data_dir0,"/Attachments/Data.Rda"))
meta_file = rbind(data.frame(patient=featurepatients, class=featuretimes, week=featureweeks, type=rep("train",length(featurepatients))),
                  data.frame(patient=v.featurepatients, class=v.featuretimes, week=v.featureweeks, type=rep("test",length(v.featurepatients))))#class=time
meta_file$id = paste0(meta_file$patient, "_", meta_file$class)
meta_file = meta_file[order(meta_file$id),]




## saving filters ---------------------------------

fcspaths = list.files(temp_fcs_dir, full.names=T)
if (length(fcspaths)==0) {
  fcspaths_ = list.files(data_dir0, recursive=T, full.names=T)
  fcspaths_ = fcspaths_[grepl("LeukocytesRdata",fcspaths_)]
  fs = llply(fcspaths_, function(fcspath_) {
    try({
      f = get(load(fcspath_))
      write.FCS(f, paste0(temp_fcs_dir, "/", gsub(".Rdata",".fcs",fileNames(fcspath_))))
    })
  }, .parallel=T)
  fcspaths = list.files(temp_fcs_dir, full.names=T)
}
fcspaths = sort(fcspaths[grepl("Unstim",fcspaths)])
fs <- read.flowSet(fcspaths)
identifier(fs) = gsub("BL","4",str_extract(gsub("Repeat_", "",fcspaths),"PTLG[0-9]+_[0-9A-Z]+"))

whn = ceiling(sqrt(length(fcspaths)))
wh = widthheight*whn

leukocytes = as(fs,Class="list")[[1]]

gating_channels <- c("CD57", "CD19", "CD4", "CD8", "IgD", "CD11c", "CD16", "CD3", "CD38", "CD27", "CD14", "CXCR5", "CCR7", "CD45RA", "CD20", "CD127", "CD33", "CD28", "CD161", "TCRgd", "CD123", "CD56", "HLADR", "CD25", 
                     "CD235ab_CD61", "CD66", "CD45", "Tbet", "CD7", "FoxP3", "CD11b")
channels.ind <- sort(Find.markers(leukocytes, gating_channels))
channels.name = leukocytes@parameters@data[channels.ind,1]
names(channels.name) = names(channels.ind)





## Mononuclear cells & Granulocytes from Leukocytes-----------------------

fslist = as(fs,Class="list")
gates_dir = paste0(temp_gates_dir, "/gates.Rdata")
if (overwrite | !file.exists(gates_dir)) {
  gates = llply(fslist, function(leukocytes) { 
    all.filters = list()
    all.gthres <- NULL # matrix for saving the gating thresholds
    
    # tryCatch({
    
    ## Mononuclear cells & Granulocytes from Leukocytes-----------------------
    if(length(deGate(leukocytes, channel = channels.ind["CD66"], all.cuts = T)) > 1){
      cd66.gate.low <- min(deGate(leukocytes, channel = channels.ind["CD66"], all.cuts = T))
    }else{
      cd66.gate.low <- deGate(leukocytes, channel = channels.ind["CD66"], use.upper = T, upper = F)
    }
    
    
    
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
    
    # all.events[4] <- nrow(mononuclear)
    # all.props[4] <- (nrow(mononuclear)/nrow(leukocytes))*100
    # 
    # all.events[5] <- nrow(granulocytes)
    # all.props[5] <- (nrow(granulocytes)/nrow(leukocytes))*100
    
    ## Saving the gating thresholds
    all.gthres["cd66.low"] <- cd66.gate.low
    all.gthres["cd66"] <- cd66.gate
    all.gthres["cd45"] <- cd45.gate
    all.filters[["mononuclear"]] = mononuclear.flowD@filter
    colnames(granulocytes.flowD@filter) = colnames(mononuclear.flowD@filter)
    all.filters[["granulocyte"]] = granulocytes.flowD@filter
    
    
    # plotDens(leukocytes, c(channels.ind["CD66"],channels.ind["CD45"]), main = "Mononuclear cells & Granulocytes", cex.lab = 2, cex.axis = 2, cex.main=2)
    # lines(mononuclear.flowD@filter, lwd = 2)
    # lines(granulocytes.flowD@filter, lwd = 2)
    
    
    # plotDens(leukocytes.temp, c(channels.ind["CD66"],channels.ind["CD45"]), main = "", cex.lab = 2, cex.axis = 2, cex.main=2); abline(h=cd45.gate); abline(v=cd66.gate)
    # lines(mononuclear.flowD.temp@filter, lwd=2)
    # lines(granulocytes.flowD@filter, lwd=2)
    
    
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
    
    # all.events[6] <- nrow(Bcells)
    # all.props[6] <- (nrow(Bcells)/nrow(mononuclear))*100
    
    # all.events[7] <- nrow(Tcells)
    # all.props[7] <- (nrow(Tcells)/nrow(mononuclear))*100
    
    # ## Saving proportions as percentage of Mononuclear cells
    # all.freqFeatures[1] <- nrow(Bcells)/nrow(mononuclear)
    # all.freqFeatures[2] <- nrow(Tcells)/nrow(mononuclear)
    
    ## Saving the gating thresholds
    all.gthres["cd3"] <- cd3.gate
    all.gthres["cd19"] <- cd19.gate
    
    # plotDens(mononuclear, c(channels.ind["CD3"],channels.ind["CD19"]), main = "B cells & T cells", cex.lab = 2, cex.axis = 2, cex.main=2); abline(h=cd19.gate, lwd=2); abline(v=cd3.gate, lwd=2)
    # lines(Bcells.flowD@filter, lwd = 2)
    # lines(Tcells.flowD@filter, lwd = 2)
    
    #############################################################################
    ## NK cells & lin- from NK.LinNeg------------------
    
    #cd14.gate.high <- deGate(NK.LinNeg, channel = channels.ind["CD4"], use.percentile = T, percentile = 0.9999)
    cd14.gate.high <- deGate(NK.LinNeg, channel = channels.ind["CD4"], use.percentile = T, percentile = 0.0001)
    cd7.gate <- deGate(NK.LinNeg, channel = channels.ind["CD7"])
    
    NKcells.flowD <- flowDensity(NK.LinNeg, channels = c('Lu175Di', 'Pr141Di'), position = c(T, T), gates = c(cd14.gate.high, cd7.gate))
    LinNeg.flowD <- flowDensity(NK.LinNeg, channels = c('Lu175Di', 'Pr141Di'), position = c(T, F), gates = c(cd14.gate.high, cd7.gate))
    
    NKcells <- getflowFrame(NKcells.flowD)
    LinNeg <- getflowFrame(LinNeg.flowD)
    
    # all.events[8] <- nrow(NKcells)
    # all.props[8] <- (nrow(NKcells)/nrow(NK.LinNeg))*100
    
    # all.events[9] <- nrow(LinNeg)
    # all.props[9] <- (nrow(LinNeg)/nrow(NK.LinNeg))*100
    
    # ## Saving proportions as percentage of Mononuclear cells
    # all.freqFeatures[3] <- nrow(NKcells)/nrow(mononuclear)
    # all.freqFeatures[4] <- nrow(LinNeg)/nrow(mononuclear)
    
    ## Saving the gating thresholds
    all.gthres["cd14.high"] <- cd14.gate.high
    all.gthres["cd7"] <- cd7.gate
    
    
    
    # min.y <-  min(exprs(NK.LinNeg)[, c('Pr141Di')])
    # max.y <- max(exprs(NK.LinNeg)[, c('Pr141Di')])+1
    # min.x <- min(exprs(NK.LinNeg)[, c('Lu175Di')])
    # max.x <-  max(exprs(NK.LinNeg)[, c('Lu175Di')])+2
    # plotDens(NK.LinNeg, c(channels.ind["CD14"],channels.ind["CD7"]), main = "NK cells & lin-", cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x, max.x), ylim = c(min.y, max.y))
    # lines(NKcells.flowD@filter, lwd = 2)
    # lines(LinNeg.flowD@filter, lwd = 2)
    
    
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
    
    # all.events[10] <- nrow(cd56low.cd16pos)
    # all.props[10] <- (nrow(cd56low.cd16pos)/nrow(NKcells))*100
    # 
    # all.events[11] <- nrow(cd56pos.cd16neg)
    # all.props[11] <- (nrow(cd56pos.cd16neg)/nrow(NKcells))*100
    # 
    # ## Saving proportions as percentage of Mononuclear cells
    # all.freqFeatures[5] <- nrow(cd56low.cd16pos)/nrow(mononuclear)
    # all.freqFeatures[6] <- nrow(cd56pos.cd16neg)/nrow(mononuclear)
    
    ## Saving the gating thresholds
    all.gthres["cd56.gate.low"] <- cd56.gate.low
    all.gthres["cd56.gate.mid"] <- cd56.gate.mid
    all.gthres["cd56.gate.high"] <- cd56.gate.high
    all.gthres["cd16.gate.low"] <- cd16.gate.low
    all.gthres["cd16.gate.mid"] <- cd16.gate.mid
    all.gthres["cd16.gate.high"] <- cd16.gate.high
    all.filters[["cd56low.cd16pos"]] = cd56low.cd16pos.flowD@filter
    all.filters[["cd56pos.cd16neg"]] = cd56pos.cd16neg.flowD@filter #might contain no cells
    
    
    
    # min.y <-  min(exprs(NKcells)[, c('Ho165Di')])-1
    # max.y <- max(exprs(NKcells)[, c('Ho165Di')])+1
    # min.x <- min(exprs(NKcells)[, c('Yb176Di')])
    # max.x <-  max(exprs(NKcells)[, c('Yb176Di')])+2
    # plotDens(NKcells, c(channels.ind["CD56"],channels.ind["CD16"]), main = "CD56low CD16+/CD56+CD16- NKcells", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    # lines(cd56low.cd16pos.flowD@filter)
    # lines(cd56pos.cd16neg.flowD@filter)
    
    
    
    
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
    
    
    # all.events[12] <- nrow(cMC)
    # all.props[12] <- (nrow(cMC)/nrow(LinNeg))*100
    # 
    # all.events[13] <- nrow(ncMC)
    # all.props[13] <- (nrow(ncMC)/nrow(LinNeg))*100
    # 
    # all.events[14] <- nrow(intMC)
    # all.props[14] <- (nrow(intMC)/nrow(LinNeg))*100
    # 
    # ## Saving proportions as percentage of Mononuclear cells
    # all.freqFeatures[7] <- nrow(cMC)/nrow(mononuclear)
    # all.freqFeatures[8] <- nrow(ncMC)/nrow(mononuclear)
    # all.freqFeatures[9] <- nrow(intMC)/nrow(mononuclear)
    
    ## Saving the gating thresholds
    all.gthres["cd14.low"] <- cd14.gate.low
    all.gthres["cd16.LinNeg"] <- cd16.gate.LinNeg
    
    # min.y <-  min(exprs(LinNeg)[, c('Ho165Di')])-1
    # max.y <- max(exprs(LinNeg)[, c('Ho165Di')])+1
    # min.x <- min(exprs(LinNeg)[, c('Lu175Di')])
    # max.x <-  max(exprs(LinNeg)[, c('Lu175Di')])+1
    # plotDens(LinNeg, c(channels.ind["CD14"],channels.ind["CD16"]), main = "cMC, ncMC, & intMC", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    # lines(cMC.flowD@filter, lwd = 2)
    # lines(ncMC.flowD@filter, lwd = 2)
    # lines(intMC.flowD@filter, lwd = 2)
    # lines(NOT.MC.flowD@filter, lwd = 2)
    # lines(cd14pos.flowD@filter, lwd =2)
    
    
    ##############################################################################################################
    ## pDCs from the NOT MC cells
    
    cd123.gate <- deGate(NOT.MC, channel = channels.ind["CD123"])
    hladr.gate.DCs <- deGate(NOT.MC, channel = channels.ind["HLADR"], use.upper = T, upper = F, tinypeak.removal = 0.9)
    
    pDC.flowD <- flowDensity(NOT.MC, channels = c('Nd148Di', 'Yb174Di'), position = c(T, T), gates = c(cd123.gate, hladr.gate.DCs))
    
    pDC <- getflowFrame(pDC.flowD)
    
    # all.events[15] <- nrow(pDC)
    # all.props[15] <- (nrow(pDC)/nrow(NOT.MC))*100
    # 
    # ## Saving proportions as percentage of Mononuclear cells
    # all.freqFeatures[10] <- nrow(pDC)/nrow(mononuclear)
    
    
    ## Saving the gating thresholds
    all.gthres["cd123"] <- cd123.gate
    all.gthres["hladr.DCs"] <- hladr.gate.DCs
    all.filters[["pDC"]] = pDC.flowD@filter
    
    # min.y <-  min(exprs(NOT.MC)[, c('Yb174Di')])
    # max.y <- max(exprs(NOT.MC)[, c('Yb174Di')])+1
    # min.x <- min(exprs(NOT.MC)[, c('Nd148Di')])
    # max.x <-  max(exprs(NOT.MC)[, c('Nd148Di')])+2
    # plotDens(NOT.MC, c(channels.ind["CD123"],channels.ind["HLADR"]), main = "pDCs", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    # lines(pDC.flowD@filter, lwd=2)
    
    
    ##############################################################################################################
    # mDCs from the NOT MC cells
    
    # cd11c.gate <- deGate(NOT.MC, channel = channels.ind["CD11c"], all.cuts = T)
    # if(length(cd11c.gate) > 1){
    #   cd11c.gate <- cd11c.gate[2]
    # }
    
    cd11c.gate <- deGate(NOT.MC, channel = channels.ind["CD11c"], use.percentile = T, percentile = 0.1)
    
    mDCs.flowD <- flowDensity(NOT.MC, channels = c('Sm147Di', 'Yb174Di'), position = c(T,T), gates = c(cd11c.gate, hladr.gate.DCs))
    
    mDC <- getflowFrame(mDCs.flowD)
    
    # all.events[16] <- nrow(mDC)
    # all.props[16] <- (nrow(mDC)/nrow(NOT.MC))*100
    # 
    # ## Saving proportions as percentage of Mononuclear cells
    # all.freqFeatures[11] <- nrow(mDC)/nrow(mononuclear)
    
    ## Saving the gating thresholds
    all.gthres["cd11c"] <- cd11c.gate
    all.filters[["mDCs"]] = mDCs.flowD@filter
    
    # min.y <-  min(exprs(NOT.MC)[, c('Yb174Di')])
    # max.y <- max(exprs(NOT.MC)[, c('Yb174Di')])+1
    # min.x <- min(exprs(NOT.MC)[, c('Sm147Di')])
    # max.x <-  max(exprs(NOT.MC)[, c('Sm147Di')])+2
    # plotDens(NOT.MC, c(channels.ind["CD11c"],channels.ind["HLADR"]), main = "mDCs", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    # lines(mDCs.flowD@filter, lwd = 2)
    
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
    
    # all.events[17] <- nrow(MMDSC)
    # all.props[17] <- (nrow(MMDSC)/nrow(cMC))*100
    # 
    # ## Saving proportions as percentage of Mononuclear cells
    # all.freqFeatures[12] <- nrow(MMDSC)/nrow(mononuclear)
    
    ## Saving the gating thresholds
    all.gthres["cd11b"] <- cd11b.gate
    all.gthres["hladr.MMDSCs"] <- hladr.gate.MMDSCs
    all.filters[["MMDSC"]] = MMDSC.flowD@filter
    
    # min.y <-  min(exprs(cMC)[, c('Yb174Di')])
    # max.y <- max(exprs(cMC)[, c('Yb174Di')])+1
    # min.x <- min(exprs(cMC)[, c('Nd144Di')])
    # max.x <-  max(exprs(cMC)[, c('Nd144Di')])+2
    # plotDens(cMC, c(channels.ind["CD11b"],channels.ind["HLADR"]), main = "M-MDSCs", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    # lines(MMDSC.flowD@filter, lwd =2)
    # lines(x=c(cd11b.gate, max.x+1), y=c(hladr.gate.MMDSCs, hladr.gate.MMDSCs), lwd = 2)
    # lines(x=c(cd11b.gate, cd11b.gate), y=c(min.y-1, hladr.gate.MMDSCs), lwd = 2)
    
    
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
    
    # all.events[18] <- nrow(cd4.Tcells)
    # all.props[18] <- (nrow(cd4.Tcells)/nrow(Tcells))*100
    # 
    # all.events[19] <- nrow(cd8.Tcells)
    # all.props[19] <- (nrow(cd8.Tcells)/nrow(Tcells))*100
    # 
    # ## Saving proportions as percentage of Mononuclear cells
    # all.freqFeatures[13] <- nrow(cd4.Tcells)/nrow(mononuclear)
    # all.freqFeatures[14] <- nrow(cd8.Tcells)/nrow(mononuclear)
    
    
    ## Saving the gating thresholds
    all.gthres["cd4"] <- cd4.gate
    all.gthres["cd8a"] <- cd8a.gate
    all.filters[["cd4t"]] = cd4.Tcells.flowD@filter
    all.filters[["cd8t"]] = cd8.Tcells.flowD@filter
    all.filters[["notcd4cd8t"]] = NOT.cd4.cd8.Tcells.flowD@filter
    
    
    # plotDens(Tcells, c(channels.ind["CD4"],channels.ind["CD8"]), main = "CD8+ Tcells / CD4+ T cells", cex.lab = 2, cex.axis = 2, cex.main=2)
    # lines(cd4.Tcells.flowD@filter, lwd = 2)
    # lines(cd8.Tcells.flowD@filter, lwd = 2)
    # lines(NOT.cd4.cd8.Tcells.flowD@filter, lwd = 2)
    
    
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
    
    # all.events[20] <- nrow(cd4.T.Naive)
    # all.props[20] <- (nrow(cd4.T.Naive)/nrow(cd4.Tcells))*100
    # 
    # all.events[21] <- nrow(cd4.T.Memory)
    # all.props[21] <- (nrow(cd4.T.Memory)/nrow(cd4.Tcells))*100
    # 
    # ## Saving proportions as percentage of Mononuclear cells
    # all.freqFeatures[15] <- nrow(cd4.T.Naive)/nrow(mononuclear)
    # all.freqFeatures[16] <- nrow(cd4.T.Memory)/nrow(mononuclear)
    
    ## Saving the gating thresholds
    all.gthres["cd45RA.high.cd4"] <- cd45RA.gate.high.cd4
    all.filters[["cd4.T.Naive"]] = cd4.T.Naive.flowD@filter
    all.filters[["cd4.T.Memory"]] = cd4.T.Memory.flowD@filter
    
    # plotDens(cd4.Tcells, c(channels.ind["CD4"],channels.ind["CD45RA"]), main = "CD4+ T naive & CD4+ T memory cells", cex.lab = 2, cex.axis = 2, cex.main=2)
    # lines(cd4.T.Naive.flowD@filter, lwd = 2)
    # lines(cd4.T.Memory.flowD@filter, lwd = 2)
    
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
    
    # all.events[22] <- nrow(cd4.Th1.Naive)
    # all.props[22] <- (nrow(cd4.Th1.Naive)/nrow(cd4.Tcells))*100
    # 
    # all.events[23] <- nrow(cd4.Th1.Memory)
    # all.props[23] <- (nrow(cd4.Th1.Memory)/nrow(cd4.Tcells))*100
    # 
    # ## Saving proportions as percentage of Mononuclear cells
    # all.freqFeatures[17] <- nrow(cd4.Th1.Naive)/nrow(mononuclear)
    # all.freqFeatures[18] <- nrow(cd4.Th1.Memory)/nrow(mononuclear)
    
    ## Saving the gating thresholds
    all.gthres["tbet.cd4t"] <- Tbet.gate.CD4T
    
    
    # min.y <-  min(exprs(cd4.Tcells)[, c('Nd143Di')])
    # max.y <- max(exprs(cd4.Tcells)[, c('Nd143Di')])+1
    # min.x <- min(exprs(cd4.Tcells)[, c('Gd160Di')])
    # max.x <- max(exprs(cd4.Tcells)[, c('Gd160Di')])+2
    # plotDens(cd4.Tcells, c(channels.ind["Tbet"],channels.ind["CD45RA"]), main = "CD4+ Th1 naive & CD4+ Th1 memory cells", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    # lines(cd4.Th1.Naive.flowD@filter, lwd=2)
    # lines(cd4.Th1.Memory.flowD@filter, lwd=2)
    
    
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
    
    # all.events[24] <- nrow(Tregs.Naive)
    # all.props[24] <- (nrow(Tregs.Naive)/nrow(cd4.T.Naive))*100
    # 
    # ## Saving proportions as percentage of Mononuclear cells
    # all.freqFeatures[19] <- nrow(Tregs.Naive)/nrow(mononuclear)
    
    
    ## Saving the gating thresholds (rotated gates)
    all.gthres["foxp3.low.naive"] <- FoxP3.gate.Naive.low
    all.gthres["foxp3.high.naive"] <- FoxP3.gate.Naive.high
    all.gthres["cd25.naive"] <- cd25.gate.Naive
    all.filters[["tregs.naive"]] = Tregs.Naive.flowD@filter
    colnames(all.filters[["tregs.naive"]]) = c("Dy162Di", "Tm169Di")
    
    # min.y <-  min(exprs(cd4.T.Naive)[, c('Tm169Di')])
    # max.y <- max(exprs(cd4.T.Naive)[, c('Tm169Di')])+2
    # min.x <- min(exprs(cd4.T.Naive)[, c('Dy162Di')])
    # max.x <-  max(exprs(cd4.T.Naive)[, c('Dy162Di')])+2
    # plotDens(cd4.T.Naive, c(channels.ind["FoxP3"],channels.ind["CD25"]), main = "Tregs Naive", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    # # abline(h=cd25.gate.Naive); abline(v=FoxP3.gate.Naive)
    # lines(Tregs.Naive.flowD@filter, lwd = 2)
    
    
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
    if (is.null(FoxP3.gate.Memory.high))
      FoxP3.gate.Memory.high = deGate(cd4.T.Memory.temp, channel = channels.ind["FoxP3"], use.upper = T, upper = T)
    cd25.gate.Memory <- deGate(cd4.T.Memory.temp, channel = channels.ind["CD25"], use.percentile = T, percentile = 0.95)
    
    Tregs.Memory.flowD.temp <- flowDensity(cd4.T.Memory.temp, channels = c('Dy162Di', 'Tm169Di'), position = c(T,T), gates = c(FoxP3.gate.Memory.low, cd25.gate.Memory))
    Tregs.Memory.flowD <- flowDensity(Tregs.Memory.flowD.temp, channels = c('Dy162Di', 'Tm169Di'), position = c(F,NA), gates = c(FoxP3.gate.Memory.high, NA))
    
    Tregs.Memory <- getflowFrame(Tregs.Memory.flowD)
    
    Tregs.Memory@exprs[,c(channels.ind["FoxP3"],channels.ind["CD25"])] <- t(t(R) %*% t(Tregs.Memory@exprs[,c(channels.ind["FoxP3"],channels.ind["CD25"])]))
    
    Tregs.Memory.flowD@filter <- rotate.data(Tregs.Memory.flowD@filter,c(channels.ind["FoxP3"],channels.ind["CD25"]),theta = pi/4)$data
    Tregs.Memory.flowD@flow.frame <- rotate.data(getflowFrame(Tregs.Memory.flowD),c(channels.ind["FoxP3"],channels.ind["CD25"]),theta = pi/4)$data
    Tregs.Memory.flowD@proportion <- (Tregs.Memory.flowD@cell.count/nrow(cd4.T.Memory))*100
    
    Tregs.Memory <- getflowFrame(Tregs.Memory.flowD)
    
    # all.events[25] <- nrow(Tregs.Memory)
    # all.props[25] <- (nrow(Tregs.Memory)/nrow(cd4.T.Memory))*100
    # 
    # ## Saving proportions as percentage of Mononuclear cells
    # all.freqFeatures[20] <- nrow(Tregs.Memory)/nrow(mononuclear)
    
    ## Saving the gating thresholds (rotated gates)
    all.gthres["foxp3.low.memory"] <- FoxP3.gate.Memory.low
    all.gthres["foxp3.high.memory"] <- FoxP3.gate.Memory.high
    all.gthres["cd25.memory"] <- cd25.gate.Memory
    all.filters[["tregs.memory"]] = Tregs.Memory.flowD@filter
    colnames(all.filters[["tregs.memory"]]) = c(channels.name["FoxP3"],channels.name["CD25"])
    
    
    # min.y <-  min(exprs(cd4.T.Memory)[, c('Tm169Di')])
    # max.y <- max(exprs(cd4.T.Memory)[, c('Tm169Di')])+2
    # min.x <- min(exprs(cd4.T.Memory)[, c('Dy162Di')])
    # max.x <-  max(exprs(cd4.T.Memory)[, c('Dy162Di')])+2
    # plotDens(cd4.T.Memory, c(channels.ind["FoxP3"],channels.ind["CD25"]), main = "Tregs Memory", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    # # abline(h=cd25.gate); abline(v=FoxP3.gate)
    # lines(Tregs.Memory.flowD@filter, lwd = 2)
    
    #############################################################################################
    ## Gamma-Delta T-cells from NOT.cd4.cd8.Tcells.flowD
    
    TCRD.gate <- deGate(NOT.cd4.cd8.Tcells, channel = channels.ind["TCRgd"])
    
    gammaDelta.Tcells.flowD <- flowDensity(NOT.cd4.cd8.Tcells, channels = c('Sm152Di', 'Er170Di'), position = c(T,NA), gates = c(TCRD.gate+0.5,NA), ellip.gate=T)
    
    gammaDelta.Tcells <- getflowFrame(gammaDelta.Tcells.flowD)
    
    # all.events[26] <- nrow(gammaDelta.Tcells)
    # all.props[26] <- (nrow(gammaDelta.Tcells)/nrow(NOT.cd4.cd8.Tcells))*100
    # 
    # ## Saving proportions as percentage of Mononuclear cells
    # all.freqFeatures[21] <- nrow(gammaDelta.Tcells)/nrow(mononuclear)
    
    ## Saving the gating thresholds 
    all.gthres["tcrd"] <- TCRD.gate
    all.filters[["gammadelta.tcell"]] = gammaDelta.Tcells.flowD@filter
    
    
    # min.y <-  min(exprs(NOT.cd4.cd8.Tcells)[, c('Er170Di')])
    # max.y <- max(exprs(NOT.cd4.cd8.Tcells)[, c('Er170Di')])+2
    # min.x <- min(exprs(NOT.cd4.cd8.Tcells)[, c('Sm152Di')])
    # max.x <-  max(exprs(NOT.cd4.cd8.Tcells)[, c('Sm152Di')])+2
    # plotDens(NOT.cd4.cd8.Tcells, c(channels.ind["TCRgd"],channels.ind["CD3"]), main = "gamma-delta T-cells", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    # lines(gammaDelta.Tcells.flowD@filter, lwd = 2)
    
    
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
    
    
    # all.events[27] <- nrow(cd8.T.Naive)
    # all.props[27] <- (nrow(cd8.T.Naive)/nrow(cd8.Tcells))*100
    # 
    # all.events[28] <- nrow(cd8.T.Memory)
    # all.props[28] <- (nrow(cd8.T.Memory)/nrow(cd8.Tcells))*100
    # 
    # 
    # ## Saving proportions as percentage of Mononuclear cells
    # all.freqFeatures[22] <- nrow(cd8.T.Naive)/nrow(mononuclear)
    # all.freqFeatures[23] <- nrow(cd8.T.Memory)/nrow(mononuclear)
    
    ## Saving the gating thresholds 
    all.gthres["cd45ra.cd8"] <- cd45RA.gate.cd8
    all.filters[["cd8.T.Naive"]] = cd8.T.Naive.flowD@filter
    all.filters[["cd8.T.Memory"]] = cd8.T.Memory.flowD@filter
    
    
    # min.y <-  min(exprs(cd8.Tcells)[, c('Nd143Di')])
    # max.y <- max(exprs(cd8.Tcells)[, c('Nd143Di')])+2
    # min.x <- min(exprs(cd8.Tcells)[, c('Nd146Di')])-2
    # max.x <-  max(exprs(cd8.Tcells)[, c('Nd146Di')])
    # plotDens(cd8.Tcells, c(channels.ind["CD8"],channels.ind["CD45RA"]), main = "CD8+ T Naive & CD8+ T Memory", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    # lines(cd8.T.Naive.flowD@filter, lwd = 2)
    # lines(cd8.T.Memory.flowD@filter, lwd = 2)
    
    #################################################################################################################
    ## CD25+CD8+ T naive & CD25+CD8+ T memory cells from CD8+ T cells
    
    Tbet.gate.CD8T <- deGate(cd8.Tcells, channel = channels.ind["Tbet"])
    
    cd25.cd8.T.Naive.flowD <- flowDensity(cd8.Tcells, channels = c('Gd160Di', 'Nd143Di'), position = c(T, T), gates = c(Tbet.gate.CD8T, cd45RA.gate.cd8))
    cd25.cd8.T.Memory.flowD <- flowDensity(cd8.Tcells, channels = c('Gd160Di', 'Nd143Di'), position = c(T, F), gates = c(Tbet.gate.CD8T, cd45RA.gate.cd8))
    
    cd25.cd8.T.Naive <- getflowFrame(cd25.cd8.T.Naive.flowD)
    cd25.cd8.T.Memory <- getflowFrame(cd25.cd8.T.Memory.flowD)
    
    # all.events[29] <- nrow(cd25.cd8.T.Naive)
    # all.props[29] <- (nrow(cd25.cd8.T.Naive)/nrow(cd8.Tcells))*100
    # 
    # all.events[30] <- nrow(cd25.cd8.T.Memory)
    # all.props[30] <- (nrow(cd25.cd8.T.Memory)/nrow(cd8.Tcells))*100
    # 
    # ## Saving proportions as percentage of Mononuclear cells
    # all.freqFeatures[24] <- nrow(cd25.cd8.T.Naive)/nrow(mononuclear)
    # all.freqFeatures[25] <- nrow(cd25.cd8.T.Memory)/nrow(mononuclear)
    
    ## Saving the gating thresholds 
    all.gthres["tbet.cd8t"] <- Tbet.gate.CD8T
    # all.gthres["cd45RA.cd8"] <- cd45RA.gate.cd8
    
    
    # min.y <-  min(exprs(cd8.Tcells)[, c('Nd143Di')])
    # max.y <- max(exprs(cd8.Tcells)[, c('Nd143Di')])+2
    # min.x <- min(exprs(cd8.Tcells)[, c('Gd160Di')])
    # max.x <-  max(exprs(cd8.Tcells)[, c('Gd160Di')])+2
    # plotDens(cd8.Tcells, c(channels.ind["Tbet"], channels.ind["CD45RA"]), main = "CD25+CD8+ T Naive & CD25+CD8+ T Memory", cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y), xlim = c(min.x, max.x))
    # lines(cd25.cd8.T.Naive.flowD@filter, lwd = 2)
    # lines(cd25.cd8.T.Memory.flowD@filter, lwd = 2)
    
    
    
    #################################################################################################################
    #################################################################################################################
    
    return(list(all.gthres=all.gthres, all.filters=all.filters))
    
    
    # },error = function(err) { return(NULL) } ) # end of tryCatch
  }, .parallel=T)  
  save(gates, file=gates_dir)
} else {
  gates = get(load(gates_dir))
}

gthres = llply(gates, function(x) x[[1]])
filters = llply(gates, function(x) x[[2]])

# initialize gating set and plotting parameters
gs <- GatingSet(fs)


# leucocyte: Mononuclear cells & Granulocytes
mf <- llply(filters, function(x) polygonGate(filterId="mononuclear", .gate=x$mononuclear))
gf <- llply(filters, function(x) polygonGate(filterId="granulocyte", .gate=x$granulocyte))
nodeIDs <- flowWorkspace::add(gs, mf)
nodeIDs <- flowWorkspace::add(gs, gf)
recompute(gs)

png(paste0(temp_gates_dir, "/01_leucocyte_mononuclearRED-granulocyteBLACK.png"), width=wh, height=wh)
par(mfrow=c(whn,whn),mar=c(5,5,5,5))
fslist = as(getData(gs), Class='list')
for (i in 1:length(fslist)) {
  f = fslist[[i]]
  plotDens(f, c(channels.ind["CD66"],channels.ind["CD45"]), main=meta_file$id[i], cex.lab = 2, cex.axis = 2, cex.main=2)
  lines(filters[[i]]$mononuclear, lwd=2, col="red")
  lines(filters[[i]]$granulocyte, lwd=2, col="black")
}
graphics.off()



## mononuclear: B cells & T cells 
btf = llply(gthres, function(x) return(quadGate("Er170Di"=x["cd3"], "Nd142Di"=x["cd19"])))
nodeIDs <- flowWorkspace::add(gs, btf, parent="mononuclear")
# bcell = CD3-CD19+
# tcell = CD3+CD19-
# nkLinNeg = CD3-CD19-
recompute(gs)

png(paste0(temp_gates_dir, "/02_mononuclear_bcellNP-tcellPN-nklinnegNN.png"), width=wh, height=wh)
par(mfrow=c(whn,whn),mar=c(5,5,5,5))
fslist = as(getData(gs, 'mononuclear'), Class='list')
for (i in 1:length(fslist)) {
  f = fslist[[i]]
  plotDens(f, c(channels.ind["CD3"],channels.ind["CD19"]), main=meta_file$id[i], cex.lab = 2, cex.axis = 2, cex.main=2)
  abline(h=gthres[[i]]["cd19"], lwd=2)
  abline(v=gthres[[i]]["cd3"], lwd=2)
}
graphics.off()


## nkLinNeg CD3-CD19-: NK cells & lin-
nkf = llply(gthres, function(x) return(quadGate("Lu175Di"=x["cd14.high"], "Pr141Di"=x["cd7"])))
nodeIDs <- flowWorkspace::add(gs, nkf, parent="CD3-CD19-")
# nk = CD14+CD7+
# linNeg = CD14+CD7-
recompute(gs)

png(paste0(temp_gates_dir, "/03_nklinneg_nkPP-linnegPN.png"), width=wh, height=wh)
par(mfrow=c(whn,whn),mar=c(5,5,5,5))
fslist = as(getData(gs, 'CD3-CD19-'), Class='list')
for (i in 1:length(fslist)) {
  f = fslist[[i]]
  plotDens(f, c(channels.ind["CD14"],channels.ind["CD7"]), main=meta_file$id[i], cex.lab = 2, cex.axis = 2, cex.main=2)
  abline(h=gthres[[i]]["cd14.high"], lwd=2)
  abline(v=gthres[[i]]["cd7"], lwd=2)
}
graphics.off()


## nk CD14+CD7+: CD56low CD16+/CD56+CD16- NKcells
cd56pcd16nf = llply(filters, function(x) return(polygonGate(filterId="cd56pos.cd16neg", .gate=x$cd56pos.cd16neg)))
cd45lcd16pf = llply(filters, function(x) return(polygonGate(filterId="cd56low.cd16pos", .gate=x$cd56low.cd16pos)))
nodeIDs <- flowWorkspace::add(gs, cd56pcd16nf, parent="CD14+CD7+")
nodeIDs <- flowWorkspace::add(gs, cd45lcd16pf, parent="CD14+CD7+")
recompute(gs)

png(paste0(temp_gates_dir, "/04_nk_NPred-PNblack.png"), width=wh, height=wh)
par(mfrow=c(whn,whn),mar=c(5,5,5,5))
fslist = as(getData(gs, 'CD14+CD7+'), Class='list')
for (i in 1:length(fslist)) {
  f = fslist[[i]]
  plotDens(f, c(channels.ind["CD56"],channels.ind["CD16"]), main=meta_file$id[i], cex.lab = 2, cex.axis = 2, cex.main=2)
  lines(filters[[i]]$cd56low.cd16pos, lwd=2, col="red")
  lines(filters[[i]]$cd56pos.cd16neg, lwd=2, col="black")
}
graphics.off()



## linNeg CD14+CD7-: cMC, ncMC, intMC, & NOT MC (DC gate)
mcf = llply(gthres, function(x) return(quadGate("Lu175Di"=x["cd14.low"], "Ho165Di"=x["cd16.LinNeg"])))
nodeIDs <- flowWorkspace::add(gs, mcf, parent="CD14+CD7-")
# cMC = CD14+CD16-
# ncMC = CD14-CD16+
# intMC = CD14+CD16+
# NOT.MC = CD14-CD16-
recompute(gs)

png(paste0(temp_gates_dir, "/05_linnegCD14+CD7-_cMCPN-ncMCNP-intMCPP-notMCNN.png"), width=wh, height=wh)
par(mfrow=c(whn,whn),mar=c(5,5,5,5))
fslist = as(getData(gs, 'CD14+CD7-'), Class='list')
for (i in 1:length(fslist)) {
  f = fslist[[i]]
  plotDens(f, c(channels.ind["CD14"],channels.ind["CD16"]), main=meta_file$id[i], cex.lab = 2, cex.axis = 2, cex.main=2)
  abline(h=gthres[[i]]["cd14.low"], lwd=2)
  abline(v=gthres[[i]]["cd16.LinNeg"], lwd=2)
}
graphics.off()


## notMC CD14-CD16-: M-MDSCs
pdcf = llply(filters, function(x) return(polygonGate(filterId="pDC", .gate=x$pDC)))
nodeIDs <- flowWorkspace::add(gs, pdcf, parent="CD14-CD16-")
recompute(gs)

png(paste0(temp_gates_dir, "/06_notMC.CD14-CD16-_pDCPP.png"), width=wh, height=wh)
par(mfrow=c(whn,whn),mar=c(5,5,5,5))
fslist = as(getData(gs, 'CD14-CD16-'), Class='list')
for (i in 1:length(fslist)) {
  f = fslist[[i]]
  plotDens(f, c(channels.ind["CD123"],channels.ind["HLADR"]), main=meta_file$id[i], cex.lab = 2, cex.axis = 2, cex.main=2)
  abline(h=gthres[[i]]["cd123"], lwd=2)
  abline(v=gthres[[i]]["hladr.DCs"], lwd=2)
}
graphics.off()


mdcf = llply(filters, function(x) return(polygonGate(filterId="mDC", .gate=x$mDCs)))
nodeIDs <- flowWorkspace::add(gs, mdcf, parent="CD14-CD16-")
recompute(gs)

png(paste0(temp_gates_dir, "/07_notMC.CD14-CD16-_mDCPP.png"), width=wh, height=wh)
par(mfrow=c(whn,whn),mar=c(5,5,5,5))
fslist = as(getData(gs, 'CD14-CD16-'), Class='list')
for (i in 1:length(fslist)) {
  f = fslist[[i]]
  plotDens(f, c(channels.ind["CD11c"],channels.ind["HLADR"]), main=meta_file$id[i], cex.lab = 2, cex.axis = 2, cex.main=2)
  abline(h=gthres[[i]]["cd11c"], lwd=2)
  abline(v=gthres[[i]]["hladr.DCs"], lwd=2)
}
graphics.off()


## cMC CD14+CD16-: mDCs
mmdscf = llply(filters, function(x) return(polygonGate(filterId="mmdsc", .gate=x$MMDSC))) # can do quad, but nah
nodeIDs <- flowWorkspace::add(gs, mmdscf, parent="CD14+CD16-")
recompute(gs)

png(paste0(temp_gates_dir, "/08_cMC.CD14+CD16-_mmdmcPN.png"), width=wh, height=wh)
par(mfrow=c(whn,whn),mar=c(5,5,5,5))
fslist = as(getData(gs, 'CD14+CD16-'), Class='list')
for (i in 1:length(fslist)) {
  f = fslist[[i]]
  plotDens(f, c(channels.ind["CD11b"],channels.ind["HLADR"]), main=meta_file$id[i], cex.lab = 2, cex.axis = 2, cex.main=2)
  abline(h=gthres[[i]]["cd11b"], lwd=2)
  abline(v=gthres[[i]]["hladr.MMDSCs"], lwd=2)
}
graphics.off()


## tcell CD3+CD19-: CD8+ Tcells / CD4+ T cells 
cd4tf = llply(filters, function(x) return(polygonGate(filterId="cd4tcell", .gate=x$cd4t)))
cd8tf = llply(filters, function(x) return(polygonGate(filterId="cd8tcell", .gate=x$cd8t)))
notcd4cd8tf = llply(filters, function(x) return(polygonGate(filterId="notcd4cd8tcell", .gate=x$notcd4cd8t)))
nodeIDs <- flowWorkspace::add(gs, cd4tf, parent="CD3+CD19-")
nodeIDs <- flowWorkspace::add(gs, cd8tf, parent="CD3+CD19-")
nodeIDs <- flowWorkspace::add(gs, notcd4cd8tf, parent="CD3+CD19-")
recompute(gs)

png(paste0(temp_gates_dir, "/09_tcellCD3+CD19-_cd4tPN-cd8tNP-cotcd4cd8tNN.png"), width=wh, height=wh)
par(mfrow=c(whn,whn),mar=c(5,5,5,5))
fslist = as(getData(gs, 'CD3+CD19-'), Class='list')
for (i in 1:length(fslist)) {
  f = fslist[[i]]
  plotDens(f, c(channels.ind["CD4"],channels.ind["CD8"]), main=meta_file$id[i], cex.lab = 2, cex.axis = 2, cex.main=2)
  lines(filters[[i]]$cd4t, lwd=2)
  lines(filters[[i]]$cd8t, lwd=2)
  lines(filters[[i]]$notcd4cd8t, lwd=2)
}
graphics.off()


## cd4tcell: CD4+ T naive & CD4+ T memory cells
cd45pf = llply(filters, function(x) return(polygonGate(filterId="cd4tnaive", .gate=x$cd4.T.Naive))) # CD4.T.Naive CD45RA+
cd45nf = llply(filters, function(x) return(polygonGate(filterId="cd4tmemory", .gate=x$cd4.T.Memory))) # CD4.T.Memory CD45RA-
nodeIDs <- flowWorkspace::add(gs, cd45pf, parent="cd4tcell")
nodeIDs <- flowWorkspace::add(gs, cd45nf, parent="cd4tcell")
recompute(gs)

png(paste0(temp_gates_dir, "/10_cd4tcell_cd4tnaive-cd4tmemory.png"), width=wh, height=wh)
par(mfrow=c(whn,whn),mar=c(5,5,5,5))
fslist = as(getData(gs, 'cd4tcell'), Class='list')
for (i in 1:length(fslist)) {
  f = fslist[[i]]
  plotDens(f, c(channels.ind["CD4"],channels.ind["CD45RA"]), main=meta_file$id[i], cex.lab = 2, cex.axis = 2, cex.main=2)
  lines(filters[[i]]$cd4.T.Naive, lwd=2)
  lines(filters[[i]]$cd4.T.Memory, lwd=2)
}
graphics.off()


tbet =  llply(gthres, function(x) return(quadGate("Gd160Di"=x["tbet.cd4t"], "Nd143Di"=x["cd45RA.high.cd4"])))
nodeIDs <- flowWorkspace::add(gs, tbet, parent="cd4tcell")
# cd4.Th1.Naive = Tbet+CD45RA+
# cd4.Th1.Memory = Tbet+CD45RA-
recompute(gs)

png(paste0(temp_gates_dir, "/11_cd4tcell_cd4th1naivePP-cd4th1memoryPN.png"), width=wh, height=wh)
par(mfrow=c(whn,whn),mar=c(5,5,5,5))
fslist = as(getData(gs, 'cd4tcell'), Class='list')
for (i in 1:length(fslist)) {
  f = fslist[[i]]
  plotDens(f, c(channels.ind["Tbet"],channels.ind["CD45RA"]), main=meta_file$id[i], cex.lab = 2, cex.axis = 2, cex.main=2)
  abline(h=gthres[[i]]["tbet.cd4t"], lwd=2)
  abline(v=gthres[[i]]["cd45RA.high.cd4"], lwd=2)
}
graphics.off()


tbetf = llply(gthres, function(x) return(quadGate("Gd160Di"=x["tbet.cd8t"], "Nd143Di"=x["cd45RA.cd8"])))
nodeIDs <- flowWorkspace::add(gs, tbetf, parent="cd8tcell")
# cd25.cd8.T.Naive = Tbet+CD45RA+
# cd25.cd8.T.Memory = Tbet+CD45RA-
recompute(gs)

png(paste0(temp_gates_dir, "/12_cd8tcell_cd25cd8tnaivePP-cd25cd8tmemoryPN.png"), width=wh, height=wh)
par(mfrow=c(whn,whn),mar=c(5,5,5,5))
fslist = as(getData(gs, 'cd8tcell'), Class='list')
for (i in 1:length(fslist)) {
  f = fslist[[i]]
  plotDens(f, c(channels.ind["Tbet"],channels.ind["CD45RA"]), main=meta_file$id[i], cex.lab = 2, cex.axis = 2, cex.main=2)
  abline(h=gthres[[i]]["tbet.cd8t"], lwd=2)
  abline(v=gthres[[i]]["cd45RA.cd8"], lwd=2)
}
graphics.off()


## Tregs Naive
tnf = llply(filters, function(x) return(polygonGate(filterId="tregs.naive", .gate=x$tregs.naive)))
nodeIDs <- flowWorkspace::add(gs, tnf, parent="cd4tnaive")
recompute(gs)

png(paste0(temp_gates_dir, "/13_cd4tnaive_tregsnaive.png"), width=wh, height=wh)
par(mfrow=c(whn,whn),mar=c(5,5,5,5))
fslist = as(getData(gs, 'cd4tnaive'), Class='list')
for (i in 1:length(fslist)) {
  f = fslist[[i]]
  plotDens(f, c(channels.ind["FoxP3"],channels.ind["CD25"]), main=meta_file$id[i], cex.lab = 2, cex.axis = 2, cex.main=2)
  lines(filters[[i]]$tregs.naive, lwd=2)
}
graphics.off()


## Tregs Memory
tmf = llply(filters, function(x) return(polygonGate(filterId="tregs.memory", .gate=x$tregs.memory)))
nodeIDs <- flowWorkspace::add(gs, tmf, parent="cd4tmemory")
recompute(gs)

png(paste0(temp_gates_dir, "/14_cd4tmemory_tregsmemory.png"), width=wh, height=wh)
par(mfrow=c(whn,whn),mar=c(5,5,5,5))
fslist = as(getData(gs, 'cd4tmemory'), Class='list')
for (i in 1:length(fslist)) {
  f = fslist[[i]]
  plotDens(f, c(channels.ind["FoxP3"],channels.ind["CD25"]), main=meta_file$id[i], cex.lab = 2, cex.axis = 2, cex.main=2)
  lines(filters[[i]]$tregs.memory, lwd=2)
}
graphics.off()


## gamma-delta T-cells
gdtf = llply(filters, function(x) return(polygonGate(filterId="gammadelta.tcell", .gate=x$gammadelta.tcell)))
nodeIDs <- flowWorkspace::add(gs, gdtf, parent="notcd4cd8tcell")
recompute(gs)

png(paste0(temp_gates_dir, "/15_notcd4cd8tcell_gammadeltatcell.png"), width=wh, height=wh)
par(mfrow=c(whn,whn),mar=c(5,5,5,5))
fslist = as(getData(gs, 'notcd4cd8tcell'), Class='list')
for (i in 1:length(fslist)) {
  f = fslist[[i]]
  plotDens(f, c(channels.ind["TCRgd"],channels.ind["CD3"]), main=meta_file$id[i], cex.lab = 2, cex.axis = 2, cex.main=2)
  lines(filters[[i]]$gammadelta.tcell, lwd=2)
}
graphics.off()


## CD8+ T Naive & CD8+ T Memory
cd8nf = llply(filters, function(x) return(polygonGate(filterId="cd8tnaive", .gate=x$cd8.T.Naive))) # CD45RA+ gate cd45ra.cd8
nodeIDs <- flowWorkspace::add(gs, cd8nf, parent="cd8tcell")
cd8mf = llply(filters, function(x) return(polygonGate(filterId="cd8tmemory", .gate=x$cd8.T.Memory)))# CD45RA- gate cd45ra.cd8
nodeIDs <- flowWorkspace::add(gs, cd8mf, parent="cd8tcell")
recompute(gs)

png(paste0(temp_gates_dir, "/16_cd8tcell_cd8tnaiveRED-cd8tmemoryBLACK.png"), width=wh, height=wh)
par(mfrow=c(whn,whn),mar=c(5,5,5,5))
fslist = as(getData(gs, 'cd8tcell'), Class='list')
for (i in 1:length(fslist)) {
  f = fslist[[i]]
  plotDens(f, c(channels.ind["CD8"],channels.ind["CD45RA"]), main=meta_file$id[i], cex.lab = 2, cex.axis = 2, cex.main=2)
  lines(filters[[i]]$cd8.T.Naive, lwd=2, col='red')
  lines(filters[[i]]$cd8.T.Memory, lwd=2, col='black')
}
graphics.off()


## plot thresholds ---------------------------------
gtsr = ceiling(sqrt(length(gthres[[1]]))) 
png(paste0(temp_gates_dir, "/gates.png"), width=gtsr*500, height=gtsr*300)
par(mfrow=c(gtsr,gtsr),mar=c(5,5,5,5))
for (i in order(names(gthres[[1]]))) {
  allthres = laply(gthres, function(x) x[i])
  plot(density(allthres, na.rm=T), main=names(gthres[[1]])[i])
}
graphics.off()



## flowtype -----------------------------------------

fslist = as(fs,Class="list")
markers = c("CD11b", "CD11c", "CD123", "CD14", "CD16", "CD19", 
            # "CD25",
            "CD3", "CD4","CD45", "CD45RA", "CD56", 
            "CD66", 
            "CD7", "CD8", 
            # "FoxP3", 
            "HLADR", "Tbet", "TCRgd")
gthresm = c("cd11b", "cd11c", "cd123", "cd14.low", "cd16.gate.mid", "cd19", 
            # "cd25...",# slanted
            "cd3", "cd4","cd45", "cd45ra.cd8", "cd56.gate.mid", 
            "cd66", #.low
            "cd7", "cd8a", 
            # "foxp3...", # slanted
            "hladr.MMDSCs", "tbet.cd8t", "tcrd")
ftl = llply(1:length(fslist), function(i) {
  flowType(Frame=fslist[[i]], 
           PropMarkers=channels.ind[markers], MarkerNames=markers, 
           MaxMarkersPerPop=6, PartitionsPerMarker=2, Methods='Thresholds', 
           Thresholds=as.list(gthres[[i]][gthresm]), 
           verbose=F, MemLimit=400)
}, .parallel=F)
ft = ldply(ftl, function(ft) ft@CellFreqs)
ftcell = unlist(lapply(ftl[[1]]@PhenoCodes, function(x){return( decodePhenotype(x, markers, ftl[[1]]@PartitionsPerMarker) )}))
meta_cell = getPhen(ftcell)
ftp = Reduce('cbind', llply(1:ncol(ft), function(xi) ft[,xi]/ft[,1]))
colnames(ft) = colnames(ftp) = ftcell
rownames(ft) = rownames(ftp) = meta_file$id

writecsv = F
for (typed in c("train","validate")) {
  typei = meta_file$type==typed
  
  meta_dir = paste(result_dir, "_",typed,"/meta", sep=""); dir.create(meta_dir, showWarnings=F, recursive=T)
  meta_file_ = meta_file[typei,-"type"]
  save(meta_file_, file=paste0(meta_dir,"/file.Rdata"))
  if (writecsv) write.csv(meta_file_, file=paste0(meta_dir,"/file.csv"), row.names=F)
  save(meta_cell, file=paste0(meta_dir,"/cell.Rdata"))
  if (writecsv) write.csv(meta_cell, file=paste0(meta_dir,"/cell.csv"), row.names=F)
  
  feat_dir = gsub("meta","feat",meta_dir); dir.create(feat_dir, showWarnings=F, recursive=T)
  feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")
  feat_file_cell_prop_dir = paste(feat_dir, "/file-cell-prop", sep="")
  ft_ = ft[typei,]
  ftp_ = ftp[typei,]
  save(ft_, file=paste0(feat_file_cell_count_dir,".Rdata"))
  if (writecsv) write.csv(ft_, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)
  save(ftp_, file=paste0(feat_file_cell_prop_dir,".Rdata"))
  if (writecsv) write.csv(ftp_, file=paste0(feat_file_cell_prop_dir,".csv"), row.names=T)
}


