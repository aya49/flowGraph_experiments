#Modified by Alice Yue 20160518

## 000_RUNME

libr("flowCore")
libr("flowDensity")
libr("flowType")
# libr("flowBin")
# libr("stringr")
# libr("snowfall")
# libr("foreach")
#install.packages("dplyr")
# libr(extrafont)
libr("stringr")
# libr("arrayMvout")
libr("pracma")
libr("quantmod")
source("code/3iTcellfunctions.R")

load("allFCS.Rdata")

UpdatedVersion <- TRUE
# Set to TRUE if you want to record MFI while gating - do this when all the gating is done (set the UpdatedVersion to TRUE), then it will go faster :)
RecordMFI <- TRUE
UpdatedMFI <- F
fixLymph <- F
fixLymph1 <- F
fixLymph2 <- F
fixNotGran <- F
fixLive <- F
fixLive2 <- F
fixTcell <- F
fixCD11c <- F
fixB2 <- F
deltwo <- c(8,14)

fixLymphL <- NULL
fixLymph1L <- NULL
fixNotGranL <- NULL
fixLiveL <- NULL
fixLive2L <- NULL
fixTcellL <- NULL
fixCD11cL <- NULL
fixB2L <- NULL

# ## Create directories
# load("uniqueGT.Rdata")
# suppressWarnings ( dir.create ( "Figures/ScatterPlots_Updated") )
# invisible(sapply(1:length(uniqueGT), function(x){
#   suppressWarnings(dir.create ( paste("Figures/ScatterPlots_Updated/", uniqueGT[x], sep="") ))}))

load("par.Rdata")
load("store.allFCS.Rdata")
load("lgl657.Rdata")

if (UpdatedVersion == TRUE) {
  # load("changed.allFCS.Rdata")
  if (file.exists("all.gthresUpdated_20160718.Rdata")) { gt <- get(load("all.gthresUpdated_20160718.Rdata")) } else { pr <- get(load("all.gthres.Rdata")) }
}
if (RecordMFI == TRUE) {
  MFI.mean.OriginalData <- get(load("MFI.mean.OriginalData.Rdata"))
  MFI.median.OriginalData <- get(load("MFI.median.OriginalData.Rdata"))
  MFI.mean.TransformedData <- get(load("MFI.mean.TransformedData.Rdata"))
  MFI.median.TransformedData <- get(load("MFI.median.TransformedData.Rdata"))
}

load("channels.ind.Rdata")
verbose_debris <- F
verbose_margin <- F
verbose_flowClean <- F
index.Clean.Error <- NULL
errorFileIndex <- NULL
MFIerrorFileIndex <- NULL
clean.chan <- as.vector(sample(channels.ind,6))


# lymphGateError <- c(11,43,63,82,83,87,108,143,155,156,202,208,210,214,230,235,239,258,260,261,273,274,298,308,319,354,362,364,396,406,414,473,487,488,491,512,517,535,542,561,583,585,667,671,682,702,703,713,720,760,784,808,814,824,827,868,900,901,917,918,929,930,950,962,963,965,1003,1006,1023,1034,1036,1083,1085,1086,1099,1100,1101,1102,1104,1141,1143,1157,1160,1167,1172,1173,1179,1195,1198,1208,1210,1237,1258,1263,1299,1329,1349,1443,1445,1454,1456,1483,1497,1518,1541,1578,1663,1694,1708,1748,1750,1752,1753,1769,1811,1823,1829,1852,1853,1861,1871,1872,1874,1891,1892)
# NKGateError <- c(309,208,273,274,487,585,760,824,827,900,1003,1085,1141,1456,1497,1913,7,130,131,240,244,268,279,423,742,921,922,1671,1732,1817)
# exDir <- dir("Figures/ScatterPlots/_exceptions", full.names=F, recursive=T)
# #exceptions
# exFiles <- NULL 
# for (i in 1:length(exDir)) {
#   exFiles <- append(exFiles, as.numeric( strsplit( errorFile2dir[i] , "_")[[1]][1] ))
# }
# 
# errorFile1 <- which(all.props[,ncol(all.props)]==0) 
# #_fix + ErrorFileIndex
# errorFile2dir <- dir("Figures/ScatterPlots/_fix", full.names=F, recursive=T)
# errorFile2 <- NULL
# for (i in 1:length(errorFile2dir)) {
#   errorFile2 <- append(errorFile2, as.numeric( strsplit( errorFile2dir[i] , "_")[[1]][1] ))
# }
# errorFile <- sort(append(errorFile1,errorFile2))

# updatedBadly <- c(154,189,362,363,370,372,386,395,574,677,793,797,877,927,929,1100,1113,1143,1155,1183,1211,1216,1249,1269,1271,1278,1282,1316,1332,1344,1372,1375,1392,1396,1407,1420,1534,1540,1542,1543,1544,1717,1724,1804,1823,1848,1858,1916)  
# updatedNewly <- sort(c(1750,1748,1368,1370,21, 1371,1158,1329,1775,965,1369,1425,1872,82,215,491,673,862,879,1100,1143,1264,1547,1720,1829))
# ev[updatedBadly,] <- get(load("all.events.Rdata"))[updatedBadly,]
# pr[updatedBadly,] <- get(load("all.props.Rdata"))[updatedBadly,]
# gt[updatedBadly,] <- get(load("all.gthres.Rdata"))[updatedBadly,]

errorFiles <- NULL
save(errorFiles, file="errorFiles.Rdata")
start <- Sys.time()

# for(i in sample(changed.allFCS, 10)){
# for(i in changed.allFCS[c(20:length(changed.allFCS))]) {
# for (i in updatedNewly) {
for(i in 1:nrow(store.allFCS)) {
  # for (i in updatedNewly) {
  # for(i in 76:nrow(store.allFCS)){for(i in sample(c(1:nrow(store.allFCS)),10) ) {
  # for(i in sample(c(1:nrow(store.allFCS)),1917) ){
  # for(i in sample(c(1:nrow(store.allFCS)),1) ){
  # for (i in which(all.props[,ncol(all.props)]==0)) {
  # 0bytes
  # for (i in c(333,334,337,709,710,712,791,841,843,844,886) {
  # for (i in which(store.allFCS$channels==18)) {
  # Weird ungated
  # for (i in c(79,204,205,275,278,324,351,538,555,598,610,741,1007,1014,1030,1521,1774,1178,1682,1788,1787,1829)) {
  # Tcell <- c(1369,845,61,221,572,702,1014,1100,1113,1425,1671,1709,1820,36,89,111,124,154,160,162,232,283,292,315,370,400,443,465,480,497,560,527,529,560,576,593,626,747,811,889,930,1004,1005,1006,1009,1057,1106,1107,1150,1234,1243,1277,1278,1279,1280,1302,1304,1305,1306,1307,1308,1310,1311,1312,1313,1344,1348,1349,1513,1638,1656,1379,1690,1712,1714,1730,1749,1818,1829,1843,1872,1873,1903,1913,1682)
  # for (i in Tcell) {
  # for (i in lymphGateError) {
  # CD5 gate too right
  # for (i in c(1595,751,753,1870,835)) {
  # for (i in c(1:1681)) {
  # for (i in errorFile) {
  # for (i in updatedNewly) {
  
  start2 <- Sys.time()
  cat(paste(i, ": Starting ", store.allFCS[i,1], " / ", store.allFCS[i,2]), "; ", sep="" )
  
  tryCatch({
    load (file =  paste("After_Clean/", store.allFCS[i,1],"/AftClean_", store.allFCS[i,2],".Rdata",sep="") )
    
    fn <- paste("Figures/ScatterPlots/_20160715/", str_pad(i, 4, pad = "0"), "_", store.allFCS[i,1], "_", store.allFCS[i,2], "smallLymph.png", sep = "" )
    png ( file = fn, width=2200, height=1800)
    par(mfrow=c(4,5),mar=(c(5, 5, 4, 2) + 0.1))
    
    cat("); Gating")
    
    
    
    ## Gating All Events. for singlets
    
    ## 1st method of plotting FSC-A_SSC-W
    if (!UpdatedVersion) {
      # gt[i,"singlets.gate"] <- deGate(f, channel = c("SSC-W"), use.percentile = T, percentile = 0.90)
      # gt[i,"singlets.gate.l"] <- deGate(f, channel = c("SSC-W"), use.percentile = T, percentile = 0.0001)
      d <- density(SPLN_L000064504_D02_009.labelled.fcsf@exprs[,"SSC-W"], na.rm=T)
      d <- smooth.spline(d$x, d$y, spar=0.4)
      d$y[which(d$y<0)] <- 0
      ma <- d$x[which.max(d$y)]
      a <- deGate(f, channel = c("SSC-W"), tinypeak.removal = .1)+5000
      v <- d$x[findValleys(d$y)]
      v <- v[which(v>ma)]
      gt[i,"singlets.gate"] <- v[which.min(v-ma)]
      b <- deGate(f, channel = c("SSC-W"), upper=F)
      if (b<ma) {
        gt[i,"singlets.gate.l"] <- b
      }
    }
    temp <- flowDensity(f, channels = c("FSC-A", "SSC-W"), position = c(NA, F), gates = c(NA, gt[i,"singlets.gate"]))
    singlets.flowD.l <- flowDensity(temp, channels = c("FSC-A", "SSC-W"), position = c(NA, T), gates = c(NA, gt[i,"singlets.gate.l"]))
#     pr[i,2] <- singlets.flowD.l@proportion <- (singlets.flowD.l@cell.count/nrow(f))*100
#     ev[i,2] <- singlets.flowD.l@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(1:16)] <- sapply(1:16, function(x){mean(f.bt@exprs[singlets.flowD.l@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(1:16)] <- sapply(1:16, function(x){median(f.bt@exprs[singlets.flowD.l@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(1:16)] <- sapply(1:16, function(x){mean(f@exprs[singlets.flowD.l@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(1:16)] <- sapply(1:16, function(x){median(f@exprs[singlets.flowD.l@index, x], na.rm = T)}) }
    
    #     ## 2nd Alternate method of plotting FSC-A_SSC-W
    #     rot<-rotate.data(f,c("FSC-A","SSC-W"),theta = -pi/2)$data
    #     singlets.high <- deGate(rot, "FSC-A",percentile=.99, use.percentile=T)
    #     singlets.low  <- deGate(rot, "FSC-A",percentile=.1,use.percentile=T)
    #     names(singlets.high) <- names(singlets.low) <- identifier(f)
    #     temp<- flowDensity(rot,channels=c("FSC-A","SSC-W"),position=c(F,NA), gates=c(singlets.high[identifier(f)],NA))
    #     temp2<- flowDensity(temp@flow.frame, channels=c("FSC-A","SSC-W"),position=c(T,NA), gates=c(singlets.low[identifier(f)],NA))
    #     temp2@filter <- rotate.data(temp2@filter,c("FSC-A","SSC-W"),theta = pi/2)$data
    #     temp2@flow.frame <- rotate.data(getflowFrame(temp2),c("FSC-A","SSC-W"),theta = pi/2)$data
    #     temp2@proportion <- (temp2@cell.count/nrow(f))*100
    #     singlets <- getflowFrame(temp2)
    #     # plotDens(f,c("FSC-A","SSC-W"),main="Ungated")
    #     # lines(temp2@filter,lwd=2)
    
    plotDens(f,c("FSC-A","SSC-W"),main="All events", xlab="01. FSC-A", ylab="06. SSC-W", xlim=c(0, 275000), ylim=c(0, 275000), devn = F)
    lines(singlets.flowD.l@filter,lwd=1)
    
    
    
    ## Gating Singlets. Plotting Live/Dead_SSC-A for live
    if ((!UpdatedVersion | fixLive) & !i%in%fixLiveL) {
      fixLiveL <- append(fixLiveL, i)
      temp <- flowDensity(singlets.flowD.l, channels=c(10,4), position=c(NA,F), gates=c(NA,40000))
      d <- density(temp@flow.frame@exprs[,10], na.rm=T)
      fp <- findpeaks(d$y, minpeakheight=.1)
      if (length(as.vector(fp))>4) {
        gt[i,"live.gate"] <- deGate(temp@flow.frame, channel = c(10))-.3
      } else if (deGate(temp@flow.frame, channel = c(10), use.percentile=T, percentile=1)<3) { #already gated
        gt[i,"live.gate"] <- deGate(temp@flow.frame, channel = c(10), use.percentile=T, percentile=1)
      } else {
        gt[i,"live.gate"] <- deGate(temp@flow.frame, channel = c(10))
      }
      #         if (gt[i,"live.gate"]>2.3 | gt[i,"live.gate"]<1.2) {
      #           temp <- flowDensity(temp, channels=c(10,4), position=c(T,NA), gates=c(1.2,NA))
      #           temp <- flowDensity(temp, channels=c(10,4), position=c(F,NA), gates=c(2.28,NA))
      #           gt[i,"live.gate"] <- deGate(temp@flow.frame, channel = c(10))-.3
      #         }
    }
    
    live.flowD <- flowDensity(singlets.flowD.l, channels = c(10,4), position = c(F,NA), gates = c(gt[i,"live.gate"],NA))
#     pr[i,3] <- live.flowD@proportion
#     ev[i,3] <- live.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(17:32)] <- sapply(1:16, function(x){mean(f.bt@exprs[live.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(17:32)] <- sapply(1:16, function(x){median(f.bt@exprs[live.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(17:32)] <- sapply(1:16, function(x){mean(f@exprs[live.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(17:32)] <- sapply(1:16, function(x){median(f@exprs[live.flowD@index, x], na.rm = T)}) }
    
    plotDens(singlets.flowD.l@flow.frame,  c(10,4), main= "Singlets", xlab="10. live", ylab="04. SSC-A", xlim=c(0, 4.5), ylim=c(0, 275000), devn = F); abline(v=gt[i,"live.gate"], lwd=2); 
    lines(live.flowD@filter, lwd=1)
    
    rm(singlets.flowD.l)
    
    
    
    ## Gating Live. Plotting FSC-A_SSC-A for Lymphocytes
    if (!UpdatedVersion) {
      gt[i,"fsc.a.gate.low"] <- deGate(live.flowD@flow.frame, channel = c(1), tinypeak.removal = 0.6, upper = F)
      if (gt[i,"fsc.a.gate.low"]>70000) {
        temp1 <- flowDensity(live.flowD, channels = c(1,4), position = c(F, NA), gates = c(gt[i,"fsc.a.gate.low"], NA))
        gt[i,"fsc.a.gate.low"] <- deGate(temp@flow.frame, channel = c(1), tinypeak.removal = 0.9)
      }
      #         if (gt[i,"fsc.a.gate.low"]<40000) {
      #           temp2 <- flowDensity(live.flowD, channels = c(1,4), position = c(T, NA), gates = c(gt[i,"fsc.a.gate.low"], NA))
      #           gt[i,"fsc.a.gate.low"] <- deGate(temp@flow.frame, channel = c(1), tinypeak.removal = 0.9)
      #         }
      #         if (gt[i,"fsc.a.gate.low"]<40000 | gt[i,"fsc.a.gate.low"]>70000) { gt[i,"fsc.a.gate.low"] <- 50000 }
      gt[i,"fsc.a.gate"] <- deGate(live.flowD@flow.frame, channel = c(1), use.percentile = T, percentile = 0.99999999)
      gt[i,"ssc.a.gate"] <- deGate(live.flowD@flow.frame, channel = c(4), use.percentile = T, percentile = 0.9999999)
    }
    if (fixLive2 & !i%in%fixLive2L) {
      fixLive2L <- append(fixLive2L, i)
      d <- density(live.flowD@flow.frame@exprs[,1], na.rm=T)
      fp <- findpeaks(d$y, minpeakheight=2e-06)
      if (length(as.vector(fp))<=4) {
        gt[i,"fsc.a.gate.low"] <- deGate(live.flowD@flow.frame, channel = c(1), tinypeak.removal = 0.6, upper = F)
      } else if (max(fp[,1])!=fp[which(fp[,2]==sort(fp[,2])[2]),1]) { # max peak isn't the second one
        gt[i,"fsc.a.gate.low"] <- deGate(live.flowD@flow.frame, channel = c(1), tinypeak.removal = 0.6, upper = F)
      } else {
        gt[i,"fsc.a.gate.low"] <- deGate(live.flowD@flow.frame, channel = c(1), tinypeak.removal = 0.6, upper = F)
      }
    }
    temp <- flowDensity(live.flowD, channels = c(1,4), position = c(T, F), gates = c(gt[i,"fsc.a.gate.low"], gt[i,"ssc.a.gate"]))
    lymph.flowD <- flowDensity(temp, channels = c(1,4), position = c(F,NA), gates = c(gt[i,"fsc.a.gate"], NA))
#     pr[i,4] <- lymph.flowD@proportion <- 100*lymph.flowD@cell.count/live.flowD@cell.count
#     ev[i,4] <- lymph.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(33:48)] <- sapply(1:16, function(x){mean(f.bt@exprs[lymph.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(33:48)] <- sapply(1:16, function(x){median(f.bt@exprs[lymph.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(33:48)] <- sapply(1:16, function(x){mean(f@exprs[lymph.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(33:48)] <- sapply(1:16, function(x){median(f@exprs[lymph.flowD@index, x], na.rm = T)}) }
    
    
    rm(live.flowD)
    
    
    
    ## Gating Lymphocytes. Plotting CD5(8)_CD11b(12) for Not/Granulocytes
    # Create CD5 gate to left (2 ovals); CD11b gate
    if (!UpdatedVersion) {
      gt[i,"CD5.gate.high"] <- deGate(lymph.flowD@flow.frame, channel = c(8), use.percentile = T, percentile = 0.99999)
      gt[i,"CD11b.gate"] <- deGate(lymph.flowD@flow.frame, channel = c(12))-.4
      gt[i,"CD11b.gate.high"] <- deGate(lymph.flowD@flow.frame, channel = c(12), use.percentile = T, percentile = 0.9999999)
      
      temp <- flowDensity(lymph.flowD, channels = c(8,12), position = c(F,T), gates = c(gt[i,"CD5.gate.high"], gt[i,"CD11b.gate"])) #upper half
      gt[i,"CD5.gate"] <- deGate(temp@flow.frame, channel = c(8), tinypeak.removal = .1)
      if (gt[i,"CD11b.gate"] < 2.3) { #if too low
        gt[i,"CD11b.gate"] <- deGate(temp@flow.frame, channel = c(12), upper=F, tinypeak.removal = .6)-.1
        if (gt[i,"CD11b.gate"] < 2.3 | gt[i,"CD11b.gate"] > 3.1) { gt[i,"CD11b.gate"] <- 2.4 }
      }
      if (gt[i,"CD11b.gate"] > 3.1) {
        temp1 <- flowDensity(temp, channels = c(8,12), position = c(T,F), gates = c(NA,gt[i,"CD11b.gate"]))
        if (temp1@cell.count>1) { gt[i,"CD11b.gate"] <- deGate(temp@flow.frame, channel = c(12), tinypeak.removal = .6)-.1 }
        if (gt[i,"CD11b.gate"] > 3.1 | gt[i,"CD11b.gate"] < 2.3) { gt[i,"CD11b.gate"] <- 2.5 }
      }
      temp <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(gt[i,"CD5.gate"], gt[i,"CD11b.gate"]))
      while ( gt[i,"CD5.gate"]>2.75 | temp@proportion<.5 | gt[i,"CD5.gate"]<2.2) {
        if (gt[i,"CD5.gate"]>2.75 | temp@proportion<.3) { #if too right
          temp1 <- flowDensity(lymph.flowD, channels = c(8,12), position = c(F,T), gates = c(gt[i,"CD5.gate"], gt[i,"CD11b.gate"]))
          gt[i,"CD5.gate"] <- deGate(temp1@flow.frame, channel=c(8), tinypeak.removal = .1)
          if (gt[i,"CD5.gate"]>2.75 | temp@proportion<.5) { gt[i,"CD5.gate"] <- gt[i,"CD5.gate"]-.35 }
        } else {
          temp1 <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(gt[i,"CD5.gate"], gt[i,"CD11b.gate"]))
          gt[i,"CD5.gate"] <- deGate(temp1@flow.frame, channel=c(8), tinypeak.removal = .1, upper=F)
          if (gt[i,"CD5.gate"]<2.2) { gt[i,"CD5.gate"] <- deGate(temp1@flow.frame, channel=c(8), tinypeak.removal = .1) }
        }
        temp <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(gt[i,"CD5.gate"], gt[i,"CD11b.gate"]))
      }
    } # else { temp <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(gt[i,"CD5.gate"], gt[i,"CD11b.gate"])) }
    
    # fix CD5 gate to right (1 oval); CD11b gate up tight
    if (fixLymph & !i%in%fixLymphL) {
      fixLymphL <- append(fixLymphL, i)
      temp <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(gt[i,"CD5.gate"], gt[i,"CD11b.gate"]))
      temp1 <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(1.5, gt[i,"CD11b.gate"]))
      temp1 <- flowDensity(temp1, channels = c(8,12), position = c(F,NA), gates = c(3,NA))
      CD5.gate <- gt[i,"CD5.gate"]
      require(quantmod)
      require(pracma)
      d <- density(temp1@flow.frame@exprs[,8], na.rm=T)
      v <- d$x[sort(findValleys(d$y))]
      p <- d$x[sort(findPeaks(d$y))]
      if (length(p)==1) {
        fp <- findpeaks(d$y, npeaks=1)
        temp1 <- flowDensity(temp1, channels = c(8,12), position = c(F,NA), gates = c(d$x[p[4]],NA))
        temp1 <- flowDensity(temp1, channels = c(8,12), position = c(T,NA), gates = c(d$x[p[3]],NA))
        d <- density(temp1@flow.frame@exprs[,8], na.rm=T)
        fp <-findpeaks(d$y, npeaks=3)
        if (is.matrix(fp)) {
          p <- d$x[sort(fp[,2])]
        } else {
          p <- d$x[sort(fp[2])]
        }
        gt[i,"CD5.gate"] <- mean(p[2],p[3])
      } else if (length(p)>=3) {
        for (k in 0:(length(v)-1)) {
          gt[i,"CD5.gate"] <- v[length(v)-k]
          temp <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(gt[i,"CD5.gate"], gt[i,"CD11b.gate"]))
          if (temp@proportion>.5 & gt[i,"CD5.gate"]>2.2) {
            break()
          }
        }
      } else if (length(p)==2) {
        if (v[1]<2.2) { #if valley low, get right half and find min between 2 peaks as gate CD5
          temp1 <- flowDensity(temp1, channels = c(8,12), position = c(T,NA), gates = c(v[1],NA))
          d <- density(temp1@flow.frame@exprs[,8], na.rm=T)
          fp <- findpeaks(d$y, npeaks=2)
          if (is.matrix(fp)) {
            dy2n <- sort(fp[,2])
          } else {
            dy2n <- sort(fp[2])
          }
          if (length(dy2n)==2) {
            gt[i,"CD5.gate"] <- d$x[ dy2n[1]+which.min(d$y[c(dy2n[1]:dy2n[2])]) ]
            temp <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(gt[i,"CD5.gate"], gt[i,"CD11b.gate"]))
            if (temp@proportion<.5) {
              temp1 <- flowDensity(temp1, channels = c(8,12), position = c(F,NA), gates = c(gt[i,"CD5.gate"],NA))
              gt[i,"CD5.gate"] <- deGate(temp1@flow.frame, channel=c(8))
            }
          } else { gt[i,"CD5.gate"] <- p[2] }
        } else {
          gt[i,"CD5.gate"] <- v[1]
          temp <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(gt[i,"CD5.gate"], gt[i,"CD11b.gate"]))
          if (temp@proportion<.5) {
            gt[i,"CD5.gate"] <- gt[i,"CD5.gate"]-.45
          }
        }
      }
      #         if (temp@proportion<.5) { #if too right
      #           gt[i,"CD5.gate"] <- deGate(temp1@flow.frame, channel=c(8), tinypeak.removal = .1)
      #           temp <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(gt[i,"CD5.gate"], gt[i,"CD11b.gate"]))
      #           if (temp@proportion<.5) {
      #             temp1 <- flowDensity(lymph.flowD, channels = c(8,12), position = c(F,T), gates = c(gt[i,"CD5.gate"], gt[i,"CD11b.gate"]))
      #             gt[i,"CD5.gate"] <- deGate(temp1@flow.frame, channel=c(8), tinypeak.removal = .1)
      #             temp <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(gt[i,"CD5.gate"], gt[i,"CD11b.gate"]))
      #             if (temp@proportion<.5) {
      #               gt[i,"CD5.gate"] <- gt[i,"CD5.gate"]-.35
      #               temp <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(gt[i,"CD5.gate"], gt[i,"CD11b.gate"]))
      #             }
      #           }
      #         } else if (gt[i,"CD5.gate"]<2.2) {
      #           temp1 <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(gt[i,"CD5.gate"], gt[i,"CD11b.gate"]))
      #           gt[i,"CD5.gate"] <- deGate(temp1@flow.frame, channel=c(8), tinypeak.removal=.1)
      #           if (gt[i,"CD5.gate"]<2.2) {
      #             temp1 <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(gt[i,"CD5.gate"], gt[i,"CD11b.gate"]))
      #             gt[i,"CD5.gate"] <- deGate(temp1@flow.frame, channel=c(8), tinypeak.removal = .1)
      #           }
      #           if (gt[i,"CD5.gate"]<2.2) { gt[i,"CD5.gate"] <- gt[i,"CD5.gate"]+.3 }
      #           temp <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(gt[i,"CD5.gate"], gt[i,"CD11b.gate"]))
      #         }
      #         # compare with old CD5.gate
      #         if (CD5.gate > 2.2 & temp2@proportion>.5) { gt[i,"CD5.gate"] <- CD5.gate } # if original not too big
      #         if (CD5.gate < 2.2 & temp@proportion<.5) { gt[i,"CD5.gate"] <- mean(c(CD5.gate, gt[i,"CD5.gate"])) } # if original too small and current too big
      
      temp <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(gt[i,"CD5.gate"], 2.2))
      CD11b.gate <- deGate(temp@flow.frame, channel=c(12), tinypeak.removal=.8, upper=F)-.15
      d <- density(temp@flow.frame@exprs[,12], na.rm=T)
      if (abs(CD11b.gate-d$x[which.max(d$y)]) <.4) {
        temp <- flowDensity(temp, channels = c(8,12), position = c(NA,F), gates = c(NA, d$x[which.max(d$y)]))
        CD11b.gate <- deGate(temp@flow.frame, channel=c(12), tinypeak.removal=.8, upper=F)-.25
      }
      if (CD11b.gate <= 3.45) { gt[i,"CD11b.gate"] <- CD11b.gate }
      
      #       # Gate Lymph strategy #1
      #         gran.flowD <- gt[i,"Ly6C.gate1"] flowDensity(temp, channels = c(8,12), position = c(F,F), gates = c(gt[i,"CD5.gate.high"], gt[i,"CD11b.gate.high"]))
      #         pr[i,5] <- gran.flowD@proportion <- 100*gran.flowD@cell.count/lymph.flowD@cell.count
      #         ev[i,5] <- gran.flowD@cell.count
      #         if (RecordMFI) { MFI.mean.OriginalData[i,c(33:48)] <- sapply(1:16, function(x){mean(f.bt@exprs[gran.flowD@index, x], na.rm = T)})
      #         MFI.median.OriginalData[i,c(33:48)] <- sapply(1:16, function(x){median(f.bt@exprs[gran.flowD@index, x], na.rm = T)})
      #         MFI.mean.TransformedData[i,c(33:48)] <- sapply(1:16, function(x){mean(f@exprs[gran.flowD@index, x], na.rm = T)})
      #         MFI.median.TransformedData[i,c(33:48)] <- sapply(1:16, function(x){median(f@exprs[gran.flowD@index, x], na.rm = T)}) }
      #         
      #         notGran.flowD <- notSubFrame(lymph.flowD@flow.frame, channels = c(8,12), position="logical", gates="missing", gran.flowD@filter)
      #         pr[i,6] <- notGran.flowD@proportion <- 100-gran.flowD@proportion
      #         ev[i,6] <- notGran.flowD@cell.count
      #         if (RecordMFI) { MFI.mean.OriginalData[i,c(65:80)] <- sapply(1:16, function(x){mean(f.bt@exprs[notGran.flowD@index, x], na.rm = T)})
      #         MFI.median.OriginalData[i,c(65:80)] <- sapply(1:16, function(x){median(f.bt@exprs[notGran.flowD@index, x], na.rm = T)})
      #         MFI.mean.TransformedData[i,c(65:80)] <- sapply(1:16, function(x){mean(f@exprs[notGran.flowD@index, x], na.rm = T)})
      #         MFI.median.TransformedData[i,c(65:80)] <- sapply(1:16, function(x){median(f@exprs[notGran.flowD@index, x], na.rm = T)}) }
      #         
      #         suppressWarnings({
      #           plotDens(lymph.flowD@flow.frame, channels = c(8,12), main = "Lymphocytes", xlab="08. CD5", ylab="12. CD11b", xlim=c(0, 4.5), ylim=c(0, 4.5), devn = F); abline(v=gt[i,"CD5.gate"], lwd=2); abline(v=gt[i,"CD5.gate.high"], lwd=2); abline(h=gt[i,"CD11b.gate"], lwd=2); abline(h=gt[i,"CD11b.gate.high"], lwd=2)
      #           lines(gran.flowD@filter, lwd=1)
      #         })
      
    } else {
      temp <- flowDensity(lymph.flowD, channels = c(8,12), position = c(T,T), gates = c(gt[i,"CD5.gate"]-.5, 2.2))
      CD11b.gate <- deGate(temp@flow.frame, channel=c(12), tinypeak.removal=.6, upper=F)-.15
      d <- density(temp@flow.frame@exprs[,12], na.rm=T)
      if (abs(CD11b.gate-d$x[which.max(d$y)]) <.4) {
        temp <- flowDensity(temp, channels = c(8,12), position = c(NA,F), gates = c(NA, d$x[which.max(d$y)]))
        CD11b.gate <- deGate(temp@flow.frame, channel=c(12), tinypeak.removal=.8, upper=F)-.25
      }
      if (CD11b.gate <= 3.45) { gt[i,"CD11b.gate"] <- CD11b.gate }
    }
    
    # Gate Lymph strategy #2 in two chunks (gating strategy #2 2016 06 29)
    
    suppressWarnings({
      plotDens(lymph.flowD@flow.frame, channels = c(8,12), main = "Lymphocytes", xlab="08. CD5", ylab="12. CD11b", xlim=c(0, 4.5), ylim=c(0, 4.5), devn = F); abline(v=gt[i,"CD5.gate.high"], lwd=2); abline(v=gt[i,"CD5.gate"], lty=2); abline(h=gt[i,"CD11b.gate"], lwd=2); abline(h=gt[i,"CD11b.gate.high"], lwd=2)
    })
    
    lymph1.flowD <- flowDensity(lymph.flowD, channels = c(8,12), position = c(NA,T), gates = c(NA, gt[i,"CD11b.gate"]))
    
    if (fixLymph1) {
      fixLymph1L <- append(fixLymph1L, i)
      # Set Ly6C.gate1
      temp <- flowDensity(lymph1.flowD, channels=c(8,9), position=c(T,NA), gates=c(gt[i,"CD5.gate"]-.3,NA))
      d <- density(temp@flow.frame@exprs[,9], na.rm=T)
      fp <- findpeaks(d$y, minpeakheight=.3, minpeakdistance=60)
      if (is.matrix(fp)) { # multiple peaks
        if (d$x[max(fp[,2])]>2) { #if there is a bump up high>2
          if (nrow(fp)>1) {
            dy2n <- sort(fp[,2])
            m <- d$x[ dy2n[length(dy2n)-1]+which.min(d$y[c(dy2n[length(dy2n)-1]:dy2n[length(dy2n)])]) ]
          } else {
            m <- deGate(temp@flow.frame, channel=c(9), upper=F, tinypeak.removal = .8)
          }
        } else {
          m <- deGate(temp@flow.frame, channel=c(9), alpha=.5)
        }
      } else { # 1 peak only
        if (d$x[fp[2]]>2) { up <- T
        m <- deGate(temp@flow.frame, channel=c(9), upper=F, tinypeak.removal = .8)
        } else {
          m <- deGate(temp@flow.frame, channel=c(9), alpha=.5)
        }
      }
      temp <- flowDensity(lymph1.flowD, channels=c(8,9), position=c(NA,T), gates=c(NA,m-.1))
      if (d$x[which.max(d$y)]>2) { # biggest bump is above 2
        temp <- flowDensity(temp, channels=c(8,9), position=c(NA,F), gates=c(NA,d$x[which.max(d$y)]-.1))
        d2 <- density(temp@flow.frame@exprs[,9], na.rm=T)
      }
      a <- deGate(temp@flow.frame, channel=c(9), upper=F, tinypeak.removal = .2)-.05
      b <- deGate(temp@flow.frame, channel=c(9), upper=F, tinypeak.removal = .8)-.05
      if (a!=b) {
        ay <- d2$y[findInterval(a, d2$x)]; if (length(ay)==0) { ay <- min(d$x)}
        by <- d2$y[findInterval(b, d2$x)]; if (length(by)==0) { by <- min(d$x)}
        if (ay>by | (a<=2 & b>2)) { a <- b }
      }
      if (a<m | d$y[findInterval(a, d$x)]> .8*max(d$y)) { a <- m }
      gt[i,"Ly6C.gate1"] <- a
      # gt[i,"Ly6C.gate1"] <- d$x[ dy2n[length(dy2n)-1]+which.min(d$y[c(dy2n[length(dy2n)-1]:dy2n[length(dy2n)])]) ]
    }
    if (fixLymph2) {
      #Set CD5.gate again
      temp <- flowDensity(lymph1.flowD, channels=c(8,9), position=c(NA,T), gates=c(NA,gt[i,"Ly6C.gate1"]))
      d <- density(temp@flow.frame@exprs[,8], na.rm=T)
      fp <- findpeaks(d$y, minpeakheight=.25, minpeakdistance=60)
      if (is.matrix(fp)) { if (nrow(fp)>1) {
        dy2n <- sort(fp[,2])
        m <- d$x[ dy2n[length(dy2n)-1]+which.min(d$y[c(dy2n[length(dy2n)-1]:dy2n[length(dy2n)])]) ]
      } else {
        m <- deGate(temp@flow.frame, channel=c(8), upper=F, tinypeak.removal=.1)
      } } else {
        m <- deGate(temp@flow.frame, channel=c(8), upper=F, tinypeak.removal=.1)
      }
      temp <- flowDensity(temp, channels=c(8,9), position=c(T,NA), gates=c(m,NA))
      d <- density(temp@flow.frame@exprs[,8], na.rm=T)
      temp <- flowDensity(temp, channels=c(8,9), position=c(F,NA), gates=c(d$x[which.max(d$y)],NA))
      gt[i,"CD5.gate"] <- deGate(temp@flow.frame, channel=c(8), upper=F, tinypeak.removal=.6)
      if (abs(gt[i,"CD5.gate"]-d$x[which.max(d$y)])<.1) { gt[i,"CD5.gate"] <- m+.1}
      
    }
    
    temp <- flowDensity(lymph1.flowD, channels=c(8,9), position=c(T,T), gates=c(gt[i,"CD5.gate"],gt[i,"Ly6C.gate1"]))
    gran.flowD <- flowDensity(temp, channels=c(8,9), position=c(F,F), gates=c(gt[i,"CD5.gate.high"], gt[i,"Ly6C.gate.high"]))
#     pr[i,5] <- gran.flowD@proportion <- gran.flowD@cell.count/lymph.flowD@cell.count
#     ev[i,5] <- gran.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(49:64)] <- sapply(1:16, function(x){mean(f.bt@exprs[gran.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(49:64)] <- sapply(1:16, function(x){median(f.bt@exprs[gran.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(49:64)] <- sapply(1:16, function(x){mean(f@exprs[gran.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(49:64)] <- sapply(1:16, function(x){median(f@exprs[gran.flowD@index, x], na.rm = T)}) }
    
    temp1 <- flowDensity(lymph.flowD, channels = c(8,12), position = c(NA,F), gates = c(NA, gt[i,"CD11b.gate"]))
    notGran.flowD <- notSubFrame(lymph1.flowD@flow.frame, channels = c(8,9), position="logical", gates="missing", gran.flowD@filter)
    notGran.flowD@index <- unique(sort(c(notGran.flowD@index, temp1@index)))
    notGran.flowD@flow.frame@exprs[temp1@index,] <- temp1@flow.frame@exprs[temp1@index,]
#     ev[i,6] <- notGran.flowD@cell.count <- notGran.flowD@cell.count + temp1@cell.count
#     pr[i,6] <- notGran.flowD@proportion <- 100*notGran.flowD@cell.count/lymph.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(65:80)] <- sapply(1:16, function(x){mean(f.bt@exprs[notGran.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(65:80)] <- sapply(1:16, function(x){median(f.bt@exprs[notGran.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(65:80)] <- sapply(1:16, function(x){mean(f@exprs[notGran.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(65:80)] <- sapply(1:16, function(x){median(f@exprs[notGran.flowD@index, x], na.rm = T)}) }
    
    suppressWarnings({
      plotDens(lymph1.flowD@flow.frame, channels=c(8,9), main = "Lymphocytes-CD11b+", xlab="08. CD5", ylab="09. Ly6C", xlim=c(0, 4.5), ylim=c(0, 4.5), devn = F); abline(v=gt[i,"CD5.gate"], lwd=2); abline(h=gt[i,"Ly6C.gate1"], lwd=2); abline(v=gt[i,"CD5.gate.high"], lwd=2); abline(h=gt[i,"Ly6C.gate.high"], lwd=2)
      lines(gran.flowD@filter, lwd=1)
      plotDens(notGran.flowD@flow.frame, channels = c(8,12), main = "Not Granulocytes", xlab="08. CD5", ylab="12. CD11b", xlim=c(0, 4.5), ylim=c(0, 4.5), devn = F); abline(v=gt[i,"CD5.gate"], lty=2); abline(v=gt[i,"CD5.gate.high"], lwd=2); abline(h=gt[i,"CD11b.gate"], lty=2); abline(h=gt[i,"CD11b.gate.high"], lwd=2)
    })
    
    rm(lymph1.flowD)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ## Gating Not Granulocytes. Plotting Ly6C(9)_CD11b(12) for Not/Monocytes
    if (fixNotGran & !i%in%fixNotGranL){
      fixNotGranL <- append(fixNotGranL, i)
      temp <- flowDensity(notGran.flowD, channels = c(9,12), position = c(F,T), gates = c(4, gt[i,"CD11b.gate"]))
      d <- density(temp@flow.frame@exprs[,9], na.rm=T)
      v <- d$x[findValleys(d$y)]
      v <- v[length(v)]
      temp <- flowDensity(temp, channels = c(9,12), position = c(T,NA), gates = c(v,NA))
      d <- density(temp@flow.frame@exprs[,9], na.rm=T)
      fp <- findpeaks(d$y, minpeakheight=.15, minpeakdistance=50)
      if (length(as.vector(fp))>4) { 
        dy2n <- sort(fp[,2])
        m <- d$x[ dy2n[length(dy2n)-1]+which.min(d$y[c(dy2n[length(dy2n)-1]:dy2n[length(dy2n)])]) ]
        temp <- flowDensity(temp, channels = c(9,12), position = c(T,NA), gates = c(m,NA))
        d <- density(temp@flow.frame@exprs[,9], na.rm=T)
        temp <- flowDensity(temp, channels = c(9,12), position = c(NA,F), gates = c(NA,d$x[which.max(d$y)]))
        gt[i,"Ly6C.gate"] <- deGate(temp@flow.frame, channel=c(9), upper=F, tinypeak.removal = .9)+.1
        if (d$x[which.max(d$y)]-d$x[findInterval(gt[i,"Ly6C.gate"], d$x)] <.075) { # if at highest density, move back to minimum density
          gt[i,"Ly6C.gate"] <- m+.1
        } 
      } else {
        gt[i,"Ly6C.gate"] <- deGate(temp@flow.frame, channel=c(9), upper=F, tinypeak.removal = .9)
        if (abs(d$x[which.max(d$y)]-d$x[findInterval(gt[i,"Ly6C.gate"], d$x)])<.1) { gt[i,"Ly6C.gate"] <- v+.1 }
      }
    }
    if (!UpdatedVersion) {
      temp <- flowDensity(notGran.flowD, channels = c(9,12), position = c(NA,T), gates = c(NA, CD11b.gate))
      temp <- flowDensity(temp, channels=c(9,12), position=c(T,NA), gates=c(deGate(temp@flow.frame, channel = c(9))))
      gt[i,"Ly6C.gate"] <- deGate(temp@flow.frame, channel=c(9), upper=F, tinypeak.removal = .9)
      temp1 <- flowDensity(notGran.flowD, channels = c(9,12), position = c(T,T), gates = c(gt[i,"Ly6C.gate"], CD11b.gate))
      if (temp1@proportion>2) { # if too left
        d <- density(temp1@flow.frame@exprs[,9], na.rm=T)
        gt[i,"Ly6C.gate"] <- d$x[which.max(d$y)]
      }
      if (temp1@proportion<.65) { gt[i,"Ly6C.gate"] <- gt[i,"Ly6C.gate"]-.25 }
    }
    
    temp <- flowDensity(notGran.flowD, channels = c(9,12), position = c(T,T), gates = c(gt[i,"Ly6C.gate"], gt[i,"CD11b.gate"]))
    mon.flowD <- flowDensity(temp, channels = c(9,12), position = c(NA, F), gates = c(NA, gt[i,"CD11b.gate.high"]))
#     pr[i,7] <- mon.flowD@proportion <- 100*mon.flowD@cell.count/notGran.flowD@cell.count
#     ev[i,7] <- mon.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(81:96)] <- sapply(1:16, function(x){mean(f.bt@exprs[mon.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(81:96)] <- sapply(1:16, function(x){median(f.bt@exprs[mon.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(81:96)] <- sapply(1:16, function(x){mean(f@exprs[mon.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(81:96)] <- sapply(1:16, function(x){median(f@exprs[mon.flowD@index, x], na.rm = T)}) }
    
    notMon.flowD <- notSubFrame(notGran.flowD@flow.frame, channels = c(9,12), position= "logical", gates = "missing", mon.flowD@filter)
#     pr[i,8] <- notMon.flowD@proportion <- 100-mon.flowD@proportion
#     ev[i,8] <- notMon.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(97:112)] <- sapply(1:16, function(x){mean(f.bt@exprs[notMon.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(97:112)] <- sapply(1:16, function(x){median(f.bt@exprs[notMon.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(97:112)] <- sapply(1:16, function(x){mean(f@exprs[notMon.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(97:112)] <- sapply(1:16, function(x){median(f@exprs[notMon.flowD@index, x], na.rm = T)}) }
    
    plotDens(notGran.flowD@flow.frame, channels = c(9,12), main = "Not Granulocytes", xlab="09. Ly6c", ylab="12. CD11b", xlim=c(0, 4.5), ylim=c(0, 4.5), devn = F); abline(v=gt[i,"Ly6C.gate"], lwd=2); abline(h=gt[i,"CD11b.gate"], lwd=2); abline(h=gt[i,"CD11b.gate.high"], lwd=2)
    lines(mon.flowD@filter, lwd=1); points(notGran.flowD@flow.frame@exprs[notMon.flowD@index,c(9,12)], col=2, cex = 0.1)
    
    rm(lymph.flowD, gran.flowD) 
    
    
    
    ## Gating Not Monocytes. Plotting CD11b(12)_SSC-H(5) for Not/Eosinophils
    if (!UpdatedVersion) {
      gt[i,"ssc.h.gate.high"] <- deGate(notMon.flowD@flow.frame, channel = c(5), use.percentile = T, percentile = .999999) 
      temp <- flowDensity(notMon.flowD, channels = c(12,5), position = c(T,F), gates = c(CD11b.gate, gt[i,"ssc.h.gate.high"]))
      temp <- flowDensity(temp, channels = c(12,5), position = c(F,T), gates = c(gt[i,"CD11b.gate.high"], deGate(notMon, channel = c(5))))
      gt[i,"ssc.h.gate"] <- deGate(temp@flow.frame, channel = c(5), tinypeak.removal = .5)
      if (gt[i,"ssc.h.gate"]>14000) { gt[i,"ssc.h.gate"] <- deGate(temp@flow.frame, channel = c(5)) }
    }
    
    temp <- flowDensity(notMon.flowD, channels = c(12,5), position = c(T,T), gates = c(gt[i,"CD11b.gate"], gt[i,"ssc.h.gate"]))
    eos.flowD <- flowDensity(temp, channels = c(12,5), position = c(F,F), gates = c(gt[i,"CD11b.gate.high"], gt[i,"ssc.h.gate.high"]))
#     pr[i,9] <-  eos.flowD@proportion <- 100*eos.flowD@cell.count/notMon.flowD@cell.count
#     ev[i,9] <- eos.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(113:128)] <- sapply(1:16, function(x){mean(f.bt@exprs[eos.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(113:128)] <- sapply(1:16, function(x){median(f.bt@exprs[eos.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(113:128)] <- sapply(1:16, function(x){mean(f@exprs[eos.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(113:128)] <- sapply(1:16, function(x){median(f@exprs[eos.flowD@index, x], na.rm = T)}) }
    
    if (exists("eos.flowD")) {
      notEos.flowD <- notSubFrame(notMon.flowD@flow.frame, channels = c(12,5), position= "logical", gates = "missing", eos.flowD@filter)
    } else { notEos.flowD <- notMon.flowD }
#     pr[i,10] <- notEos.flowD@proportion <- 100-eos.flowD@proportion
#     ev[i,10] <- notEos.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(129:144)] <- sapply(1:16, function(x){mean(f.bt@exprs[notEos.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(129:144)] <- sapply(1:16, function(x){median(f.bt@exprs[notEos.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(129:144)] <- sapply(1:16, function(x){mean(f@exprs[notEos.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(129:144)] <- sapply(1:16, function(x){median(f@exprs[notEos.flowD@index, x], na.rm = T)}) }
    
    plotDens(notMon.flowD@flow.frame, channels = c(12,5), main = "Not Monocytes", xlab="12. CD11b", ylab="05. SSC-H", xlim=c(0, 4.5), ylim=c(0, 200000), devn = F); abline(v=gt[i,"CD11b.gate"], lwd=2); abline(v=gt[i,"CD11b.gate.high"], lwd=2); abline(h=gt[i,"ssc.h.gate"], lwd=2); abline(h=gt[i,"ssc.h.gate.high"], lwd=2)
    lines(eos.flowD@filter, lwd=1); points(notMon.flowD@flow.frame@exprs[notEos.flowD@index,c(12,5)], col=2, cex = 0.1)
    
    rm(notGran.flowD, mon.flowD)
    
    
    
    ## Gating Not Eosinophils. Plotting CD161(13)_CD19(16) for CD161+/-
    
    if (!UpdatedVersion) {
      gt[i,"CD161.gate"] <- deGate(notEos.flowD@flow.frame, channel = c(13))+.2
      gt[i,"CD161.gate.high"] <- deGate(notEos.flowD@flow.frame, channel = c(13), use.percentile = T, percentile = 0.999999)
      temp <- flowDensity(notEos.flowD, channels = c(13,16), position = c(T,NA), gates = c(gt[i,"CD161.gate"], NA))
      temp <- flowDensity(temp, channels = c(13,16), position = c(F,F), gates = c(gt[i,"CD161.gate.high"], deGate(temp@flow.frame, channel = c(16))))
      gt[i,"CD19.gate"] <- deGate(temp@flow.frame, channel = c(16), tinypeak.removal =.9)+.15
      
      temp <- flowDensity(notEos.flowD, channels = c(13,16), position = c(NA, F), gates = c(NA, gt[i,"CD19.gate"]))
      gt[i,"CD161.gate"] <- deGate(temp@flow.frame, channel=c(13))-.15
      if (gt[i,"CD161.gate"]>2) {gt[i,"CD161.gate"] <- gt[i,"CD161.gate"]-.2}
    }
    
    temp <- flowDensity(notEos.flowD, channels = c(13,16), position = c(T, NA), gates = c(gt[i,"CD161.gate"], NA))
    CD161.flowD <- flowDensity(temp, channels = c(13,16), position = c(F, F), gates = c(gt[i,"CD161.gate.high"], gt[i,"CD19.gate"]))
#     pr[i,11] <-  CD161.flowD@proportion <- 100*CD161.flowD@cell.count/notEos.flowD@cell.count
#     ev[i,11] <- CD161.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(145:160)] <- sapply(1:16, function(x){mean(f.bt@exprs[CD161.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(145:160)] <- sapply(1:16, function(x){median(f.bt@exprs[CD161.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(145:160)] <- sapply(1:16, function(x){mean(f@exprs[CD161.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(145:160)] <- sapply(1:16, function(x){median(f@exprs[CD161.flowD@index, x], na.rm = T)}) }
    
    notCD161.flowD <- notSubFrame(notEos.flowD@flow.frame, channels = c(13,16), position= "logical", gates = "missing", CD161.flowD@filter)
#     pr[i,22] <- notCD161.flowD@proportion <- 100-notCD161.flowD@cell.count
#     ev[i,22] <- notCD161.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(321:336)] <- sapply(1:16, function(x){mean(f.bt@exprs[notCD161.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(321:336)] <- sapply(1:16, function(x){median(f.bt@exprs[notCD161.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(321:336)] <- sapply(1:16, function(x){mean(f@exprs[notCD161.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(321:336)] <- sapply(1:16, function(x){median(f@exprs[notCD161.flowD@index, x], na.rm = T)}) }
    
    plotDens(notEos.flowD@flow.frame, channels = c(13,16), main = "Not Eosinophils", xlab="13. CD161", ylab="16. CD19", xlim=c(0, 4.5), ylim=c(0, 4.5), devn = F); abline(v=gt[i,"CD161.gate"], lwd=2); abline(v=gt[i,"CD161.gate.high"], lwd=2); abline(h=gt[i,"CD19.gate"], lwd=2)
    lines(CD161.flowD@filter, lwd=1); points(notEos.flowD@flow.frame@exprs[CD161.flowD@index,c(13,16)], col=2, cex = 0.1)
    
    
    rm(notMon.flowD, eos.flowD)
    
    
    
    ## Gating CD161+. Plotting CD5(8)_CD11b(13) for NK T-/Cells
    if (!UpdatedVersion) {
      # d <- density(CD161@exprs[,8])
      # t <- CD5.gate2
      # if (CD5.gate2 < d$x[which.max(d$y)] | CD5.gate2-d$x[which.max(d$y)]<.5) {
      #   CD5.gate2 <- d$x[which.max(d$y)]+.3
      # }
    }
    temp <- flowDensity(CD161.flowD, channels = c(8,13), position = c(NA,T), gates = c(NA, gt[i,"CD161.gate"]))
    NK.flowD <- flowDensity(temp, channels = c(8,13), position = c(F,F), gates = c(gt[i,"CD5.gate2"], gt[i,"CD161.gate.high"]))
    NK.flowD@proportion <- 100*NK.flowD@cell.count/CD161.flowD@cell.count
#     pr[i,12] <- NK.flowD@proportion <- 100*NK.flowD@cell.count/CD161.flowD@cell.count
#     ev[i,12] <- NK.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(161:176)] <- sapply(1:16, function(x){mean(f.bt@exprs[NK.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(161:176)] <- sapply(1:16, function(x){median(f.bt@exprs[NK.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(161:176)] <- sapply(1:16, function(x){mean(f@exprs[NK.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(161:176)] <- sapply(1:16, function(x){median(f@exprs[NK.flowD@index, x], na.rm = T)}) }
    
    temp <- flowDensity(CD161.flowD, channels = c(8,13), position = c(T,T), gates = c(gt[i,"CD5.gate2"], gt[i,"CD161.gate"]))
    NKT.flowD <- flowDensity(temp, channels = c(8,13), position = c(F,F), gates = c(gt[i,"CD5.gate.high"], gt[i,"CD161.gate.high"]))
#     pr[i,17] <- NKT.flowD@proportion <- 100*NKT.flowD@cell.count/CD161.flowD@cell.count
#     ev[i,17] <- NKT.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(241:256)] <- sapply(1:16, function(x){mean(f.bt@exprs[NKT.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(241:256)] <- sapply(1:16, function(x){median(f.bt@exprs[NKT.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(241:256)] <- sapply(1:16, function(x){mean(f@exprs[NKT.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(241:256)] <- sapply(1:16, function(x){median(f@exprs[NKT.flowD@index, x], na.rm = T)}) }
    
    plotDens(CD161.flowD@flow.frame, channels = c(8,13), main = "CD161+", xlab="08. CD5 #2~", ylab="13. CD161", xlim=c(0, 4.5), ylim=c(0, 4.5), devn = F); abline(v=gt[i,"CD5.gate"], lty=2); abline(v=gt[i,"CD5.gate2"], lwd=2); abline(v=gt[i,"CD5.gate.high"], lwd=2); abline(h=gt[i,"CD161.gate"], lwd=2); abline(h=gt[i,"CD161.gate.high"], lwd=2)
    lines(NK.flowD@filter, lwd=1); points(CD161.flowD@flow.frame@exprs[NK.flowD@index,c(8,13)], col=2, cex = 0.1)
    lines(NKT.flowD@filter, lwd=1); points(CD161.flowD@flow.frame@exprs[NKT.flowD@index,c(8,13)], col=2, cex = 0.1)
    
    rm(notEos.flowD)
    
    
    
    ## Gating NK Cells. Plotting CD11b(12)_Ly6C(9)
    if (!UpdatedVersion) {
      gt[i,"Ly6C.gate2"] <- deGate(NK.flowD@flow.frame, channel=c(9))-.1
      temp <- flowDensity(NK.flowD, channels=c(9,12), position=c(F,NA), gates = c(gt[i,"Ly6C.gate2"],NA)) # left half
      gt[i,"CD11b.gate2"] <- deGate(temp@flow.frame, channel=c(12))-.1 #use left half for horizontal gate
      if (flowDensity(temp, channels= c(9,12), position=c(NA,T), gates=c(NA,gt[i,"CD11b.gate2"]))@proportion < 30) { # if left top half too small
        temp <- flowDensity(temp, channels= c(9,12), position=c(NA,F), gates=c(gt[i,"CD11b.gate2"],NA)) # cut bottom left half
        d <- density(temp@flow.frame@exprs[,12], na.rm=T)
        gt[i,"CD11b.gate2"] <- d$x[which.max(d$y)]+.2
      }
      temp <- flowDensity(NK.flowD, channels=c(9,12), position=c(NA,T), gates = c(NA,gt[i,"CD11b.gate2"])) # right half
      gt[i,"Ly6C.gate2"] <- deGate(temp@flow.frame, channel=c(9))+.2
    }
    NKimmatLy6C.flowD <- flowDensity(NK.flowD, channels = c(9,12), position = c(T,F), gates = c(gt[i,"Ly6C.gate2"],gt[i,"CD11b.gate2"]))
#     pr[i,13] <- NKimmatLy6C.flowD@proportion
#     ev[i,13] <- NKimmatLy6C.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(177:192)] <- sapply(1:16, function(x){mean(f.bt@exprs[NKimmatLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(177:192)] <- sapply(1:16, function(x){median(f.bt@exprs[NKimmatLy6C.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(177:192)] <- sapply(1:16, function(x){mean(f@exprs[NKimmatLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(177:192)] <- sapply(1:16, function(x){median(f@exprs[NKimmatLy6C.flowD@index, x], na.rm = T)}) }
    
    NKmatLy6C.flowD <- flowDensity(NK.flowD, channels = c(9,12), position = c(T,T), gates = c(gt[i,"Ly6C.gate2"],gt[i,"CD11b.gate2"]))
#     pr[i,14] <- NKmatLy6C.flowD@proportion
#     ev[i,14] <- NKmatLy6C.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(193:208)] <- sapply(1:16, function(x){mean(f.bt@exprs[NKmatLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(193:208)] <- sapply(1:16, function(x){median(f.bt@exprs[NKmatLy6C.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(193:208)] <- sapply(1:16, function(x){mean(f@exprs[NKmatLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(193:208)] <- sapply(1:16, function(x){median(f@exprs[NKmatLy6C.flowD@index, x], na.rm = T)}) }
    
    NKimmatNotLy6C.flowD <- flowDensity(NK.flowD, channels = c(9,12), position = c(F,F), gates = c(gt[i,"Ly6C.gate2"],gt[i,"CD11b.gate2"]))
#     pr[i,15] <- NKimmatNotLy6C.flowD@proportion
#     ev[i,15] <- NKimmatNotLy6C.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(209:224)] <- sapply(1:16, function(x){mean(f.bt@exprs[NKimmatNotLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(209:224)] <- sapply(1:16, function(x){median(f.bt@exprs[NKimmatNotLy6C.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(209:224)] <- sapply(1:16, function(x){mean(f@exprs[NKimmatNotLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(209:224)] <- sapply(1:16, function(x){median(f@exprs[NKimmatNotLy6C.flowD@index, x], na.rm = T)}) }
    
    NKmatNotLy6C.flowD <- flowDensity(NK.flowD, channels = c(9,12), position = c(F,T), gates = c(gt[i,"Ly6C.gate2"],gt[i,"CD11b.gate2"]))
#     pr[i,16] <- NKmatNotLy6C.flowD@proportion
#     ev[i,16] <- NKmatNotLy6C.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(225:240)] <- sapply(1:16, function(x){mean(f.bt@exprs[NKmatNotLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(225:240)] <- sapply(1:16, function(x){median(f.bt@exprs[NKmatNotLy6C.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(225:240)] <- sapply(1:16, function(x){mean(f@exprs[NKmatNotLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(225:240)] <- sapply(1:16, function(x){median(f@exprs[NKmatNotLy6C.flowD@index, x], na.rm = T)}) }
    
    plotDens(NK.flowD@flow.frame, channels = c(9,12), main = "NK-cell", ylab="12. CD11b #2", xlab="09. Ly6c #2", xlim=c(0, 4.5), ylim=c(0, 4.5), devn = F); abline(h=gt[i,"CD11b.gate"], lty=2); abline(h=gt[i,"CD11b.gate2"], lwd=2); abline(v=gt[i,"Ly6C.gate"], lty=2); abline(v=gt[i,"Ly6C.gate2"], lwd=2)
    lines(NKimmatLy6C.flowD@filter, lwd=1)
    lines(NKmatNotLy6C.flowD@filter, lwd=1)
    
    rm(NKimmatNotLy6C.flowD, NKimmatLy6C.flowD, NKmatNotLy6C.flowD, NKmatLy6C.flowD)
    
    
    
    ## Gating NK T-cells. Plotting CD11b(12)_Ly6C(9)
    if (!UpdatedVersion) {
      temp <- flowDensity(NKT.flowD, channels=c(9,12), position=c(NA,F), gates=c(NA,gt[i,"CD11b.gate2"]))
      gt[i,"Ly6C.gate3"] <- deGate(temp@flow.frame, channel=c(9))+.2
    }
    NKTnotCD11bLy6C.flowD <- flowDensity(NKT.flowD, channels = c(9,12), position = c(T,F), gates = c(gt[i,"Ly6C.gate3"],gt[i,"CD11b.gate2"]))
#     pr[i,18] <- NKTnotCD11bLy6C.flowD@proportion
#     ev[i,18] <- NKTnotCD11bLy6C.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(257:272)] <- sapply(1:16, function(x){mean(f.bt@exprs[NKTnotCD11bLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(257:272)] <- sapply(1:16, function(x){median(f.bt@exprs[NKTnotCD11bLy6C.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(257:272)] <- sapply(1:16, function(x){mean(f@exprs[NKTnotCD11bLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(257:272)] <- sapply(1:16, function(x){median(f@exprs[NKTnotCD11bLy6C.flowD@index, x], na.rm = T)}) }
    
    NKTCD11bLy6C.flowD <- flowDensity(NKT.flowD, channels = c(9,12), position = c(T,T), gates = c(gt[i,"Ly6C.gate3"],gt[i,"CD11b.gate2"]))
#     pr[i,19] <- NKTCD11bLy6C.flowD@proportion
#     ev[i,19] <- NKTCD11bLy6C.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(273:288)] <- sapply(1:16, function(x){mean(f.bt@exprs[NKTCD11bLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(273:288)] <- sapply(1:16, function(x){median(f.bt@exprs[NKTCD11bLy6C.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(273:288)] <- sapply(1:16, function(x){mean(f@exprs[NKTCD11bLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(273:288)] <- sapply(1:16, function(x){median(f@exprs[NKTCD11bLy6C.flowD@index, x], na.rm = T)}) }
    
    NKTnotCD11bNotLy6C.flowD <- flowDensity(NKT.flowD, channels = c(9,12), position = c(F,F), gates = c(gt[i,"Ly6C.gate3"],gt[i,"CD11b.gate2"]))
#     pr[i,20] <- NKTnotCD11bNotLy6C.flowD@proportion
#     ev[i,20] <- NKTnotCD11bNotLy6C.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(289:304)] <- sapply(1:16, function(x){mean(f.bt@exprs[NKTnotCD11bNotLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(289:304)] <- sapply(1:16, function(x){median(f.bt@exprs[NKTnotCD11bNotLy6C.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(289:304)] <- sapply(1:16, function(x){mean(f@exprs[NKTnotCD11bNotLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(289:304)] <- sapply(1:16, function(x){median(f@exprs[NKTnotCD11bNotLy6C.flowD@index, x], na.rm = T)}) }
    
    NKTCD11bNotLy6C.flowD <- flowDensity(NKT.flowD, channels = c(9,12), position = c(F,T), gates = c(gt[i,"Ly6C.gate3"],gt[i,"CD11b.gate2"]))
#     pr[i,21] <- NKTCD11bNotLy6C.flowD@proportion
#     ev[i,21] <- NKTCD11bNotLy6C.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(305:320)] <- sapply(1:16, function(x){mean(f.bt@exprs[NKTCD11bNotLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(305:320)] <- sapply(1:16, function(x){median(f.bt@exprs[NKTCD11bNotLy6C.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(305:320)] <- sapply(1:16, function(x){mean(f@exprs[NKTCD11bNotLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(305:320)] <- sapply(1:16, function(x){median(f@exprs[NKTCD11bNotLy6C.flowD@index, x], na.rm = T)}) }
    
    plotDens(NKT.flowD@flow.frame, channels = c(9,12), main = "NKT-cell", ylab="12. CD11b #2", xlab="09. Ly6c #3", xlim=c(0, 4.5), ylim=c(0, 4.5), devn = F); abline(h=gt[i,"CD11b.gate"], lty=2); abline(h=gt[i,"CD11b.gate2"], lwd=2); abline(v=gt[i,"Ly6C.gate"], lty=2); abline(v=gt[i,"Ly6C.gate2"], lty=2); abline(v=gt[i,"Ly6C.gate3"], lwd=2)
    lines(NKTCD11bLy6C.flowD@filter, lwd=1)
    lines(NKTnotCD11bNotLy6C.flowD@filter, lwd=1)
    
    rm(CD161.flowD, NK.flowD, NKT.flowD, NKTCD11bNotLy6C.flowD, NKTnotCD11bNotLy6C.flowD, NKTCD11bLy6C.flowD, NKTnotCD11bLy6C.flowD)
    
    
    
    ## Gating not CD161+. Plotting mhcii(7)_CD5(8) for Not/T-cells
    if (!UpdatedVersion) {
      gt[i,"mhcii.gate"] <- deGate(notCD161.flowD@flow.frame, channel=c(7))
      gt[i,"mhcii.gate.low"] <- deGate(notCD161.flowD@flow.frame, channel=c(7), upper=F)
      
      temps <- notCD161.flowD
      temps@flow.frame <- rotate.data(notCD161.flowD@flow.frame,c(7,8), theta=-pi/8)$data
      gt[i,"mhcii.gate.slant"] <- deGate(temps@flow.frame, channel=c(7))
      
      temp <- flowDensity(temps, channels = c(7,8), position = c(F,NA), gates = c(gt[i,"mhcii.gate.slant"], NA))
      temp@filter <- rotate.data(temp@filter,c(7,8),theta = pi/8)$data; temp@flow.frame <- rotate.data(temp@flow.frame,c(7,8),theta = pi/8)$data
      
      gt[i,"CD5.gate2"] <- deGate(temp@flow.frame, channel=c(8))+.15
      if (gt[i,"CD5.gate2"]>2.5) { gt[i,"CD5.gate2"] <- deGate(temp@flow.frame, channel=c(8), upper=F, tinypeak.removal = .5)  }
      if (gt[i,"CD5.gate2"]>3) {gt[i,"CD5.gate2"] <- 1.7}
    } else {
      temp <- flowDensity(rotate.data(notCD161.flowD@flow.frame,c(7,8), theta=-pi/8)$data, channels = c(7,8), position = c(F,NA), gates = c(gt[i,"mhcii.gate.slant"], NA))
      temp@filter <- rotate.data(temp@filter,c(7,8),theta = pi/8)$data; temp@flow.frame <- rotate.data(temp@flow.frame,c(7,8),theta = pi/8)$data
    }
    Tcell.flowD <- flowDensity(temp, channels = c(7,8), position = c(T,T), gates = c(gt[i,"mhcii.gate.low"], gt[i,"CD5.gate2"]))
#     pr[i,23] <- Tcell.flowD@proportion <- 100*Tcell.flowD@cell.count/notCD161.flowD@cell.count
#     ev[i,23] <- Tcell.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(337:352)] <- sapply(1:16, function(x){mean(f.bt@exprs[Tcell.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(337:352)] <- sapply(1:16, function(x){median(f.bt@exprs[Tcell.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(337:352)] <- sapply(1:16, function(x){mean(f@exprs[Tcell.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(337:352)] <- sapply(1:16, function(x){median(f@exprs[Tcell.flowD@index, x], na.rm = T)}) }
    
    notTcell.flowD <- notSubFrame(notCD161.flowD@flow.frame, channels = c(7,8), position= "logical", gates = "missing", Tcell.flowD@filter)
#     pr[i,25] <- notTcell.flowD@proportion <- 100-Tcell.flowD@proportion
#     ev[i,25] <- notTcell.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(369:384)] <- sapply(1:16, function(x){mean(f.bt@exprs[notTcell.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(369:384)] <- sapply(1:16, function(x){median(f.bt@exprs[notTcell.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(369:384)] <- sapply(1:16, function(x){mean(f@exprs[notTcell.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(369:384)] <- sapply(1:16, function(x){median(f@exprs[notTcell.flowD@index, x], na.rm = T)}) }
    
    plotDens(notCD161.flowD@flow.frame, channels = c(7,8), main = "Not CD161+", xlab="07. MHCII", ylab="08. CD5 #2", xlim=c(0, 4.5), ylim=c(0, 4.5), devn = F); abline(v=gt[i,"mhcii.gate.low"], lwd=2); abline(v=gt[i,"mhcii.gate"], lty=2); abline(h=gt[i,"CD5.gate"], lty=2); abline(h=gt[i,"CD5.gate2"], lwd=2)
    lines(Tcell.flowD@filter, lwd=1); points(notCD161.flowD@flow.frame@exprs[Tcell.flowD@index,c(7,9)], col=2, cex = 0.1)
    
    
    
    ## Gating T-cells. Plotting Ly6C(9) for Ly6C+ T-cells
    if (fixTcell) {
      d <- density(Tcell.flowD@flow.frame@exprs[,9], na.rm=T)
      tryCatch({fp <- findpeaks(d$y, npeaks=3, minpeakheight=.15, minpeakdistance=50)}, error=function(err) {fp <- findpeaks(d$y, minpeakheight=.15, minpeakdistance=50)})
      if (is.matrix(fp)) { dy2n <- sort(fp[,2]) } else { dy2n <- fp[2] }
      if (length(dy2n)==3) {
        gt[i,"Ly6C.gate4"] <- d$x[ dy2n[1]+which.min(d$y[c(dy2n[1]:dy2n[2])]) ]-.15
      } else if (length(dy2n)==2) {
        m <- d$x[ dy2n[1]+which.min(d$y[c(dy2n[1]:dy2n[2])]) ]
        temp <- flowDensity(Tcell.flowD, channels=c(9,8), position=c(F,NA), gates=c(m,NA))
        temp <- flowDensity(temp, channels=c(9,8), position=c(T,NA), gates=c(d$x[dy2n[1]],NA))
        a <- deGate(temp@flow.frame, channel = c(9), alpha=.5)
        if (a<m & abs(a-d$x[dy2n[1]])>.25) {
          gt[i,"Ly6C.gate4"] <- a-.1
        } else { gt[i,"Ly6C.gate4"] <- m-.1}
      } else {
        temp <- flowDensity(Tcell.flowD, channels = c(9,8), position = c(T,NA), gates = c(d$x[dy2n[1]], NA))
        temp <- flowDensity(temp, channels = c(9,8), position = c(F,NA), gates = c(2.6,NA))
        a <- deGate(temp@flow.frame, channel = c(9), alpha=.5)-.1
        if (a>2.6) {
          temp <- flowDensity(temp, channels = c(9,8), position = c(F,NA), gates = c(a,NA))
          a <- deGate(temp@flow.frame, channel = c(9), tinypeak.removal=.1)-.1
        }
        b <- deGate(Tcell.flowD@flow.frame, channel = c(9), upper=F)
        if (b<d$x[dy2n[1]]) {
          b <- d$x[dy2n[1]]+(d$x[dy2n[1]]-b)
        }
        if (abs(a-1.95)< abs(b-1.95)) {
          gt[i,"Ly6C.gate4"] <- a
        } else { gt[i,"Ly6C.gate4"] <- b }
      }
    }
    if (!UpdatedVersion) {
      d <- density(Tcell.flowD@flow.frame@exprs[,9])
      fp <- findpeaks(d$y, npeaks=3)
      if (nrow(fp)>1) {
        dy2n <- sort(fp[,2])
      } else {
        dy2n <- sort(fp[2])
      }
      gt[i,"Ly6C.gate4"]  <- d$x[ dy2n[1]+which.min(d$y[dy2n[1]:dy2n[2]])]-.15
      #         gt[i,"Ly6C.gate.high"] <- deGate(Tcell.flowD@flow.frame, channel = c(9), use.percentile = T, percentile = .999999)
      #         temp <- flowDensity(Tcell.flowD, channels = c(9,8), position = c(T,NA), gates = c(gt[i,"Ly6C.gate2"],NA))
      #         if (temp@proportion<45 & temp@proportion>20) {
      #           gt[i,"Ly6C.gate4"] <- gt[i,"Ly6C.gate2"]
      #         } else {
      #           de <- deGate(Tcell.flowD@flow.frame, channel = c(9), tinypeak.removal = .6)-.1
      #           temp <- flowDensity(Tcell.flowD, channels = c(9,8), position = c(T,NA), gates = c(de,NA))
      #           if (temp@proportion<40 & temp@proportion>20) {
      #             gt[i,"Ly6C.gate4"] <- de
      #           } else {
      #             gt[i,"Ly6C.gate4"] <- mean(c(de-.2,gt[i,"Ly6C.gate2"]))
      #           }
      #         }
    }
    temp <- flowDensity(Tcell.flowD@flow.frame, channels = c(9,8), position = c(T,NA), gates = c(gt[i,"Ly6C.gate4"],NA))
    TcellLy6C.flowD <- flowDensity(temp, channels = c(9,8), position = c(F,NA), gates = c(gt[i,"Ly6C.gate.high"],NA))
#     pr[i,24] <- TcellLy6C.flowD@proportion <- 100*TcellLy6C.flowD@cell.count/Tcell.flowD@cell.count
#     ev[i,24] <- TcellLy6C.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(353:368)] <- sapply(1:16, function(x){mean(f.bt@exprs[TcellLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(353:368)] <- sapply(1:16, function(x){median(f.bt@exprs[TcellLy6C.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(353:368)] <- sapply(1:16, function(x){mean(f@exprs[TcellLy6C.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(353:368)] <- sapply(1:16, function(x){median(f@exprs[TcellLy6C.flowD@index, x], na.rm = T)}) }
    
    plot(density(Tcell.flowD@flow.frame@exprs[,9], na.rm=T), main = "T-cell", xlab="09. Ly6C #4", xlim=c(0, 4.5), ylab="Count/Density"); abline(v=gt[i,"Ly6C.gate.high"], lwd=2); abline(v=gt[i,"Ly6C.gate"], lty=2); abline(v=gt[i,"Ly6C.gate4"], lwd=2)
    
    rm(TcellLy6C.flowD)
    
    
    
    ## Gating Not T-cells. Plotting CD19(16)_CD11c(14) for cDC/B-cell
    
    temps <- notTcell.flowD
    temps@flow.frame <- rotate.data(notTcell.flowD@flow.frame, c(16,14), theta=pi/4)$data
    if (!UpdatedVersion) {
      gt[i,"CD19.gate"] <- deGate(notTcell.flowD@flow.frame, channel = c(16)) #, use.percentile = T, percentile = 0.9999999)
      d0 <- density(notTcells.flowD@flow.frame@exprs[,14], na.rm=T) # help to find max
      gt[i,"CD11c.gate.slant.low"] <- deGate(notTcells.flowD@flow.frame, channel=c(14), tinypeak.removal = .2)
      temp <- flowDensity(notTcells.flowD, channels=c(16,14), position=c(NA,T), gates=c(NA,gt[i,"CD11c.gate.slant.low"]))
      gt[i,"CD11c.gate.slant.high"] <- deGate(temp@flow.frame, channel=c(14))-.25
      gt[i,"CD11c.gate"] <- deGate(notTcell.flowD@flow.frame, channel=c(14))+.2
      gt[i,"CD11c.gate.high"] <- deGate(notTcell.flowD@flow.frame, channel=c(14), use.percentile = T, percentile = 0.9999999)
      temp <- flowDensity(notTcell.flowD, channels=c(16,14), position=c(T,NA), gates=c(gt[i,"CD19.gate"],NA))
      gt[i,"CD19.gate2"] <- deGate(temp@flow.frame, channel=c(16), upper=F, tinypeak.removal = .4)
      if (abs(gt[i,"CD19.gate"]-gt[i,"CD19.gate2"]) <.8) {gt[i,"CD19.gate2"] <- gt[i,"CD19.gate"]}
    }
    
    temp <- flowDensity(temps, channels = c(16,14), position = c(NA,T), gates = c(NA, gt[i,"CD11c.gate.slant.high"]))
    temp@filter <- rotate.data(temp@filter,c(16,14),theta = -pi/4)$data; temp@flow.frame <- rotate.data(temp@flow.frame,c(16,14),theta = -pi/4)$data
    temp <- flowDensity(temp, channels = c(16,14), position = c(NA,T), gates = c(NA, gt[i,"CD11c.gate"]))
    cDC.flowD <- flowDensity(temp, channels = c(16,14), position = c(F,F), gates = c(gt[i,"CD19.gate"], gt[i,"CD11c.gate.high"]))
#     pr[i,26] <- cDC.flowD@proportion <- 100*cDC.flowD@cell.count/notTcell.flowD@cell.count
#     ev[i,26] <- cDC.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(385:400)] <- sapply(1:16, function(x){mean(f.bt@exprs[cDC.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(385:400)] <- sapply(1:16, function(x){median(f.bt@exprs[cDC.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(385:400)] <- sapply(1:16, function(x){mean(f@exprs[cDC.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(385:400)] <- sapply(1:16, function(x){median(f@exprs[cDC.flowD@index, x], na.rm = T)}) }
    
    temp <- flowDensity(temps, channels = c(16,14), position = c(NA,F), gates = c(NA, gt[i,"CD11c.gate.slant.low"]))
    temp@filter <- rotate.data(temp@filter,c(16,14),theta = -pi/4)$data; temp@flow.frame <- rotate.data(temp@flow.frame,c(16,14),theta = -pi/4)$data
    Bcell.flowD <- flowDensity(temp, channels = c(16,14), position = c(T,F), gates = c(gt[i,"CD19.gate2"], gt[i,"CD11c.gate"]))
#     pr[i,29] <- Bcell.flowD@proportion <- 100*Bcell.flowD@cell.count/notTcell.flowD@cell.count
#     ev[i,29] <- Bcell.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(433:448)] <- sapply(1:16, function(x){mean(f.bt@exprs[Bcell.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(433:448)] <- sapply(1:16, function(x){median(f.bt@exprs[Bcell.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(433:448)] <- sapply(1:16, function(x){mean(f@exprs[Bcell.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(433:448)] <- sapply(1:16, function(x){median(f@exprs[Bcell.flowD@index, x], na.rm = T)}) }
    
    plotDens(notTcell.flowD@flow.frame, channels = c(16,14), main = "Not T-cell", xlab="16. CD19", ylab="14. CD11c", xlim=c(0, 4.5), ylim=c(0, 4.5), devn = F); abline(v=gt[i,"CD19.gate"], lwd=2); abline(v=gt[i,"CD19.gate2"], lwd=2); abline(h=gt[i,"CD11c.gate"] , lwd=2); abline(h=gt[i,"CD11c.gate.high"], lwd=2)
    lines(cDC.flowD@filter, lwd=1)
    lines(Bcell.flowD@filter, lwd=1)
    
    rm(notCD161.flowD, Tcell.flowD)
    
    
    
    ## Gating cDC. Plotting CD11b(12)_mhcii(7)
    
    if (!UpdatedVersion) {
      gt[i,"CD11b.gate.lowlow"] <- deGate(cDC.flowD@flow.frame, channel = c(12), use.percentile = T, percentile = 0.000000000000000000000001)
      gt[i,"CD11b.gate3"] <- deGate(cDC.flowD@flow.frame, channel = c(12))
    }
    temp <- flowDensity(cDC.flowD, channels = c(12,7), position = c(T,T), gates = c(gt[i,"CD11b.gate.lowlow"], gt[i,"mhcii.gate"]))
    cDCCD8.flowD <- flowDensity(temp, channels = c(12,7), position = c(F,NA), gates = c(gt[i,"CD11b.gate3"], NA))
#     pr[i,27] <- cDCCD8.flowD@proportion <- 100*cDCCD8.flowD@cell.count/cDC.flowD@cell.count
#     ev[i,27] <- cDCCD8.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(401:416)] <- sapply(1:16, function(x){mean(f.bt@exprs[cDCCD8.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(401:416)] <- sapply(1:16, function(x){median(f.bt@exprs[cDCCD8.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(401:416)] <- sapply(1:16, function(x){mean(f@exprs[cDCCD8.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(401:416)] <- sapply(1:16, function(x){median(f@exprs[cDCCD8.flowD@index, x], na.rm = T)}) }
    
    temp <- flowDensity(cDC.flowD, channels = c(12,7), position = c(T,T), gates = c(gt[i,"CD11b.gate3"], gt[i,"mhcii.gate"]))
    cDCCD11b.flowD <- flowDensity(temp, channels = c(12,7), position = c(F,NA), gates = c(gt[i,"CD11b.gate.high"], NA))
#     pr[i,28] <- cDCCD11b.flowD@proportion <- 100*cDCCD11b.flowD@cell.count/cDC.flowD@cell.count
#     ev[i,28] <- cDCCD11b.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(417:432)] <- sapply(1:16, function(x){mean(f.bt@exprs[cDCCD11b.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(417:432)] <- sapply(1:16, function(x){median(f.bt@exprs[cDCCD11b.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(417:432)] <- sapply(1:16, function(x){mean(f@exprs[cDCCD11b.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(417:432)] <- sapply(1:16, function(x){median(f@exprs[cDCCD11b.flowD@index, x], na.rm = T)}) }
    
    plotDens(cDC.flowD@flow.frame, channels = c(12,7), main = "cDC", xlab="12. CD11b #3", ylab="07. MHCII", xlim=c(0, 4.5), ylim=c(0, 4.5), devn = F); abline(v=gt[i,"CD11b.gate.lowlow"], lwd=2); abline(v=gt[i,"CD11b.gate2"], lty=2); abline(v=gt[i,"CD11b.gate3"], lwd=2); abline(v=gt[i,"CD11b.gate"], lty=2); abline(v=gt[i,"CD11b.gate.high"], lwd=2); abline(h=gt[i,"mhcii.gate"], lwd=2)
    lines(cDCCD8.flowD@filter, lwd=1)
    lines(cDCCD11b.flowD@filter, lwd=1)
    
    rm(cDCCD8.flowD, cDCCD11b.flowD)
    
    
    
    ## Gating B-cell. Plotting CD5(8)_CD21(15) for B1/2 B-cells
    temps <- Bcell.flowD
    temps@flow.frame <- rotate.data(Bcell.flowD@flow.frame, c(8,15), theta=-pi/16)$data
    if (!UpdatedVersion) {
      gt[i,"CD5.gate.slant"] <- deGate(temps@flow.frame, channel=c(8))+.4
    }
    temp <- flowDensity(temps, channels = c(8,15), position = c(T,NA), gates = c(gt[i,"CD5.gate.slant"], NA))
    temp@filter <- rotate.data(temp@filter,c(8,15),theta = pi/16)$data; temp@flow.frame <- rotate.data(temp@flow.frame,c(8,15),theta = pi/16)$data
    BcellB1.flowD <- flowDensity(temp, channels = c(8,15), position = c(F,F), gates = c(gt[i,"CD5.gate.high"], gt[i,"CD21.gate.high"]))
#     pr[i,30] <- BcellB1.flowD@proportion <- 100*BcellB1.flowD@cell.count/Bcell.flowD@cell.count
#     ev[i,30] <- BcellB1.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(449:464)] <- sapply(1:16, function(x){mean(f.bt@exprs[BcellB1.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(449:464)] <- sapply(1:16, function(x){median(f.bt@exprs[BcellB1.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(449:464)] <- sapply(1:16, function(x){mean(f@exprs[BcellB1.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(449:464)] <- sapply(1:16, function(x){median(f@exprs[BcellB1.flowD@index, x], na.rm = T)}) }
    
    BcellB2.flowD <- notSubFrame(Bcell.flowD@flow.frame, channels=c(8,15), position= "logical", gates = "missing", BcellB1.flowD@filter)
#     pr[i,31] <- BcellB2.flowD@proportion <- 100*BcellB2.flowD@cell.count/Bcell.flowD@cell.count
#     ev[i,31] <- BcellB2.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(465:480)] <- sapply(1:16, function(x){mean(f.bt@exprs[BcellB2.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(465:480)] <- sapply(1:16, function(x){median(f.bt@exprs[BcellB2.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(465:480)] <- sapply(1:16, function(x){mean(f@exprs[BcellB2.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(465:480)] <- sapply(1:16, function(x){median(f@exprs[BcellB2.flowD@index, x], na.rm = T)}) }
    
    plotDens(Bcell.flowD@flow.frame, channels = c(8,15), main = "B-cell", xlab="08. CD5 #2", ylab="15. CD21", xlim=c(0, 4.5), ylim=c(0, 4.5), devn = F); abline(v=gt[i,"CD5.gate2"], lty=2); abline(v=gt[i,"CD5.gate"], lty=2); abline(v=gt[i,"CD5.gate.high"], lwd=2); abline(h=gt[i,"CD21.gate.high"], lwd=2)
    lines(BcellB1.flowD@filter, lwd=1)
    
    rm(notTcell.flowD, cDC.flowD)
    
    
    
    ## Gating B2 B-cells. Plotting CD23(11)_CD21(15)
    temps <- temps2 <- BcellB2.flowD
    temps@flow.frame <- rotate.data(BcellB2.flowD@flow.frame,c(11,15), theta=-pi/8)$data
    temps2@flow.frame <- rotate.data(BcellB2.flowD@flow.frame,c(11,15), theta=pi/8)$data
    if (fixB2 & !i%in%fixB2L) {
      fixB2L <- append(fixB2L,i)
      d <- density(temps@flow.frame@exprs[,15], na.rm=T)
      fp <- findpeaks(d$y, minpeakheight=.5)
      if (is.matrix(fp)) { if (nrow(fp)>1) {
        ma <- d$x[ fp[ which.max(fp[,1]), 2] ]
      } else { ma <- d$x[as.vector(fp)[2]]
      } } else { ma <- d$x[fp[2]] }
      temp <- flowDensity(temps, channels = c(11,15), position = c(NA,F), gates = c(NA, ma))
      gt[i,"CD21.gate.slant.low"] <- deGate(temp@flow.frame, channel=c(15), upper=F, tinypeak.removal=.6)+.15
      
      d <- density(temps2@flow.frame@exprs[,15], na.rm=T)
      fp <- findpeaks(d$y, minpeakheight=.5)
      if (is.matrix(fp)) { if (nrow(fp)>1) {
        ma <- d$x[ fp[ which.max(fp[,1]), 2] ]
      } else { ma <- d$x[as.vector(fp)[2]]
      } } else { ma <- d$x[fp[2]] }
      gt[i,"CD21.gate.slant.high"] <- deGate(temps2@flow.frame, channel=c(15))
      if (gt[i,"CD21.gate.slant.high"]<ma) {
        temp <- flowDensity(temps2, channels = c(11,15), position = c(NA,T), gates = c(NA, ma-.2))
        gt[i,"CD21.gate.slant.high"] <- deGate(temp@flow.frame, channel=c(15))
      }
      
    }
    if (!UpdatedVersion) {
      gt[i,"CD23.gate"] <- deGate(BcellB2, channel = c(11))
      gt[i,"CD23.gate.lowlow"] <- deGate(BcellB2, channel = c(11), upper = F, use.percentile = T, percentile = .0001)
      gt[i,"CD23.gate.low"] <- deGate(BcellB2, channel = c(11), upper = F, use.percentile = T, percentile = .01)
      gt[i,"CD23.gate.high"] <- deGate(BcellB2, channel = c(11), use.percentile = T, percentile = 0.9999999)
      d1 <- density(temps@flow.frame@exprs[,15], na.rm=T)
      d2 <- density(temps@flow.frame@exprs[,11], na.rm=T)
      gt[i,"CD23.gate.slant"] <- d2$x[which.max(d2$y)]
      if (gt[i,"CD23.gate.slant"]<1.5) {gt[i,"CD23.gate.slant"] <- 1.8}
      gt[i,"CD21.gate.slant.low"] <- deGate(temps@flow.frame, channel=c(15))+.5
      if (gt[i,"CD21.gate.slant.low"] >= d1$x[which.max(d1$y)]) { #if too high
        gt[i,"CD21.gate.slant.low"] <- deGate(temps@flow.frame, channel=c(15), use.percentile = T, percentile = .11)-.2
        if (gt[i,"CD21.gate.slant.low"]<2) {gt[i,"CD21.gate.slant.low"] <- 2.2}
      }
      # plotDens(BcellB2s, channels = c(11,15), main = "B-cell B2", devn = F); abline(v=gt[i,"CD23.gate.slant"], lwd=2); abline(h=gt[i,"CD21.gate.slant.low"], lwd=2)
      
      gt[i,"CD21.gate.slant.high"] <- deGate(temps2@flow.frame, channel=c(15))+.25
      # plotDens(BcellB2s2, channels=c(11,15), main = "B-cell B2", devn = F); abline(h=gt[i,"CD21.gate.slant.high"], lwd=2)
      
      d <- density(temps2@flow.frame@exprs[,15], na.rm=T)
      gt[i,"CD21.gate"] <- d$x[which.max(d$y)]
      gt[i,"CD21.gate.low"] <- deGate(BcellB2.flowD@flow.frame, channel=c(15), use.percentile = T, percentile = 0.000001)
    }
    
    temp <- flowDensity(temps2, channels = c(11,15), position = c(NA,T), gates = c(NA, gt[i,"CD21.gate.slant.high"]))
    temp@filter <- rotate.data(temp@filter,c(11,15),theta = -pi/8)$data; temp@flow.frame <- rotate.data(temp@flow.frame,c(11,15),theta = -pi/8)$data
    temp <- flowDensity(temp, channels = c(11,15), position = c(T,NA), gates = c(gt[i,"CD23.gate.lowlow"], NA))
    MZB.flowD <- flowDensity(temp, channels = c(11,15), position = c(F,NA), gates = c(gt[i,"CD23.gate.high"], NA))
#     pr[i,32] <- MZB.flowD@proportion <- 100*MZB.flowD@cell.count/BcellB2.flowD@cell.count
#     ev[i,32] <- MZB.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(481:496)] <- sapply(1:16, function(x){mean(f.bt@exprs[MZB.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(481:496)] <- sapply(1:16, function(x){median(f.bt@exprs[MZB.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(481:496)] <- sapply(1:16, function(x){mean(f@exprs[MZB.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(481:496)] <- sapply(1:16, function(x){median(f@exprs[MZB.flowD@index, x], na.rm = T)}) }
    
    temp <- flowDensity(temps, channels = c(11,15), position = c(F,F), gates = c(gt[i,"CD23.gate.slant"], gt[i,"CD21.gate.slant.low"]))
    temp@filter <- rotate.data(temp@filter,c(11,15),theta = pi/8)$data; temp@flow.frame <- rotate.data(temp@flow.frame,c(11,15),theta = pi/8)$data
    temp@filter <- rotate.data(temp@filter,c(11,15),theta = pi/8)$data; temp@flow.frame <- rotate.data(temp@flow.frame,c(11,15),theta = pi/8)$data
    temp <- flowDensity(temp, channels = c(11,15), position = c(NA,F), gates = c(NA, gt[i,"CD21.gate.slant.high"]))
    temp@filter <- rotate.data(temp@filter,c(11,15),theta = -pi/8)$data; temp@flow.frame <- rotate.data(temp@flow.frame,c(11,15),theta = -pi/8)$data
    preB.flowD <- flowDensity(temp, channels = c(11,15), position = c(T,T), gates = c(gt[i,"CD23.gate.low"], gt[i,"CD21.gate.low"]))
#     pr[i,33] <- preB.flowD@proportion <- 100*preB.flowD@cell.count/BcellB2.flowD@cell.count
#     ev[i,33] <- preB.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(497:512)] <- sapply(1:16, function(x){mean(f.bt@exprs[preB.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(497:512)] <- sapply(1:16, function(x){median(f.bt@exprs[preB.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(497:512)] <- sapply(1:16, function(x){mean(f@exprs[preB.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(497:512)] <- sapply(1:16, function(x){median(f@exprs[preB.flowD@index, x], na.rm = T)}) }
    
    temp <- flowDensity(temps, channels = c(11,15), position = c(NA,T), gates = c(NA, gt[i,"CD21.gate.slant.low"]))
    temp@filter <- rotate.data(temp@filter,c(11,15),theta = (pi/8))$data; temp@flow.frame <- rotate.data(temp@flow.frame,c(11,15),theta = (pi/8))$data
    temp@filter <- rotate.data(temp@filter,c(11,15),theta = (pi/8))$data; temp@flow.frame <- rotate.data(temp@flow.frame,c(11,15),theta = (pi/8))$data
    temp <- flowDensity(temp, channels = c(11,15), position = c(NA,F), gates = c(NA, gt[i,"CD21.gate.slant.high"]))
    temp@filter <- rotate.data(temp@filter,c(11,15),theta = -pi/8)$data; temp@flow.frame <- rotate.data(temp@flow.frame,c(11,15),theta = -pi/8)$data
    folB.flowD <- flowDensity(temp, channels = c(11,15), position = c(F,NA), gates = c(gt[i,"CD23.gate.high"], NA))
#     pr[i,34] <- folB.flowD@proportion <- 100*folB.flowD@cell.count/BcellB2.flowD@cell.count
#     ev[i,34] <- folB.flowD@cell.count
#     if (RecordMFI) { MFI.mean.OriginalData[i,c(513:528)] <- sapply(1:16, function(x){mean(f.bt@exprs[folB.flowD@index, x], na.rm = T)})
#     MFI.median.OriginalData[i,c(513:528)] <- sapply(1:16, function(x){median(f.bt@exprs[folB.flowD@index, x], na.rm = T)})
#     MFI.mean.TransformedData[i,c(513:528)] <- sapply(1:16, function(x){mean(f@exprs[folB.flowD@index, x], na.rm = T)})
#     MFI.median.TransformedData[i,c(513:528)] <- sapply(1:16, function(x){median(f@exprs[folB.flowD@index, x], na.rm = T)}) }
#     
    plotDens(BcellB2.flowD@flow.frame, channels = c(11,15), main = "B-cell B2", xlab="11. CD23", ylab="15. CD21", xlim=c(0, 4.5), ylim=c(0, 4.5), devn = F); abline(v=gt[i,"CD23.gate.lowlow"], lwd=2); abline(v=gt[i,"CD23.gate.low"], lwd=2); abline(v=gt[i,"CD23.gate.high"], lwd=2); abline(v=gt[i,"CD23.gate"], lty=2); abline(h=gt[i,"CD21.gate.low"], lwd=2); abline(h=gt[i,"CD21.gate"], lty=2) 
    lines(preB.flowD@filter, lwd=1)
    lines(folB.flowD@filter, lwd=1)
    lines(MZB.flowD@filter, lwd=1)
    
    rm(Bcell.flowD, BcellB1.flowD, BcellB2.flowD, MZB.flowD, preB.flowD, folB.flowD)
    
    
    graphics.off()
    cat("; Time is: ",TimeOutput(start2),"\n",sep="")
    
  }, error = function(err) {
    cat(paste("ERROR in Gating:  ",err));
    cat("; Time is: ",TimeOutput(start2),"\n",sep="")
    load("errorFiles.Rdata"); errorFiles <- append(errorFiles, i); save(errorFiles, file="errorFiles.Rdata")
    graphics.off()
  })
}
cat("Total time is: ",TimeOutput(start),"\n",sep="")

# ungated FCS files
ungated <- NULL
for (i in 1:col(gt)) {
  a <- which(gt[,i]<=0)
  if (length(a>0)) {
    gt[a,i] <- mean(gt[-a,i])
    ungated <- append(ungated, a)
  }
}
ungated <- unique(ungated) #REGATE!!!
save(ungated, file="ungated.Rdata")

Events_Proportions_TableUpdated <- NULL
cname <- NULL
for(i in 1:ncol(pr)){
  Events_Proportions_TableUpdated <- cbind(Events_Proportions_TableUpdated, ev[,i], pr[,i])
  cname <- c(cname, colnames(ev)[i], colnames(pr)[i])
}
Events_Proportions_TableUpdated <- Events_Proportions_TableUpdated[, -c(69,70)]
colnames(Events_Proportions_TableUpdated) <- cname

Tts_Proportions_Table <- cbind(store.allFCS, Events_Proportions_TableUpdated)
Tts_Proportions_Table <- Tts_Proportions_Table[-which(Tts_Proportions_Table[,6]==0),]

write.csv(Events_Proportions_TableUpdated, file =  paste("Events_Proportions_TableUpdated_", gsub(" ", "_", gsub("[-:]","", as.character(Sys.time()))), ".csv",sep=""))
write.csv(Tts_Proportions_Table, file =  paste("Tts_Proportions_Table_", gsub(" ", "_", gsub("[-:]","", as.character(Sys.time()))), ".csv", sep=""))
if (UpdatedVersion) {
  write.csv(pr, file =  paste("all.propsUpdated.csv",sep="")); save(pr, file="all.propsUpdated.Rdata")
  write.csv(ev, file =  paste("all.eventsUpdated.csv",sep="")); save(ev, file="all.eventsUpdated.Rdata")
  write.csv(gt, file =  paste("all.gthresUpdated_20160721.csv",sep="")); save(gt, file="all.gthresUpdated_20160721.Rdata")
} else {
  write.csv(pr, file =  paste("all.props.csv",sep="")); save(pr, file="all.props.Rdata")
  write.csv(ev, file =  paste("all.events.csv",sep="")); save(ev, file="all.events.Rdata")
  write.csv(gt, file =  paste("all.gthres.csv",sep="")); save(gt, file="all.gthres.Rdata")
}
if (RecordMFI) {
  write.csv(cbind(store.allFCS, MFI.mean.OriginalData), file="MFI.mean.OriginalData.csv", row.names=F); save(MFI.mean.OriginalData, file="MFI.mean.OriginalData.Rdata")
  write.csv(cbind(store.allFCS, MFI.median.OriginalData), file="MFI.median.OriginalData.csv", row.names=F); save(MFI.median.OriginalData, file="MFI.median.OriginalData.Rdata")
  write.csv(cbind(store.allFCS, MFI.mean.TransformedData), file="MFI.mean.TransformedData.csv", row.names=F); save(MFI.mean.TransformedData, file="MFI.mean.TransformedData.Rdata")
  write.csv(cbind(store.allFCS, MFI.median.TransformedData), file="MFI.median.TransformedData.csv", row.names=F); save(MFI.median.TransformedData, file="MFI.median.TransformedData.Rdata")
}
