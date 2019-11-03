## input: gates + fcm files,  meta data paths
## output: feat_file_cell_count, meta_file
## process: 
## - takes gates + fcm files (values randomly generated via normal distribution; see different pos_ for how positive controls are changed
## - compiles flowtype vectors together to create cell count matrix
## - reformats meta file (meta info for fcm files)
## - creates cell meta file (meta info for fcm cell populations)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## input directories
input_dir = "/mnt/f/Brinkman group/current/Alice/gating_projects/pregnancy" # make controls based off of pregnancy data set

## ouput directories
# see for loop below



## libraries
source("source/_func.R")
libr(c("flowCore", "flowType", "flowDensity", "flowViz",
       "CytoML", "flowWorkspace",
       "pracma", "tools", "MASS", "KernSmooth", #"fitdistrplus",
       "foreach", "doMC", "plyr", "stringr")) # too much, but will use later on

## cores
no_cores = detectCores()-1
registerDoMC(no_cores)


## options
writecsv = F
nsample = 1000 # # of samples to create per data set; don't change this, i hard coded in cytodx and save fcm below on which files to use
nctrl = .5 # % of control samples
markern = 4 # # of markers
maxmarker = 6

normean = 300000 # 301234.7 for pregnancy
normsd = 0 # 64734.05 for pregnancy; if >0, remember to save as count not countAdj & run 01_feat_normalizeadj

lastlsdp = .1 # when generating only last layer cell proportions, use lastlsdp*normean/(2^markern) as sd


start = Sys.time()


# define number of cells in each fcs file
# for pregnancy data set, mean=301234.7, sd=64734.05
# ncells = floor(rnorm(nsample,300000,65000)) # number of cells for each sample
ncells = rnorm(nsample,normean,normsd) # number of cells for each sample

## meta/file
meta_file = data.frame(id=paste0("a",1:nsample), class="exp")
meta_file$class[1:(nctrl*nsample)] = "control"
meta_file$train = rep(F,nrow(meta_file))
meta_file$train[(nctrl*nsample+1):(nctrl*nsample+(1-nctrl)*nsample/2)] = meta_file$train[1:(nctrl*nsample/2)] = T



## feat/file-count-count: flowtype

# load sample fcs file
f = read.FCS(paste0(input_dir,"/samplefcs.fcs"))

# markers
gthresm = markers = LETTERS[1:markern]

# marker indices in f@exprs
ci = c(1:markern); names(ci) = markers 

# marker thresholds
cvd = rnorm(ncells[1],2,1)
p50 = quantile(cvd, .5)
p75 = quantile(cvd, .75)
p25 = quantile(cvd, .25)
thress0 = llply(markers, function(x) p50); names(thress0) = markers
thress1 = thress2 = thress4 = thress0
thress1[[gthresm[1]]] = thress2[gthresm[1:2]] = p25
thress4[gthresm[1:4]] = quantile(cvd, .501)


# start = Sys.time()
# for (ds in c(paste0("pos",1:26),paste0("ctrl",0:9))) {
for (ds in c(paste0("pos",1:32))) {
  # for (ds in c("pos5")) {
  # clear/load memory
  
  # start
  start2 = Sys.time()
  
  
  # ouput directories
  result_dir = paste0(root, "/result/",ds)
  fcs_dir = paste0(result_dir,"/fcs"); dir.create(fcs_dir, showWarnings=F, recursive=T)
  meta_dir = paste(result_dir, "/meta", sep=""); dir.create(meta_dir, showWarnings=F, recursive=T)
  feat_dir = paste(result_dir, "/feat", sep=""); dir.create(feat_dir, showWarnings=F)
  feat_file_cell_count_dir = paste(feat_dir,"/file-cell-countAdj",sep="")
  
  # make cell names
  f@exprs = matrix(rnorm(ncells[1]*length(markers),2,1),nrow=ncells[1])
  a = flowType(Frame=f, PropMarkers=ci, MarkerNames=markers, 
               MaxMarkersPerPop=min(markern,maxmarker), PartitionsPerMarker=2, 
               Thresholds=thress0, 
               Methods='Thresholds', verbose=F, MemLimit=60)
  ftcell = unlist(lapply(a@PhenoCodes, function(x){return( decodePhenotype(x, markers, a@PartitionsPerMarker) )}))
  ftcell_ = str_count(ftcell,"[+|-]")
  lastlcp = ftcell[ftcell_==length(markers)]
  lastlcpm = llply(str_extract_all(lastlcp,"[A-Z][+|-]"), function(x) 
    grepl("[+]",x) )
  lastlallposi = which(sapply(lastlcpm, function(x) all(x)))
  lastlallnegi = which(sapply(lastlcpm, function(x) all(!x)))
  
  
  # flowtype
  ftl = llply(loopInd(1:nsample,no_cores), function(ii) {
    llply(ii, function(i) {
      # v1 randomized matrix
      f@exprs = matrix(rnorm(ncells[i]*length(markers),2,1), nrow=ncells[i])
      # v2 randomized matrix: randomly define last layer pops
      lastlm = normean/length(lastlcp)
      lastlsd = .1*lastlm
      lastln = round(rnorm(length(lastlcp),lastlm,lastlsd))
      lastlncs = cumsum(lastln)
      
      fex = matrix(p25, nrow=sum(lastln), ncol=length(markers))
      fex[1:lastlncs[1], lastlcpm[[1]]] = p75
      for (j in 2:length(lastlcp)) 
        fex[(lastlncs[j-1]+1):lastlncs[j], lastlcpm[[j]]] = p75
      
      if (!grepl("ctrl",ds))
        if (as.numeric(gsub("pos","",ds))>6 ) f@exprs = fex
      
      
      
      thress = thress0
      if (i>(nsample*nctrl) & grepl("pos",ds)) {
        # make base graph for plot
        if (i == nsample*nctrl+1 & grepl("pos",ds)) {
          ft = flowType(Frame=f, PropMarkers=ci, MarkerNames=markers,
                        MaxMarkersPerPop=min(markern,maxmarker), PartitionsPerMarker=2,
                        Thresholds=thress,
                        Methods='Thresholds', verbose=F, MemLimit=60)
          ftcell = unlist(lapply(ft@PhenoCodes, function(x)
            decodePhenotype(x, ft@MarkerNames, rep(2,length(ft@MarkerNames))) ))
          ftv0 = ftv0_ = ft@CellFreqs
          ftv0 = round(ftv0/ftv0[1],3)
        }
        
        ap = f@exprs[,1]>thress[[1]]
        bp = f@exprs[,2]>thress[[2]]
        cp = f@exprs[,3]>thress[[3]]
        dp = f@exprs[,4]>thress[[4]]
        # ep = f@exprs[,5]>thress[[5]]
        double = ap & bp
        triple = ap & bp & cp
        quad   = ap & bp & cp & dp
        # quint = ap & bp & cp & dp & ep
        
        # change f values
        if (ds=="pos1") { # A+ > .75; A- > .25
          thress = thress1 
        } 
        else if (ds=="pos2") { # A+, B+ > .75; A-, B- > .25
          thress = thress2 
        } 
        else if (ds=="pos3") { # A-B+ > A+B+ x1.5
          tm = sum(double)/2
          f@exprs[sample(which(bp & !ap),tm),1] = p75 # 
          f@exprs[sample(which(ap & !bp),tm),1] = p25 # 
          f@exprs[sample(which(!ap & !bp),tm),1] = p25 # 
        }
        else if (ds=="pos4") { # A-B+ > A+B+ x1.5; D+c- > D+c+ x1.5
          tm = round(sum(double)/2)
          f@exprs[sample(which(bp & !ap),tm),1] = p75 # 
          f@exprs[sample(which(ap & !bp),tm),1] = p25 # 
          f@exprs[sample(which(!ap & !bp),tm),1] = p25 # 
          
          f@exprs[sample(which(dp & !cp),tm),3] = p75 # 
          f@exprs[sample(which(cp & !dp),tm),3] = p25 # 
          f@exprs[sample(which(!dp & !cp),tm),3] = p25 # 
        }
        else if (ds=="pos5") { # A-B+C+ > A+B+C+ x1.5
          tm = sum(triple)/2
          f@exprs[sample(which(bp & cp & !ap),tm),1] = p75 # bc
          f@exprs[sample(which(ap & cp & !bp),tm),1] = p25 # ac
          f@exprs[sample(which(ap & bp & !cp),tm),1] = p25 # ab
          f@exprs[sample(which(!ap & !bp & !cp),tm),1] = p75 # a
          
          # tm = sum(triple)/2/3
          # f@exprs[which(bp & cp & !ap)[1:tm],1] = p75 # bc
          # f@exprs[which(ap & cp & !bp)[1:tm],1] = p25 # ac
          # f@exprs[which(ap & bp & !cp)[1:tm],1] = p25 # ab
          # f@exprs[which(!ap & !bp & !cp)[1:tm],1] = p75 # a
          # 
          # f@exprs[which(bp & cp & !ap)[(tm+1):(tm*2)],2] = p25 
          # f@exprs[which(ap & cp & !bp)[(tm+1):(tm*2)],2] = p75 
          # f@exprs[which(ap & bp & !cp)[(tm+1):(tm*2)],2] = p25 
          # f@exprs[which(!ap & !bp & !cp)[(tm+1):(tm*2)],2] = p75 
          # 
          # f@exprs[which(bp & cp & !ap)[(tm*2+1):(tm*3)],3] = p25 
          # f@exprs[which(ap & cp & !bp)[(tm*2+1):(tm*3)],3] = p25 
          # f@exprs[which(ap & bp & !cp)[(tm*2+1):(tm*3)],3] = p75 
          # f@exprs[which(!ap & !bp & !cp)[(tm*2+1):(tm*3)],3] = p75 
        }
        else if (ds=="pos6") { # A-B+C+ > A+B+C+ x1.5; b+D+c- > b+D+c+ x1.5
          tm = sum(triple)/2
          f@exprs[sample(which(bp & cp & !ap),tm),1] = p75 # bc
          f@exprs[sample(which(ap & cp & !bp),tm),1] = p25 # ac
          f@exprs[sample(which(ap & bp & !cp),tm),1] = p25 # ab
          f@exprs[sample(which(!ap & !bp & !cp),tm),1] = p75 # a
          
          f@exprs[sample(which(bp & cp & !dp),tm),3] = p75 # bc
          f@exprs[sample(which(bp & dp & !cp),tm),3] = p25 # ac
          f@exprs[sample(which(dp & cp & !bp),tm),3] = p25 # ab
          f@exprs[sample(which(!bp & !cp & !dp),tm),3] = p75 # a
          
          # tn = sum(quad)/2
          # f@exprs[sample(which(!ap & bp & cp & dp),tn),1] = p75 # bcd
          # f@exprs[sample(which(ap & !bp & cp & dp),tn),1] = p25 # acd
          # f@exprs[sample(which(ap & bp & !cp & dp),tn),1] = p25 # abd
          # f@exprs[sample(which(ap & bp & cp & !dp),tn),1] = p25 # abc
          # f@exprs[sample(which(!ap & bp & !cp & !dp),tn),1] = p75 # ab
          # f@exprs[sample(which(!ap & !bp & cp & !dp),tn),1] = p75 # ac
          # f@exprs[sample(which(!ap & !bp & !cp & dp),tn),1] = p75 # ad
          # f@exprs[sample(which(ap & !bp & !cp & !dp),tn),1] = p25 # a
        }
        else if (ds=="pos7") { # A+B+C+D+ > x2
          f@exprs = rbind(f@exprs, f@exprs[ap & bp & cp & dp,])
          
          # tm = sum(triple)/2/3
          # f@exprs[sample(which(bp & cp & !ap),tm),1] = p75 # bc
          # f@exprs[sample(which(!bp & cp & ap),tm),2] = p75 # bc
          # f@exprs[sample(which(bp & !cp & ap),tm),3] = p75 # bc
        } 
        else if (ds=="pos8") { # A-B-C-D-E- > x2
          f@exprs = rbind(f@exprs, f@exprs[!ap & !bp & !cp & !dp,])
          
          # tn = sum(quint)/2
          # f@exprs[sample(which(!ap & bp & cp & dp & ep),tn),1] = p75
          # f@exprs[sample(which(ap & !bp & cp & dp & ep),tn),1] = p25
          # f@exprs[sample(which(ap & bp & !cp & dp & ep),tn),1] = p25
          # f@exprs[sample(which(ap & bp & cp & !dp & ep),tn),1] = p25
          # f@exprs[sample(which(ap & bp & cp & dp & !ep),tn),1] = p25
          # 
          # f@exprs[sample(which(!ap & bp & cp & !dp & !ep),tn),1] = p75
          # f@exprs[sample(which(!ap & bp & !cp & dp & !ep),tn),1] = p75
          # f@exprs[sample(which(!ap & bp & !cp & !dp & ep),tn),1] = p75
          # f@exprs[sample(which(!ap & !bp & cp & dp & !ep),tn),1] = p75
          # f@exprs[sample(which(!ap & !bp & cp & !dp & ep),tn),1] = p75
          # f@exprs[sample(which(!ap & !bp & !cp & dp & ep),tn),1] = p75
          # 
          # f@exprs[sample(which(ap & bp & !cp & !dp & !ep),tn),1] = p25
          # f@exprs[sample(which(ap & !bp & cp & !dp & !ep),tn),1] = p25
          # f@exprs[sample(which(ap & !bp & !cp & dp & !ep),tn),1] = p25
          # f@exprs[sample(which(ap & !bp & !cp & !dp & ep),tn),1] = p25
          # 
          # f@exprs[sample(which(!ap & !bp & !cp & !dp & !ep),tn),1] = p75
        } 
        else if (ds=="pos9") { # same as above but both
          f@exprs = rbind(f@exprs, f@exprs[ap & bp & cp & dp,])
          f@exprs = rbind(f@exprs, f@exprs[!ap & !bp & !cp & !dp,])
        } 
        else if (ds=="pos10") { # A+, B+, c+
          f@exprs = Reduce(rbind,list(f@exprs, f@exprs[sample(which(ap),sum(double)),], f@exprs[sample(which(bp),sum(double)),], f@exprs[sample(which(cp),sum(double)),]))
        } 
        else if (ds=="pos11") { # 1.5x A+
          tripleind = which(fex[,1]==p75)
          tm = sum(double)
          f@exprs = rbind(fex,fex[sample(tripleind,tm),])
        } 
        else if (ds=="pos12") { # 1.5x A+, B+
          tripleind = which(fex[,1]==p75)
          tm = sum(double)
          f@exprs = rbind(fex,fex[sample(tripleind,tm),])
          tripleind = which(fex[,2]==p75)
          f@exprs = rbind(f@exprs,fex[sample(tripleind,tm),])
        } 
        else if (ds=="pos13") { #1.5x A+B+
          tripleind = which(fex[,1]==p75 & fex[,2]==p75)
          tm = sum(double)/2
          f@exprs = rbind(fex,fex[sample(tripleind,tm),])
        } 
        else if (ds=="pos14") { #1.5x A+B+, D+c+
          tm = sum(double)/2
          tripleind = which(fex[,1]==p75 & fex[,2]==p75)
          f@exprs = rbind(fex,fex[sample(tripleind,tm),])
          tripleind = which(fex[,4]==p75 & fex[,3]==p75)
          f@exprs = rbind(f@exprs,fex[sample(tripleind,tm),])
        } 
        else if (ds=="pos15") { #1.5x A+B+C+
          tm = sum(triple)/2
          tripleind = which(fex[,1]==p75 & fex[,2]==p75 & fex[,3]==p75)
          f@exprs = rbind(fex,fex[sample(tripleind,tm),])
        } 
        else if (ds=="pos16") { #1.5x A+B+C+, C+D+b+
          tm = sum(triple)/2
          tripleind = which(fex[,1]==p75 & fex[,2]==p75 & fex[,3]==p75)
          f@exprs = rbind(fex,fex[sample(tripleind,tm),])
          tripleind = which(fex[,3]==p75 & fex[,4]==p75 & fex[,2]==p75)
          f@exprs = rbind(f@exprs,fex[sample(tripleind,tm),])
        } 
        else if (ds=="pos17") { # A+B+C+D+E+ > x2
          f@exprs = rbind(fex, matrix(p75,nrow=lastln[lastlallposi], ncol=length(markers)))
        }
        else if (ds=="pos18") { # A-B-C-D-E- > x2
          f@exprs = rbind(fex, matrix(p25,nrow=lastln[lastlallnegi], ncol=length(markers)))
        }
        else if (ds=="pos19") { # same as above but both
          f@exprs = rbind(fex, matrix(p75,nrow=lastln[lastlallposi], ncol=length(markers)))
          f@exprs = rbind(f@exprs, matrix(p25,nrow=lastln[lastlallnegi], ncol=length(markers)))
        } 
        else if (ds=="pos20") { # A+, B+, c+
          tripleind = which(fex[,1]==p75)
          tm = sum(double)
          f@exprs = rbind(fex,fex[sample(tripleind,tm),])
          tripleind = which(fex[,2]==p75)
          f@exprs = rbind(f@exprs,fex[sample(tripleind,tm),])
          tripleind = which(fex[,3]==p75)
          f@exprs = rbind(f@exprs,fex[sample(tripleind,tm),])
        } 
        else if (ds=="pos21") { # 1.5x A+
          tripleind = which(ap)
          tm = sum(double)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
        } 
        else if (ds=="pos22") { # 1.5x A+, B+
          tm = sum(double)
          tripleind = which(ap)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
          tripleind = which(bp)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
        } 
        else if (ds=="pos23") { #1.5x A+B+
          tm = sum(double)/2
          tripleind = which(ap & bp)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
        } 
        else if (ds=="pos24") { #1.5x A+B+, C+D+
          tm = sum(double)/2
          tripleind = which(ap & bp)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
          tripleind = which(dp & cp)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
        } 
        else if (ds=="pos25") {
          tm = sum(triple)/2
          tripleind = which(ap & bp & cp)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
        }
        else if (ds=="pos26") { #1.5x A+B+C+, B+C+D+
          tm = sum(triple)/2
          tripleind = which(ap & bp & cp)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
          tripleind = which(cp & dp & bp)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
        } 
        else if (ds=="pos27") { #1.5x A+B+; decrease other ones accordingly
          tm = .33*nrow(f@exprs) - sum(ap & bp) #sum(double)/3
          f@exprs[c(sample(which(!ap & !bp),tm/3) , sample(which(ap & !bp),tm/3) , sample(which(!ap & bp),tm/3)),c(1,2)] = p75
        } 
        else if (ds=="pos28") { #1.5x A+B+, C+D+, like above
          tm = sum(double)/2
          f@exprs[c(sample(which(!ap & !bp),tm/3) , sample(which(ap & !bp),tm/3) , sample(which(!ap & bp),tm/3)),c(1,2)] = p75
        }
        else if (ds=="pos29") { #1.5x A+, A+B+C+
          tm = sum(triple)/2
          tripleind = which(ap & bp & cp)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
          tm = sum(ap)*1.5-sum(ap)
          tripleind = which(f@exprs[,1]>thress[[1]] & f@exprs[,2]>thress[[2]])
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
        }
        else if (ds=="pos30") { #1.5x A+, A+B+C+
          tm = sum(ap)/2
          tripleind = which(ap)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
          tm = sum(triple)/2
          tripleind = which(cp & ap & bp)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
        }
        else if (ds=="pos31") { #1.5x A+, B+C+D+
          tm = sum(ap)/2
          tripleind = which(ap)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
          tm = sum(triple)/2
          tripleind = which(cp & dp & bp)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
        }
        else if (ds=="pos32") {
          tm = sum(triple)/2
          tripleind = which(!ap & !bp & !cp)
          f@exprs = rbind(f@exprs,f@exprs[sample(tripleind,tm),])
        }
        
        if (i == nsample*nctrl+1 & grepl("pos",ds)) {
          # if (ds%in%c("pos1","pos2","pos9")) la=1
          # if (ds%in%c("pos3","pos5","pos6","pos7","pos8")) la=3
          # if (ds%in%c("pos4")) la=4
          
          ft = flowType(Frame=f, PropMarkers=ci, MarkerNames=markers,
                        MaxMarkersPerPop=min(markern,maxmarker), PartitionsPerMarker=2,
                        Thresholds=thress,
                        Methods='Thresholds', verbose=F, MemLimit=60)
          ftcell = unlist(lapply(ft@PhenoCodes, function(x)
            decodePhenotype(x, ft@MarkerNames, rep(2,length(ft@MarkerNames))) ))
          ftv = ftv_ = ft@CellFreqs
          ftv = round(ftv/ftv[1],3)
          names(ftv) = ftcell
          a = getPhenCP(cp=ftcell,no_cores=no_cores)
          al = layout_gr(a$gr$e,a$gr$v)
          # alp = layout_gr(a$grp$e,a$grp$v)
          al = gpdf(al)
          al$v$color = ifelse(ftv_>ftv0_,"increase","decrease")
          vind_ = !grepl("-",al$v$name)
          # alp = gpdf(alp)
          # alp$v$color = ifelse(ftv_[vind_]>ftv0_[vind_],"increase","decrease") #ftv[vind_]
          al$v$size=1
          al$v$sizeb=1
          
          al$v$label = paste0(al$v$name,":",ftv)
          al$v$v_ind = abs(ftv_-ftv0_)/ftv0_ >.05 #& !grepl("[-]",al$v$name)
          al$v$label_ind = al$v$v_ind & !grepl("[-]",al$v$name)
          al$e$e_ind = al$e[,1]%in%al$v$name[al$v$v_ind] & al$e[,2]%in%al$v$name[al$v$v_ind]
          
          gp = gggraph(al)

          ggsave(paste0(meta_dir,"/all_sig.png"), plot=gp, scale = 1, width =9, height =11, units = "in", dpi = 300, limitsize = TRUE)
          
          # alp$v$label = paste0(alp$v$name,":",ftv[vind_])
          # vind = abs(ftv[vind_]-ftv0[vind_])/ftv0[vind_] >.05
          # gp = gggraph(alp, v_ind=vind, 
          #              e_ind=alp$e[,1]%in%alp$v$name[vind] & alp$e[,2]%in%alp$v$name[vind])
          # gp = gp +
          #   geom_label_repel(
          #     # data=alp$v[str_count(alp$v$name,"[-+]")==la & vind,],
          #     data=alp$v[vind & !grepl("[-]",al$v$name)[vind_],],
          #     aes(x=x,y=y,label=label, color=color),
          #     nudge_x = -.1, direction = "y", hjust = 1, segment.size = 0.2)
          # 
          # ggsave(paste0(meta_dir,"/all_sigpos.png"), plot=gp, scale = 1, width =9, height =11, units = "in", dpi = 300, limitsize = TRUE)
        }
      }
      fe = f@exprs
      colnames(fe) = LETTERS[1:ncol(fe)]
      if (i%in%c(1:5,251:255,501:505,751:755) & grepl("pos",ds)) save(fe,file=paste0(fcs_dir,"/a",i,".Rdata"))
      flowType(Frame=f, PropMarkers=ci, MarkerNames=markers, 
               MaxMarkersPerPop=min(markern,maxmarker), PartitionsPerMarker=2, 
               Thresholds=thress, 
               Methods='Thresholds', verbose=F, MemLimit=60)@CellFreqs
    })
  }, .parallel=T)
  ft = Reduce(rbind,unlist(ftl,recursive=F))
  colnames(ft) = ftcell
  rownames(ft) = meta_file$id
  
  save(meta_file, file=paste0(meta_dir,"/file.Rdata"))
  if (writecsv) write.csv(meta_file, file=paste0(meta_dir,"/file.csv"), row.names=F)
  
  # compile ---------------------------
  m0 = as.matrix(ft)
  save(m0, file=paste0(feat_file_cell_count_dir,".Rdata"))
  if (writecsv) write.csv(m0, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)
  
  time_output(start2, ds)
  rm(list=c("ftl","ft","m0")); gc()
  
}
time_output(start)

