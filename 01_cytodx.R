# i'm going to apply cytoDX on all the data sets! for comparison with lnpropexpect and flowtype cell counts
# i make the fcs for controls and positive controls in the data script... so i'll do it there!

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr(c("flowCore", "flowType", "flowDensity", "flowViz",
       "CytoML", "flowWorkspace", "CytoDx",
       "foreach", "doMC", "plyr", "stringr")) # too much, but will use later on

## cores
no_cores = detectCores()-1
registerDoMC(no_cores)



start = Sys.time()


fcss = mfs = NULL
l_ply(c("bodenmiller","flowcap6","flowcap6_ctrl","flowcap6_pos","genentech",paste0("ctrl",0:9),paste0("pos",1:5),"pregnancy"), function(data) { try({
  start1 = Sys.time()
  print(data)
  
  ## load processed fcs
  if (grepl("flowcap6",data)) {
    # flowcap-ii
    if (!exists("fcfcss")) {
      fcs_dirs = list.files("/mnt/f/Brinkman group/current/Alice/flowCAP-II/data/FCS", full.names=T, pattern=".FCS")
      sm = get(load(paste0(root, "/result/","flowcap6","/meta/file.Rdata")))
      fcsid = as.numeric(gsub(".Rdata","",fileNames(fcs_dirs)))
      fcs_dirs = fcs_dirs[match(sm$id, fcsid)] # only get tube 6
      fcfcss = llply(fcs_dirs, read.FCS)
      names(fcfcss) = sm$id
    } 
    fcss = fcfcss
    sm0 = get(load(paste0(root, "/result/",data,"/meta/file.Rdata")))
    markers = str_extract(fcss[[1]]@parameters@data[,2],"[A-Za-z0-9]+")
    if (grepl("flowcap6_pos",data)) {
      new_ids = sm0$id[sm0$class=="mix"]
      # match samples for mixing
      normali = as.character(sm0$id[sm0$class=="control"])
      amli = as.character(sm0$id[sm0$class=="aml"])
      for (i in new_ids) {
        normind = sample(normali, 5)
        amlind = sample(amli, 5)
        f = fcss[[sample(1:length(fcss),1)]]
        inn = nrow(f@exprs)/10
        f@exprs = Reduce(rbind,llply(append(normind, amlind), function(x) 
          fcss[[x]]@exprs[sample(1:nrow(fcss[[x]]@exprs),inn),] ))
        fcss[[as.character(i)]] = f
      }
    } else if (grepl("flowcap6_ctrl",data)) {
      fcss = fcss[as.character(sm0$id)]
    }
    
  } else if (grepl("pregnancy",data)) {
    # pregnancy
    gs = load_gs("/mnt/f/Brinkman group/current/Alice/gating_projects/pregnancy/gs")
    sm0 = get(load("/mnt/f/Brinkman group/current/Alice/flowtype_metric/result/pregnancy/meta/file.Rdata"))
    fcss = as(gs_pop_get_data(gs),Class="list")
    fcss = llply(fcss, function(f) {
      f@exprs = f@exprs[,c(23,50,40,45,20,11,18,51,14,16,21,35,27)]
      f
    })
    names(fcss) = gsub("Repeat_","",names(fcss))
    names(fcss) = gsub("BL","4",names(fcss))
    names(fcss) = str_extract(names(fcss),"PTLG[0-9]+_[0-9]")
    fcss = fcss[sm0$id]
    markers = c(# "CD11b", "CD11c", 
      "CD123", "CD14", "CD16", 
      # "CD19", 
      # "CD25",
      "CD3", "CD4","CD45", "CD45RA", "CD56", 
      "CD66", 
      "CD7", "CD8", 
      # "FoxP3", 
      # "HLADR", 
      "Tbet", "TCRgd")
    
  } else if (data=="bodenmiller") {
    fcs_dirs = list.files("/mnt/f/Brinkman group/current/Alice/gating_projects/HDCytoData_Bodenmiller/fcs", full.names=T, pattern=".fcs")
    sm0 = get(load(paste0(root, "/result/bodenmiller/meta/file.Rdata")))
    fcss = get(load("/mnt/f/Brinkman group/current/Alice/gating_projects/HDCytoData_Bodenmiller/fcs.Rdata"))
    fcss = llply(fcss, function(f) {
      f@exprs = f@exprs[,c(3,9, 11, 12, 21, 29, 33)]
      f
    })
    fcss = fcss[sm0$id]
    markers = c("CD3","CD4","CD20","CD33","CD14","IgM","CD7") # HLA-DR
    
  } else if (data=="genentech") {
    gs = load_gs("/mnt/f/Brinkman group/current/Alice/gating_projects/genetch/Tube_003gs")
    fcss = as(gs_pop_get_data(gs, "Myeloid"),Class="list")
    names(fcss) = gsub("%|.fcs","",names(fcss))
    sm0 = get(load(paste0(root, "/result/genentech/meta/file.Rdata")))
    fcss = fcss[sm0$id]
    markers = get(load("/mnt/f/Brinkman group/current/Alice/gating_projects/genetch/flowType/Tube_003/MarkerNames_myeloid.Rdata"))
    fcss = llply(fcss, function(f) {
      mi = match(markers,f@parameters@data$desc)
      f@exprs = f@exprs[,mi]
      f
    })
    
  } else { # ctrl/pos
    sm0 = get(load(paste0(root, "/result/",data,"/meta/file.Rdata")))
    fcss = llply(sm0$id, function(x) get(load(paste0(root, "/result/",data,"/fcs/",x,".Rdata"))),.parallel=T)
    names(fcss) = sm0$id
  }
  
  
  for (uc in unique(sm0$class)) {
    if (uc=="control") next()
    sm = sm0[sm0$class%in%c(uc,"control"),]
    fcs = fcss[as.character(sm$id)]
    
    c_dir = paste0(root,"/result/",data,"/cytodx/",uc)
    dir.create(c_dir,showWarnings=F,recursive=T)
    
    ## collate fcs data
    train_data = Reduce(rbind,llply(as.character(sm$id[sm$train]), function(x) {
      if (typeof(fcs[[x]])=="S4") {
        a = fcs[[x]]@exprs
      } else {
        a = fcs[[x]]
      }
      data.frame(a,x,sm$class[as.character(sm$id)==x])
    }))
    test_data = Reduce(rbind,llply(as.character(sm$id[!sm$train]), function(x) {
      if (typeof(fcs[[x]])=="S4") {
        a = fcs[[x]]@exprs
      } else {
        a = fcs[[x]]
      }
      as = rep(x,nrow(a))
      ay = rep(sm$class[as.character(sm$id)==x],nrow(a))
      data.frame(a,as,ay)
    }))
    colnames(train_data) = colnames(test_data) = c(markers,"xSample","y")
    if (sum(is.na(train_data))>0) {
      b = train_data[,1:(ncol(train_data)-2)]
      b = b-min(b[!is.na(b)])
      b[is.na(b)] = 0
      train_data[,1:(ncol(train_data)-2)] = b
    }
    if (sum(is.na(test_data))>0) {
      b = test_data[,1:(ncol(test_data)-2)]
      b = b-min(b[!is.na(b)])
      b[is.na(b)] = 0
      test_data[,1:(ncol(test_data)-2)] = b
    }
    
    ## build models
    start1 = Sys.time()
    x = as.matrix(train_data[,markers])
    x = x-min(x[!is.na(x)])
    x[is.na(x)] = 0
    fit1 = CytoDx.fit(x=x,
                      y=train_data$y,
                      xSample = train_data$xSample,
                      reg=F, family="binomial", parallelCore=no_cores)
    # Use decision tree to find the cell subsets that are associated; uses anova to compare means
    png(paste0(c_dir,"/tree1.png"), width=1000, height=800)
    tg1 = treeGate(P=fit1$train.Data.cell$y.Pred.s0, x=train_data[,markers])
    graphics.off()
    save(fit1,file=paste0(c_dir,"/fit1.Rdata"))
    save(tg1,file=paste0(c_dir,"/tree1.Rdata"))
    
    start1 = Sys.time()
    fit1 = CytoDx.fit(x=as.matrix(test_data[,markers]),
                      y=test_data$y,
                      xSample = test_data$xSample,
                      reg=F, family="binomial", parallelCore=no_cores)
    # Use decision tree to find the cell subsets that are associated; uses anova to compare means
    png(paste0(c_dir,"/tree2.png"), width=1000, height=800)
    tg1 = treeGate(P=fit1$train.Data.cell$y.Pred.s0, x=test_data[,markers])
    graphics.off()
    save(fit1,file=paste0(c_dir,"/fit2.Rdata"))
    save(tg1,file=paste0(c_dir,"/tree2.Rdata"))
  }
  
  rm(list=c("tg1","fit1","test_data","train_data")); gc()
  
  time_output(start1)
})},.parallel=T)
time_output(start)

# # tg1t = tg1$splits[tg1$splits[,1]>0,]
# tg1t = tg1$frame
# yval = NULL
# yval[""] = .5
# cp = ""
# # yval: fitted value of the response at the node
# for (i in 1:nrow(tg1t)) {
#   if (tg1t[i,1]=="<leaf>") {
#     lasti = substr(cp,nchar(cp),nchar(cp))
#     if (!grepl("[+|-]",lasti)) {
#       cp1 = paste0(cp,"-")
#       yval[cp1] = tg1t[i,1]
#     }
#   }
# }




