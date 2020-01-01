## input: count adjusted matrix (flowcap data set only) (must run 00_data_flowcap.R before running this, every time)
## output: count adjusted matrix + frankenstein samples
## process: 
## - add a 3rd category of samples as a linear combination (weight = 1/(2*matchsamples)) of randomly sampled (5) control & (5) aml samples -- sampling done once, combinations same for all flowcap panels
## - one frankenstein sample for every control sample is made

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)


## libraries
source("source/_func.R")
libr(c("stringr",
       "foreach","doMC","plyr"))


## cores
no_cores = detectCores()-1
registerDoMC(no_cores)


## options
writecsv = F
options(stringsAsFactors=FALSE)
options(na.rm=T)

matchsamples = 5 #pos: number of samples from each normal and aml to mix
amlprop = 1/2 #ctrl: proportion of control sample to make as "aml"



start = Sys.time()

randomindp = randomindc = NULL
result_dirs =list.dirs(paste0(root,"/result"),full.names=T,recursive=F)
for (result_dir in result_dirs) {
  if (!grepl("flowcap",result_dir) | grepl("pos|ctrl",result_dir)) next
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta")
  meta_file_dir = paste(meta_dir, "/file", sep="")
  feat_dir = paste(result_dir, "/feat", sep="")
  feat_file_cell_countAdj_dir = paste(feat_dir, "/file-cell-countAdj", sep="")
  
  
  for (ds in c("pos","ctrl")) {
    
    ## output directories
    result_dir_ = paste0(result_dir,"_",ds)
    meta_dir_ = paste0(result_dir_,"/meta"); dir.create(meta_dir_, recursive=T, showWarnings=F)
    meta_file_dir_ = paste(meta_dir_, "/file", sep="")
    feat_dir_ = paste0(result_dir_,"/feat"); dir.create(feat_dir_, recursive=T, showWarnings=F)
    feat_file_cell_countAdj_dir_ = paste(feat_dir_, "/file-cell-countAdj", sep="")
    feat_file_cell_prop_dir_ = paste0(feat_dir_, "/file-cell-prop")
    
    
    ## prepare data
    file.copy(file.path(meta_dir,list.files(meta_dir)), meta_dir_, overwrite=T)
    meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
    feat_file_cell_countAdj0 = get(load(paste0(feat_file_cell_countAdj_dir,".Rdata")))
    
    if (ds=="pos") {
      ## match samples for mixing
      normali = which(meta_file0$class=="control")
      amli = which(meta_file0$class=="aml")
      if (is.null(randomindp)) {
        for (i in 1:length(amli)) {
          randomindp[[i]] = list()
          randomindp[[i]]$normal = sample(normali, matchsamples)
          randomindp[[i]]$aml = sample(amli, matchsamples)
        }
      }
      weight = 1/(2*matchsamples)
      
      feat_file_cell_countAdj2 = as.matrix(Reduce(rbind,llply(randomindp,function(ri)
        weight*colSums(feat_file_cell_countAdj0[append(ri$normal,ri$aml),]) )))
      rownames(feat_file_cell_countAdj2) = c((1+nrow(feat_file_cell_countAdj0)):(nrow(feat_file_cell_countAdj2)+nrow(feat_file_cell_countAdj0)))
      m0 = rbind(feat_file_cell_countAdj0,feat_file_cell_countAdj2)
      p0 = m0/m0[,colnames(m0)==""]
      
      
      ## meta/file
      meta_file = rbind(meta_file0, data.frame(
        class=rep("mix", length(randomindp)), 
        id=max(meta_file0$id)+c(1:length(randomindp)), 
        train=append(rep(F,floor(length(randomindp)/2)), rep(T,ceiling(length(randomindp)/2))),
        subject=rep(0, length(randomindp)),
      tube=rep(meta_file0$tube[1], length(randomindp))
        ))
      
    } else if (ds=="ctrl") {
      ## extract controls
      m0 = feat_file_cell_countAdj0[meta_file0$class=="control",]
      p0 = m0/m0[,colnames(m0)==""]
      
      ## meta/file
      meta_file = meta_file0[meta_file0$class=="control",]
      if (is.null(randomindc)) 
        randomindc = sample(1:nrow(meta_file), floor(amlprop*nrow(meta_file)))
      meta_file$class[randomindc] = "exp"
    }
    
    ## save
    save(p0, file=paste0(feat_file_cell_prop_dir_,".Rdata"))
    save(m0, file=paste0(feat_file_cell_countAdj_dir_,".Rdata"))
    save(meta_file, file=paste0(meta_file_dir_,".Rdata"))
  }
}

time_output(start)
