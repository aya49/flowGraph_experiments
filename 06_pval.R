## input: features 
## output: p values and their corelation between train/test

set.seed(10)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr(c("stringr","colorspace", "Matrix", "plyr",
       "vegan", # libr(proxy)
       "foreach","doMC",
       "kernlab"))



#Setup Cores
no_cores = 6 # detectCores()-5
registerDoMC(no_cores)



## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

# cvn-fold cross validation

overwrite = T #overwrite distances?
writecsv = F

readcsv = F

pthres = .1 # p value sig threshold for t test
good_count = 5
# minfold = 5 # minimum number of samples in each fold/class
good_sample = minfold = 5 # each class must have more thatn good_sample amount of samples or else prune
# countThres = 1000 #insignificant if count under

id_col = "id"
target_col = "class"
order_cols = NULL

control = "control"


for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)[-16]) {
  # result_dir = paste0(root, "/result/impc_panel1_sanger-spleen") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta") # meta files directory
  meta_file_dir = paste(meta_dir, "/file", sep="") #meta for rows (sample)
  meta_cell_dir = paste(meta_dir, "/cell", sep="") #meta for rows (sample)
  feat_dir = paste(result_dir, "/feat", sep="") #feature files directory
  
  
  ## output directories
  pval_dir = paste0(result_dir,"/pval"); suppressWarnings(dir.create (pval_dir))
  
  
  
  #data paths
  feat_count = "feat_cell_countAdj"
  feat_types = list.files(path=feat_dir,pattern=".Rdata")
  # feat_types = feat_types[!feat_types=="file-cell-count.Rdata"]
  feat_types = gsub(".Rdata","",feat_types)
  feat_types = feat_types[!grepl("KO|Max",feat_types)]
  # feat_types = feat_types[!grepl("Freqp_orig",feat_types)]
  
  meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
  

  start = Sys.time()
  
  mc0 = Matrix(get(load(paste0(feat_dir,"/", feat_count,".Rdata"))))
  
  #load different features matrix and calculate distance matrix
  # for (feat_type in feat_types_) {
  foldsip0 = llply(feat_types, function(feat_type) {
    tryCatch({
    cat("\n", feat_type, " ",sep="")
    start2 = Sys.time()
    
    ## upload and prep feature matrix
    m0 = Matrix(get(load(paste0(feat_dir,"/", feat_type,".Rdata"))))
    iscell = grepl("_",colnames(m0))[10]
    
    sm = meta_file[match(rownames(m0),meta_file$id),]
    # good_sample
    tcl = table(sm$class); 
    delrow = rownames(m0)%in%sm$id[sm$class%in%names(tcl)[tcl<=minfold]]
    sm = sm[!delrow,]
    m = m0[!delrow,]
    if (all(sm$class==sm$class[1])) return(F)
    
    # n-fold cv indices
    controli = sm$class=="control"
    foldsip = list()
    if ("type"%in%colnames(sm)) {
      for (uc in unique(sm$class[!controli])) {
        uci = sm$class==uc
        foldsip[[uc]][[1]]$test = sm$id[sm$type=="test" & uci]
        foldsip[[uc]][[1]]$train = sm$id[sm$type=="train" & uci]
      }
    } else if (min(table(sm$class))<(minfold*3)) { # only two folds? then still train/test
      for (uc in unique(sm$class[!controli])) {
        uci = which(sm$class==uc)
        te_ind = sample(uci,floor(nrow(sm)/2))
        foldsip[[uc]][[1]]$test = sm$id[sample(uci,floor(nrow(sm)/2))]
        foldsip[[uc]][[1]]$train = sm$id[uci[!uci%in%te_ind]]
      }
    } else {
      for (uc in unique(sm$class[!controli])) {
        uci = which(sm$class==uc)
        cvn = min(10,length(uci)/minfold)
        foldsip_test = split(uci, cut(seq_along(uci), cvn, labels=F))
        foldsip[[uc]] = llply(foldsip_test, function(x) 
          list(test=sm$id[x], train=sm$id[uci[!uci%in%x]]) )
      }
    }
    
    # calculate p value and plot; plot needs clean up, right now plots same number of plots per row, will want to fix so that plot only up until a nmber of folds per row so it's a unique class per row
    for (uc in names(foldsip)) {
      for (fold in 1:length(foldsip[[uc]])) {
        foldsip[[uc]][[fold]]$controln = sum(controli)
        
        tri = foldsip[[uc]][[fold]]$train
        tei = foldsip[[uc]][[fold]]$test
        
        pv_tr = sapply(1:ncol(m), function(j) 
          t.test(m[controli,j], m[tri,j])$p.value)
        pv_te = sapply(1:ncol(m), function(j) 
          t.test(m[controli,j], m[tei,j])$p.value)
        
        names(pv_tr) = names(pv_te) = colnames(m)
        foldsip[[uc]][[fold]]$p_train = pv_tr
        foldsip[[uc]][[fold]]$p_test = pv_te
        
        pv_trl = log(pv_tr)
        pv_tel = log(pv_te)
        
        foldsip[[uc]][[fold]]$pcorr = cor(pv_trl,pv_tel, method="pearson")
        foldsip[[uc]][[fold]]$pcorrp = cor.test(pv_trl,pv_tel, method="pearson")$p.value
      }
    }
    
    return(foldsip)
    time_output(start2, feat_type)
    }, error = function(err) { cat(paste("ERROR:  ",err)); return(T) })
  }, .parallel=T)
  save(foldsip0, file=paste0(pval_dir,"/pval.Rdata"))
  
  
  for (feat_type in names(foldsip0)) {
    nh = length(foldsip0[[feat_type]])
    nw = max(sapply(foldsip0[[feat_type]],length))
    png(paste0(pval_dir,"/",feat_type,".png"),
        width=nw*500,height=nh*500)
    par(mfrow=c(nh,nw))
    
    for (uc in names(foldsip0[[feat_type]])) {
      for (fold in 1:length(foldsip0[[feat_type]][[uc]])) {
        fsip = foldsip0[[feat_type]][[uc]][[fold]]
        
        pv_trl = log(fsip$p_train)
        pv_tel = log(fsip$p_test)

        phens = names(pv_trl)
        if(grepl("_",phens[10])) phens = sapply(str_split(names(pv_trl),"_"), function(x) x[1])
        mcm = colMeans(mc0[c(which(rownames(mc0)%in%meta_file$id[meta_file$class=="control"]),fsip$train,fsip$test), phens])
        mcm_ = log(mcm); mcm_ = 2*mcm_/max(mcm_) # for plot dot size
        
        dcol = rep("black",length(pv_trl))
        dcol[fsip$p_train<pthres] = "blue"
        dcol[fsip$p_test<pthres] = "red"
        dcol[fsip$p_train<pthres & fsip$p_test<pthres] = "purple"
        
        plot(pv_trl, pv_tel, pch=16, cex=mcm_, col=dcol,
             main=paste0("ln(train) vs ln(test) ",feat_type," feature pvalues
             \nsize=2*ln(meancount)/max(ln(meancount))
\nsig train=blue, test=red, both=purple
            \nclass=",uc," (ntrain=",length(fsip$train)," ntest=",length(fsip$test)," nctrl=",sum(fsip$controln),")
            \n pearson corr=",round(fsip$pcorr,3)," pvalue=",round(fsip$pcorrp,3)))
      }
    }
    graphics.off()
  }
  time_output(start)
}


