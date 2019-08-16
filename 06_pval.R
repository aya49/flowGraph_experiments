## input: features 
## output: p values and their corelation between train/test sets; p value qq plots, recall/precision, number of significant features

set.seed(10)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)


## libraries
source("source/_func.R")
libr(c("stringr","Matrix", "plyr",
       "metap",
       "foreach","doMC"))


## cores
no_cores = detectCores()-1
registerDoMC(no_cores)



## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

# cvn-fold cross validation

# overwrite = F #overwrite?

good_sample = minfold = 5 # each class must have more thatn good_sample amount of samples or else prune



adjust = c("BY", "none")#,"fisher","none","BH","bonferroni") #pvalue adjustment; "lanc",
pthres = .05 # p value sig threshold for t test
good_count = 5



start = Sys.time()

result_dirs = list.dirs(paste0(root, "/result"), full.names=T, recursive=F)
for (result_dir in result_dirs) {
  print(result_dir)
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta")
  meta_file_dir = paste(meta_dir, "/file", sep="")
  feat_dir = paste(result_dir, "/feat", sep="") # feature files
  
  
  ## output directories
  sum_p_dir = paste0(result_dir,"/feat_pval"); dir.create(sum_p_dir, recursive=T, showWarnings=F)
  sum_m_dir = paste0(result_dir,"/feat_mean"); dir.create(sum_m_dir, showWarnings=F)
  
  
  
  start1 = Sys.time()
  
  
  # feature paths
  feat_types = list.files(path=feat_dir,pattern=".Rdata")
  feat_types = gsub(".Rdata","",feat_types)
  
  # load files
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
  
  feats = llply(feat_types, function(feat_type) {
    ## upload data and prep feature matrix + meta ------------
    m0 = get(load(paste0(feat_dir,"/", feat_type,".Rdata")))
    sm = meta_file[match(rownames(m0),meta_file$id),]
    
    # delete classes with too few samples
    tcl = table(sm$class)
    inrow = rownames(m0) %in% 
      sm$id[sm$class %in% names(tcl)[tcl>minfold]]
    if (all(!inrow)) return(F)
    sm = sm[inrow,]
    if (all(sm$class==sm$class[1])) return(F)
    m = m0[inrow,]
    
    ## save mean values for each class ----------------------
    mname = paste0(sum_m_dir,"/",feat_type,".Rdata")
    fmean = llply(unique(sm$class), function(uc) {
      a = apply(m[sm$class=="control",],2,mean_geo)
      names(a) = colnames(m)
      return(a)
    })
    names(fmean) = unique(sm$class)
    save(fmean, file=mname)
    return(list(sm=sm,m=m,mname=mname))
  }, .parallel=T)
  names(feats) = feat_types
  
  
  # calculate p values --------------------------------
  # for (feat_type in feat_types_) {
  for (feat_type in feat_types) {
    start2 = Sys.time()
    cat("\n", feat_type, " ",sep="")

    m = feats[[feat_type]]$m
    sm = feats[[feat_type]]$sm
    # m[is.infinite(m)] = min(m[!is.infinite(m)])
    # m[is.nan(m) | is.na(m)] = 0
    
    controli = sm$class=="control"
    # controln = sum(controli)

    # foldsip = foldres = NULL # save original
    for (uc in unique(sm$class[!controli])) {
      uci = sm$class==uc
      
      ## make test/train(x10) indices
      if ("train"%in%colnames(sm)) {
        tri = list(c = which(controli & sm$train),
                   e = which(uci & sm$train))
        tei = list(c = which(controli & !sm$train),
                   e = which(uci & !sm$train))
      } else {
        wu = which(uci)
        wui = sample(1:length(wu),length(wu)/2)
        wuii = setdiff(1:length(wu),wui)
        tri = list(c = which(controli),
                   e = wui)
        tei = list(c = which(controli),
                   e = wuii)
      }
      
      ## calc pvalues
      cm = m[controli,,drop=F]
      em = m[uci,,drop=F]
      cmtrs = m[tri$c,,drop=F]
      cmtes = m[tei$c,,drop=F]
      emtrs = m[tri$e,,drop=F]
      emtes = m[tei$e,,drop=F]
      
      pvs = llply(loopInd(1:ncol(m),no_cores), function(jj) {
        laply(jj, function(j) {
          # cm = m[controli,j];
          cm_ = cm[,j]
          em_ = em[,j]
          cmtr = cmtrs[,j]
          cmte = cmtes[,j]
          emtr = emtrs[,j]
          emte = emtes[,j]
          return(c(
            tryCatch({ 
              t.test(em_,cm_)$p.value }, error=function(e) {1}),
            tryCatch({ 
              wilcox.test(em_,cm_)$p.value}, warning=function(w) {1}),
            
            tryCatch({ 
              t.test(emtr,cmtr)$p.value }, error=function(e) {1}),
            tryCatch({ 
              wilcox.test(emtr,cmtr)$p.value},warning=function(w) {1}),
            
            tryCatch({ 
              t.test(emte,cmte)$p.value }, error=function(e) {1}),
            tryCatch({ 
              wilcox.test(emte,cmte)$p.value}, warning=function(w) {1})
          ))
        })
      },.parallel=T)
      pvs = Reduce(rbind,pvs); rownames(pvs) = colnames(m)
      pvs[is.na(pvs)|is.nan(pvs)] = 1
      pvs = list(t=list(all=pvs[,1],train=pvs[,3],test=pvs[,5]),
                 wilcox=list(all=pvs[,2],train=pvs[,4],test=pvs[,6]))

      for (ptype in names(pvs)) {
        for (adj in adjust) {
          pvs_ = NULL
          for (tretype in names(pvs[[ptype]])) {
            # pv = sapply(pvs, function(x) x[[ptype]][[tretype]]); 
            # pv[is.nan(pv) | is.na(pv)] = 1
            pv = pvs[[ptype]][[tretype]]

            ## adjust or combo p values
            if (adj%in%c("lanc","fisher")) {
              # if (grepl("group",feat_type)) next()
              
              phens = names(pv)
              if(grepl("_",phens[10])) phens = sapply(str_split(names(pv),"_"), function(x) x[1])
              nonp = gsub("[-]|[+]","",phens)
              allplus = which(grepl("[-]|[+]",phens) & !duplicated(nonp))
              if (length(allplus)==length(phens)) next
              # nonpu = unique(nonp)
              groupi = match(nonp, nonp)
              groups = llply(allplus, function(i) {
                ii = which(groupi==groupi[i])
                if (length(ii)>1) return(ii)
                return(NULL)
              })
              names(groups) = gsub("-","+",phens[allplus])
              groups = plyr::compact(groups)
              
              pv_ = a = rep(1,length(groups)); names(a) = names(groups)
              ap = sapply(groups, function(x) any(pv[x]<1))
              
              # don't use, not sure degrees of freedom...
              if (adj=="lanc") {
                if (any(ap)) pv_[ap] = 
                    laply(groups[ap],function(x) invchisq(pv[x],2)$p)
              }
              if (adj=="fisher") {
                if (any(ap)) pv_[ap] = 
                    laply(groups[ap],function(x) sumlog(pv_[x])$p)
              }
              names(pv_) = names(groups)
            } else {
              pv_ = p.adjust(pv, method=adj)
            }
            pvs_[[tretype]] = pv_
            # foldsip[[uc]]$p[[ptype]][[adj]][[tretype]] = pv_
          } # tretype
          pname = paste0(sum_p_dir,"/",feat_type,"_",uc,"_",ptype,"-",adj,".Rdata")
          save(pvs_,file=pname)
        } # adj
      } # ptype
    } #uc
    time_output(start2)
  }#, .parallel=F)) # feat_types
  
  time_output(start1)
} # result

time_output(start) # 10min


