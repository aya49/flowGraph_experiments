## input: features 
## output: p values and their corelation between train/test sets; p value qq plots, recall/precision, number of significant features

set.seed(10)

## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

## libraries
source("source/_func.R")
libr(c("stringr","colorspace", "Matrix", "plyr",
       "lattice", "gridExtra", # libr(proxy)
       "metap",
       "igraph",
       "foreach","doMC",
       "kernlab"))

graph_dir = paste0(root,"/pval_graphs.Rdata")



#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)



## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

# cvn-fold cross validation

overwrite = T #overwrite?
writecsv = F

readcsv = F

adjust = c("BY", "none")#,"fisher","none","BH","bonferroni") #pvalue adjustment; "lanc",
pthres = .05 # p value sig threshold for t test
good_count = 5
# minfold = 5 # minimum number of samples in each fold/class
good_sample = minfold = 5 # each class must have more thatn good_sample amount of samples or else prune
cvn0 = 10 # cvn-fold cross validation
testn = 1/2 #proportion of samples to make into test sample if none specified
# countThres = 1000 #insignificant if count under

id_col = "id"
target_col = "class"
order_cols = NULL

control = "control"
feat_count = "file-cell-countAdj"

## calcuate p values!

# format of p values:
# tri = foldsip[[uc]]$indices$train
# tei = foldsip[[uc]]$indices$test
# pv_tr_ = foldsip[[uc]]$p[[ptype]][[adj]]$train
# pv_tr2_ = foldsip[[uc]]$p[[ptype]][[adj]]$p_train2
# pv_te_ = foldsip[[uc]]$p[[ptype]][[adj]]$p_test
# pv_all_ = foldsip[[uc]]$p[[ptype]][[adj]]$p_all

start = Sys.time()
table = pvals = NULL
for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)){#[c(2,6,7,8)]) {
  if (grepl("pregnancy",result_dir)) next
  
  print(result_dir)
  data = fileNames(result_dir)
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta")
  meta_file_dir = paste(meta_dir, "/file", sep="") #meta for rows (sample)
  feat_dir = paste(result_dir, "/feat", sep="") #feature files directory
  meta_cell_childpn_names_dir = paste(meta_dir, "/cell_childpn_names",sep="") #specifies a phenotypes children and splits them into +/- (only for when both -/+ exists)
  # meta_cell_parent_names_dir = paste(meta_dir, "/cell_parent_names",sep="") #specifies a phenotypes parents
  
  ## output directories
  # unlink(paste0(result_dir,"/pval"))
  pvalsource_dir = paste0(result_dir,"/pval/src"); suppressWarnings(dir.create (pvalsource_dir, recursive=T))
  
  pvalgr_dir = paste0(result_dir,"/pval/graph")
  dir.create(pvalgr_dir, recursive=T, showWarnings=F)
  
  # make base graph
  meta_dir = paste0(result_dir,"/meta")
  meta_cell_childpn_names_dir = paste(meta_dir, "/cell_childpn_names",sep="") #specifies a phenotypes children and splits them into +/- (only for when both -/+ exists)
  meta_cell_childpn_names = get(load(paste0(meta_cell_childpn_names_dir, ".Rdata")))
  meta_cell_childpn_names_ = ldply(names(meta_cell_childpn_names), function(x) {
    to = unlist(meta_cell_childpn_names[[x]])
    data.frame(from=rep(x,length(to)), to=to) 
  })
  # meta_cell_parent_names = get(load(paste0(meta_cell_parent_names_dir, ".Rdata")))
  gr_e00 = as.matrix(meta_cell_childpn_names_)
  gr_v00 = append("",unique(as.vector(gr_e00)))
  gr_e00 = rbind(as.matrix(ldply(gr_v00[str_count(gr_v00,"[+-]")==1], function(x) c("",x))),gr_e00)
  gr_e00 = as.data.frame(gr_e00); colnames(gr_e00) = c("from","to")
  
  
  #data paths
  feat_types = list.files(path=feat_dir,pattern=".Rdata")
  feat_types = feat_types[!feat_types=="file-cell-count.Rdata"]
  feat_types = gsub(".Rdata","",feat_types)
  # feat_types = feat_types[!grepl("KO|Max",feat_types)]
  
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
  
  
  start1 = Sys.time()
  
  #load different features matrix and calculate distance matrix
  # for (feat_type in feat_types_) {
  result = llply(feat_types, function(feat_type) #for (feat_type in feat_types) 
  {
    # tryCatch({
    if (!file.exists(paste0(pvalsource_dir,"/",feat_type,".Rdata")) | overwrite) {
      start2 = Sys.time()
      cat("\n", feat_type, " ",sep="")
      
      ## upload and prep feature matrix + meta
      m0 = as.matrix(Matrix(get(load(paste0(feat_dir,"/", feat_type,".Rdata")))))
      sm = meta_file[match(rownames(m0),meta_file$id),]
      # good_sample
      tcl = table(sm$class); 
      delrow =rownames(m0)%in%sm$id[sm$class%in%names(tcl)[tcl<=minfold]]
      sm = sm[!delrow,]
      if (all(sm$class==sm$class[1])) return(F)
      m = m0[!delrow,]
      
      m[is.infinite(m)] = min(m[!is.infinite(m)])
      m[is.nan(m) | is.na(m)] = 0
      
      controli = sm$class=="control"
      controln = sum(controli)
      
      # make base graph
      all_sig_mc = colMeans(m[sm$class=="control",])
      etf = grepl("_",colnames(m)[1])
      if (etf) {
        # elist = as.data.frame(Reduce("rbind",str_split(all_sig_,"_")))
        # match0e = match(data.frame(t(elist)), data.frame(t(gr_e0)))
        # gr_e$p[match0e] = pv_all_
        # gr_e$mean_uc[match0e] = all_sig_me
        # gr_e$mean_ctrl[match0e] = all_sig_mc
        gr_e0 = as.data.frame(Reduce("rbind",str_split(colnames(m),"_")))
        colnames(gr_e0) = c("from","to")
        gr_e0$mean_ct = all_sig_mc
        gr_v0 = data.frame(name=append("",unique(as.vector(gr_e0))))
      } else {
        # gr = gr0 - setdiff(gr_v0, names(all_sig))
        gr_e0 = gr_e00[gr_e00[,1]%in%colnames(m) & gr_e00[,2]%in%colnames(m),]
        gr_v0 = data.frame(name=colnames(m), mean_ct=all_sig_mc)
      }
      
      
      
      foldsip = foldres = NULL # save original
      for (uc in unique(sm$class[!controli])) {
        
        
        ## save grph
        all_sig_me = colMeans(m[sm$class==uc,])
        gr_e = gr_e0
        gr_v = gr_v0
        if (etf) {
          gr_e$mean_uc = all_sig_me
        } else {
          gr_v$mean_uc = all_sig_me
        }
        a = list(e=gr_e,v=gr_v)
        save(a, file=paste0(pvalgr_dir,"/",feat_type,"_",uc,".Rdata"))
        
        
        ## make test/train(x10) indices
        uci = sm$class==uc
        if (sum(uci)<minfold) next
        
        if ("type"%in%colnames(sm)) {
          foldsip[[uc]]$indices$test = as.character(sm$id[sm$type=="test" & uci])
          foldsip[[uc]]$indices$train = as.character(sm$id[sm$type=="train" & uci])
        } else {
          testni = sample(which(uci), max(minfold,testn*sum(uci)))
          foldsip[[uc]]$indices$test = as.character(sm$id[testni])
          foldsip[[uc]]$indices$train = as.character(sm$id[!c(1:nrow(sm))%in%testni & uci])
        }
        
        # n-fold cv train indices; if only two folds? then still train/test; number of folds rounded
        trind = sample(foldsip[[uc]]$indices$train)
        if (length(trind)>(minfold*3)) {
          cvn = floor(min(cvn0,length(trind)/minfold))
          foldsip[[uc]]$indices$train = split(trind, cut(seq_along(trind), cvn, labels=F))
        }
        
        tri = foldsip[[uc]]$indices$train
        tei = foldsip[[uc]]$indices$test
        
        ## calculate p value
        # pvs = llply(1:ncol(m), function(j) {
        pvs = llply(loopInd(1:ncol(m),no_cores), function(jj) {
          llply(jj, function(j) {
            cm = m[controli,j]; 
            cmm = m[!controli,j]; 
            trm = m[unlist(tri),j]; 
            trm2 = NULL; if (is.list(tri)) trm2 = lapply(tri, function(trii) m[trii,j])
            tem = m[tei,j]
            cmtf = !all(cm==cm[1])
            
            t = wilcox = list()
            t$test = t$train = t$train2 = t$all = 
              wilcox$test = wilcox$train = wilcox$train2 = wilcox$all = 1
            if (!all(tem==tem[1]) & cmtf) {
              t$test = t.test(cm, tem)$p.value
              wilcox$test = wilcox.test(cm, tem)$p.value
            }
            if (!all(trm==trm[1]) & cmtf) {
              t$train = t.test(cm, trm)$p.value
              wilcox$train = wilcox.test(cm, trm)$p.value
              if (!is.null(trm2)) {
                t$train2 = median(sapply(trm2, function(trm2i) {
                  if (!all(trm2i==trm2i[1])) return(1)
                  return(t.test(cm, trm2i)$p.value)
                }))
                wilcox$train2 = median(sapply(trm2, function(trm2i) {
                  if (!all(trm2i==trm2i[1])) return(1)
                  return(wilcox.test(cm, trm2i)$p.value)
                }))
              }
            }
            if (!all(cmm==cmm[1]) & cmtf) {
              t$all = t.test(cm, cmm)$p.value
              wilcox$all = wilcox.test(cm, cmm)$p.value
            }
            return(list(t=t, wilcox=wilcox))
          })
        },.parallel=F)
        pvs = unlist(pvs, recursive=F)
        # })
        
        for (ptype in names(pvs[[1]])) {
          for (adj in adjust) {
            for (tretype in names(pvs[[1]][[ptype]])) {
              pv = sapply(pvs, function(x) x[[ptype]][[tretype]]); 
              pv[is.nan(pv) | is.na(pv)] = 1
              names(pv) = colnames(m)
              
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
              
              foldsip[[uc]]$p[[ptype]][[adj]][[tretype]] = pv_
            } # tretype
            
            pv_tr_ = foldsip[[uc]]$p[[ptype]][[adj]]$train
            pv_tr2_ = foldsip[[uc]]$p[[ptype]][[adj]]$train2
            pv_te_ = foldsip[[uc]]$p[[ptype]][[adj]]$test
            pv_all_ = foldsip[[uc]]$p[[ptype]][[adj]]$all
            
            ## calculate correlations between p values test & train/2
            pv_trl = -log(pv_tr_)
            pv_trl2 = -log(pv_tr2_)
            pv_tel = -log(pv_te_)
            # pv_alll
            pcorr = cor(pv_trl,pv_tel, method="spearman")
            pcorrp = cor.test(pv_trl,pv_tel, method="spearman")$p.value
            pcorr2 = cor(pv_trl2,pv_tel, method="spearman")
            pcorrp2 = cor.test(pv_trl2,pv_tel, method="spearman")$p.value
            
            # get sig
            tr_sig = pv_tr_<pthres
            tr2_sig = pv_tr2_<pthres
            te_sig = pv_te_<pthres
            all_sig = pv_all_<pthres
            
            overlap = sum(tr_sig & te_sig)
            rec = overlap/sum(te_sig)
            prec = overlap/sum(tr_sig)
            overlap2 = sum(tr2_sig & te_sig)
            rec2 = overlap2/sum(te_sig)
            prec2 = overlap2/sum(tr2_sig)
            
            # # calculate connectedness
            # all_sig_ = names(all_sig)[all_sig]
            # gr = NULL
            # if (length(all_sig_)>2) {
            #   # make edge list & graph
            #   etf = grepl("_",all_sig_[1])
            #   if (etf) {
            #     elist = as.data.frame(Reduce("rbind",str_split(all_sig_,"_")))
            #     colnames(elist) = c("from","to")
            #     gr = graph_from_data_frame(elist)
            #   } else {
            #     gr = gr0 - setdiff(gr0_v, all_sig_)
            #     
            #     # all_sig_ = all_sig_[order(sapply(all_sig_, str_count, "[+-]"))]
            #     # elist = ldply(all_sig_, function(node) {
            #     #     parent = meta_cell_parent_names[[node]]
            #     #     if (!is.null(parent)) {
            #     #       parent = parent[parent%in%all_sig_]
            #     #       if (length(parent)>0) 
            #     #         return( ldply(parent, function(x) data.frame(from=x, to=node)) )
            #     #     }
            #     #     return(NULL)
            #     # })
            #     # if (nrow(elist)>0) { # if there are edges
            #     #   if (length(V(gr)[[]])<all_sig_) { # if there is unconnected nodes
            #     #     all_sig_rest = all_sig_[!all_sig_%in%]
            #     #     gr = graph_from_data_frame(d=elist, vertices=all_sig_rest, directed=T)
            #     #   } else {
            #     #     gr = graph_from_edgelist(as.matrix(elist))
            #     #   }
            #     # } else { # no edges
            #     # }
            #     
            #   }
            #   
            #   # # graph stats
            #   # con = components(gr) # membership (cluster id/feat), csize (cluster sizes), no (of clusers)	
            #   # 
            #   # cl_no = con$no
            #   # out_no = sum(con$csize<3)
            #   # out_no2 = sum(con$csize<3)
            # }
            # foldsip[[uc]]$p[[ptype]][[adj]]$graph = gr
            
            
            # save table
            foldres = 
              rbind(
                foldres, data.frame(
                  data=data, feature=feat_type, class=uc, sig_test=ptype, adjust.combine=adj, p_thres=pthres, 
                  n_train_folds=ifelse(is.list(tri),1,length(tri)), n_samples_train=length(unlist(tri)), n_samples_test=length(tei), n_samples_control=sum(controln),
                  
                  m_test_sig=sum(te_sig), m_test=length(te_sig), 
                  m_test_perc=sum(te_sig)/length(te_sig),
                  m_train_sig=sum(tr_sig), m_train_sig2=sum(tr2_sig), m_train=length(tr_sig), 
                  m_train_perc=sum(tr_sig)/length(tr_sig), m_train_perc2=sum(tr2_sig)/length(tr2_sig),
                  m_overlap_sig=overlap, m_overlap_sig2=overlap2, 
                  
                  corr_spear=pcorr, corr_spear_p=pcorrp,
                  corr_spear2=pcorr2, corr_spear_p2=pcorrp2,
                  recall=rec, precision=prec, f=2*((prec*rec)/(prec+rec)),
                  recall2=rec2, precision2=prec2, f2=2*((prec2*rec2)/(prec2+rec2))
                ))
            
            
            # ## make graph
            # all_sig_ = names(all_sig)[all_sig]
            # 
            # gr_e = gr_e0
            # # gr0 = graph_from_edgelist(gr_e0)
            # gr = NULL
            # if (length(all_sig_)>2) {
            #   # make edge list & graph
            #   if (etf) {
            #     gr_e$p = pv_all_
            #     gr_e$mean_uc = all_sig_me
            #     gr_e$mean_ctrl = all_sig_mc
            #     gr = graph_from_data_frame(gr_e,directed=T)
            #   } else {
            #     # gr = gr0 - setdiff(gr_v0, names(all_sig))
            #     gr_e$p = ifelse(gr_e0[,1]%in%all_sig_ & gr_e0[,2]%in%all_sig_,0,1)
            #     gr_v = data.frame(
            #       name=names(all_sig), p=pv_all_,
            #       mean_uc=all_sig_me, mean_ctrl=all_sig_mc)
            # 
            #     gr = graph_from_data_frame(gr_e,directed=T,vertices=gr_v)
            #   }
            # }
            # if (is.null(gr)) next
            # save(gr, file=paste0(pvalgr_dir,"/",feat_type,"_",uc,"_",ptype,"_",adj,".Rdata"))
            
          } # adj
        } # ptype
        
      }
      save(foldsip, file=paste0(pvalsource_dir,"/",feat_type,".Rdata"))
      save(foldres, file=paste0(pvalsource_dir,"/",feat_type,"_table.Rdata"))
      time_output(start2)
    } # if
    foldsip = get(load(paste0(pvalsource_dir,"/",feat_type,".Rdata")))
    foldres = get(load(paste0(pvalsource_dir,"/",feat_type,"_table.Rdata")))
    
    return(list(foldsip=foldsip, foldres=foldres))
  }, .parallel=T)
  pvals[[data]] = llply(result, function(x) x$foldsip)
  names(pvals[[data]]) = feat_types
  table = rbind(table, ldply(result, function(x) x$foldres))
  
  time_output(start1)
} # result
save(table, file=paste0(root,"/pval_table.Rdata"))
save(pvals, file=paste0(root,"/pval.Rdata"))
time_output(start)


