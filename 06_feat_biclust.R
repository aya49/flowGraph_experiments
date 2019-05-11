## input: features 
## output: biclusters + some nice plots (if number of clusters/k required, the number of ground truth classes were provided)


## root directory
root = "/mnt/f/Brinkman group/current/Alice/flowtype_metrics"
setwd(root)


## libraries
source("source/_funcAlice.R")
source("source/_funcdist.R")
source("source/_bayesianbiclustering.R")
libr(c("biclust", "NMF","fabia","GrNMF", #devtools::install_github("jstjohn/GrNMF")
       "pheatmap",
       "foreach","doMC",
       "stringr", "Matrix", "plyr",
       "tcltk"))

## setup Cores for parallel processing (parallelized for each feature)
no_cores = 6#detectCores()-3
registerDoMC(no_cores)








## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

plotsize = 500

overwrite = T #overwrite biclust?
# writecsv = F

good_count = 1 #trim matrix; only keep col/rows that meet criteria for more than 3 elements
good_sample = 3 #trim matrix; only keep rows that are a part of a class with more than 3 samples

k = 6 # until what layer of cells to use
countThres = 1000 #a cell is insignificant if count under cell CountThres so delete -- only for matrices that have cell populations as column names
target_cols = c("class","gender") #the interested column in meta_file
control = "control" #control value in target_col column
id_col = "id" #the column in meta_file matching rownames in feature matrices
order_cols = NULL #if matrix rows should be ordered by a certain column
# split_col = "tube" # if certain rows in matrices should be analyzed in isolation, split matrix by this column in meta_file

bcmethods = c("plaid","CC","bimax","BBbinary","nmf","CC","GrNMF","fabia") #biclustering methods; GrNMF-<weight of graph regularization>
bmethodspar = list(#knn=c(1:6), 
  plaid=NA,
  CC=NA,
  bimax=NA,
  BBbinary=NA,
  nmf=c("nsNMF","lee","brunet"),
  GrNMF=c(0,1,5,10),
  fabia=NA)

#,"quest", "CC", "spectral", "Xmotifs", have to change this manually in function...
onlysigBB = T #only extract significant or all biclusters from BB-binary B2PS biclustering?
pval_thres = .05
Kr = 6; Kc = 20 #number of row and column biclusters for binary bayesian biclustering
pthres = .01 #pthreshold for choosing biclusters out of all bayesian biclusters
min_iter = 100 #BB-binary
sig_biclust_thres = .025 # * max contribution: threshold at whifeature matrices
qthres = .15 # quantile of nmf type methods factors; how large does factor effect need to be for bicluster to be significant

plot_size_bar = c(700,700)
plot_size = c(300,300)

for (result_dir in list.dirs(paste0(root, "/result"), full.names=T, recursive=F)) {
  # result_dir = paste0(root, "/result/impc_panel1_sanger-spleen") # data sets: flowcap_panel1-7, impc_panel1_sanger-spleen
  
  ## input directories
  meta_dir = paste0(result_dir,"/meta") # meta files directory
  meta_file_dir = paste(meta_dir, "/file", sep="") #meta for rows (sample)
  meta_cell_child_names_dir = paste(meta_dir, "/cell_childpn_names", sep="") #lists children cell population of each cell population; only used for graph regularized non matrix factorization
  meta_cell_dir = paste(meta_dir, "/cell", sep="")
  feat_dir = paste(result_dir, "/feat", sep="") #feature files directory
  
  ## output directories
  biclust_dir = paste(result_dir,  "/clust/feat", sep="") # path to save biclusterings: row clusterings, row labels, column colusterings
  biclust_plot_dir = paste(biclust_dir,  "/plot", sep="")
  dir.create (biclust_plot_dir,showWarnings=F, recursive=T)
  
  
  #feature matrix paths
  feat_types = list.files(path=feat_dir,pattern=".Rdata")
  feat_types = gsub(".Rdata","",feat_types)
  feat_count = "file-cell-countAdj" # cell count features used to trim matrix
  #feat_types = list.files(path=result_dir,pattern=glob2rx("matrix*.Rdata"))
  
  
  
  start = Sys.time()
  
  # read meta_file and count matrix (count matrix only to trim feature matrix of small cell populations -- only for matrices that have cell populations as column names)
  mc = Matrix(get(load(paste0(feat_dir,"/", feat_count,".Rdata"))))
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
  
  # make adjacency matrix ccm
  meta_cell_child_names = get(load(paste0(meta_cell_child_names_dir,".Rdata")))
  meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
  ccm = Matrix(0,nrow=nrow(meta_cell), ncol=nrow(meta_cell),sparse=T, dimnames=list(meta_cell$phenotype, meta_cell$phenotype))
  for (parent in names(meta_cell_child_names)) {
    tochild = meta_cell$phenotype %in% 
      unlist(meta_cell_child_names[names(meta_cell_child_names)==parent])
    ccm[parent,tochild] = ccm[tochild,parent] = 1
  }
  
  ## for each feature
  a = llply(feat_types, function(feat_type) {
    tryCatch({
      cat("\n", feat_type, " ",sep="")
      start2 = Sys.time()
      
      ## upload and prep feature matrix
      m = mbinary = m0 = as.matrix(Matrix(get(load(paste0(feat_dir,"/", feat_type,".Rdata")))))
      mbinary[mbinary != 0] = 1 #make matrix binary (for p values TRIM only)
      # sm = meta_file[match(rownames(m0),meta_file[,id_col]),]
      m = m[apply(m, 1, function(x) any(x>0)), apply(m, 2, function(x) any(x>0))]
      sm = meta_file[match(rownames(m),meta_file[,id_col]),]
      
      colhascell = !grepl("_",colnames(m)[1])
      
      # for (target_col in target_cols) {
      #   if (!target_col%in%colnames(meta_file)) next
      # las0 = sm$label
      
      #trim matrix
      # if (colhascell) mm = trimMatrix(m0,TRIM=T, mc=mc, sampleMeta=meta_file[match(colnames(m0),meta_file[,id_col]),], sampleMeta_to_m1_col=id_col, target_col=target_col, control=control, order_cols=order_cols, colsplitlen=NULL, k=k, countThres=countThres, goodcount=good_count, good_sample=good_sample)
      # if (is.null(mm)) next
      # m = mm$m
      # sm = mm$sm
      
      #list out class labels
      # class = sm[,target_col]; 
      # class_unique = unique(class)
      # if (length(class_unique)<2) next
      
      # dirname = paste0(target_col, "_",
      #                  ifelse(length(class_unique)<10, 
      #                         paste0(paste0(class_unique, ".",sapply(class_unique, function(x) sum(sm[,target_col]==x))), collapse="-"), length(class_unique)))
      # biclust_dir1 = paste0(biclust_dir,"/",dirname); dir.create(biclust_dir1, showWarnings=F)
      # biclust_plot_dir1 = paste0(biclust_plot_dir,"/",dirname); dir.create(biclust_plot_dir1, showWarnings=F)
      
      bcname1 = paste0(biclust_dir, "/", feat_type)
      # bcname2 = paste0(biclust_plot_dir1, "/", feat_type)
      
      # ## split up analysis of feature matrix rows by split_col
      # if (is.null(split_col)) {
      #   split_ind = list(all = 1:nrow(sm))
      # } else {
      #   split_ids = unique(sm[,split_col])
      #   split_ids = split_ids[!is.na(split_ids)]
      #   split_ind = lapply(split_ids, function(split_id) which(sm[,split_col]==split_id) )
      #   names(split_ind) = split_ids
      # }
      
      # for (tube in names(split_ind)) {
      #   m = m_ordered[split_ind[[tube]],]
      #   if (!sum(m_ordered!=0)>0) next
      #   sm = sm_split = sm[split_ind[[tube]],]
      #   if (length(unique(sm[,target_col]))<=1) next
      
      ## for each biclustering method
      # # plot
      # plotn = length(unlist(bmethodspar))
      # plotn2 = ceiling(sqrt(plotn))
      # png(bcname2, width=plotn2*plotsize, height=plotn2*plotsize)
      # par(mfrow=c(plotn2,plotn2))
      
      clusts = list()
      for (bcmethod in names(bmethodspar)) { #requires binary matrix
        for (par in bmethodspar[[bcmethod]]) {
          bcb = bc = NULL
          cat("\n ",bcmethod," ",par, " ",sep="")
          start2 = Sys.time()
          if (bcmethod=="BBbinary" & !grepl("pval[A-z]*TRIM",feat_type)) next
          
          # # where to save biclustering
          # bcname0 = paste("/",bcmethod, "_", feat_type, "_splitby-",split_col,".", tube, "_layer", str_pad(k, 2, pad = "0"), "_countThres-", countThres, sep = "")
          # bcname = paste0(biclust_dir, bcname0)
          
          # bicluster
          # if (overwrite | !file.exists(paste0(bcname,".Rdata"))) {
          if (bcmethod == "plaid") bc = biclust(as.matrix(m), method=BCPlaid(), row.release=.3,col.release=.7, back.fit=10, verbose = F)
          if (bcmethod == "CC") bc = biclust(as.matrix(m), method=BCCC(), number=Kr)
          if (bcmethod == "Xmotifs") bc = biclust(as.matrix(m), method=BCXmotifs(), number=Kr, ns=50, nd=500, alpha=10)
          if (bcmethod == "spectral") bc = biclust(as.matrix(m), method=BCSpectral(), numberOfEigenvalues=10)
          if (bcmethod == "bimax") bc = biclust(as.matrix(m), method=BCBimax(),number=Kr)
          if (bcmethod == "quest") bc = biclust(as.matrix(m), method=BCQuest(), number=Kr, ns=50)
          if (!is.null(bc)) if (bc@Number==0) next
          
          bcb = NULL
          if (bcmethod=="GrNMF" & colhascell) {
            #build binary relation graph between features
            cpind = match(colnames(m),colnames(ccm))
            medge = ccm[cpind,cpind]
            # bicluster
            grnmfm = t(as.matrix(m-min(m)))
            grnmfe = as.matrix(medge)
            bcb = grnmf(grnmfm, grnmfe, k=Kr, lambda_multiple=par, n_iter=max(ncol(m),1000), converge=1e-06, dynamic_lambda=T)
            
            if (is.null(bcb)) next
            bcb$rowxfactor = bcb$V
            bcb$factorxcol = t(bcb$U)
            if (all(is.na(bcb$rowxfactor))) next
          }
          
          # unused
          if (bcmethod=="fabia") {
            bcb0 = fabia(as.matrix(abs(m)), p=Kr,alpha=0.01,cyc=max(ncol(m),1000),spl=0,spz=0.5,non_negative=0,random=1.0,center=2,norm=1,scale=0.0,lap=1.0,nL=0,lL=0,bL=0)
            if (is.null(bcb0)) next
            
            bcb = list(source=bcb0)
            bcb$rowxfactor = bcb0@L; if(all(bcb$rowxfactor==0)) next
            bcb$factorxcol = bcb0@Z
            
            plot(bcb0)
          }
          
          if (grepl("nmf",bcmethod)) { tryCatch({ 
            bcb = nmf(abs(as.matrix(m)), rank=Kr, method=par)#, nrun=10, method=list("lee", "brunet", "nsNMF"))
            if (is.null(bcb)) next
            
            bcb$rowxfactor = basis(bcb)
            bcb$factorxcol = coef(bcb)
          }, error = function(err) { cat(paste("nmf error:  ",err)); bcb = NULL }) }
          
          # adjust format of biclustering to match output of biclust()
          if (grepl("nmf|GrNMF|fabia",bcmethod)) {
            # threshold to determine significant bicluster
            rthres = quantile(bcb$rowxfactor,qthres)
            cthres = quantile(bcb$factorxcol,qthres)
            
            # get the framework of biclust output to put nmf output into
            bc = biclust(array(0,dim=c(2,2)), method=BCPlaid()) #get the framwork
            bc@RowxNumber = array(F, dim=dim(bcb$rowxfactor)) 
            for(ri in 1:nrow(bcb$rowxfactor)) {
              mi = which.max(bcb$rowxfactor[ri,])
              if (max(bcb$rowxfactor[ri,])>rthres) bc@RowxNumber[ri,mi] = T
            }
            bc@NumberxCol = array(F, dim=dim(bcb$factorxcol)) 
            for(ri in 1:ncol(bcb$factorxcol)) {
              mi = which.max(bcb$factorxcol[,ri])
              if (max(bcb$factorxcol[,ri])>cthres) bc@NumberxCol[mi,ri] = T
              
              bc@Number = nrow(bcb$factorxcol)
              bc@info = list(rowxfactor=bcb$rowxfactor,factoxcol=bcb$factorxcol)
            } 
            
            ## plot x 4: factor heatplot + contribution of row/col
            bcname3 = paste0(biclust_plot_dir, "/", feat_type,"_", bcmethod, ifelse(!is.na(par),paste0("-", par),""), "_factors.png")
            png(bcname3, width=plotsize, height=plotsize)
            par(mfrow=c(2,2))
            
            rowxfactor = bc@info$rowxfactor
            plot(sort(rowxfactor[,1]),type="l", main="contribution of rows to each factor\n each line = factor" )
            for (fi in 2:ncol(rowxfactor)) lines(sort(rowxfactor[,fi]))
            
            aheatmap(rowxfactor, main="row x factor")
            
            factorxcol = bc@info$factoxcol
            plot(sort(factorxcol[1,]),type="l", main="contribution of cols to each factor\n each line = factor" )
            for (fi in 2:nrow(factorxcol)) lines(sort(factorxcol[fi,]))
            aheatmap(factorxcol, main="factor x col")
            
            graphics.off()
          }
          
          
          if (bcmethod == "BBbinary") {
            # create binary matrix as input into BB-binary
            bcb = B2PS(as.matrix(mbinary), sideData=NULL, Kt=Kr, Kp=Kc, iterations=max(ncol(mbinary)/20,min_iter), alpha_p = 1, alpha_t = 1, alpha_e = 1, alpha_sd = 1)
            # THETA is trasnposed!
            bcb$theta = t(bcb$Theta[,,2])
            
            # get significant biclusters only
            theta = bcb$theta
            theta = theta[sort(unique(bcb$transcript.clusters)),sort(unique(bcb$patient.clusters))]
            bcb$p = array(1,dim = c(nrow(bcb$theta), ncol(bcb$theta)))
            for (i in unique(bcb$transcript.clusters)) {
              for (j in unique(bcb$patient.clusters)) {
                bcb$p[i,j] = t.test.single(as.vector(theta),bcb$theta[i,j])
                if (is.na(bcb$p[i,j])) bcb$p[i,j] = 1
              }
            }
            # adjust format of biclustering to match output of biclust()
            bc = biclust(array(0,dim=c(2,2)), method=BCPlaid())
            bc@info = bcb
            if (onlysigBB) {
              bc@info$cid = cid = which(bcb$p<pval_thres, arr.ind=T)
              if (nrow(cid) == 0) {
                if (length(theta)>=6 & !all(bcb$p==1)) {
                  bc@info$cid = cid = which(bcb$p<=sort(bcb$p)[4], arr.ind=T)
                } else {
                  # save(NULL, file=paste0(bcname,".Rdata"))
                  next
                }
              }
            } else {
              bc@info$cid = cid = which(bcb$p<=1, arr.ind=T)
            }
            
            # cidmatrix = data.frame(row=rep(0,nrow(cid)), col=rep(0,nrow(cid)))
            # cidmatrix[cid] = 1
            # bc@info$sigcluster = cidmatrix
            bc@Number = nrow(cid)
            bc@RowxNumber = array(F, dim=c(nrow(m),nrow(cid))); colnames(bc@RowxNumber) = rep(1,ncol(bc@RowxNumber))
            bc@NumberxCol = array(F, dim=c(nrow(cid), ncol(m))); rownames(bc@NumberxCol) = rep(1,nrow(bc@NumberxCol))
            for (i in 1:nrow(cid)) { #column = cluster
              rowclust = cid[i,1]
              colclust = cid[i,2]
              colnames(bc@RowxNumber)[i] = rownames(bc@NumberxCol)[i] = bcb$p[rowclust,colclust]
              rows = bcb$transcript.clusters==rowclust
              cols = bcb$patient.clusters==colclust
              bc@RowxNumber[rows, i] = T
              bc@NumberxCol[i, cols] = T
              # cat("\n", sum(bcb$transcript.clusters==rowclust), ", ", sum(bcb$patient.clusters==colclust), sep="")
            }
            rowclust = bc@info$transcript.clusters
            colclust = bc@info$patient.clusters
          }
          
          if (nrow(bc@RowxNumber) != nrow(m) | ncol(bc@RowxNumber) != bc@Number) bc@RowxNumber = t(bc@RowxNumber)
          if (ncol(bc@NumberxCol) != ncol(m) | nrow(bc@NumberxCol) != bc@Number) bc@NumberxCol = t(bc@NumberxCol)
          
          ## get & save clustering & labels!
          
          if (bcmethod != "BBbinary") {
            rowclust = rowxcluster_to_cluster(bc@RowxNumber)
            colclust = clusterxcol_to_cluster(bc@NumberxCol)
          } 
          # rowlabel = sm[,target_col]
          
          if (length(unique(rowclust))==1) next
          
          # names(rowclust) = names(rowlabel) = rownames(m)
          names(rowclust) = rownames(m)
          names(colclust) = colnames(m)
          clusts[[paste0(bcmethod,".",par)]] = list()
          clusts[[paste0(bcmethod,".",par)]][["x"]] = list(x=rowclust, y=colclust, bc=bc)
          # f1 = f.measure.comembership(rowlabel,rowclust)
          
          # save biclustering
          # bc0 = list(source=bc,rowclust=rowclust,colclust=colclust)
          # bc0 = list(source=bc,rowclust=rowclust,colclust=colclust,rowlabel=rowlabel)
          # save(bc0, file=paste0(bcname,".Rdata"))
          # write.csv(rowclust, file=paste0(bcname,"_rowclust.csv"))
          # write.csv(colclust, file=paste0(bcname,"_colclust.csv"))
          # write.csv(rowlabel, file=paste0(bcname,"_rowlabel.csv"))
          
          # #save row/col as csv
          # try({
          #   bcgene = bc@RowxNumber; rownames(bcgene) = sm[,target_cols[1]]
          #   bcgene = bcgene[apply(bcgene, 1, function(x) all(!x)),]
          #   write.csv(bcgene,file=paste0(bcname,"_row.csv"))
          # })
          # try ({
          #   bccol = bc@NumberxCol
          #   if (ncol(bccol)==bc@Number) bccol = t(bccol)
          #   colnames(bccol) = colnames(m)
          #   bccol = bccol[,apply(bccol, 2, function(x) all(!x))]
          #   write.csv(bccol,file=paste0(bcname,"_col.csv"))
          # })
          
          
          ## plot ---------------------------------
          
          ## pretty heatmap
          # prepare col/row annotation for pretty heatmap
          if (bcmethod=="BBbinary") {
            colcluster_temp = bc@info$patient.clusters              
            colcluster_temp[!colcluster_temp%in%bc$colclust] = 0
            rowcluster_temp = bc@info$transcript.clusters
            rowcluster_temp[!rowcluster_temp%in%bc$rowclust] = 0
          } else {
            colcluster_temp = unlist(apply(bc@NumberxCol, 2, function(x) 
              ifelse(sum(x)==0, 0, max(which(x))) ))
            rowcluster_temp = unlist(apply(bc@RowxNumber, 1, function(x)
              ifelse(sum(x)==0, 0, max(which(x))) ))
          }
          # plot
          for (target_col in target_cols) {
            tryCatch ({
              if (!target_col%in%colnames(sm)) next
              
              #list out class labels
              class = sm[,target_col];
              class_unique = unique(class)
              if (length(class_unique)<2) next
              
              # plot path
              dirname = paste0(target_col, "_",
                               ifelse(length(class_unique)<10,
                                      paste0(paste0(class_unique, ".",sapply(class_unique, function(x) sum(sm[,target_col]==x))), collapse="-"), length(class_unique)))
              # biclust_dir1 = paste0(biclust_dir,"/",dirname); dir.create(biclust_dir1, showWarnings=F)
              biclust_plot_dir1 = paste0(biclust_plot_dir,"/",dirname); dir.create(biclust_plot_dir1, showWarnings=F)
              # bcname1 = paste0(biclust_dir, "/", feat_type)
              bcname2 = paste0(biclust_plot_dir1, "/", feat_type, "_", bcmethod, ifelse(!is.na(par),paste0("-", par),""))
              
              
              # prep data
              col_annot = data.frame(colcluster=factor(colcluster_temp))
              row_annot = data.frame(rowcluster=factor(rowcluster_temp), class=factor(class))
              row_order = order(class, row_annot$rowcluster, decreasing=T)
              col_order = order(col_annot$colcluster, decreasing=T)
              annotation_row = row_annot[row_order,]
              annotation_col = data.frame(col_annot[col_order,])
              
              mh = as.matrix(m)[row_order,col_order]
              rownames(mh) = 1:nrow(mh)
              rownames(annotation_col) = colnames(mh)
              rownames(annotation_row) = 1:nrow(mh)
              # rownames(col_annot)=names(colclust)
              
              # pheatmap(mh, main=paste0(target_col),#,"; precision, recall, f measure comembership:\n",paste0(c(f1$p,f1$r,f1$f_comember),collapse=", ")),
              #          annotation_row=annotation_row, annotation_col=annotation_col,
              #          show_rownames = F, 
              #          show_colnames=F,
              #          cluster_cols = T, cluster_rows = T,
              #          # cellwidth = 3, cellheight = 3, 
              #          # fontsize = 3, 
              #          filename = paste0(bcname2, "_pheatmap.png"))
              
              pheatmap(mh, main=paste0(target_col),#,"; precision, recall, f measure comembership:\n",paste0(c(f1$p,f1$r,f1$f_comember),collapse=", ")),
                       annotation_row=annotation_row, annotation_col=annotation_col,
                       show_rownames = F, 
                       show_colnames=F,
                       cluster_cols = F, cluster_rows = F,
                       # cellwidth = 3, cellheight = 3, 
                       fontsize = 3, 
                       filename = paste0(bcname2, "_pheatmap-sorted.png"),
                       width=3, height=2) #inches
              
            }, error = function(err) { cat(paste("pheatmap error:  ",err)); return(T) })
            
            ## plot row clusters against different sample target_cols
            tryCatch({
              
              # png(paste0(clust_path_plot, "_stats_",target_col,".png", sep=""), height=plot_size[1]*rowcolno, width=plot_size[2]*rowcolno)
              # par(mar=c(10,5,3,2), mfrow=rep(rowcolno,2))
              target_col_valuesL = list()
              target_col_topics = c()
              for (BCi in 1:bc@Number) {
                target_col_values = sm[bc@RowxNumber[,BCi],target_col]
                target_col_valuesL[[BCi]] = target_col_valuesT = table(target_col_values)
                target_col_topics = union(target_col_topics,names(target_col_valuesT))
                # try({
                #   barplot(target_col_valuesT, las=2, xlab=target_col, ylab="# of samples in bicluster with target_col on x axis", main=paste0(length(target_col_values)," samples in bicluster ", BCi))
                # })
              }
              # graphics.off()
              
              target_col_topics = sort(target_col_topics)
              target_col_valuesM = sapply(target_col_valuesL, function(x) {
                sapply(target_col_topics, function(y) {
                  if (y %in% names(x)) return (x[y])
                  return (0)
                })
              })
              target_col_valuesM = t(target_col_valuesM)
              if (all(colnames(target_col_valuesM)==colnames(target_col_valuesM)[1])) target_col_valuesM = t(target_col_valuesM)
              colnames(target_col_valuesM) = target_col_topics
              rownames(target_col_valuesM) = c(1:nrow(target_col_valuesM))
              # avm = target_col_valuesM/20
              
              png(paste0(bcname2, "_vs.png", sep=""), height=500, width=500)
              par(mar=c(10,5,5,8), xpd=TRUE)
              
              colour = rainbow(ncol(target_col_valuesM))
              barplot(t(target_col_valuesM), xlab="bicluster",ylab="# of samples", col=colour, main=paste0("# samples in biclusters with target_col of ",target_col,"\n",target_col))#,"; precision, recall, f measure comembership:\n",paste0(c(f1$p,f1$r,f1$f_comember),collapse=", ")))
              legend("topright",legend=colnames(target_col_valuesM),fill=colour,inset=c(-.2,0))
              # plot(row(avm), col(avm),
              #   cex=avm,
              #   xlim=c(0.5,nrow(avm)+0.5), ylim=c(0.5,ncol(avm)+0.5),
              #   axes=FALSE, ann=FALSE
              # )
              # text(row(avm), col(avm), paste0(target_col_valuesM, "\nsamples"), col="brown", pos=3)
              # axis(1,at=1:nrow(avm),labels=rownames(avm),cex.axis=0.8)
              # axis(2,at=1:ncol(avm),labels=colnames(avm),cex.axis=0.8)
              # title(xlab="bicluster",ylab=target_col)
              # box()
              graphics.off()
            }, error = function(err) { cat(paste("target_col plot error:  ",err)); return(T) })
            
          }
          # # Specifying clustering from distance matrix
          # drows = dist(test, method = "minkowski")
          # dcols = dist(t(test), method = "minkowski")
          # pheatmap(test, clustering_distance_rows = drows, clustering_distance_cols = dcols)
          
          
          
          # png(paste0(clust_path_plot, "_bar.png", sep=""), height=plot_size_bar[1], width=plot_size_bar[2])
          # par(mar=c(2,6,3,2))
          # try({
          #   plotclust(bc,as.matrix(m))
          # })
          # graphics.off()
          
          # png(paste0(clust_path_plot, "_bubble.png", sep=""), height=plot_size_bar[1], width=plot_size_bar[2])
          # par(mar=c(20,20,20,20))
          # try({
          #   bubbleplot(as.matrix(m),bc)
          # })
          # graphics.off()
          
          
          ## plot heatmaps
          tryCatch ({
            bcname2 = paste0(biclust_plot_dir, "/", feat_type, "_", bcmethod, ifelse(!is.na(par),paste0("-", par),""))
            
            png(paste0(bcname2,"_heatmap0.png"), height=plot_size_bar[1], width=plot_size_bar[2])
            par(mar=c(5,3,6,5))
            try({
              heatmapBC(as.matrix(m),bc, order=T, local=T, outside=T)
            })
            graphics.off()
            
            rowcolno = ceiling(sqrt(bc@Number))
            png(paste0(bcname2,"_heatmap.png"), height=plot_size[1]*rowcolno, width=plot_size[2]*rowcolno)
            par(mar=c(25,50,20,50))
            par(mfrow=rep(rowcolno,2))
            for (BCi in 1:bc@Number) {
              drawHeatmap(as.matrix(m),bc,BCi,plotAll=T)
            }
            graphics.off()
          }, error = function(err) { cat(paste("heatmap error:  ",err)); return(T) })
          
          
          
          
        } # par
        
        # #temporary
        # if (nrow(bc@RowxNumber) != nrow(m) | ncol(bc@RowxNumber) != bc@Number) bc@RowxNumber = t(bc@RowxNumber)
        # if (ncol(bc@NumberxCol) != ncol(m) | nrow(bc@NumberxCol) != bc@Number) bc@NumberxCol = t(bc@NumberxCol)
        # if (is.null(bc0$rowclust) & bcmethod!="BB-binary") bc0$rowclust = rowxcluster_to_cluster(bc@RowxNumber)
        # try({
        #   # get hard clustering; put samples into larger cluster
        #   la = as.numeric(factor(sm[,target_col]))
        #   if (bcmethod=="BB-binary") {
        #     cl = bc0$rowclust = bc@info$transcript.clusters
        #   } else {
        #     cl = bc0$rowclust
        #   }
        #   bc@info$f1 = f.measure.comembership(la,cl)
        # })
        # save(bc, file=paste0(bcname,".Rdata")); if (writecsv) write.csv(as.matrix(bc), file=paste0(checkm(bc,bcname),".Rdata"))
        
        # rowlabel = as.numeric(factor(sm[,target_col]))
        
        # }
        
      } # bcmethod
      save(clusts, file=paste0(bcname1,".Rdata"))
      # } #target_col
      time_output(start2)
      return(F)
    }, error = function(err) { cat(paste("ERROR:  ",err)); return(T) })
    return(F)
  }, .parallel=T)
  
  
  
  time_output(start)
  
}


