## Input: scores --> Output: tables and plots
# aya43@sfu.ca 20170419

#Directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)
options(include.rownames=F)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="") #metadata for cell populations (phenotype)
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="") #metadata for FCM files
# sampleMetaTrain_dir = paste0("attachments/AMLTraining.csv") #which FCM files are testing/training files
matrix_dir = paste(result_dir, "/matrix", sep="") #feature matrices
rw_dir = paste(result_dir,  "/rw", sep="")
dist_dir = paste(result_dir, "/dist", sep="")
plot_dir = paste(result_dir, "/plots", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dir[i])) }
plot_dist_dir = paste(plot_dir, "/dist", sep=""); for(i in 1:length(plot_dist_dir)) { suppressWarnings(dir.create(plot_dist_dir[i])) }
plot_dist_source_dir = paste0(plot_dist_dir, "/source"); for(i in 1:length(plot_dist_source_dir)) { suppressWarnings(dir.create(plot_dist_source_dir[i])) }

#Output
plot_dist_eval_dir = paste(plot_dir, "/dist_eval", sep=""); for(i in 1:length(plot_dist_eval_dir)) { suppressWarnings(dir.create(plot_dist_eval_dir[i])) }
dist_score_dir = paste(dist_dir, "/dist_score", sep=""); for (i in 1:length(dist_score_dir)) { suppressWarnings(dir.create(dist_score_dir[i])) }
cl_score_result_dir = paste0(dist_score_dir,"/score_cl_list")
dist_clusterf_dir = paste(dist_score_dir, "/dist_cluster_f", sep=""); for (i in 1:length(dist_clusterf_dir)) { suppressWarnings(dir.create(dist_clusterf_dir[i])) }
dist_clustercl_dir = paste(dist_score_dir, "/dist_cluster_cl", sep=""); for (i in 1:length(dist_clustercl_dir)) { suppressWarnings(dir.create(dist_clustercl_dir[i])) }
dist_kern_dir = paste(dist_dir, "/dist_kern", sep=""); for (i in 1:length(dist_kern_dir)) { suppressWarnings(dir.create(dist_kern_dir[i])) }
dist_score_all_dir = paste0(dist_score_dir,"/dist_score_all")
dist_comp_dir = paste(result_dir, "/dist_compare", sep="")
plot_compare_dir = paste(plot_dir, "/dist_compare", sep=""); for(i in 1:length(plot_compare_dir)) { suppressWarnings(dir.create(plot_compare_dir[i])) }
dist_comp_dist_dir = paste(result_dir, "/dist_compare_dist", sep="")
tables_dir = "tablesD.txt" #open this file up, too many escape characters, just do manually: replace \{ with {, \} with }, backslashes with \

#Libraries/Functions
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
libr("stringr")
libr("foreach")
libr("doMC")
libr("doSNOW")
libr("lubridate") #if there are date variables
libr("FastKNN")
libr("fpc")
libr("cluster")
libr("mclust")
libr("igraph")
libr("kernlab")
libr("densitycut") #devtools::install_bitbucket("jerry00/densitycut_dev")
libr("vegan")
libr("Rtsne")
libr("fpc")
libr("PerfMeas")
libr("clues")
libr("clusteval")
libr("xtable")

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)















#Options for script
width = 700; height = 700 # of each plot in png

cltypes = c("distmatrix")#,"knn","kmed","lv","spec1","spec","hc","dc1","dc")#,"spec","dc") #"rw1" (produces one cluster...)
cltypesclass = c("knn") # which is classification
cltypesorigm = c("spec","dc") #uses original matrix
sctypes_ivdist = c("NCA", "silmed") #metrics for distmatrix only

cltypes0 = c("distmatrix", "spec1","knn","kmed","hc","lv") #for plot and table making
cltypespar0 = list(distmatrix="none", knn=c(1), kmed="none", lv=c(.1), spec="none", spec1="none", hc=c("ward.D2"), rw1=c(.05)) #for plot and table making
writeMeta = c("feature","featend","layer","dist") #meta data to include in tables

matrix_type_features = c("CountAdj","Prop","LogFold","Pval","rw_Pval","rw_CountAdj","Child_entropy","Parent_entropy","Child_prop","Child_pnratio","Parent_contrib","Parent_effort","Freqp0.5")
ignoredist = ".csv|simmatrix"

dis = c("euclidean", "maximum", "manhattan", "canberra", "binary", "bray", "kulczynski", "minkowski", "morisita", "horn", "binomial", "rbf","poly","vanilla","tanh","laplace","bessel","anova","spline")

interested = c("gene") #sampleMeta columns to plot
splitby = c("none") #match with above
nophenodevcol = c("")
phenodevfeat = "contrib|effort|pval|logfold"
plotcols = c("feature","layer","dist") #split columns to plot
plotcols0 = c("Features","Layers","Distance Metrics") #for display

evscore = "f_comember_co"
evclassscore = "f_comember_co_1"
ivscore = c("silmed","pearsongamma")
ivdscore = c("NCA")
originalscore = c("f_comember_co","f_comember_co_1","silmed","pearsongamma","NCA")
displayscore = c("F Measure","F Measure","Median Silhouette Index","Pearson Gamma","NCA")

doDistPlot = F #plot distances against each other
overwritedist = F
doRtsne = T
doPCA = T

plotpc = 5 #number of pca pc to plot

quantilnorm = T #quantile norm distances for comparison?



#Prepare data
distmfile1 = list.files(dist_clusterf_dir, recursive=F, full.names=T, pattern=".Rdata$")
fs = file.size(distmfile1)
distmfile10 = distmfile1 = distmfile1[order(fs)]; fs = sort(fs)
distmfile1names = fileNames(distmfile1)
#head(distmfile1names[fs<30000])
distMeta1 = distMetafun(distmfile1,c(dis,"other"),features=matrix_type_features)
# distMeta1$size = fs

distmfile2 = list.files(dist_dir, recursive=F, full.names=T, pattern=".Rdata$")
distMeta2 = distMetafun(distmfile2,dis=c(dis,"other"),features=matrix_type_features) #for comparing dist matrices














## collating scores
start=Sys.time()

f = foreach(di=distmfile1) %dopar% { get(load(di)) }
cm = foreach(di=distmfile1) %dopar% { get(load(gsub("dist_cluster_f","dist_cluster_cl",di))) }
names(f) = names(cm) = fileNames(distmfile1)

result = foreach(i = names(f)) %dopar% { #dist
  b = NULL
  a = NULL
  aind = 1
  for (j in names(f[[i]])) { #col
    for (k in names(f[[i]][[j]])) { #dind
      for (l in names(f[[i]][[j]][[k]])) { #cltype
        for (m in names(f[[i]][[j]][[k]][[l]])) { #par
          if (l=="knn") { clusters = 2
          } else { clusters = length(unique(cm[[i]][[j]][[k]][[l]]$clt[,m])) }
          b = rbind(b, c(i,j,k,l,m,clusters))
          aa = sapply(f[[i]][[j]][[k]][[l]][[m]], function(x) { #some are lists, so need to get them all to be vector
            if (!is.na(x) & !is.na(x)) return(x)
            return(NA)
          })
          names(aa) = names(f[[i]][[j]][[k]][[l]][[m]])
          a[[aind]] = aa[!is.na(aa)]
          aind = aind+1
        }
      }
    }
  }
  return(list(rn=b, score=a))
}

TimeOutput(start)











## collate scores and their attributes ------------------------------------------

scores = foreach(ii=1:length(result)) %dopar% { result[[ii]]$score }
scores01 = scores = unlist(scores, recursive=F)
rn0 = foreach(ii=1:length(result),.combine="rbind") %dopar% { result[[ii]]$rn }
colnames(rn0) = c("fpath","colnam","dindname","cltype","par","clusters")

nrow(rn0)==length(scores) #should be true

#prepare attributes table
rn01 = rn = cbind(distMeta1[match(rn0[,1],fileNames(distMeta1$path)),],rn0[,-1])

# #temporary cut, to make it go faster, collate only needed stuff for thesis
# nopropfreq = !grepl("Prop|Freqp",rn$type,ignore.case=F)
# nonorm = ((grepl("effort|contrib|prop|pnratio",rn$type) & rn$norm=="none") | (!grepl("effort|contrib|prop|pnratio",rn$type) & rn$norm=="cellpop"))
# norand = rn$rand==0
# rn = rn[nonorm & nopropfreq & norand,]
# scores = scores[nonorm & nopropfreq & norand]

#collate scores into a table (NA for missing)
metricnames = Reduce('union',lapply(scores, function(x) names(x)))
# score0 = foreach(sc=scores, .combine='rbind') %dopar% {
#   a = rep(NA,length(metricnames))
#   a[match(names(sc),metricnames)] = sc
#   return(a)
# }
score0 = list()
for(sc in scores) {
  a = rep(NA,length(metricnames))
  a[match(names(sc),metricnames)] = sc
  score0[[length(score0)+1]] = a
}
score0 = Reduce("rbind",score0)
# score00 = Reduce("rbind",score0[grepl("rw_",rn01$path)])
rownames(score0) = rn[,1]
colnames(score0) = metricnames

score1 = cbind(rn[,-1],score0)
rownames(score1) = rn[,1]
for (cltype in score1$cltype) {
  if (any(grepl("none",unique(score1[score1$cltype==cltype,"par"]))) & length(unique(score1[score1$cltype==cltype,"par"]))>1) {
    score1[score1$cltype==cltype,"par"] = "none"
    print(unique(score1[score1$cltype==cltype,"par"]))
  }
}

save(score1,file=paste0(dist_score_all_dir,".Rdata"))
write.csv(score1,file=paste0(dist_score_all_dir,".csv"))

TimeOutput(start)










## trim and order score matrix -------------------------------------------------

score1 = get(load(paste0(dist_score_all_dir,".Rdata")))
scorecol = 18 #column where score starts, can also just look at the score only matrix...

#temporary cut, to make it go faster, collate only needed stuff for thesis
nopropfreq = !grepl("Prop|Freqp",score1$type,ignore.case=F)
nonorm = ((grepl("effort|contrib|prop|pnratio",score1$type) & score1$norm=="none") | (!grepl("effort|contrib|prop|pnratio",score1$type) & score1$norm=="cellpop"))
norand = score1$rand==0

score2 = score1[nopropfreq & nonorm & norand & !score1$weighted,]
score2$featend = sapply(score2$featend, function(x) {
  if (grepl("TRIM",x)) return("TRIM")
  return("")
})

# score2 = score2[order(unlist(score2[,"feature"]),score2[,"layer"],score2[,"dist"],score2[,"featend"]),]
score2 = score2[!duplicated(sapply(1:nrow(score2), function(x) paste0(score2[x,c("layer","dist","feature","featend","colnam","dindname","cltype","par","weighted","weightedorig")],collapse=" "))),]











## make latex score tables ------------------------------------------------------
#splits tables into trimmed and untrimmed, regardless of where those trimming p values come from (i.e. prop/countadj)

sink(tables_dir)
colfunc = colorRampPalette(c("red","black"))
colours = colfunc(4) #shade scores relative to the same interested col, variable and score

scorelist = list()

#for each interested column
for (colnam in unique(score2$colnam)) { 
  scorelist[[colnam]] = list()
  score3 = score2[score2$colnam==colnam,]
  if (str_extract(colnam,"interested-[a-z]+")%in%nophenodevcol) score3 = score3[!grepl(phenodevfeat,score3$feature,ignore.case=T) & !grepl("TRIM",score3$featend),]
  
  #for each variable type in interested column
  # dindname = ifelse(str_extract(colnam,"interested-[a-z]+")%in%nophenodevcol,"all","avg")
  for (dindname in unique(score3$dindname)) {
    score4 = score3[score3$dindname==dindname,]
    
    #for each classification/clustering algorithm
    for (cltype in cltypes0) {
      score5 = score4[score4$cltype==cltype,]
      
      #for each classification/clustering parameter
      for (parr in cltypespar0[cltype]) {
        if (!"none"%in%cltypespar0[cltype]) score5 = score5[score5$par==parr,]
        
        scores = NULL
        if (cltype%in%c("distmatrix")) { scores = cbind(scores,score5[,c(ivdscore)]); colnames(scores) = ivdscore
        } else if (cltype%in%cltypesclass) { scores = cbind(scores,score5[,paste0(evscore,"_1")]); colnames(scores) = paste0(evscore,"_1")
        } else { scores = score5[,c(evscore,ivscore)]; colnames(scores) = c(evscore,ivscore) }
        if (is.null(dim(scores))) next
        
        #meta data columns
        scoresrn = score5[,match(writeMeta,colnames(score5))]
        
        #colour scores
        for (ncols in 1:ncol(scores)) {
          one = scores[,ncols] > quantile(score4[,colnames(scores)[ncols]],.95,na.rm=T)
          two = scores[,ncols] > quantile(score4[,colnames(scores)[ncols]],.8,na.rm=T)
          three = scores[,ncols] > quantile(score4[,colnames(scores)[ncols]],.5,na.rm=T)
          scores[,ncols] = as.character(signif(as.numeric(scores[,ncols]),4))
          scores[one & !is.na(scores[,ncols]),ncols] = paste0("\\textcolor{one}{", signif(as.numeric(scores[one & !is.na(scores[,ncols]),ncols]),4), "}")
          scores[!one & two & !is.na(scores[,ncols]),ncols] = paste0("\\textcolor{two}{", signif(as.numeric(scores[!one & two & !is.na(scores[,ncols]),ncols]),4), "}")
          scores[!one & !two & three & !is.na(scores[,ncols]),ncols] = paste0("\\textcolor{three}{", signif(as.numeric(scores[!one & !two & three & !is.na(scores[,ncols]),ncols]),4), "}")
        }
        
        #print un/trimmed tables
        print(paste("\n",cltype, colnam, dindname))
        tabless = cbind(scoresrn,scores)
        tabless = tabless[order(tabless[,"feature"],tabless[,"layer"],tabless[,"dist"],tabless[,"featend"]),]
        
        if (any(grepl("TRIM",tabless$featend))) {
          tables1 = tabless[grepl("TRIM",tabless$featend) & !duplicated(sapply(1:nrow(tabless), function(x) paste0(tabless[x,c(1:4)],collapse=" "))),]
          rownames(tables1) = NULL
          print(xtable(tables1),include.rownames=F)#, sanitize.text.function = function(x) x)
        } 
        tables2 = tabless[grepl("TRIM",tabless$featend) & !duplicated(sapply(1:nrow(tabless), function(x) paste0(tabless[x,c(1:4)],collapse=" "))),]
        rownames(tables2) = NULL
        print(xtable(tables2),include.rownames=F)
      }
    }
  }
}
sink()













## plot each cltype, colnam, column = dindname & par -----------------------------------
# score2 = score1[nopropfreq & nonorm & norand & !score1$weighted,]

start = Sys.time()

score20 = score2[score2$cltype%in%cltypes0,]

#set colours for features
colours = rainbow(length(unique(score20$feature)))
names(colours) = unique(score20$feature)
colours["Child_prop"] = "yellow2"

#for each interested column
errorsplot = foreach (colnam = unique(score20$colnam)) %dopar% { #unique(score2$cltype)) %dopar% {
  tryCatch({
    score3 = score20[score20$colnam==colnam,]
    if (str_extract(colnam,"interested-[a-z]+")%in%nophenodevcol) {
      score3 = score3[!grepl("contrib|effort|pval|logfold",score3$feature,ignore.case=T) & !grepl("TRIM",score3$featend),]
    } else if (sum(!grepl("all",score3$dindname))>0) { score3 = score3[!grepl("all",score3$dindname),] }
    
    #for each score
    for (score in scorecol:ncol(score3)) {
      
      #for each classification/clustering
      for (cltype in cltypes0) {
        
        #scorename to use as column name
        scorename = colnames(score3)[score]
        if (cltype%in%c("distmatrix") & !scorename%in%c(ivscore,ivdscore)) { next
        } else if (cltype%in%cltypesclass & !scorename%in%evclassscore) { next
        } else if (!cltype%in%c("distmatrix",cltypesclass) & !scorename%in%c(evscore,ivscore)) { next }
        
        #scorename to use as as display
        scorename0 = displayscore[match(scorename,displayscore)]
        
        score4 = score3[score3$cltype==cltype,] 
        
        #plot scores on same scale; set min and max score values on plot scale
        if (sum(is.na(score4[,score]))>0) score4[is.na(score4[,score]),score] = min(score4[!is.na(score4[,score]),score])
        if (sum(score4[,score]==Inf)>0) score4[score4[,score]==Inf,score] = max(score4[!score4[,score]==Inf | !score4[,score]==-Inf,score])
        if (sum(score4[,score]==-Inf)>0) score4[score4[,score]==-Inf,score] = min(score4[!score4[,score]==Inf | !score4[,score]==-Inf,score])
        if (quantile(score4[score4[,score]!=0,score]/max(score4[score4[,score]!=0,score]),.9,na.rm=T)<.01) { next }
        
        ncolplot = length(unique(unique(score4$par)))*length(unique(score4$dindname))
        
        #for each parameter
        pars = cltypespar0 #unique(score4$par)
        for (parr in pars) {
          score5 = score4
          if (parr!="none") score5 = score4[score4$par==parr,]
          
          ## plot
          # pngname = paste0(plot_compare_dir,"_thesis/", scorename,"__",cltype,"__",colnam,"__par-",parr,".png")
          pngname = paste0(plot_compare_dir,"/", scorename,"__",cltype,"__",colnam,"__par-",parr,".png")
          png(filename=pngname, height=height*length(plotcols), width=width)#*(length(unique(score4$dindname))))
          par(mfcol=c(length(plotcols),length(pars)), mar=c(13,6,3,15), xpd=T)
          
          # for (dindname in ifelse(str_extract(colnam,"interested-[a-z]+")%in%nophenodevcol,"all","avg")) {
          for (dindname in sort(unique(score4$dindname))) {
            score6 = score5[score5$dindname==dindname,]
            if (dindname=="avg") score6 = score5[!score5$dindname%in%c("all","avg"),]
            main1 = paste0(colnam,"-",dindname)
            
            y_score = score6[,score]
            y_scoresize = .1 + (y_score-min(y_score))/max(y_score-min(y_score))
            
            #for each plot
            for (splitcols in plotcols) {
              #y label
              splitcols0 = plotcols0[match(splitcols,plotcols0)]
              
              #x label
              x_splitnames = sort(unique(score6[,splitcols]))
              if (length(grep("rbf",x_splitnames))>0) x_splitnames[grep("rbf",x_splitnames)] = "Gaussian Kernel" #for display
              if (is.na(as.numeric(x_splitnames[1]))) x_splitnames = sapply(x_splitnames, function(x) paste(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)), sep=""))
              
              #values
              x_split0 = as.integer(factor(score6[,splitcols]))
              x_split = jitter(x_split0, factor=1)
              xy_boxplot = lapply(sort(unique(x_split0)), function(x) y_score[x_split0==x])
              
              #plot
              boxplot(xy_boxplot, lwd = 1, outline=F, ylim=c(min(score3[!is.na(score3[,score]) & score3[,score]!=Inf,score]),max(score3[!is.na(score3[,score]) & score3[,score]!=Inf,score])),
                      main = splitcols0, #paste0(main1,"\n",cltype,"_par-", parr,"  ;;  ", splitcols, " vs ", scorename), 
                      xaxt ='n', yaxt='n', ylab=scorename0, cex.lab=2, cex.main=2, font.lab=2) #,xaxt ='n', xlab=testcol
              axis(1, at=1:length(x_splitnames), labels=x_splitnames, las=2, cex.axis = 1.5)
              axis(2, cex.axis=1.5)
              points(x_split,y_score,col=colours[score6$feature], pch=16,cex=(y_scoresize+1))
              if (sum(score6$weighted)>0) points(x_split[score6$weighted],y_score[score6$weighted],col="black",cex=(y_scoresize[score6$weighted]+2))
              if (sum(score6$weightedorig)>0) points(x_split[score6$weightedorig],y_score[score6$weightedorig],col="blue",cex=(y_scoresize[score6$weightedorig]+2))
              if (sum(grepl("TRIM",score6$featend))>0) points(x_split[grepl("TRIM",score6$featend)],y_score[grepl("TRIM",score6$featend)],col="brown",cex=(y_scoresize[grepl("TRIM",score6$featend)]+2))
              if (sum(grepl("TRIM",score6$featend))>0) {
                legend("topright",
                       legend=c(unique(score6$feature),"Trimmed"),
                       col=c(colours[unique(score6$feature)],"brown"),
                       pch =c(rep(16,length(unique(score6$feature))), 1),
                       cex=1.5,
                       inset=c(-.25,0))
              } else {
                legend("topright",
                       legend=c(unique(score6$feature)),
                       col=c(colours[unique(score6$feature)]),
                       pch =c(rep(16,length(unique(score6$feature)))),
                       cex=1.5,
                       inset=c(-.25,0))
              }
            }
          }
          graphics.off()
        }
      }
    }
    return(F)
  }, error = function(e) {
    cat(paste("ERROR:  ",e)); graphics.off(); return(T)
  })
}

TimeOutput(start)













## Plot distances ----------------------------------------------------------

if (doDistPlot) {
  start = Sys.time()
  
  if (overwritedist | !file.exists(paste0(dist_comp_dir,".Rdata"))) {
    cat("samples in common ")
    
    samples0 = foreach(i=distMeta2$path) %dopar% {
      d0 = get(load(distMeta2$path[i]))
      return(rownames(as.matrix(d0)))
    }
    samples1 = Reduce('intersect',samples0)
    
    TimeOutput(start)
    
    cat(", collating distances ")
    d = d0 = foreach(i=distMeta2$path, .combine="rbind") %dopar% {
      d0 = as.matrix(get(load(distMeta2$path[i])))
      ind = match(samples1,rownames(d0))
      d0 = d0[ind,ind]
      d = unlist(as.list(as.dist(d0)))
      return(d)
    }
    rownames(d) = rownames(d0) = distMeta2$path
    save(d0,file=paste0(dist_comp_dir,".Rdata"))
  } else {
    d0 = get(load(paste0(dist_comp_dir,".Rdata")))
  }
  
  d1 = t(apply(d0,1,function(x) exp(log(x,getlv(max(x),100)))/100 ))
  
  TimeOutput(start)
  
  
  cat(", normalizing ")
  if (quantilnorm) {
    colramp = colorRampPalette(c(3,"white",2))(20)
    # x11(); par(mfrow=c(2,1))
    # plot(density(log(d0+1)[,1]),col=colramp[1],lwd=3)
    d = normalize.quantiles(as.matrix(d0))
    # plot(density(log(d+1)[,1]),col=colramp[1],lwd=3)
  }
  TimeOutput(start)
  
  
  cat(", distances ")
  if (overwritedist | !file.exists(paste0(dist_comp_dist_dir,".Rdata"))) {
    did = dist(d)
    did1 = dist(d1)
  } else {
    did0 = get(load(paste0(dist_comp_dist_dir,".Rdata")))
    did = did0$did
    did1 = did0$did1
  }
  TimeOutput(start)
  
  
  
  
  
  
  
  
  
  
  
  
  if (doRtsne) {
    
    cat(", Rtsne ")
    tsne = Rtsne(as.matrix(did),is_distance=T, perplexity = 1)
    tsne1 = Rtsne(as.matrix(did1),is_distance=T, perplexity = 1)
    
    TimeOutput(start)
    
    
    cat(", Plotting Rtsne ")
    dm2 = distMeta2[match(rownames(distMeta2)),]
    dm2[dm2[,"weighted"],"dist"] = "weighted manhattan"
    sc2 = score1[match(fileNames(distMeta2$path),fileNames(rownames(score1))),]
    plotscores = c(ivdscore,ivscore) #plot scores for distance matrix
    
    colour = rainbow(max(length(unique(dm2[dm2,"feature"])),length(unique(dm2[dm2,"layer"]))))
    
    png(filename=paste0(plot_compare_dir,"/dist__",paste(plotscores,collapse="-"),".png"), width=width*length(plotscores), height=height)
    par(mfcol=c(length(plotscores),1), mar=c(10,5,5,5))
    
    #for each score (column)
    for (plotscore in plotscores) {
      silm0 = sc2[,plotscore]
      silm0[is.na(silm0)|silm0==Inf] = min(silm0[!is.na(silm0)& !silm0==Inf])
      
      for (ii in c(1:2)) {
        if (ii==1) {
          x = tsne
          main1 = paste0("size = ",plotscore,"\ndistance measures; quantnormalized")
        } else if (ii==2) {
          x = tsne1
          main1 = paste0("size = ",plotscore,"distance measures; 0-1 normalized")
        }
        
        plot(x$Y,t='n',main=main1, xlim=c(min(x$Y[,1])-.3*(max(x$Y[,1])-min(x$Y[,1])),max(x$Y[,1])))
        points(x$Y[,1],x$Y[,2],col=colour[as.integer(factor(dm2[,"feature"]))],pch=16,cex=exp(silm0/max(silm0)))
        points(x$Y[,1],x$Y[,2],col=colour[as.integer(factor(dm2[,"layer"]))],cex=exp((silm0/max(silm0))+.5))
        if (sum(dm2$weighted)>0) points(x$Y[dm2$weighted,1],x$Y[dm2$weighted,2],col="black",cex=exp((silm0/max(silm0))[dm2$weighted]+1))
        if (sum(dm2$weightedorig)>0) points(x$Y[dm2$weightedorig,1],x$Y[dm2$weightedorig,2],col="blue",cex=exp((silm0/max(silm0))[dm2$weightedorig]+1))
        
        legend("topleft",
               legend=c(paste0("feature (1st ring): ",unique(dm2[,"feature"])),
                        paste0("layer (2nd ring): ",unique(dm2[,"layer"])), "weighted (3rd ring)"),
               col=c(colour[as.integer(factor(unique(dm2[,"feature"])))],
                     colour[as.integer(factor(unique(dm2[,"layer"])))], "black"),
               pch =c(rep(16,length(unique(dm2[,"feature"]))),
                      rep(1,length(unique(dm2[,"layer"]))), 1) )
      }
    }
    graphics.off()
  }
  
  
  
  
  
  
  
  
  
  
  
  
  if(doPCA) {
    
    cat(", PCA-ing")
    
    pc = prcomp(d)
    pc1 = prcomp(d1)
    
    
    cat(", Plotting PCA... ")
    dm2 = distMeta2[match(rownames(distMeta2)),]
    dm2[dm2[,"weighted"],"dist"] = "weighted manhattan"
    sc2 = score1[match(fileNames(distMeta2$path),fileNames(rownames(score1))),]
    
    colour = rainbow(max(length(unique(dm2[dm2,"feature"])),length(unique(dm2[dm2,"layer"]))))
    cp1 = c(15:18,7:14)[1:length(unique(dm2[,"featend"]))]
    cp2 = c(0,1,5,2,6)[length(unique(dm2[,"dist"]))] #empty inside
    
    #for each score (column)
    for (plotscore in plotscores) {
      silm0 = sc2[,plotscore]
      silm0[is.na(silm0)|silm0==Inf] = min(silm0[!is.na(silm0)& !silm0==Inf])
      
      for (ii in c(1:2)) {
        if (ii==1) {
          x = pc
          main1 = paste0("size = ",plotscore,"\ndistance measures; quantnormalized")
        } else if (ii==2) {
          x = pc1
          main1 = paste0("size = ",plotscore,"distance measures; 0-1 normalized")
        }
        
        png(filename=paste0(plot_compare_dir,"/distcomparePCA_",plotscore,"_",ii,".png"), width=600, height=600*(1+plotpc))
        par(mfrow=c(plotpc+1,1))
        plot(pc$sdev^2/sum(pc$sdev^2), xlim = c(0, 15), type = "b", pch = 16, xlab = "principal components", ylab = "variance explained")
        legend("topleft",
               legend=c(paste0("feature (1st ring): ",unique(dm2[,"feature"])),
                        paste0("feattype: ",unique(dm2[,"featend"])),
                        paste0("layer (2nd ring): ",unique(dm2[,"layer"])),
                        paste0("distmeasure: ",unique(dm2[,"dist"])), "weighted (3rd ring)"),
               col=c(colour[as.integer(factor(unique(dm2[,"feature"])))],
                     rep("black",length(unique(dm2[,"featend"]))), 
                     colour[as.integer(factor(unique(dm2[,"layer"])))],
                     rep("black",length(unique(dm2[,"dist"]))), "black"),
               pch =c(rep(16,length(unique(dm2[,"feature"]))),
                      cp1[as.integer(factor(unique(dm2[,"featend"])))],
                      rep(1,length(unique(dm2[,"layer"]))),
                      cp2[as.integer(factor(unique(dm2[,"dist"])))], 1) )
        
        for (i in 1:plotpc) {
          plot(pc$x[,i], pc$x[,i+1], t='n', main = "distance measures\nsize = exp of avg med silhouette index of KO genes", xlab = paste0("PC_",i), ylab = paste0("PC_",i+1),xlim=c(min(pc$x[,i])-50,max(pc$x[,i])))
          points(pc$x[,1],pc$x[,2],col=colour[as.integer(factor(dm2[,"feature"]))],pch=cp1[as.integer(factor(dm2[,"featend"]))],cex=exp(1.5*silm0/max(silm0)))
          points(pc$x[,1],pc$x[,2],col=colour[as.integer(factor(dm2[,"layer"]))],pch=cp2[as.integer(factor(dm2[,"dist"]))],cex=exp((1.5*silm0/max(silm0))+.5))
          points(pc$x[dm2$weighted,1],pc$x[dm2$weighted,2],col="black",cex=exp((2*silm0/max(silm0))[dm2$weighted]+1))
          points(pc$x[dm2$weightedorig,1],pc$x[dm2$weightedorig,2],col="blue",cex=exp((2*silm0/max(silm0))[dm2$weightedorig]+1))
          points(0, 0, pch = 3, cex = 4, lwd = 4)
        }
        graphics.off()
      }
    }
    
    TimeOutput(start)
  }
  
}