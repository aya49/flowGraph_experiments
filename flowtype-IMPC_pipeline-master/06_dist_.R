## Input: original features --> Output: distance matrices
# aya43@sfu.ca 20161220

#Directory
root = "~/projects/IMPC"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#,"Sanger_MLN","CIPHE","TCP","H")

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir,  "/", panelL, "/", centreL, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir,  "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir,  "/", panelL, "/", centreL, "/matrix", sep="")

#Output
dist_dir = paste(result_dir,  "/", panelL, "/", centreL, "/dist", sep=""); for (i in 1:length(dist_dir)) { suppressWarnings(dir.create (dist_dir[i])) }
dist_type_dir = NULL
phenoChild_dir = paste(result_dir,  "/", panelL, "/", centreL, "/phenoChild.Rdata",sep="")
phenoChild_ind_dir = paste(result_dir,  "/", panelL, "/", centreL, "/phenoChild_ind.Rdata",sep="")
phenoChildpn_dir = paste(result_dir,  "/", panelL, "/", centreL, "/phenoChildpn.Rdata",sep="")
phenoChildpn_ind_dir = paste(result_dir,  "/", panelL, "/", centreL, "/phenoChildpn_ind.Rdata",sep="")

#Libraries/Functions
source("~/projects/IMPC/code/_funcAlice.R")
libr("stringr")
libr("colorspace")
libr("vegan") # libr(proxy)
libr("foreach")
libr("doMC")
libr("lubridate") #if there are date variables
libr("kernlab")

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)








#Options for script
# sampletimes = 5
overwrite = F #overwrite distances?
writecsv = F

cellCountThres = c(200) #insignificant if count under

dis = c("manhattan")#, "euclidean") #distances metrics to use #, "binomial", "cao", "jaccard","mahalanobis", "gower","altGower","morisita","horn", "manhattan", "maximum", "binary", "minkowski", "raup", "mountford"
disnoavg = c("euclidean", "manhattan", "canberra", "bray", "kulczynski", "morisita", "horn", "binomial")  #dis measures that don't average over number of features
disindist = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski") #use dist function when possible, instead of vegdist
disnoneg = c("canberra") #dis measures that can't handle negative values
# disinkernel = c("rbf","poly","vanilla","tanh","laplace","bessel","anova","spline") #kernels

normalize = c("none","cellpop", "layer") # by none (for child matrices only), cell pop, layer

matrix_count = c("CountAdj")



















start = Sys.time()


#for (ci in 4:1) {
for (ci in 1:length(paste0(panelL,centreL))) {
  start1 = Sys.time()
  centre = paste0(panelL," ",centreL)[ci]
  cat("\n",centre,sep="")
  
  #Prepare data
  matrix_type = list.files(path=paste(result_dir,  "/", panelL, "/", centreL, sep="")[ci],pattern=glob2rx("matrix*.Rdata"))
  matrix_type = matrix_type[!matrix_type=="matrixCount.Rdata" & !matrix_type=="matrixProp.Rdata"]
  matrix_type = gsub("matrix|.Rdata","",matrix_type)
  matrix_type = matrix_type[!grepl("leaf",matrix_type)]
  matrix_type = matrix_type[!grepl("KO|Max",matrix_type)]
  matrix_type = matrix_type[!grepl("Count_sample",matrix_type)]
  matrix_type = matrix_type[grepl("CountAdj",matrix_type)]
  
  # matrix_type = c("CountAdj", "Parent_entropy", "Child_entropy", 
  #                 paste0("Child_entropyTRIM_CountAdj_",c(1:sampletimes)), "Child_entropyTRIM_CountAdj", "Child_entropyTRIM_Prop", 
  #                 paste0("Parent_entropyTRIM_CountAdj_",c(1:sampletimes)), "Parent_entropyTRIM_CountAdj", "Parent_entropyTRIM_Prop", 
  #                 paste0("LogFoldTRIM_CountAdj_",c(1:sampletimes)), "LogFoldTRIM_CountAdj", "LogFoldTRIM_Prop", 
  #                 paste0("PvalTRIM_CountAdj_",c(1:sampletimes)), "PvalTRIM_CountAdj", "PvalTRIM_Prop",
  #                 paste0("LogFold_CountAdj_",c(1:sampletimes)), "LogFold_CountAdj", "LogFold_Prop", 
  #                 paste0("Pval_CountAdj_",c(1:sampletimes)), "Pval_CountAdj", "Pval_Prop",
  #                 "Child_pnratio", "Child_prop",
  #                 paste0("Parent_effortTRIM_CountAdj_",c(1:sampletimes)), "Parent_effortTRIM_CountAdj", "Parent_effortTRIM_Prop", 
  #                 paste0("Parent_contribTRIM_CountAdj_",c(1:sampletimes)), "Parent_contribTRIM_CountAdj", "Parent_contribTRIM_Prop", 
  #                 paste0("Child_pnratioTRIM_CountAdj_",c(1:sampletimes)), "Child_pnratioTRIM_CountAdj", "Child_pnratioTRIM_Prop", 
  #                 paste0("Child_propTRIM_CountAdj_",c(1:sampletimes)), "Child_propTRIM_CountAdj", "Child_propTRIM_Prop",
  #                 paste0("Parent_effort_CountAdj_",c(1:sampletimes)), "Parent_effort_CountAdj", "Parent_effort_Prop", 
  #                 paste0("Parent_contrib_CountAdj_",c(1:sampletimes)), "Parent_contrib_CountAdj", "Parent_contrib_Prop")

  
  
  
  
  
  
  m0 = get(load(paste0(matrix_dir[ci], matrix_count,".Rdata")))
  phenoMeta = get(load(phenoMeta_dir[ci]))
  sampleMeta = get(load(sampleMeta_dir[ci]))
  
  k0 = c(max(phenoMeta$phenolevel)) #1,4, # how many layers to consider i.e. k=max(phenolevel) only
  
  
  #convert to csv
  # m_paths = paste0(matrix_dir[ci],matrix_type,".Rdata")
  # foreach(mp = m_paths) %dopar% {
  #   mm = get(load(mp))
  #   if (is.null(dim(mm))) mm = Reduce('cbind',mm)
  #   write.table(as.matrix(mm),sep=",",row.names=F,col.names=F,file=gsub(".Rdata",".csv",mp))
  # }
  # fl = gsub("result/","/home/ayue/projects/IMPC/result/",gsub(".Rdata",".csv",m_paths))
  # write.table(fl[!grepl("Prop",fl)],file=paste0(paste0(result_dir,  "/", panelL, "/", centreL)[ci] ,"/featlist.csv"),sep=",",row.names=F,col.names=F)
  
  fl = read.csv(paste0(paste0(result_dir,  "/", panelL, "/", centreL)[ci] ,"/featlist.csv"),header=F)[,1]
  result = foreach(fli = fl) %dopar% {
    mm = get(load(gsub(".csv",".Rdata",fli)))
    if (is.null(dim(mm))) {
      lg = as.numeric(factor(sampleMeta$gene[match(rownames(mm[[1]]),sampleMeta$fileName)]))
      lge = as.numeric(factor(sampleMeta$gender[match(rownames(mm[[1]]),sampleMeta$fileName)]))
    } else {
      lg = as.numeric(factor(sampleMeta$gene[match(rownames(mm),sampleMeta$fileName)]))
      lge = as.numeric(factor(sampleMeta$gender[match(rownames(mm),sampleMeta$fileName)]))
    }
    write.table(lg,file=gsub(".csv","_Gene.txt",fli),row.names=F,col.names=F,quote=F,sep=",")
    write.table(lge,file=gsub(".csv","_Gender.txt",fli),row.names=F,col.names=F,quote=F,sep=",")
  }

  
  
  
  
  
  
  
  
  
  
  
  
  
  #load different features matrix and calculate distance matrix
  # for (mcp in matrix_type_) {
  a = foreach(mcp=matrix_type) %dopar% {
    tryCatch({
      cat("\n", mcp, " ",sep="")
      start2 = Sys.time()
      
      #start a list of phenotypes included in each distance matrix calculation such that none are repeated
      leavePhenotype = list()
      doneAll = F
      
      #load feature matrix
      mresult = Loadintermatrices(paste0(matrix_dir[ci], mcp,".Rdata"))
      mml0 = mresult$mml
      mmlname = names(mml0)
      pt = mpi = mresult$pt
      gt = mresult$gt
      
      #get to-delete low count phenotype indices; CountAdj should be first one
      for (countThres in cellCountThres) { cat("\ncountThres: ",countThres," >",sep="")
        
        #get to-delete high no of marker phenotypes
        for (k in k0) { cat(" level",k," ",sep="")
          
          #trim feature matrix
          mmlresult = trimMatrices(mml0,m0,pt,gt,phenoMeta,leavePhenotype,doneAll, countThres,k)
          if (is.null(mmlresult)) return(NULL)
          m = mmlresult$mml[[1]]
          pm = mmlresult$pm
          leavePhenotype = mmlresult$leavePhenotype
          doneAll = mmlresult$doneAll
          
          if (is.null(m)) return(NULL)
          if (is.null(dim(m))) {
            a = Reduce('cbind',m); if (all(a==0)) next
          } else { if (all(m==0)) next }
          
          #for every distance type
          a = match(disnoneg,dis)
          loop.ind = 1:length(dis); if (sum(!is.na(a))>0) loop.ind = loop.ind[-a]
          # foreach(i=loop.ind) %dopar% { #for each phenotype
          for (i in loop.ind) {
            cat(", ", length(dis)-i+1, ":", dis[i], " ", sep="")
            
            #assume after splitting dist filename by "_", distance is second element
            dname = paste(dist_dir[ci], "/", mcp, "_", dis[i], "_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_normalize-", sep = "" )
            
            ## calculate distances
            if (is.null(dim(m))) { #edge feature
              
              if ("none"%in%normalize) { #none
                dname0 = paste0(dname, "none.Rdata", sep="")
                if (overwrite | !file.exists(dname0)) {
                  if (dis[i]%in%disindist) { d = dist(a, method=dis[i])
                  } else { d = vegdist(a, method=dis[i]) }
                  
                  save(d, file=dname0)
                  if (writecsv) write.csv(as.matrix(d), file=gsub(".Rdata",".csv",checkm(d,dname0)))
                }
              } 
              if (("cellpop"%in%normalize | "layer"%in%normalize)) {
                if ("cellpop"%in%normalize) { #cell pop
                  dc = matrix(0,nrow=nrow(m[[1]]),ncol=nrow(m[[1]]))
                  dname1 = paste0(dname, "cellpop.Rdata", sep="")
                }
                if ("layer"%in%normalize) { #layer
                  dl = lapply(unique(pm$phenolevel), function(y) return(matrix(0, nrow=nrow(m[[1]]), ncol=nrow(m[[1]]))) )
                  names(dl) = as.character(unique(pm$phenolevel))
                  dname2 = paste0(dname, "layer.Rdata", sep="")
                }
                if (overwrite | !file.exists(dname1) | !file.exists(dname2)) {
                  for (x in 1:length(m)) {
                    if (dis[i]%in%disindist) { dx = as.matrix(dist(m[[x]],method=dis[i]))
                    } else { dx = as.matrix(vegdist(m[[x]],method=dis[i])) }
                    if (dis[i]%in%disindist) dx = dx/ncol(m[[x]])
                    if ("cellpop"%in%normalize) dc = dc + dx
                    if ("layer"%in%normalize) dl[[as.character(pm$phenolevel[x])]] = dl[[as.character(pm$phenolevel[x])]] + dx #1 indexing (layer starts at 0)
                  }
                  if ("cellpop"%in%normalize) {
                    save(dc, file=dname1)
                    if (writecsv) write.csv(as.matrix(dc), file=gsub(".Rdata",".csv",checkm(d,dname1)))
                  }
                  if ("layer"%in%normalize) {
                    dl = lapply(unique(pm$phenolevel), function(x) {
                      a = dl[[as.character(x)]]
                      a[a>0] = a[a>0]/sum(pm$phenolevel==x)
                      return(a)
                    })
                    a = Reduce('+',dl); a[a>0] = a[a>0]/length(dl) #get mean of all matrices
                    dl = as.dist(a)
                    
                    save(dl, file=dname2)
                    if (writecsv) write.csv(as.matrix(dl), file=gsub(".Rdata",".csv",checkm(d,dname2)))
                  }
                }
              }
              
            } else { #node feature
              if ("cellpop"%in%normalize) { #cell pop
                dname1 = paste0(dname, "cellpop.Rdata")
                if (overwrite | !file.exists(dname1)) {
                  if (dis[i]%in%disindist) { d = dist(m, method=dis[i])
                  } else { d = vegdist(m, method=dis[i]) }
                  
                  save(d, file=dname1)
                  if (writecsv) write.csv(as.matrix(d), file=gsub(".Rdata",".csv",checkm(d,dname1)))
                }
              }
              if ("layer"%in%normalize) { #layer
                dname2 = paste0(dname, "layer.Rdata", sep="")
                if (overwrite | !file.exists(dname2)) {
                  layers = unique(pm$phenolevel)
                  #if (sum(layers==0)>0) layers = layers[-which(layers==0)]
                  if (dis[i]%in%disindist) { dd = lapply(layers, function(x) return(as.matrix(dist(m[,which(pm$phenolevel==x)], method=dis[i]))) )
                  } else { dd = lapply(layers, function(x) return(as.matrix(vegdist(m[,which(pm$phenolevel==x)], method=dis[i]))) ) }
                  names(dd) = as.character(layers)
                  
                  #dis measures that don't average over number of features need extra processing
                  if (dis[i]%in%disnoavg) dd = lapply(layers, function(x) {
                    a = dd[[as.character(x)]]
                    a[a>0] = a[a>0]/sum(pm$phenolevel==x)
                    return(a)
                  })
                  a = Reduce('+',dd); a[a>0] = a[a>0]/length(dd) #get mean of all matrices
                  dl = as.dist(a)
                  
                  save(dl,file=dname2)
                  if (writecsv) write.csv(as.matrix(d), file=paste0(checkm(d,dname2),".csv"))
                }
              }
            }
            
          } #dis
        } #layer
      } #countThres
      
      TimeOutput(start2)
    }, error = function(err) { cat(paste("ERROR:  ",err)) })
  }
  TimeOutput(start1)
  
}
TimeOutput(start)















#delete distance matrices with all 0
distmfile = list.files(dist_dir[ci], recursive=F, full.names=T, pattern=".Rdata$")
a = foreach (i = 1:length(distmfile),.combine="c") %dopar% {
  a = F
  d = get(load(distmfile[i]))
  if (any(is.na(as.matrix(d)))) { rm(d); return(T) }
  if (all(as.matrix(d)==as.matrix(d)[1,1])) { rm(d); return(T) }
  rm(d)
  return(a)
}
file.remove(distmfile[a])




