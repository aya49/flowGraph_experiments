# plotcluster: combine different features and discriminant plot
# aya43@sfu.ca 20161220

#root = "~/projects/IMPC"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#,"Sanger_MLN","CIPHE","TCP","H")

#Input
phenoMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrix", sep="")
# mcp_types = c("/original", "/child_pn","/child_prop", "/trim") 
# matrix_type1 = c("CountAdj") #, "Prop", "Count",  #countAdj comes first, to set which columns are deleted based on cellcountthres
# matrix_type2 = c("Child_prop")
# matrix_type3 = c("Child_pnratio", "Child_entropy")
# matrix_type4 = c("Child_pnratio", "Child_entropy")
matrix_count = c("CountAdj")
matrix_type = c("Child_entropyTRIM", "Parent_entropyTRIM", "LogFoldTRIM_CountAdj", "PvalTRIM_CountAdj", "Parent_effortTRIM", "Parent_contribTRIM", "Child_pnratioTRIM", "Child_propTRIM")#, "LogFold_CountAdj", "CountAdj", "Pval_CountAdj", "Parent_entropy", "Child_entropy", "Parent_effort", "Parent_contrib", "Child_pnratio", "Child_prop") 
phenoChild_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChild.Rdata",sep="")
phenoChild_ind_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChild_ind.Rdata",sep="")
phenoChildpn_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChildpn.Rdata",sep="")
phenoChildpn_ind_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChildpn_ind.Rdata",sep="")

cellCountThres = c(200,200,500,500,500) #insignificant if count under
pclustermethods = c("bc","vbc","adc","awc","arc","anc") #"dc","mvdc","nc","wnc", 
pclustermethods2 = c("bc","vbc","adc","awc","arc","anc") #assume two clusters

interested = c("gene", "gender", "colony", "strain", "fur", "specimen", "sample") #sampleMeta columns to plot


#Output
plot_dir = paste(result_dir, "/", panelL, "/", centreL, "/plots", sep="")
pc_dir = paste(plot_dir, "/plotcluster", sep=""); for(i in 1:length(pc_dir)) { suppressWarnings(dir.create(pc_dir[i])) }

source("~/projects/IMPC/code/_funcAlice.R")
libr("stringr")
libr("colorspace")
libr("fpc") # libr(proxy)
libr("foreach")
libr("doMC")
libr("lubridate") #if there are date variables
libr("Matrix")


no_cores = 2#detectCores()-1
registerDoMC(no_cores)




width = 700; height = 700 # of each plot in png
rowplot = 3
maxplot = rowplot*9

target_col="gene"


start = Sys.time()

for (ci in 1:length(paste0(panelL,centreL))) {
  start1 = Sys.time()
  centre = paste0(panelL," ",centreL)[ci]
  cat("\n",centre,sep="")
  control = c("+_+","+_Y"); if (!grepl("Sanger",centre)) control = "WildType"
  
  phenoMeta = get(load(phenoMeta_dir[ci]))
  sampleMeta = get(load(sampleMeta_dir[ci]))
  #k0 = c(1,3,max(phenoMeta$phenolevel)) # how many markers to consider i.e. k=max(phenolevel) only
  k0 = c(max(phenoMeta$phenolevel)) # how many markers to consider i.e. k=max(phenolevel) only
  m0 = get(load(paste0(matrix_dir[ci], matrix_count,".Rdata")))
  interestedCols = which(colnames(sampleMeta)%in%interested)
  
  for (mcp in matrix_type) {
    cat("\n", mcp, " ",sep="")
    start2 = Sys.time()
    leavePhenotype = list()
    doneAll = F
    if (!file.exists(paste0(matrix_dir[ci], mcp,".Rdata"))) {cat("doesn't exist"); next}
    mm = get(load(paste0(matrix_dir[ci], mcp,".Rdata")))
    if (!is.null(dim(mm))) { mpi = colnames(mm)
    } else { mpi = names(mm) }
    
    #get to-delete low count phenotype indices; CountAdj should be first one
    for (countThres in cellCountThres[ci]) {
      cat("\ncountThres: ",countThres," > ",sep="")
      lowCountpt = colnames(m0)[apply(m0,2,function(x) all(x<=countThres))]
      if (!length(lowCountpt)>0) lowCountpt = c()
      
      #get to-delete high no of marker phenotypes
      for (k in k0) {
        cat("level",k," ",sep="")
        highLevelpt = phenoMeta[phenoMeta$phenolevel>k,"phenotype"]
        if (!length(highLevelpt)>0) highLevelpt = c()
        
        
        # cat("length of countthres: ", length(lowCountInd))
        # cat("\nlength of highlevelind: ", length(highLevelInd))
        ## Load & fix cell count/countAdj/proportion matrix -----------------------------------------------------
        dpi0 = union(lowCountpt,highLevelpt)
        # cat("\nlength of dpi: ", length(dpi0), "\n")
        
        m = mm
        
        lpi = setdiff(mpi,dpi0) #phenotypes to leave in analysis
        
        #check if indices already calculated for on this matrix, if yes, skip, else trim and do
        if (Position(function(x) identical(x, lpi), leavePhenotype, nomatch = 0) > 0) {cat("; child skipped", sep=""); next}
        leavePhenotype[[paste0(k="k.",k,"_countThres.",countThres)]] = lpi
        
        if (!doneAll & length(lpi)==length(mpi)) doneAll = T
        
        if (!is.null(dim(mm))) { m = m[,match(lpi,mpi)]
        } else {
          ml = m[match(lpi,mpi)]
          m = foreach(k=1:length(m),.combine="cbind") %dopar% { return(m[[k]]) }
        }
        
        pm = phenoMeta[match(lpi,phenoMeta$phenotype),]
        sm = sampleMeta[match(rownames(m),sampleMeta$fileName),]
        
        foreach(pcm=pclustermethods) %dopar% {
          for (coli in interestedCols) {
            t = as.integer(factor(sm[,coli]))
            tname = sort(unique(sm[,coli]))
            
            title = paste0(mcp, "_", pcm, "_col-",colnames(sm)[coli],"_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres,"_-_",colnames(sm)[coli])
            pcname = paste0(pc_dir[ci],"/plotcluster_", title, ".png")
            
            if (colnames(sm)[coli]==target_col) {
              t[t%in%which(tname%in%control)] = 0
              tname = c(control[control%in%tname],tname[!tname%in%control])
            } 
            if (pcm%in%pclustermethods2) {
              
              if (length(unique(t))>2) {
                numplots0 = length(unique(t)) #actual number of plots
                numpng = 1
                if (min(t)==0) numplots0 = numplots0-1 #don't make one for control if it's the control column we're looking at
                if (numplots0>maxplot) { #set up multiple png's, too many plots on one png may not work...
                  numpng = ceiling(numplots0/maxplot)
                }
                
                for (j in 1:numpng) {
                  j0 = j
                  numplots = maxplot
                  if (numpng==1) { numplots = numplots0; j0=0
                  } else if (j==numpng) { numplots = numplots0-((j-1)*maxplot) }
                  colplot = ceiling(numplots/rowplot) # all genes + WT only, >3 samples
                  
                  png(file=pcname , width=width*colplot, height=height*rowplot)
                  par(mar=(c(5,5,5,5) + 0.1), mfrow=c(rowplot,colplot))
                  for (ploti in unique(sort(t))[(((j-1)*maxplot)+1):min(maxplot*j,numplots0)]) {
                    if (ploti==0) next
                    if (min(t)==0) {
                      plotcluster(m[t%in%c(0,ploti),],t[t%in%c(0,ploti)]==ploti, main=paste0(title,"\n(1) ",tname[ploti], "[",sum(t%in%ploti)," samples], control"))
                    } else {
                      plotcluster(m,t==ploti, main=paste0(title,"\n(1) ",tname[ploti], ", other"))
                    }
                  }
                }
              }
            } else {
              png(file=pcname , width=width, height=height*2)
              par(mfrow=c(2,1))
              plotcluster(m,t,main=paste0(title,"\n"))
              plot.new()
              legend("topleft",legend=tname)
            }
            graphics.off()
          }
          # face <- rFace(300,dMoNo=2,dNoEy=0)
          # grface <- as.integer(attr(face,"grouping"))
          # plotcluster(face,grface==2)
          # plotcluster(face,grface, method=pclustermethods[1])
          # plotcluster(face,grface, method=pclustermethods[2])
          # plotcluster(face,grface, method=pclustermethods[3])
          # plotcluster(face,grface, method=pclustermethods[4])
        }
      }
    }
    TimeOutput(start2)
  }
  TimeOutput(start1)
}

TimeOutput(start)


