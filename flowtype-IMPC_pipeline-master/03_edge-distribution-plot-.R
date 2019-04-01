## Input: matrixChild_prop matrices --> Output: plots of distribution of edge+
# aya43@sfu.ca 20171203

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
phenoMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrix", sep="")
matrixChild_prop_dir = paste(matrix_dir, "Child_prop",sep="")

#Output
plot_dir = paste(result_dir, "/", panelL, "/", centreL, "/plots", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dir[i])) }
plot_edge_dir = paste0(plot_dir,"/edge_distribution"); for(i in 1:length(plot_edge_dir)) { suppressWarnings(dir.create (plot_edge_dir)) }


#Libraries/Functions
source("~/projects/IMPC/code/_funcAlice.R")
libr("foreach")
libr("doMC")
libr("FSA") #set bin widths for hist() function

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)









#Options for script

allorsig = c("","TRIM_") # choose cell population from significant only or all
matrix_type = c("CountAdjPEER","CountAdj")
adjustL = c("none","BH","BY","bonferroni") #"BH" #pvalue adjustment; "none" if don't adjust
target_col = "gene"
controlL = c("[+]_[+]|[+]_Y","[+]_[+]|[+]_Y","WildType","WildType","WildType") #control value in target_col column
overwriteprob = F

binwidth = .005 #edge values from 0-1; specify a histograph bin width

sampleno = 10 # sample sampleno of cells in matrix to plot
plotno = 40 # number of plots to make



start = Sys.time()

for (ci in 1:length(paste0(panelL,centreL))) {
  centre = paste0(panelL," ",centreL)[ci]
  sampleMeta = get(load(sampleMeta_dir[ci]))
  phenoMeta = get(load(phenoMeta_dir[ci]))
  
  control = controlL[ci]
  
  WTfilename = sampleMeta$fileName[grepl(control,sampleMeta$gene)]
  
  for (aos in allorsig) {
    matrix_type0 = matrix_type
    adjustL0 = adjustL
    if (aos=="") {
      matrix_type0 = adjustL0 = ""
    }
    
    for (mcp in matrix_type0) {
      for (adj in adjustL0) {
        prob_dr = paste0(gsub("matrix","prob_matrix",matrixChild_prop_dir[ci]),aos,mcp,adj,".Rdata")
        if (overwriteprob | !file.exists(prob_dr)) {
          if (!file.exists(paste0(matrixChild_prop_dir,aos,mcp,adj,".Rdata"))) next()
          mp = get(load(paste0(matrixChild_prop_dir,aos,mcp,adj,".Rdata"))) #child prop matrix
          
          # loop.ind = loopInd(names(mp),no_cores)
          result = foreach (parenti = 1:length(mp)) %dopar% {
            parent = parent_temp = names(mp)[parenti]
            parent_temp = gsub("[+]","\\[+\\]|",parent_temp)
            parent_temp = gsub("[-]","\\[-\\]|",parent_temp)
            if (grepl("[|]$",parent_temp)) parent_temp = substr(parent_temp,1,nchar(parent_temp)-1)
            child_addedmarkers = gsub(parent_temp,"",colnames(mp[[parenti]]))
            prob_addedmarkers = grepl("[+]",child_addedmarkers)
            
            #extract edges
            mpp0 = mp[[parenti]][,prob_addedmarkers]
            
            #adjust
            if (sum(prob_addedmarkers)==0) return(NULL)
            if (sum(prob_addedmarkers)==1) {
              mpp0 = matrix(mpp0,ncol=1)
              colnames(mpp0) = colnames(mp[[parenti]])[prob_addedmarkers]
              rownames(mpp0) = rownames(mp[[parenti]])
            }
            
            #unify rows
            rowind = match(sampleMeta$fileName,rownames(mpp0))
            mpp = rbind(mpp0,rep(NA,ncol(mpp0)))
            rowind[is.na(rowind)] = nrow(mpp)
            mpp = mpp[rowind,]
            
            #adjut again
            if (sum(prob_addedmarkers)==1) {
              mpp = matrix(mpp,ncol=1)
              colnames(mpp) = colnames(mp[[parenti]])[prob_addedmarkers]
              rownames(mpp) = sampleMeta$fileName
            }
            colnames(mpp) = paste0(parent,"_",colnames(mpp))
            
            return(log(mpp)) # because i exp the original child_prop matrix to make numbers bigger so that distance looks better
          }
          names(result) = names(mp)
          mp_prob0 = remove_null(result)
          
          mp_prob = Reduce("cbind",mp_prob0)
          rownames(mp_prob) = sampleMeta$fileName
          save(mp_prob, file=prob_dr)
        } else {
          mp_prob = get(load(prob_dr)) #child prop matrix
        }
        
        
        
        
        
        ## plot ------------------------------------------
        cellpops0 = sample(1:ncol(mp_prob),min(ncol(mp_prob),sampleno*plotno))
        cellpops = loopInd(cellpops0,plotno)
        for (i in 1:plotno) {
          cellpop = cellpops[[i]]
          plotvars = lapply(cellpop, function(x){
            y = cbind(mp_prob[,x],rep(colnames(mp_prob)[x],nrow(mp_prob)))
            return(y)
          })
          plotvars = Reduce("rbind",plotvars)
          plotvars = plotvars[!apply(plotvars,1,function(x) any(is.na(x))),]
          plotvarsWT = plotvars[rownames(plotvars)%in%WTfilename,]
          plotvarsKO = plotvars[!rownames(plotvars)%in%WTfilename,]
          dens_plot(val=as.numeric(plotvarsWT[,1]),group=plotvarsWT[,2],multipleoverlay="identity",title=paste("distribution of prob_m (",sampleno," cell pops sampled)"),filename=paste0(root,"/",plot_edge_dir[ci],"/edge_distribution_",aos,mcp,adj,gsub("[:]",".",as.character(Sys.time())),"_WT.html"), binwidth=binwidth)
          dens_plot(val=as.numeric(plotvarsKO[,1]),group=plotvarsKO[,2],multipleoverlay="identity",title=paste("distribution of prob_m (",sampleno," cell pops sampled)"),filename=paste0(root,"/",plot_edge_dir[ci],"/edge_distribution_",aos,mcp,adj,gsub("[:]",".",as.character(Sys.time())),"_KO.html"), binwidth=binwidth)
          
          # #alternative density plotting code
          # j0 = unique(plotvarsWT[,2])
          # j = 1
          # val = as.numeric(plotvarsWT[plotvarsWT[,2]==j0[j],1])
          # plot(density(val), ylim=c(0,1500))
          # hist(val, w=binwidth, add=T)
          # colour = gg_color_hue(length(unique(plotvarsWT[,2]))-1)
          # for (j in 2:length(j0)) {
          #   val = as.numeric(plotvarsWT[plotvarsWT[,2]==j0[j],1])
          #   lines(density(as.numeric(plotvarsWT[plotvarsWT[,2]==j0[j],1])),col=colour[j-1])
          #   hist(val, w=binwidth, col=colour[j-1], add=T)
          # }
          
        }

      }
    }
  }
  
}
TimeOutput(start)









