## Input: original features --> Output: bicluster & plots
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
sampleMeta_dir = paste(result_dir,  "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir,  "/", panelL, "/", centreL, "/matrix", sep="")

#Output
plot_dir = paste(result_dir,  "/", panelL, "/", centreL, "/plots", sep="")

#Libraries/Functions
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_func.R")
libr("foreach")
libr("doMC")
libr("stringr")

#Setup Cores
no_cores = 2#detectCores()-3
registerDoMC(no_cores)








#Options for script
plot_size = c(300,600)

matrix_type = c( "Pval_CountAdjBH"                           ,"Pval_CountAdjbonferroni"                   ,"Pval_CountAdjBY"                           ,"Pval_CountAdjPEERBH"                      
                 ,"Pval_CountAdjPEERbonferroni"               ,"Pval_CountAdjPEERBY"                       ,"Pval_CountAdjPEER"                         ,"Pval_CountAdj")
matrix_count = c("CountAdj") #(needed if sample x cell pop matrices are used) count matrix, to get rid of low cell count cells
no_files = c(1:5)
pthresholds = c(.025,.02,.01,.005,.0025,.001)














start = Sys.time()

for (ci in 1:length(paste0(panelL,centreL))) {
  start1 = Sys.time()
  centre = paste0(panelL," ",centreL)[ci]
  cat("\n",centre,sep="")
  
  mc = get(load(paste0(matrix_dir[ci], matrix_count,".Rdata")))
  sampleMeta = get(load(sampleMeta_dir[ci]))
  
  fe = foreach(mcp=matrix_type) %dopar% {
    tryCatch({
      cat("\n", mcp, " ",sep="")
      start2 = Sys.time()
      
      ## upload and prep feature matrix
      m0 = get(load(paste0(matrix_dir[ci], mcp,".Rdata")))
      m0cn = colnames(m0)
      m0rn = rownames(m0)
      
      #get feature layers if feature names represent cell types
      colsplitlen = cell_type_layers(m0cn)
      k0 = c(1:max(colsplitlen))
      
      #get to-delete low count phenotype indices; CountAdj should be first one
      for (countThres in cellCountThres) { cat("\ncountThres: ",countThres," >",sep="")
        colIndexC = rep(F,ncol(m0))
        if(!is.null(colsplitlen)) {
          mcordercol = match(colnames(m0),colnames(mc))
          mcordercol = mcordercol[!is.na(mcordercol)]
          mc0 = mc[,mcordercol]
          colIndexC = apply(mc0, 2, function(x) any(x>=countThres))
        } 
        
        dname = paste(plot_dir[ci], "/pval-amount_", mcp, "_countThres-", countThres,"_KOsamplecountover-",good_sample, ".png", sep = "")

        
        
        layer_by_p = list()
        layer_by_p[["col"]] = NULL #layer (x no_files) x p value thresholds
        layer_by_p[["row"]] = NULL #layer x p value thresholds
        layer_by_p[["p"]] = NULL #layer x p value thresholds
        
        #get to-delete high no of marker phenotypes
        for (k in k0) { cat(" level",k," ",sep="")
          colIndexL = colsplitlen <= k
          
          #trim matrix based on cell count and layer
          m = m0[,colIndexC & colIndexL]
          colIndex0 = apply(m, 2, function(x) sum(x!=0)>0 )
          rowIndex0 = apply(m, 1, function(x) sum(x!=0)>0 )
          m = m[rowIndex0,colIndex0]
          
          #check matrix; if matrix is a list, merge; if matrix is empty, skip
          if (is.null(m)) next
          if (all(m==0)) next 
          
          colIndexP = sapply(-log(pthresholds), function(p) {
            sig_per_col = apply(m, 2, function(x) sum(abs(x)>p) )
            nocol_with_signofiles = sapply(no_files, function(no_file) sum(sig_per_col>=no_file))
            return(nocol_with_signofiles)
          })
          rowIndexP = NULL
          totalP = NULL
          for (pthres in pthresholds) {
            p = -log(pthres)
            rowIndexP = append(rowIndexP,sum(apply(m, 1, function(x) sum(abs(x)>p) )>0))
            totalP = append(totalP,sum(abs(m)>p))
          }
          
          colnames(colIndexP) = names(rowIndexP) = names(totalP) = pthresholds
          layer_by_p[["col"]] = rbind(layer_by_p[["col"]],colIndexP)
          layer_by_p[["row"]] = rbind(layer_by_p[["row"]],rowIndexP)
          layer_by_p[["p"]] = rbind(layer_by_p[["p"]],totalP)
        } #layer
        
        
        png(dname, height=plot_size[1]*(length(no_files)+2), width=plot_size[2])
        par(mfrow=c(length(no_files)+2,1), mar=c(5,10,5,10))
        
        
        #col
        colour = rainbow(length(k0))
        for (fi in 1:length(no_files)) {
          values = layer_by_p[["col"]][((length(k0)-1)*length(no_files))+fi,]
          plot(values, type="l", col=colour[length(k0)], xaxt="n", xlab="p value threshold", ylab=paste0("# of features in ", no_files[fi], " or more samples"), ylim=c(min(layer_by_p[["col"]]),max(layer_by_p[["col"]])), cex.lab=1.5, cex.axis=1.5)
          text(1:length(pthresholds),values,values, col=colour[length(k0)], pos=3, cex=1.5)
          axis(side=1, at=c(1:length(pthresholds)), labels=pthresholds, cex.axis=1.5)
          for (ki in (length(k0)-1):1) {
            values = layer_by_p[["col"]][((ki-1)*length(no_files))+fi,]
            lines(values, col=colour[ki])
            text(1:length(pthresholds), values, values, col=colour[ki], pos=3, cex=1.5)
          }
          legend("topright",legend=k0, fill=colour, title="up till layer:", cex=1.5,inset=c(-.21,0))
        }
        
        #row
        values = layer_by_p[["row"]][length(k0),]
        plot(values, type="l", col=colour[length(k0)], xaxt="n", xlab="p value threshold", ylab=paste0("# of samples with at least one significant feature"), ylim=c(min(layer_by_p[["row"]]),max(layer_by_p[["row"]])), cex.lab=1.5, cex.axis=1.5)
        text(c(1:length(pthresholds)),values,values, col=colour[length(k0)], pos=3, cex=1.5)
        axis(side=1, at=c(1:length(pthresholds)), labels=pthresholds, cex.axis=1.5)
        for (ki in (length(k0)-1):1) {
          values = layer_by_p[["row"]][ki,]
          lines(values, col=colour[ki])
          text(1:length(pthresholds),values,values, col=colour[ki], pos=3, cex=1.5)
        }
        legend("topright",legend=k0, fill=colour, title="up till layer:", cex=1.5)
        
        #total
        values = layer_by_p[["p"]][length(k0),]
        plot(values, type="l", col=colour[length(k0)], xaxt="n", xlab="p value threshold", ylab=paste0("# of matrix cells with significant value"), ylim=c(min(layer_by_p[["p"]]),max(layer_by_p[["p"]])), cex=2, cex.lab=1.5, cex.axis=1.5)
        text(1:length(pthresholds),values,values, col=colour[length(k0)], pos=3, cex=1.5)
        axis(side=1, at=c(1:length(pthresholds)), labels=pthresholds, cex.axis=1.5)
        for (ki in (length(k0)-1):1) {
          values = layer_by_p[["p"]][ki,]
          lines(values, col=colour[ki])
          text(1:length(pthresholds),values, values, col=colour[ki], pos=3, cex=1.5)
        }
        legend("topright",legend=k0, fill=colour, title="up till layer:", cex=1.5)
        
        graphics.off()
      } #countThres
      TimeOutput(start2)
    }, error = function(err) { cat(paste("ERROR:  ",err)); return(T) })
    return(F)
  }
  TimeOutput(start1)
}
TimeOutput(start)




