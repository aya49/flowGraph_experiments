# aya43@sfu.ca 20170131
# creates list of children for each non-leaf node and creates children/parent list[[matrices]] (-/+ are only for phenotypes where both -,+ data exists)

root = "~/projects/IMPC"
result_dir = "result"
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#, "Sanger_MLN","CIPHE","TCP","H")

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoMeta", sep="")
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")

matrix_weight_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixMaxCountAdj",sep="")
matrixCountAdj_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixCountAdj.Rdata", sep="") #the entropy and child proportions are based on original count matrix
matrixProp_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixProp.Rdata", sep="") #the entropy and child proportions are based on original count matrix
matrixLogFold_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixLogFold",sep="")
matrixKOCountAdj_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixKOCountAdj",sep="")
matrixPvalTRIM_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixPvalTRIM",sep="")

phenoChild_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChild.Rdata",sep="")
phenoChild_ind_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChild_ind.Rdata",sep="")
phenoChildpn_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChildpn.Rdata",sep="")
phenoChildpn_ind_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoChildpn_ind.Rdata",sep="")
phenoParent_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoParent.Rdata",sep="")
phenoParent_ind_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoParent_ind.Rdata",sep="")
phenoParentpn_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoParentpn.Rdata",sep="")
phenoParentpn_ind_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoParentpn_ind.Rdata",sep="")

#Output
matrixChild_pnratio_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixChild_pnratio",sep="")
matrixChild_prop_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixChild_prop",sep="")
matrixChild_entropy_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixChild_entropy",sep="")
matrixParent_entropy_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixParent_entropy",sep="")
matrixParent_contrib_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixParent_contrib",sep="")
matrixParent_effort_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrixParent_effort",sep="")
matrixFreqp_dir = paste(result_dir, "/matrixFreqp",sep="")

#Libraries/Functions
libr(Matrix)
libr(arules)
libr(stringr)
libr(entropy)
libr(foreach)
libr(doMC)
source("~/projects/IMPC/code/_funcAlice.R")

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)






#Options for script
pheno_type = c("CountAdj")#, "Prop")
conf = .5 #only get frequent itemsets above confidence

trimOnly = F







start = Sys.time()

for (ci in length(paste0(panelL,centreL)):1) {
  centre = paste0(panelL," ",centreL)[ci]
  cat("\n",centre,sep="")
  
  start1 = Sys.time()
  

  #Prepare data
  m = get(load(matrixCountAdj_dir[ci]))
  phenoMeta = get(load(paste0(phenoMeta_dir[ci],".Rdata")))
  sampleMeta = get(load(sampleMeta_dir[ci]))
  
  phenoChild = get(load(phenoChild_dir[ci]))
  phenoChild_ind = get(load(phenoChild_ind_dir[ci]))
  phenoChildpn = get(load(phenoChildpn_dir[ci]))
  phenoChildpn_ind = get(load(phenoChildpn_ind_dir[ci]))
  
  mpt = list()
  mw = list()
  mlf = list()
  mkc = list()
  orig = list()
  for (pcpi in 1:length(pheno_type)) {
    pcp = pheno_type[pcpi]
    mpt[[pcpi]] = get(load(paste0(matrixPvalTRIM_dir[ci], "_",pcp,".Rdata")))
    mw[[pcpi]] = get(load(paste0(matrix_weight_dir[ci], "_",pcp,".Rdata")))
    mlf[[pcpi]] = get(load(paste0(matrixLogFold_dir[ci], "_",pcp,".Rdata")))
    mkc[[pcpi]] = get(load(paste0(matrixKOCountAdj_dir[ci], "_",pcp,".Rdata")))
    if (!is.numeric(as.integer(substrRight(pcp,1)))) orig[[pcpi]] = get(load(paste0(result_dir[ci],"/matrix",pcp,".Rdata")))
  }
  phenoParent = get(load(phenoParent_dir[ci]))
  phenoParent_ind = get(load(phenoParent_ind_dir[ci]))
  phenoParentpn = get(load(phenoParentpn_dir[ci]))
  phenoParentpn_ind = get(load(phenoParentpn_ind_dir[ci]))
  
  
  
  
  
  
  
  
  
  
  ## CountAdj trim ---------------------------------------------------------------
  cat("trim count matrix ")
  
  #trim based on pvalue trim
  for (pcpi in 1:length(pheno_type)) {
    pcp = pheno_type[pcpi]
    if (is.numeric(as.integer(substrRight(pcp,1)))) next
    
    cat(" ",pheno_type[pcpi],sep="")
    matrix_orig = orig[[pcpi]]
    lpi = intersect(colnames(matrix_orig), colnames(mpt[[pcpi]]))
    
    matrix_origTRIM = matrix_orig[!is.na(match(rownames(matrix_orig),rownames(mpt[[pcpi]]))),match(lpi,colnames(matrix_orig))]
    matrix_origTRIM[as.matrix(mpt[[pcpi]][,!is.na(match(colnames(mpt[[pcpi]]),colnames(matrix_origTRIM)))])==0] = 0
    matrix_origTRIM = Matrix(matrix_origTRIM, sparse=T)
    
    save(matrix_origTRIM, file=paste0(result_dir[ci],"/matrix",pcp,"TRIM_",pcp,".Rdata"))
    write.csv(as.matrix(matrix_origTRIM), file=paste0(result_dir[ci],"/matrix",pcp,"TRIM_",pcp,".csv"))
  }
  
  
  
  
  
  
  
  
  ## freqpat ----------------------------------------------------------------------
  cat(", frequent pattern")
  
  for (pcpi in 1:length(pheno_type)) {
    pcp = pheno_type[pcpi]
    
    #prepare pvalue and weight matrices
    mpt0 = mpt[[pcpi]]
    mw0 = mw[[pcpi]]
    mw0 = mw0[match(rownames(mpt0),rownames(mw0)),match(colnames(mpt0),colnames(mw0))]
    mw0[which(mpt0==0)] = 0
    mwl = mw0/max(mw0)#exp(log(mw0,getlv(max(mw0),100)))/100
    
    #sort markers (all)
    mwlcn0 = strsplit( gsub("([+-])","\\1~",colnames(mw[[pcpi]])), "~" )
    mwlcn0 = lapply(mwlcn0,function(x) sort(x))
    names(mwlcn0) = sapply(mwlcn0, function(x) paste0(x,collapse=""))
    
    #sort markers (matrix weights only)
    mwlcn = strsplit( gsub("([+-])","\\1~",colnames(mw0)), "~" )
    mwlcn = lapply(mwlcn,function(x) sort(x))
    names(mwlcn) = sapply(mwlcn, function(x) paste0(x,collapse=""))
    
    #make transactions for each sample
    trl = lapply(1:nrow(mwl),function(i) as(mwlcn[which(mwl[i,]>0)],"transactions") )
    
    #no matter +/- change; score of a marker set is confidence (~<1) + lift (~1)
    matrixFreqp = foreach(i = 1:nrow(mwl), .combine="rbind") %dopar% {
      tryCatch({
        #arules
        tr = trl[[i]]
        weight = mwl[i,mpt0[i,]!=0]
        transactionInfo(tr) = as.data.frame(weight)
        is = weclat(tr, parameter=list(support=min(weight)), control=list(verbose=F))
        tryCatch({ ru = ruleInduction(is,confidence=conf)
        }, error = function(err) { ru = ruleInduction(is,tr,confidence=conf) })
        ru1 = as.data.frame(inspect(ru))[,-2]
        if (!nrow(ru1)>0) return(rep(0,length(mwlcn0)))
        
        #sort (rule subset) markers
        ru1set = paste0(gsub("[{,}]","",ru1[,1]),gsub("[{,}]","",ru1[,2]))
        ru1set = strsplit( gsub("([+-])","\\1~",ru1set), "~" )
        ru1set = lapply(ru1set,function(x) sort(x))
        
        #no of duplicates & score for each (rule subset) markers
        ru1dupno = rep(1,length(ru1set))
        ru1score = rep(0,length(ru1set))
        names(ru1dupno) = names(ru1score) = sapply(ru1set, function(x) paste0(x,collapse=""))
        
        #avg & delete score of duplicated marker sets
        ru1setdup = duplicated(ru1set)
        ru1score[!ru1setdup] = rowSums(ru1[which(!ru1setdup),c(3,4)])
        while (sum(ru1setdup)>0) {
          ru1setdup1 = sapply(names(ru1score), function(x) identical(names(ru1score)[[which(ru1setdup)[1]]],x))
          ru1dupno[!ru1setdup & ru1setdup1] = sum(ru1setdup1)
          ru1score[!ru1setdup & ru1setdup1] = sum(ru1[ru1setdup1,c(3,4)])/sum(ru1setdup1)
          ru1setdup[ru1setdup1] = F
        }
        ru1set = ru1set[ru1score>0] 
        ru1dupno = ru1dupno[ru1score>0] 
        ru1score = ru1score[ru1score>0] 
        
        mi = rep(0,length(mwlcn0))
        ind = match(names(mwlcn0),names(ru1score))
        mi[!is.na(ind)] = exp(ru1score[ind[!is.na(ind)]])
        return(mi)
      }, error = function(err) { return(rep(0,length(mwlcn0)))})
    }
    rownames(matrixFreqp) = rownames(mwl)
    colnames(matrixFreqp) = colnames(mw[[pcpi]])
    matrixFreqp = matrixFreqp[apply(matrixFreqp,1,function(x) any(x>0)),apply(matrixFreqp,2,function(x) any(x>0))]
    
    save(matrixFreqp, file=paste0(matrixFreqp_dir[ci],conf,"_orig_",pcp,".Rdata"))
    write.csv(matrixFreqp, file=paste0(matrixFreqp_dir[ci],conf,"_orig_",pcp,".csv"))
  }
  
  TimeOutput(start)
  
  
  
  
  
  
  
  
  
  ## Childentropy ---------------------------------------------------------------
  cat(", childentropy")
  
  if (!trimOnly) {
    matrixChild_entropy = foreach(i = 1:length(phenoChildpn), .combine='cbind') %dopar% { #for each phenotype
      en = rep(0,nrow(m))
      
      parent = m[,phenoChildpn_ind[i]]
      children = m[,unlist(phenoChildpn[[i]])]
      if (length(phenoChildpn[[i]])==1) children = matrix(children,ncol=1)
      
      no_child = length(phenoChildpn[[i]][[1]])
      non0parents = which(parent>0)
      en[non0parents] = sapply(non0parents, function(x) return(entropy(children[x,]/parent[x])/no_child)) #average entropy over # of markers added
      names(en) = rownames(m)
      
      return(entropy=en) #ratio = +child_prop / -child_prop
    }
    colnames(matrixChild_entropy) = phenoMeta$phenotype[phenoChildpn_ind]
    
    save(matrixChild_entropy, file=paste0(matrixChild_entropy_dir[ci],".Rdata"))
    write.csv(matrixChild_entropy, file=paste0(matrixChild_entropy_dir[ci],".csv"))
  } else {
    load(paste0(matrixChild_entropy_dir[ci],".Rdata"))
    colnames(matrixChild_entropy) = phenoMeta$phenotype[phenoChildpn_ind]
    
    save(matrixChild_entropy, file=paste0(matrixChild_entropy_dir[ci],".Rdata"))
    write.csv(matrixChild_entropy, file=paste0(matrixChild_entropy_dir[ci],".csv"))
  }
  
  #trim based on pvalue trim
  for (pcpi in 1:length(pheno_type)) { 
    pcp = pheno_type[pcpi]
    cat(" ",pheno_type[pcpi],sep="")
    
    lpi = intersect(colnames(matrixChild_entropy), colnames(mpt[[pcpi]]))
    matrixChild_entropyTRIM = matrixChild_entropy[!is.na(match(rownames(matrixChild_entropy),rownames(mpt[[pcpi]]))),match(lpi,colnames(matrixChild_entropy))]
    matrixChild_entropyTRIM[as.matrix(mpt[[pcpi]][,!is.na(match(colnames(mpt[[pcpi]]),colnames(matrixChild_entropyTRIM)))])==0] = 0
    matrixChild_entropyTRIM = Matrix(matrixChild_entropyTRIM, sparse=T)
    
    save(matrixChild_entropyTRIM, file=paste0(matrixChild_entropy_dir[ci],"TRIM_",pcp,".Rdata"))
    write.csv(as.matrix(matrixChild_entropyTRIM), file=paste0(matrixChild_entropy_dir[ci],"TRIM_",pcp,".csv"))
  }
  
  TimeOutput(start)
  
  
  
  
  
  
  
  ## Parententropy ---------------------------------------------------------------
  
  
  cat(", parententropy")
  if (!trimOnly) {
    matrixParent_entropy = foreach(i = 1:length(phenoParentpn), .combine='cbind') %dopar% { #for each phenotype
      en = rep(0,nrow(m))
      
      childr = m[,phenoParentpn_ind[i]]
      parent = m[,phenoParentpn[[i]]]
      if (length(phenoParentpn[[i]])==1) parent = matrix(parent,ncol=1)
      
      no_parent = length(phenoParentpn[[i]])
      non0childrs = which(childr>0)
      en[non0childrs] = sapply(non0childrs, function(x) return(entropy(childr[x]/parent[x,])/no_parent)) #average entropy over # of markers added
      names(en) = rownames(m)
      return(en) #ratio = +child_prop / -child_prop
    }
    colnames(matrixParent_entropy) = phenoMeta$phenotype[phenoParentpn_ind]
    
    save(matrixParent_entropy, file=paste0(matrixParent_entropy_dir[ci],".Rdata"))
    write.csv(matrixParent_entropy, file=paste0(matrixParent_entropy_dir[ci],".csv"))
  } else {
    load(paste0(matrixParent_entropy_dir[ci],".Rdata"))
    colnames(matrixParent_entropy) = phenoMeta$phenotype[phenoParentpn_ind]
    
    save(matrixParent_entropy, file=paste0(matrixParent_entropy_dir[ci],".Rdata"))
    write.csv(matrixParent_entropy, file=paste0(matrixParent_entropy_dir[ci],".csv"))
  }
  
  #trim based on pvalue trim
  for (pcpi in 1:length(pheno_type)) { 
    pcp = pheno_type[pcpi]
    
    lpi = intersect(colnames(matrixParent_entropy), colnames(mpt[[pcpi]]))
    matrixParent_entropyTRIM = matrixParent_entropy[!is.na(match(rownames(matrixParent_entropy),rownames(mpt[[pcpi]]))), match(lpi,colnames(matrixParent_entropy))]
    matrixParent_entropyTRIM[as.matrix(mpt[[pcpi]][,!is.na(match(colnames(mpt[[pcpi]]),colnames(matrixParent_entropyTRIM)))])==0] = 0
    
    save(matrixParent_entropyTRIM, file=paste0(matrixParent_entropy_dir[ci],"TRIM_",pcp,".Rdata"))
    write.csv(as.matrix(matrixParent_entropyTRIM), file=paste0(matrixParent_entropy_dir[ci],"TRIM_",pcp,".csv"))
  }
  
  TimeOutput(start)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ## Childprop ---------------------------------------------------------------
  cat("; childprop")
  
  if (!trimOnly) {
    childprop = foreach(i = 1:length(phenoChild)) %dopar% { #for each phenotype
      #for (i in 1:length(phenoChild)) {
      parent = m[,phenoChild_ind[i]]
      
      # Child proportion matrix
      children = m[,unlist(phenoChild[[i]])]
      children[which(children<1)] = 1
      parent0 = parent
      parent0[which(parent<1)] = 1
      childprop = exp(children/parent0)
      rownames(childprop) = rownames(m)
      return(childprop)
    }
    names(childprop) = phenoMeta$phenotype[phenoChild_ind]
    
    save(childprop, file=paste0(matrixChild_prop_dir[ci],".Rdata"))
    childpropREDUCE = foreach(k=1:length(childprop),.combine="cbind") %dopar% { return(childprop[[k]]) }
    write.csv(childpropREDUCE, file=paste0(matrixChild_prop_dir[ci],".csv"))
  } else {
    load(paste0(matrixChild_prop_dir[ci],".Rdata"))
    names(childprop) = phenoMeta$phenotype[phenoChild_ind]
    
    save(childprop, file=paste0(matrixChild_prop_dir[ci],".Rdata"))
    childpropREDUCE = foreach(k=1:length(childprop),.combine="cbind") %dopar% { return(childprop[[k]]) }
    write.csv(childpropREDUCE, file=paste0(matrixChild_prop_dir[ci],".csv"))
  }
  
  #trim based on pvalue (phenotype)
  for (pcpi in 1:length(pheno_type)) { 
    pcp = pheno_type[pcpi]
    cat(" ",pheno_type[pcpi],sep="")
    childpropTRIM = trimlistmatrix(mpt[[pcpi]],refdel=0, childprop,targetdel=0, ncore=no_cores)
    
    save(childpropTRIM, file=paste0(matrixChild_prop_dir[ci],"TRIM_",pcp,".Rdata"))
    childpropTRIMREDUCE = foreach(k=1:length(childpropTRIM),.combine="cbind") %dopar% { return(childpropTRIM[[k]]) }
    write.csv(as.matrix(childpropTRIMREDUCE), file=paste0(matrixChild_prop_dir[ci],"TRIM_",pcp,".csv"))
  }
  TimeOutput(start)
  
  
  
  
  
  
  
  
  
  
  
  ## Childpnratio ---------------------------------------------------------------
  cat(", childpnratio")
  
  if (!trimOnly) {
    pnratio = foreach(i = 1:length(phenoChildpn)) %dopar% { #for each phenotype
      parent = m[,phenoChildpn_ind[i]]
      
      # P/N ratio matrix
      pnratio = NULL
      
      pos = m[,phenoChildpn_ind[i]]-m[,phenoChildpn[[i]][[1]]]
      neg = m[,phenoChildpn[[i]][[1]]]
      pnratio = pos/neg ##get rid of 0, Inf
      pnratio[pos==0 & neg==0] = 1
      pnratio[pos==0 & neg!=0] = 1/neg[pos==0 & neg!=0]
      pnratio[pos!=0 & neg==0] = pos[pos!=0 & neg==0]
      pnratio = matrix(log(pnratio),nrow=nrow(m))
      
      if (sum(parent==0)>0) {
        pnratio[which(parent==0),] = rep(0,length(phenoChildpn[[i]][[1]]))
      } 
      rownames(pnratio) = rownames(m)
      colnames(pnratio) = colnames(pos)
      
      return(pnratio) #ratio = +child_prop / -child_prop
    }
    names(pnratio) = phenoMeta$phenotype[phenoChildpn_ind]
    
    save(pnratio, file=paste0(matrixChild_pnratio_dir[ci],".Rdata"))
    pnratioREDUCE = foreach(k=1:length(pnratio),.combine="cbind") %dopar% { return(pnratio[[k]]) }
    write.csv(pnratioREDUCE, file=paste0(matrixChild_pnratio_dir[ci],".csv"))
  } else {
    load(paste0(matrixChild_pnratio_dir[ci],".Rdata"))
    names(pnratio) = phenoMeta$phenotype[phenoChildpn_ind]
    
    save(pnratio, file=paste0(matrixChild_pnratio_dir[ci],".Rdata"))
    pnratioREDUCE = foreach(k=1:length(pnratio),.combine="cbind") %dopar% { return(pnratio[[k]]) }
    write.csv(pnratioREDUCE, file=paste0(matrixChild_pnratio_dir[ci],".csv"))
  }
  
  #trim based on pvalue trim
  for (pcpi in 1:length(pheno_type)) { 
    pcp = pheno_type[pcpi]
    cat(" ",pheno_type[pcpi],sep="")
    pnratioTRIM = trimlistmatrix(mpt[[pcpi]],refdel=0, pnratio,targetdel=0, ncore=no_cores)
    
    save(pnratioTRIM, file=paste0(matrixChild_pnratio_dir[ci],"TRIM_",pcp,".Rdata"))
    save(matrixChild_entropyTRIM, file=paste0(matrixChild_entropy_dir[ci],"TRIM_",pcp,".Rdata"))
    pnratioTRIMREDUCE = foreach(k=1:length(pnratioTRIM),.combine="cbind") %dopar% { return(pnratioTRIM[[k]]) }
  }
  
  TimeOutput(start)
  
  
  
  
  
  
  
  
  
  
  
  ## Parent Contribution ---------------------------------------------------------------
  cat(", parentcontribution")
  for (pcpi in 1:length(pheno_type)) { 
    pcp = pheno_type[pcpi]
    cat(" ",pcp)
    if (trimOnly) {
      #build original values of wt and ko
      matrixWTCountAdjFULL = foreach(i=1:ncol(mlf[[pcpi]]),.combine="cbind") %dopar% { return(mkc[[pcpi]][,i]/exp(mlf[[pcpi]][,i])) }
      colnames(matrixWTCountAdjFULL) = colnames(mkc[[pcpi]])
      contrib = foreach(i = 1:length(phenoParent)) %dopar% { #for each phenotype
        change = mkc[[pcpi]][,phenoParent_ind[i]] - matrixWTCountAdjFULL[,phenoParent_ind[i]]
        koparents = mkc[[pcpi]][,phenoParent[[i]]]
        wtparents = matrixWTCountAdjFULL[,phenoParent[[i]]]
        if (length(phenoParent[[i]])==1) {
          koparents = matrix(koparents,ncol=1)
          wtparents = matrix(wtparents,ncol=1)
        }
        changeparents = koparents - wtparents
        if (length(which(changeparents==0))>0) changeparents[changeparents==0] = 1
        
        contrib = change / changeparents
        
        rownames(contrib) = rownames(mkc[[pcpi]])
        colnames(contrib) = colnames(mkc[[pcpi]])[phenoParent[[i]]]
        return(contrib)
      }
      names(contrib) <- phenoMeta$phenotype[phenoParent_ind]
      
      save(contrib, file=paste0(matrixParent_contrib_dir[ci],"_",pcp,".Rdata"))
      contribREDUCE = foreach(k=1:length(contrib),.combine="cbind") %dopar% { return(contrib[[k]]) }
      write.csv(contribREDUCE, file=paste0(matrixParent_contrib_dir[ci],"_",pcp,".csv"))
      
    } else {
      load(paste0(matrixParent_contrib_dir[ci],"_",pcp,".Rdata"))
      names(contrib) <- phenoMeta$phenotype[phenoParent_ind]
      
      save(contrib, file=paste0(matrixParent_contrib_dir[ci],"_",pcp,".Rdata"))
      contribREDUCE = foreach(k=1:length(contrib),.combine="cbind") %dopar% { return(contrib[[k]]) }
      write.csv(contribREDUCE, file=paste0(matrixParent_contrib_dir[ci],"_",pcp,".csv"))
    } 
    
    #trim based on pvalue trim
    cat(" ",pheno_type[pcpi],sep="")
    contribTRIM = trimlistmatrix(mpt[[pcpi]],refdel=0, contrib,targetdel=0, ncore=no_cores)
    
    save(contribTRIM, file=paste0(matrixParent_contrib_dir[ci],"TRIM_",pcp,".Rdata"))
    a = foreach(k=1:length(contribTRIM),.combine="cbind") %dopar% { return(contribTRIM[[k]]) }; rownames(a) = sampleMeta$gene[match(rownames(a),sampleMeta$fileName)]
    write.table(as.matrix(a), file=paste0(matrixParent_contrib_dir[ci],"TRIM_",pcp,".csv"), sep=",", col.names=F)
  }
  TimeOutput(start)
  
  
  
  
  
  
  
  
  
  
  ## Parent effort -------------------------------------------------------
  cat(", parenteffort")
  for (pcpi in 1:length(pheno_type)) { 
    pcp = pheno_type[pcpi]
    cat(" ",pcp)
    if (trimOnly) {
      if (!exists("matrixWTCountAdjFULL")) {
        matrixWTCountAdjFULL = foreach(i=1:ncol(mlf[[pcpi]]),.combine="cbind") %dopar% { return(mkc[[pcpi]][,i]/exp(mlf[[pcpi]][,i])) }
        colnames(matrixWTCountAdjFULL) = colnames(mkc[[pcpi]])
      }
      effort = foreach(i = 1:length(phenoParent)) %dopar% {
        effort = mlf[[pcpi]][,phenoParent_ind[i]] - mlf[[pcpi]][,phenoParent[[i]]]
        if (length(phenoParent[[i]])==1) effort = matrix(effort,ncol=1)
        rownames(effort) = rownames(mkc[[pcpi]])
        colnames(effort) = colnames(mkc[[pcpi]])[phenoParent[[i]]]
        return(effort)
      }
      names(effort) <- phenoMeta$phenotype[phenoParent_ind]
      
      save(effort, file=paste0(matrixParent_effort_dir[ci],"_",pcp,".Rdata"))
      effortREDUCE = foreach(k=1:length(effort),.combine="cbind") %dopar% { return(effort[[k]]) }
      write.csv(effortREDUCE, file=paste0(matrixParent_effort_dir[ci],"_",pcp,".csv"))
    } else {
      load(paste0(matrixParent_effort_dir[ci],"_",pcp,".Rdata"))
      names(effort) <- phenoMeta$phenotype[phenoParent_ind]
      
      save(effort, file=paste0(matrixParent_effort_dir[ci],"_",pcp,".Rdata"))
      effortREDUCE = foreach(k=1:length(effort),.combine="cbind") %dopar% { return(effort[[k]]) }
      write.csv(effortREDUCE, file=paste0(matrixParent_effort_dir[ci],"_",pcp,".csv"))
    }
    
    #trim based on pvalue trim
    cat(" ",pheno_type[pcpi],sep="")
    effortTRIM = trimlistmatrix(mpt[[pcpi]],refdel=0, effort,targetdel=0, ncore=no_cores)
    
    save(effortTRIM, file=paste0(matrixParent_effort_dir[ci],"TRIM_",pcp,".Rdata"))
    b = foreach(k=1:length(effortTRIM),.combine="cbind") %dopar% { return(effortTRIM[[k]]) }; rownames(b) = sampleMeta$gene[match(rownames(b),sampleMeta$fileName)]
    write.table(as.matrix(b), file=paste0(matrixParent_effort_dir[ci],"TRIM_",pcp,".csv"), sep=",", col.names=F)
  }
  
  TimeOutput(start)
  
}

TimeOutput(start)




