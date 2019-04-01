# Uses different distance/linear measures to calculate distance & plot samples
# aya43@sfu.ca 20161220

root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/matrix", sep="")
# mcp_types = c("/original", "/child_pn","/child_prop", "/trim") 
# matrix_type1 = c("CountAdj") #, "Prop", "Count",  #countAdj comes first, to set which columns are deleted based on cellcountthres
# matrix_type2 = c("Child_prop")
# matrix_type3 = c("Child_pnratio", "Child_entropy")
# matrix_type4 = c("Child_pnratio", "Child_entropy")

# node only
# matrix_type = c("Child_entropy", "Child_pnratio", "Child_prop")
# matrix_weights = c(.3,.3,.)
# weight_matrix = c("CountAdj")

# with pvalued features
matrix_count = c("CountAdj")
matrix_type0 = list(c(273,471,460,505))
matrix_weights0 = list(c(.25,.25,.25,.25))
weight_matrix = c("MaxCountAdj_CountAdj")
# logweight = T #log weight first?
# absweight = T

phenoChild_dir = paste(result_dir, "/phenoChild.Rdata",sep="")
phenoChild_ind_dir = paste(result_dir, "/phenoChild_ind.Rdata",sep="")
phenoChildpn_dir = paste(result_dir, "/phenoChildpn.Rdata",sep="")
phenoChildpn_ind_dir = paste(result_dir, "/phenoChildpn_ind.Rdata",sep="")
cellCountThres = c(200) #insignificant if count under

ignoredist = ".csv"


#Output
dist_dir = paste(result_dir, "/dist", sep=""); for (i in 1:length(dist_dir)) { suppressWarnings(dir.create (dist_dir[i])) }
dist_type_dir = NULL
# for(i in 1:length(mcp_types)) { for(j in 1:length(dist_dir)) { dist_type_dir = append(dist_type_dir, paste0(dist_dir[j], mcp_types[i]))} }
# for (i in 1:length(dist_type_dir)) { suppressWarnings(dir.create (dist_type_dir[i])) }

libr(stringr)
libr(colorspace)
libr(vegan) # libr(proxy)
libr(foreach)
libr(doMC)
libr(lubridate) #if there are date variables
libr(pracma) #for dot product
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")


no_cores = 5#detectCores()-1
registerDoMC(no_cores)



start = Sys.time()


phenoMeta = get(load(phenoMeta_dir))
sampleMeta = get(load(sampleMeta_dir))

start2 = Sys.time()
centre = paste0(panelL,centreL)
cat("\n",paste0(panelL," ",centreL))

sampleMeta0 = get(load(sampleMeta_dir))
distmfile = list.files(dist_dir, recursive=F, full.names=T, pattern=".Rdata")
delpvalind = grepl(ignoredist,distmfile,ignore.case=T)
if (sum(delpvalind)>0) distmfile = distmfile[!delpvalind]
distmfilenames = fileNames(distmfile)

# morecombos = ""
# while (!grepl("y",morecombos,ignore.case=T) | !grepl("n",morecombos,ignore.case=T)) {
#   prompt <- "do you want to input more linear combos? (type \'yes\' or \'no\'): "
#   morecombos <- readline(prompt)
# }
# distmfilenames #list distance files
# if (grepl("y",morecombos,ignore.case=T)) {
#   input = ""
#   inputno = length(matrix_type0)+1
#   while(!grepl("x",input,ignore.case=T))
#     prompt <- "\ntype \'x\' to exit, or type your combos (weight #) (e.g. 0.5 10 0.5 20 0.5 2): "
#   input = strsplit(readline(prompt), " ")[[1]]
#   if (!grepl("x",input,ignore.case=T)) {
#     input = as.numeric(input)
#     matrix_type0[[inputno]] = input[seq(from=2,to=length(input),by=2)]
#     matrix_weights0[[inputno]] = input[seq(from=1,to=length(input),by=2)]
#     cat(matrix_weights0[[inputno]],distmfilenames[matrix_type0[[inputno]]])
#     prompt <- "\nplease confirm (\'yes\' or \'no\')"
#     input = readline(prompt)
#     if (grepl("n",input,ignore.case=T)) {
#       matrix_type0[[inputno]] = NULL
#       matrix_weights0[[inputno]] = NULL
#     }
#   }
# }

for (mi in 1:length(matrix_type0)) {
  leavePhenotype = list()
  doneAll = F
  
  start2 = Sys.time()
  matrix_type. = distmfile[matrix_type0[[mi]]]
  matrix_weights = matrix_weights0[[mi]]
  cat("\n"); cat(paste(matrix_weights,matrix_type., sep=" x "),sep=",  ")
  
  #load different matrices
  mresult = Loadintermatrices(matrix_type.)
  mml = mresult$mml
  mmlname = names(mml)
  pt = mresult$pt
  gt = mresult$gt
  
  # Calculate distance (normalize by cell population) --------------------------------------------------------
  
  if (doneAll) dname = paste(dist_dir, "/linear_", paste0(matrix_type.,"-",matrix_weights, collapse="_"), "_FULL_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_normalize-cellpop", sep="")
  dname = paste(dist_dir, "/linear_", paste0(matrix_type.,"-",matrix_weights, collapse="_"), "_layer-", str_pad(k, 2, pad = "0"), "_countThres-", countThres, "_normalize-cellpop", sep = "" )
  if (file.exists(dname)) { cat(" exists, skipped."); next }
  
  #weigh each feature matrix; scalendev = 0 (uhh...), T (scaling is done by dividing the centred columns by their SD), F (none)
  for (s in 1:3) {
    # for (l in 1:2) {
    # if (l==1 & logweight) dfinal = dl
    # if (l==2 & absweight) dfinal = d
    
    # if (s==1) { scalendev=0; dn = dfinal
    if (s==1) { scalendev=0; dn = mml
    } else {
      if (s==2) scalendev=T
      if (s==3) scalendev=F
      dn = foreach(i=1:length(mml)) %dopar% { 
        a = scale(as.vector(mml[[i]]), scale=scalendev)
        a = a-min(a)
        return(matrix(a,nrow=nrow(mml[[i]])))
      }
    }
    d0 = foreach(i=1:length(dn)) %dopar% { return(matrix_weights[i]*dn[[i]]) }
    d1 = Reduce("+",d0)
    colnames(d1) = rownames(d1) = gt
    
    # if (l==1 & logweight) {
    #   save(d1, file=paste0(dname,"_logweight_scaleSD-",scalendev,".Rdata"))
    #   write.csv(d1, file=paste0(dname,"_logweight_scaleSD-",scalendev,".csv"))
    # }
    # if (l==2 & absweight) {
    #   save(d1, file=paste0(dname,"_absweight_scaleSD-",scalendev,".Rdata"))
    #   write.csv(d1, file=paste0(dname,"_absweight_scaleSD-",scalendev,".csv"))
    # }
    save(d1, file=paste0(dname,"_scaleSD-",scalendev,".Rdata"))
    write.csv(d1, file=paste0(dname,"_scaleSD-",scalendev,".csv"))
  }
  TimeOutput(start2)
}



TimeOutput(start)



#require(RDRToolbox)
#libr(rgl)
## Isomap ##



# libr(dbscan)
# oc = optics(d, eps=10, minPts=3) # no results (gravitate towards each other)
