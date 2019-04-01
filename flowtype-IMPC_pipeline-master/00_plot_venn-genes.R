## Alice Yue 20180219
## Input: gene lists from each centre
## Output: venn diagram of overlapping genes & its number of files

root = "~/projects/IMPC"
setwd(root)

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)


data_dir = "/mnt/f/FCS data/IMPC/IMPC-Summary/flowData_Spreadsheets/"
centre1 = c("BCM", "CIPHE", "GMC", "H", "TCP")
centre = c(centre1, "Sanger")
centre_sheets_dir = c(paste0(centre1, "/", centre1, "_Panel1-flowData.csv"), "Sanger/T-cell_Panel/T-cell_SPLEEN-flowData.csv")

source("~/projects/IMPC/code/_funcAlice.R")
libr("VennDiagram")

genes = list()
for (i in 1:length(centre_sheets_dir)) {
  sheet_dir = centre_sheets_dir[i]
  sheet = read.csv(paste0(data_dir,sheet_dir), header=T,sep=',',row.names=1)
  cat("\n\n",sheet_dir,": ",colnames(sheet))
  gene = as.character(sheet[,"Genotype"])
  gene[gene == "+_+" | gene == "+_Y"] = "Wildtype"
  gene = gene[!grepl("failed",gene)]
  gene = gsub("(b)","",gene, fixed=T)
  gene = gsub("(e)","",gene, fixed=T)
  gene = gsub(" ","",gene, fixed=T)
  gene = gsub("(tm1)","",gene, fixed=T)
  gene = gsub(" KO","",gene, fixed=T)
  gene = gsub("(tm1.1)","",gene, fixed=T)
  gene = strsplit(gene,"_")
  gene = sapply(gene, function(x) x[1])
  
  # gene = toupper(gene)
  # gene = sapply(gene, function(x) {
  #   xx = x[1]
  #   if (length(x)>1) {
  #     if (x[2]=="X" | x[2]=="+") {
  #       xx = paste0(x[1],"_X")
  #     }
  #   }
  #   return (xx)
  # })
  
  # gene = strsplit(genes, "[(]")
  # gene = sapply(gene, function(x) x[1])
  # gene = strsplit(genes, "\\s+", fixed=T)
  # gene = sapply(gene, function(x) x[1])
  
  
  genes[[centre[i]]] = gene
  # cat(paste(sort(unique(genes[[centre[i]]])), set='\n'))
}

venn.plot = venn.diagram(
  genes[c(2:6)],
  paste0("gene_venn_2-6_plain.tiff")
)



