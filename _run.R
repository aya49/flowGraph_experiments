root = "/mnt/f/Brinkman group/current/Alice/flowtype_metric"
setwd(root)

source("source/00_data_bodenmiller.R")
source("source/00_data_flowcap.R")
source("source/00_data_pregnancy.R")
source("source/00_data_genetech.R")
source("source/00_data_pos.ctrl.R")
source("source/01_feat_normalizeadj.R")
source("source/02_data_flowcap2.R")
source("source/03_feat_childparent.R")
# source("source/05_stat_checkmatrix.R")
source("source/05_feat_pregnancy.R")
source("source/06_pval.R")
source("source/07_pval_stats.R")
source("source/08_pval_plot.R")



