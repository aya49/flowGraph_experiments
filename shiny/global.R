# IMPC/IMPC-Results/3i/Panel_BM-cell/semisupervised/. Input folders are flowType and ParentPopulation. plus Settingup_argument.RData I mentioned.

# rstudioapi::openProject("/mnt/f/Brinkman group/current/Alice/flowGraph") # CHANGE THIS
# runApp('/mnt/f/Brinkman group/current/Alice/flowGraph/flowType3/shiny')
# install.packages("/mnt/f/Brinkman group/current/Alice/flowtype_metric/src/flowTypeFilter/flowTypeFilterC_1.0_R_x86_64-pc-linux-gnu.tar.gz",repos=NULL, type="source")
setwd("/mnt/f/Brinkman group/current/Alice/flowtype_metric/src/shiny")
source("_func.R")
# source("moduleChangeTheme.R")

start_dir <<- c(data="/mnt/f/FCS data/Tobias Kollmann")
# rslt_dir <<- c(data="/mnt/f/Brinkman group/current/Alice/flowtype_metric/result")
# ft_dir <<- "/mnt/f/FCS data/Tobias Kollmann/flowType_results/Bcell_panel"
# fcs_dir <<- "/mnt/f/FCS data/Tobias Kollmann/flowType_data/Bcell_panel"
fg_dir <<- "/mnt/f/Brinkman group/current/Alice/flowtype_metric/result/epichipc_bcell/fg.Rdata"
fg <<- NULL
# fg <<- get(load(fg_dir))

# library(flowType)
# library(flowCore)
library(flowDensity)
library(flowTypeFilterC)
library(flowGraph)
library(flowWorkspace)
library(shinyFiles)
library(shiny)
library(ggiraph)
library(ggplot2)
library(visNetwork)
library(igraph)
library(DT)

library(shiny)
library(shinydashboard)
library(shinyFiles)
# library(shinyBS)
library(shinyWidgets)
library(bsplus)
library(ggiraph)
# library(ShinySky) #devtools::install_github("AnalytixWare/ShinySky") # search box
library(shinycssloaders) # loading spinners

# runApp('/mnt/f/Brinkman group/current/Alice/flowGraph/flowType3/shiny')

options(shiny.maxRequestSize=100000*1024^2) # file upload size limit 10gb


