# Shiny APP for flowGraph visulaization
# Written by: Mehrnoush Malek (MM), Alice Yue (AY)
# Revision: March 2020, AY

# IMPC/IMPC-Results/3i/Panel_BM-cell/semisupervised/. Input folders are flowType and ParentPopulation. plus Settingup_argument.RData I mentioned.

# runApp('/mnt/f/Brinkman group/current/Alice/flowGraph/flowType3/shiny')
setwd("/mnt/f/Brinkman group/current/Alice/flowtype_metric/src/flowType3/shiny")

library(shiny)
library(shinydashboard)
library(flowTypeFilter)
library(shinythemes)
library(RchyOptimyx)
library(shinyFiles)
library(shinyWidgets)
library(flowType)
library(rdist)
library(visNetwork)
library(ggiraph)
library(ggplot2)
library(flowCore)

# # get an example flowtype file
# ft <- get(load("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/semisupervised/flowType/BM 01092015_L000096981_001.labelled.fcs_FT.Rdata"))
# 
# # get max layers
# max.p <- max(unlist(lapply(ft@PhenoCodes, function(x) 
#     length(which(as.numeric(unlist(strsplit(x,split="")))!=0)))))


ui <- shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(),
    titlePanel("something nice", title="tab name"),
    shinydashboard::dashboardSidebar( 
        # shinyFiles::shinyDirButton("ft_folder", label="folder with flowType files (.Rdata)", "Select"),
        # shinyFiles::shinyDirButton("fcs_folder", label="folder with processed fcs files (.Rdata/.fcs)", "Select"),
        # shinyFiles::shinyDirButton("fg_folder", label="folder where flowGraph object is saved", "Select"),
        # textInput("root_name", label="name of root node", value ="Live"),
        # shiny::actionButton("go", "Start the analysis")
    ),
    shinydashboard::dashboardBody(
        shiny::fluidRow(
            column(3, DT::dataTableOutput("table_summary")),
            column(9, ggiraph::girafeOutput("plot_hierarchy"))
        ),
        shiny::fluidRow(
            column(3, ggiraph::girafeOutput('plot_box',height=800)),
            column(9, plotOutput('plot_dens',width=1000,height=900)))
    )
)




