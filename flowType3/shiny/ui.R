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


ui <- dashboardPage(
    dashboardHeader(title="testing!"),
    dashboardSidebar(
        sidebarMenu(
            shinydashboard::dashboardSidebar( 
                shinyFiles::shinyDirButton("ft_dir", label="flowType folder", "Select"),
                shinyFiles::shinyDirButton("fcs_dir", label="fcs sample folder", "Select"),
                shinyFiles::shinyDirButton("fg_dir", label="flowGraph folder", "Select")
                # shiny::radioButtons(inputId="random", label="Plot random samples?",
                #                     choices=c("YES","NO (select a sample)"), selected="YES"),
                # shinyFiles::shinyFilesButton("fcs", multiple=F, label="Choose a sample ", "Upload"),
                # shiny::radioButtons(inputId="calc", label="Load cell counts?",
                #                     choices=c("YES","NO"), selected="YES"),
                # shiny::radioButtons(inputId="method1", label="Statistical test",
                #                     choices=c("p-value","Wilcoxon"), selected="p-value"),
                # shiny::radioButtons(inputId="method2", label="P-value correction method",
                #                     choices=c("BH","Bonferroni","fdr","NULL"), selected="BH"),
                # shiny::sliderInput("slider1", "Cell count", min=2, max=1000, step=100, value=30),
                # shiny::sliderInput("slider2", "Max layer", min=2, max=max.p, step=1, value=min(max.p, 6)),
                # shiny::actionButton("go", "Start the analysis")
            )
        )
    ),
    dashboardBody(
        fluidRow(
            column(3, DT::dataTableOutput("table_summary")),
            column(9, ggiraph::girafeOutput("plot_hierarchy"))
        ),
        fluidRow(
            column(3, ggiraph::girafeOutput('plot_box',height=800)),
            column(9, plotOutput('plot_dens',width=1000,height=900)))
    )
)




