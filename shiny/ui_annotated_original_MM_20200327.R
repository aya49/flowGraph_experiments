#Shiny APP for flowType/Rchy visulaization and checking outliers
#Written by Mehrnoush Malek (MM)
#Revise: March 2020, MM

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

# get an example flowtype file
ft <- get(load("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/semisupervised/flowType/BM 01092015_L000096981_001.labelled.fcs_FT.Rdata"))

# get max layers
max.p <- max(unlist(lapply(ft@PhenoCodes, function(x) 
    length(which(as.numeric(unlist(strsplit(x,split="")))!=0)))))


ui <- shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(),
    titlePanel("something nice", title="tab name"),
    shinydashboard::dashboardSidebar( 
        shinyFiles::shinyDirButton("result", label="flowType folder", "Select"),
        shinyFiles::shinyDirButton("files", label="fcs sample folder", "Select"),
        shinyFiles::shinyFilesButton("arguments", multiple=F, label="Choose arguments", "Upload"),
        shiny::radioButtons(inputId="random", label="Plot random samples?",
                     choices=c("YES","NO (select a sample)"), selected="YES"),
        shinyFiles::shinyFilesButton("fcs", multiple=F, label="Choose a sample ", "Upload"),
        shiny::radioButtons(inputId="calc", label="Load cell counts?",
                     choices=c("YES","NO"), selected="YES"),
        shiny::radioButtons(inputId="method1", label="Statistical test",
                     choices=c("p-value","Wilcoxon"), selected="p-value"),
        shiny::radioButtons(inputId="method2", label="P-value correction method",
                     choices=c("BH","Bonferroni","fdr","NULL"), selected="BH"),
        shiny::sliderInput("slider1", "Cell count", min=2, max=1000, step=100, value=30),
        shiny::sliderInput("slider2", "Max layer", min=2, max=max.p, step=1, value=min(max.p, 6)),
        shiny::actionButton("go", "Start the analysis")
    ),
    shinydashboard::dashboardBody(
        shiny::fluidRow(column(7, visNetworkOutput('plotgraph',height=800))),
        shiny::fluidRow(column(12, ggiraphOutput('boxplot',width=1000,height=900))),
        shiny::fluidRow(plotOutput('plotdens', width=1250,height=600))
    )
)




