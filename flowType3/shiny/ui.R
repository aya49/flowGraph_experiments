#Shiny APP for flowType/Rchy visulaization and checking outliers
#Written by Mehrnoush Malek (MM)
#Revise: March 2020, MM

# IMPC/IMPC-Results/3i/Panel_BM-cell/semisupervised/. Input folders are flowType and ParentPopulation. plus Settingup_argument.RData I mentioned.

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
tmp <- load("~/data/semisupervised/flowType/BM 01092015_L000096981_001.labelled.fcs_FT.Rdata")
ft <- get(tmp)
rm(tmp)
max.p <- max(unlist(lapply(ft@PhenoCodes, function(x) length(which(as.numeric(unlist(strsplit(x,split = "")))!=0)))))
ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar( shinyDirButton("result",label="flowType folder", "Select"),
                    shinyDirButton("files",label="sample folder", "Select"),
                    radioButtons(inputId="calc", label="Load cell counts?",
                                 choices=c("YES","NO"),selected = "YES"),
                    radioButtons(inputId="method1", label="Method?",
                                 choices=c("p-value","Wilcoxon"),selected = "p-value"),
                    radioButtons(inputId="method2", label="Correction?",
                                 choices=c("BH","Bonferroni","fdr","NULL"),selected = "BH"),
                    sliderInput("slider2","Population level", min=2, max=max.p, step =1, value=6),
                  sliderInput("slider1","Cell count", min=2, max=1000, step =100, value=30),
                  radioButtons(inputId="random", label="Plot random samples?",
                               choices=c("YES","NO (select a sample)"),selected = "YES"),
                  shinyFilesButton("fcs",multiple=F, label = "Choose a sample ", "Upload"),
                  actionButton("go", "Start the analysis")),
  dashboardBody(
    fluidRow(column(7, visNetworkOutput('plotgraph',height = 800))),
    fluidRow(column(12, ggiraphOutput('boxplot',width = 1000,height = 900))),
    fluidRow(plotOutput('plotdens',width = 1250,height = 600)))
)




