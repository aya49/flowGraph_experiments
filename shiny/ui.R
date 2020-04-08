# Shiny APP for flowGraph visulaization
# Written by: Mehrnoush Malek (MM), Alice Yue (AY)
# Revision: March 2020, AY

# IMPC/IMPC-Results/3i/Panel_BM-cell/semisupervised/. Input folders are flowType and ParentPopulation. plus Settingup_argument.RData I mentioned.

# runApp('/mnt/f/Brinkman group/current/Alice/flowGraph/flowType3/shiny')
setwd("/mnt/f/Brinkman group/current/Alice/flowtype_metric/src/flowType3/shiny")

library(shiny)
library(shinydashboard)
# library(shinythemes)
library(shinyFiles)
library(DT)
# library(shinyWidgets)
library(ggiraph)

options(shiny.maxRequestSize=10000*1024^2) # file upload size limit 10gb

ui <- dashboardPage(
    dashboardHeader(title="testing!"),

    dashboardSidebar(
        sidebarMenu(
            # tags$script('$( "#fg" ).on( "click", function() { this.value = null; });'),
            shinydashboard::dashboardSidebar(
                shiny::fileInput(
                    inputId="fg", label="select flowGraph file (.Rdata/.rds)",
                    buttonLabel="browse", multiple=FALSE,
                    accept = c("Rdata/RDS", "Rdata/RDS", "*")),

                shinyFiles::shinyDirButton("ft_dir", label="flowType folder", "select"),
                shinyFiles::shinyDirButton("fcs_dir", label="fcs sample folder", "select"),
                # shinyFiles::shinyDirButton("fg_dir", label="flowGraph folder", "Select"),

                # Horizontal line ----
                tags$hr(),

                # root name
                shiny::textInput( # renderText({ input$root_name })
                    inputId="root_name", label="name of root phenotype",
                    value="live", width=NULL, placeholder="live"),
                shiny::sliderInput(
                    inputId="p_thres",
                    label="cell hierarchy plot: p-value threshold",
                    min=.01, max=1, step=.01, value=.01),
                shiny::sliderInput(
                    inputId="max_nodes",
                    label="cell hierarchy plot: # of phenotypes to show (max)",
                    min=2, max=100, step=1, value=20),
                shiny::radioButtons(
                    inputId="dotplot", label="boxplot: include dotplot?",
                    choices=c("yes", "no"), selected=c("yes")),
                shiny::radioButtons(
                    inputId="outlier", label="boxplot: include outliers?",
                    choices=c("yes", "no"), selected=c("yes")),
                shiny::radioButtons(
                    inputId="paired", label="boxplot: should samples be paired?",
                    choices=c("yes", "no"), selected=c("no"))

                # shiny::actionButton('go', 'start analysis')
            )
        )
    ),
    dashboardBody(
        fluidRow(
            column(9, DT::dataTableOutput("table_summary")),
            column(3, shiny::verbatimTextOutput('table_row'))
        ),
        fluidRow(
            column(9, ggiraph::girafeOutput("plot_hierarchy")),
            column(3, ggiraph::girafeOutput('plot_box',height=800))
        ),
        fluidRow(
            column(12, shiny::plotOutput('plot_dens',width=1000,height=900))
        )
    )
)




