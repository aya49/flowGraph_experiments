# Shiny APP for flowGraph visulaization
# Written by: Mehrnoush Malek (MM), Alice Yue (AY)
# Revision: March 2020, AY

# header <- shinydashboard::dashboardHeader(
# )

# sidebar <- shinydashboard::dashboardSidebar(
    # width=350,
    # uiChangeThemeDropdown(dropDownLabel="change theme", defaultTheme="grey_dark"),

    # menuItem("top", href="#top", icon=icon("dashboard")),
    # menuItem("summary statistics table", href="#table_bg", icon=icon("table")),
    # menuItem("cell hierarchy", href="#ch_bg", icon=icon("sitemap")),
    # menuItem("boxplots", href="#box_bg", icon=icon("boxes")),
    # menuItem("scatterplots", href="#dens_bg", icon=icon("chart"))

    # # root name
    # shiny::textInput( # renderText({ input$root_name })
    #     inputId="root_name", label="name of root phenotype",
    #     value="live", width=NULL, placeholder="live"),
    # shiny::checkboxInput(inputId="show_bgedges", value=FALSE,
    #                      label="cell hierarchy plot: show background edges"),

    # shiny::actionButton('go', 'start analysis')
# )


body <- shinydashboard::dashboardBody(
    # Custom theme
    # uiChangeThemeOutput(),
    title="testing!",
    # shinyDashboardThemes(theme=dbtheme),

    tags$script(HTML("$('body').addClass('fixed');")), # fix header and sidebar

    shinydashboard::box(
        id="top", title="flowGraph", status="primary",
        width=12,
        shiny::uiOutput("about")
    ),

    shinydashboard::box(
        title="1. select flowGraph file.",
        status="primary", solidHeader=TRUE, collapsible=TRUE,
        width=3, height='auto',
        shiny::fileInput(
            inputId="fg", label="select a flowGraph file (.Rdata/.rds)",
            buttonLabel="browse", multiple=FALSE,
            accept=c("Rdata/RDS", "Rdata/RDS", "*")
        )
    ),

    shinydashboard::box(
        id="table_bg", title="2. select a summary statistic from table.",
        status="primary", solidHeader=TRUE, collapsible=TRUE,
        width=9,
        shinydashboard::box(
            title="options",
            status="primary",
            width=4,
            shiny::sliderInput(
                inputId="p_thres",
                label="p-value threshold",
                min=.01, max=1, step=.01, value=.01
            ),
            shiny::selectInput(
                inputId="effect_size",
                label="minimum effect size",
                choices=c("negligible", "small", "medium", "large"),
                selected="medium"
            ),
            shiny::sliderInput(
                inputId="adjust0",
                label="(SpecEnr only) max % of 0's allowed",
                min=.1, max=1, step=.1, value=.5
            )
        ),
        shinydashboard::box(
            title="summary statistics table",
            status="primary",
            width=8,
            shiny::verbatimTextOutput('table_row'),
            DT::dataTableOutput("table_summary") %>%
                shinycssloaders::withSpinner(color="#0dc5c1")
        )
    ),

    shinydashboard::tabBox(
        id="ch_bg", title="3. select a cell population from cell hierarchy.",
        # status="primary", solidHeader=TRUE, collapsible=TRUE,
        width=12, height='auto', side="right",
        shiny::tabPanel(
            height='auto',
            "cell hierarchy plot",
            shiny::sliderInput(
                inputId="max_nodes",
                label="# of phenotypes to show (max)",
                min=2, max=100, step=1, value=20),
            visNetwork::visNetworkOutput("plot_visnet", height='600') %>%
                shinycssloaders::withSpinner(color="#0dc5c1")
        ),
        shiny::tabPanel(
            height='auto',
            "other summary plots",
            shinydashboard::box(
                title="QQ plot",
                status="primary",
                width=6, height='auto',
                shiny::plotOutput('plot_qq') %>%
                    shinycssloaders::withSpinner(color="#0dc5c1")
            ),
            shinydashboard::box(
                title="p-value vs feature difference plot",
                status="primary",
                width=6, height='auto',
                shiny::plotOutput("plot_pVSdifference") %>%
                    shinycssloaders::withSpinner(color="#0dc5c1")
            )
        )
    ),

    shinydashboard::box(
        id="box_bg", title="boxplot options",
        status="primary",
        width=12, height='auto',
        shiny::column(
            4, shiny::checkboxInput(inputId="dotplot", value=TRUE,
                                    label="boxplot: include dotplot")
        ),
        shiny::column(
            4, shiny::checkboxInput(inputId="outlier", value=TRUE,
                                    label="boxplot: include outliers")
        ),
        shiny::column(
            4, shiny::checkboxInput(inputId="paired", value=FALSE,
                                    label="boxplot: pair samples")
        )
    ),
    shinydashboard::box(
        title="boxplot",
        status="primary", solidHeader=TRUE, collapsible=TRUE,
        width=12, height='auto',
        shinydashboard::box(
            status="primary",
            width=4, height='auto',
            shiny::plotOutput('plot_box') %>%
                shinycssloaders::withSpinner(color="#0dc5c1")
        ),
        shinydashboard::box(
            status="primary",
            width=4, height='auto',
            shiny::plotOutput('plot_box1') %>%
                shinycssloaders::withSpinner(color="#0dc5c1")
        ),
        shinydashboard::box(
            status="primary",
            width=4, height='auto',
            shiny::plotOutput('plot_box2') %>%
                shinycssloaders::withSpinner(color="#0dc5c1")
        )
    ),

    shinydashboard::box(
        id="dens_bg", title="scatterplot FCS & flowType files.",
        status="primary",
        width=12, height='auto',

        shiny::fluidRow(
            shiny::column(
                1,
                shinyFiles::shinyDirButton(id="fcs_dir", "select", "select")
            ),
            shiny::column(
                11,
                shiny::verbatimTextOutput("fcs_dir_out")
            )
        ),
        shiny::fluidRow(
            shiny::column(
                1,
                shinyFiles::shinyDirButton(id="ft_dir", "select", "select")
            ),
            shiny::column(
                11,
                shiny::verbatimTextOutput("ft_dir_out")
            )
        )
    ),
    shinydashboard::box(
        title = "scatterplot(s)",
        status="primary", solidHeader=TRUE, collapsible=TRUE,
        width=12, height='auto',
        shinydashboard::box(
            status="primary",
            width=6, height='auto',
            shiny::plotOutput('plot_dens1') %>%
                shinycssloaders::withSpinner(color="#0dc5c1")
        ),
        shinydashboard::box(
            status="primary",
            width=6, height='auto',
            shiny::plotOutput('plot_dens2') %>%
                shinycssloaders::withSpinner(color="#0dc5c1")
        )
    )
)


dashboardPage(dashboardHeader(title="testing!"),
              shinydashboard::dashboardSidebar(disable=TRUE), body)
