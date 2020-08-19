body <- shinydashboard::dashboardBody(


    #### 0 about ####
    shiny::fluidRow(
        shiny::column(
            7,
            shinydashboard::box(
                background="black",
                # collapsible=TRUE,
                width=NULL,
                shiny::uiOutput("about_ui")
            ),

            #### 1 load flowGraph file, options, summary stats data frame ####
            shinydashboard::box(
                title="1. load flowGraph file",
                status="primary", solidHeader=TRUE,
                # collapsible=TRUE,
                width=NULL,
                shiny::fluidRow(
                    shiny::column(
                        width=6,
                        shiny::fileInput(
                            inputId="fg",
                            label="select a flowGraph file (.Rdata/.rds)",
                            buttonLabel="browse", multiple=FALSE,
                            accept=c("Rdata/RDS", "Rdata/RDS", "*")
                        ),
                        shiny::uiOutput("fg_descript_ui"),
                            shiny::selectInput(
                                inputId="xsigfrom",
                                label="SpecEnr sig (from)",
                                choices=c("all", "insig", "sig"),
                                selected="all"
                            ),
                            shiny::selectInput(
                                inputId="xsigto",
                                label="SpecEnr sig (to)",
                                choices=c("all", "insig", "sig"),
                                selected="all"
                            ),
                            shiny::selectInput(
                                inputId="ysigfrom",
                                label="raw sig (from)",
                                choices=c("all", "insig", "sig"),
                                selected="all"
                            ),
                            shiny::selectInput(
                                inputId="ysigto",
                                label="raw sig (to)",
                                choices=c("all", "insig", "sig"),
                                selected="all"
                            )
                    ),
                    shiny::column(
                        width=6,
                        shiny::uiOutput("ui_layer"),
                        shiny::uiOutput("ui_count"),
                        shiny::sliderInput(
                            inputId="p_t1",
                            label="p-value threshold (SpecEnr)",
                            min=.01, max=1, step=.01, value=.05),
                        shiny::sliderInput(
                            inputId="p_t2",
                            label="p-value threshold (raw)",
                            min=.01, max=1, step=.01, value=.05),
                        shiny::selectInput(
                            inputId="effect_size",
                            label="minimum effect size",
                            choices=c("negligible", "small", "medium", "large"),
                            selected="medium"
                        ),
                        shiny::sliderInput(
                            inputId="adjust0",
                            label="(SpecEnr only) max % of 0's allowed",
                            min=.1, max=1, step=.1, value=.4
                        )
                    ),
                    # fg data frame
                    shiny::column(
                        width=12,
                        DT::dataTableOutput("table_summary")
                    ),
                    # go! make hierarchy and summary plots
                    shiny::column(
                        width=12,
                        shiny::actionButton(
                            inputId="go",
                            label="generate/refresh plots"
                        ),
                        tags$style(".swal-modal {width: 80%;}"),
                        style="text-align: center"
                    )
                )
            )
        ),
        shiny::column(
            5,
            #### 2 summary plots ####
            shinydashboard::box(
                title="2. plots",
                status="primary", solidHeader=TRUE,
                # collapsible=TRUE,
                width=NULL,

                shinydashboard::box(
                    title="SpecEnr vs edge prop difference", status="primary",
                    width=NULL,
                    shiny::plotOutput("plot_SpecEnrVSpropdiff", height=350) %>%
                        shinycssloaders::withSpinner(
                            type=6, size=.5, proxy.height="250px")
                ),
                shinydashboard::box(
                    title="raw vs edge prop difference", status="primary",
                    width=NULL,
                    shiny::plotOutput("plot_rawVSpropdiff", height=350) %>%
                        shinycssloaders::withSpinner(
                            type=6, size=.5, proxy.height="250px")
                ),
                shinydashboard::box(
                    title="edge SpecEnr unadjusted p-value vs raw",
                    status="primary",
                    width=NULL,
                    shiny::plotOutput("plot_edgeSpecEnrVSraw", height=350) %>%
                        shinycssloaders::withSpinner(
                            type=6, size=.5, proxy.height="250px")
                )
            )
        )
    ),
    shiny::fluidRow(
        shinydashboard::box(
            title="SpecEnr vs raw p-value difference on edge",
            status="primary",
            width=12,
                shiny::column(
                    6,
                    actionButton("reset", label = "Reset selection"),
                    ggiraph::ggiraphOutput(
                        "plot_SpecEnrVSraw", width="100%", height="700px") %>%
                        shinycssloaders::withSpinner(
                            type=6, size=.5, proxy.height="250px")
                ),
                shiny::column(
                    6,
                    shiny::tableOutput("datatab")
                )
        ),
        shinydashboard::tabBox(
            title="SpecEnr filtered ^",
            # status="primary",
            width=6,
            shiny::tabPanel(
                "sig>sig",
                    ggiraph::ggiraphOutput(
                        "plot_SpecEnrVSraw1", width="100%", height="700px") %>%
                        shinycssloaders::withSpinner(
                            type=6, size=.5, proxy.height="250px")
            ),
            shiny::tabPanel(
                "sig>insig",
                ggiraph::ggiraphOutput(
                    "plot_SpecEnrVSraw2", width="100%", height="700px") %>%
                    shinycssloaders::withSpinner(
                        type=6, size=.5, proxy.height="250px")
            ),
            shiny::tabPanel(
                "insig>sig",
                ggiraph::ggiraphOutput(
                    "plot_SpecEnrVSraw3", width="100%", height="700px") %>%
                    shinycssloaders::withSpinner(
                        type=6, size=.5, proxy.height="250px")
            ),
            shiny::tabPanel(
                "insig>insig",
                ggiraph::ggiraphOutput(
                    "plot_SpecEnrVSraw4", width="100%", height="700px") %>%
                    shinycssloaders::withSpinner(
                        type=6, size=.5, proxy.height="250px")
            )
        ),
        shinydashboard::tabBox(
            title="raw filtered ^",
            # status="primary",
            width=6,
            shiny::tabPanel(
                "sig>sig",
                ggiraph::ggiraphOutput(
                    "plot_SpecEnrVSraw5", width="100%", height="700px") %>%
                    shinycssloaders::withSpinner(
                        type=6, size=.5, proxy.height="250px")
            ),
            shiny::tabPanel(
                "sig>insig",
                ggiraph::ggiraphOutput(
                    "plot_SpecEnrVSraw6", width="100%", height="700px") %>%
                    shinycssloaders::withSpinner(
                        type=6, size=.5, proxy.height="250px")
            ),
            shiny::tabPanel(
                "insig>sig",
                ggiraph::ggiraphOutput(
                    "plot_SpecEnrVSraw7", width="100%", height="700px") %>%
                    shinycssloaders::withSpinner(
                        type=6, size=.5, proxy.height="250px")
            ),
            shiny::tabPanel(
                "insig>insig",
                ggiraph::ggiraphOutput(
                    "plot_SpecEnrVSraw8", width="100%", height="700px") %>%
                    shinycssloaders::withSpinner(
                        type=6, size=.5, proxy.height="250px")
            )

        )
    ),
    shiny::fluidRow(
        shiny::column(
            3,
            shinydashboard::box(
                title="3. cell population information",
                status="primary", solidHeader=TRUE,
                # collapsible=TRUE,
                width=NULL,
                shiny::fluidRow(
                    shiny::column(
                        width=12,
                        shiny::textInput(
                            inputId="type_cpop", label="enter a cell opulation (not edge)",
                            value=""
                        ),
                        style="text-align: center"
                    ),
                    shiny::column(
                        width=12,
                        shiny::actionButton(inputId='go_cpop', label="go"),
                        style="text-align: center"
                    ),
                    shiny::column(
                        width=12,
                        shiny::uiOutput("cpop_info_ui") %>%
                            shinycssloaders::withSpinner(
                                type=6, size=.5, proxy.height="250px")
                    )
                )
            )
        ),
        shiny::column(
            9,
            shinydashboard::box(
                title="boxplots",
                status="primary",
                # collapsible=TRUE,
                width=NULL,
                shiny::fluidRow(
                    shiny::column(
                        4, shiny::checkboxInput(
                            inputId="dotplot", value=TRUE, label="dotplot")),
                    shiny::column(
                        4, shiny::checkboxInput(
                            inputId="outlier", value=TRUE, label="outliers")),
                    shiny::column(
                        4, shiny::checkboxInput(
                            inputId="paired", value=FALSE, label="pair samples"))
                ),
                shiny::fluidRow(
                    shiny::column(
                        width=4,
                        shiny::plotOutput("plot_box")
                    ),
                    shiny::column(
                        width=4,
                        shiny::plotOutput("plot_box1")
                    ),
                    shiny::column(
                        width=4,
                        shiny::plotOutput("plot_box2")
                    )
                )
            )
        )
    )
)


shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(disable=TRUE),
    shinydashboard::dashboardSidebar(disable=TRUE), body)

