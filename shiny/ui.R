body <- shinydashboard::dashboardBody(
    tags$head(
        tags$link(rel="stylesheet", type="text/css", href="www/css/drawer.css")
    ),






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
                        shiny::fluidRow(
                            shiny::column(
                                width=3,
                                shinyFiles::shinyDirButton(id="fcs_dir","select","select")
                            ),
                            shiny::column(
                                width=9,
                                shiny::verbatimTextOutput("fcs_dir_out")
                            )
                        ),
                        shiny::fluidRow(
                            shiny::column(
                                width=3,
                                shinyFiles::shinyDirButton(id="ft_dir","select","select")
                            ),
                            shiny::column(
                                width=9,
                                shiny::verbatimTextOutput("ft_dir_out")
                            )
                        )
                    ),
                    shiny::column(
                        width=6,
                        shiny::textInput( # renderText({ input$root_name })
                            inputId="root_name", label="name of root phenotype",
                            value="live", width=NULL, placeholder="live"),
                        shiny::sliderInput(
                            inputId="p_thres",
                            label="p-value threshold",
                            min=.01, max=1, step=.01, value=.01),
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
                title="2. summary statistic plots",
                status="primary", solidHeader=TRUE,
                # collapsible=TRUE,
                width=NULL,

                shinydashboard::box(
                    title="QQ plot", status="primary",
                    width=NULL,
                    shiny::uiOutput("qq_opts_ui"),
                    shiny::plotOutput("plot_qq", height=350) %>%
                        shinycssloaders::withSpinner(
                            type=6, size=.5, proxy.height="250px")
                ),
                shinydashboard::box(
                    title="p-value vs mean feat difference plot",
                    status="primary",
                    width=NULL,
                    shiny::plotOutput("plot_pVSdifference", height=350) %>%
                        shinycssloaders::withSpinner(
                            type=6, size=.5, proxy.height="250px")
                )
            )
        )
    ),
    shiny::fluidRow(
        shiny::column(
            8,
            shinydashboard::box(
                title="cell hierarchy network",
                status="primary",
                width=NULL,
                visNetwork::visNetworkOutput(
                    "graph", width="100%", height="100vh") %>%
                    shinycssloaders::withSpinner(
                        type=6, size=.5, proxy.height="250px")
            )
        ),
        shiny::column(
            4,
            shinydashboard::box(
                width=NULL, status="primary",
                shiny::uiOutput('vn_slider_ui')
            ),
            shinydashboard::box(
                title="3. cell population information",
                status="primary", solidHeader=TRUE,
                # collapsible=TRUE,
                width=NULL,
                shiny::fluidRow(
                    shiny::column(
                        width=12,
                        shiny::uiOutput('type_cpop_ui'),
                        style="text-align: center"
                    ),
                    shiny::column(
                        width=12,
                        shiny::uiOutput('go_cpop_ui'),
                        style="text-align: center; padding-bottom: 10px" # padding-top: 100px;
                    ),
                    shiny::column(
                        width=12,
                        shiny::uiOutput("cpop_info_ui") %>%
                            shinycssloaders::withSpinner(
                                type=6, size=.5, proxy.height="250px")
                    ),
                    shiny::column(
                        width=12,
                        shiny::uiOutput(outputId="click_sp_ui"),
                        style="padding-top: 5px; padding-bottom: 15px; text-align: center"
                    )
                ),
                shinydashboard::box(
                    title="boxplots",
                    status="primary",
                    # collapsible=TRUE,
                    width=NULL,
                    shiny::uiOutput("box_opts_ui"),
                    shiny::uiOutput("box_plots_ui")
                    # bsplus::bs_carousel(id="box_carosel", use_indicators=TRUE) %>%
                    #     bsplus::bs_append(shiny::plotOutput(outputId="plot_box")) %>%
                    #     bsplus::bs_append(shiny::plotOutput(outputId="plot_box1")) %>%
                    #     bsplus::bs_append(shiny::plotOutput(outputId="plot_box2"))

                    # ,
                    # shiny::column(
                    #     width=12,
                    #     shiny::uiOutput("box_click_ui"),
                    #     style="padding-top: 5px; padding-bottom: 5px; text-align: center; border-bottom: solid silver 1px"
                    # )
                )
            )

        )
    )

)


shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(disable=TRUE),
    shinydashboard::dashboardSidebar(disable=TRUE), body)

