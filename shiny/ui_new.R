material_card <- function(..., header=NULL, bgcolor="white") {
    div(class="card", header,
        div(class="card-content", ...,
            style=sprintf("background-color: %s", bgcolor)
        )
    )
}

body <- shinydashboard::dashboardBody(
    # title="testing!",

    # tags$script(HTML("$('body').addClass('fixed');")), # fix header and sidebar


    # htmlTemplate(
    #     filename <- "www/index.html",
    #     cell_hierarchy_ui <- function(id) {
    #         ns <- shiny::NS(id)
    #

    # tags$head(tags$script(HTML(button_enter))),

    material_card(

        # top box
        shiny::fluidRow(

            # left column
            shiny::column(
                width=3,

                # load fg file
                shiny::fluidRow(
                    shiny::column(
                        width=12,
                        shiny::fileInput(
                            inputId="fg",
                            label="select a flowGraph file (.Rdata/.rds)",
                            buttonLabel="browse", multiple=FALSE,
                            accept=c("Rdata/RDS", "Rdata/RDS", "*")
                        ),
                        style="font-size: 8pt; border-bottom: solid silver 1px"
                    ),
                    shiny::column(
                        width=12,
                        shiny::uiOutput("fg_descript_ui")
                    )
                ),
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
                ),
                style="border-right: solid silver 1px"
            ),

            # middle column
            shiny::column(
                width=6,
                shiny::fluidRow(

                    # fg data frame
                    shiny::column(
                        width=12,
                        DT::dataTableOutput("table_summary"),
                        style="border-bottom: solid silver 1px"
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
                ),
                style="font-size: 8pt"
            ),

            # right column (options)
            shiny::column(
                width=3,
                shiny::fluidRow(
                    shiny::column(
                        width=12,
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
                    )
                ),
                style="font-size: 8pt; border-right: solid silver 1px"
            )
        )
    ),

    material_card(

        # main box
        shiny::fluidRow(

            # left column
            shiny::column(
                width=3,
                shiny::fluidRow(
                    # shiny::column(
                    #     width=12,
                    #     shiny::uiOutput("sum_stat_ui"),
                    #     style="border-bottom: solid silver 1px"
                    # ),
                    shiny::column(
                        width=12,
                        tags$h4("QQ plot")
                    ),
                    shiny::column(
                        width=12,
                        shiny::uiOutput("qq_opts_ui")
                    ),
                    shiny::column(
                        width=12,
                        shiny::plotOutput("plot_qq", height=250) %>%
                            shinycssloaders::withSpinner(
                                type=6, size=.5, proxy.height="250px"),
                        style="border-bottom: solid silver 1px"
                    ),
                    shiny::column(
                        width=12,
                        tags$h4("p-value vs mean feat difference plot")
                    ),
                    shiny::column(
                        width=12,
                        shiny::plotOutput("plot_pVSdifference", height=250) %>%
                            shinycssloaders::withSpinner(
                                type=6, size=.5, proxy.height="250px")
                    )
                ),
                style="border-right: solid silver 1px; border-top: solid silver 1px"
            ),

            # center column
            shiny::column(
                width=6,
                shiny::fluidRow(
                    shiny::column(
                        width=12,
                        shiny::uiOutput('vn_slider_ui'),
                        style="text-align: center"
                    ),
                    # visnet cell hierarchy plot
                    shiny::column(
                        width=12,
                        shiny::div(
                            visNetwork::visNetworkOutput(
                                "graph", width="100%") %>%
                                shinycssloaders::withSpinner(
                                    type=6, size=.5, proxy.height="250px"),
                            style="margin-top: 10px"
                        )
                    ),
                    style="border-top: solid silver 1px"
                )
            ),

            # right column
            shiny::column(
                width=3,
                shiny::fluidRow(
                    shiny::column(
                        width=12,
                        shiny::uiOutput('type_cpop_ui'),
                        style="text-align: center"
                    ),
                    shiny::column(
                        width=12,
                        shiny::uiOutput('go_cpop_ui'),
                        style="text-align: center; border-bottom: solid silver 1px; padding-bottom: 10px" # padding-top: 100px;
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
                        style="padding-top: 5px; padding-bottom: 5px; text-align: center; border-bottom: solid silver 1px"
                    )
                ),
                shiny::fluidRow(
                    shiny::column(
                        width=12,
                        tags$h4("Boxplots")
                    ),
                    shiny::column(
                        width=12,
                        shiny::uiOutput("box_opts_ui")
                    ),
                    shiny::column(
                        width=12,
                        shiny::uiOutput("box_plots_ui")
                    )
                    # ,
                    # shiny::column(
                    #     width=12,
                    #     shiny::uiOutput("box_click_ui"),
                    #     style="padding-top: 5px; padding-bottom: 5px; text-align: center; border-bottom: solid silver 1px"
                    # )
                ),
                style="border-left: solid silver 1px; border-top: solid silver 1px"
            )
        )
    )
)

shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(disable=TRUE),
    shinydashboard::dashboardSidebar(disable=TRUE), body)


