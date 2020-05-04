material_card <- function(style_more="", ...) {
    shiny::div(..., style=paste0(
        "background-color: white; ",
        "width: calc(100vw - 40px); ",
        "height: calc(100vh - 40px); ",
        "padding: 10px; ", #inside
        # "margin: 10px; ", #outside
        "border: 1px solid silver; ",
        "box-shadow: 0 0 18px #888888; ",
        "font-size: 8pt; ", style_more)
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

    # top box
    material_card(
        style_more="margin: 10px",
        shiny::fluidRow(

            # left column
            shiny::column(
                width=7,

                # load fg file
                shiny::fluidRow(
                    shiny::column(
                        width=12,
                        shiny::uiOutput("about_ui")
                    ),
                    style="padding: 10px"
                ),
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
                                width=4,
                                shinyFiles::shinyDirButton(id="fcs_dir","select","select")
                            ),
                            shiny::column(
                                width=8,
                                shiny::verbatimTextOutput("fcs_dir_out")
                            )
                        ),
                        shiny::fluidRow(
                            shiny::column(
                                width=4,
                                shinyFiles::shinyDirButton(id="ft_dir","select","select")
                            ),
                            shiny::column(
                                width=8,
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
                    style="border-bottom: solid silver 1px; padding: 10px"
                ),
                shiny::column(
                    width=12,
                    shiny::uiOutput("dt_ui"),
                    style="font-size: 8pt; margin-bottom: 15px; margin-top: 10px; "
                ),
                # go! make hierarchy and summary plots
                shiny::column(
                    width=12,
                    shiny::actionButton(
                        inputId="go",
                        label="generate/refresh plots"
                    ),
                    tags$style(".swal-modal {width: 80%;}"),
                    style="text-align: center; margin-bottom: 15px; padding-top: 10px"
                )
            ),

            # right column
            shiny::column(
                width=5,
                shiny::fluidRow(
                    # shiny::column(
                    #     width=12,
                    #     shiny::uiOutput("sum_stat_ui"),
                    #     style="border-bottom: solid silver 1px"
                    # ),
                    shiny::column(
                        width=12,
                        tags$h4("QQ plot"),
                        style="margin-top: 15px;"
                    ),
                    shiny::column(
                        width=12,
                        shiny::uiOutput("qq_opts_ui")
                    ),
                    shiny::column(
                        width=12,
                        shiny::plotOutput("plot_qq", height=350) %>%
                            shinycssloaders::withSpinner(
                                type=6, size=.5, proxy.height="250px"),
                        style="border-bottom: solid silver 1px; margin-bottom: 20px; "
                    ),
                    shiny::column(
                        width=12,
                        tags$h4("p-value vs mean feat difference plot"),
                        style="margin-top: 5px;"
                    ),
                    shiny::column(
                        width=12,
                        shiny::plotOutput("plot_pVSdifference", height=350) %>%
                            shinycssloaders::withSpinner(
                                type=6, size=.5, proxy.height="250px"),
                        style="margin-bottom: 20px;"
                    )
                ),
                style="border-left: solid silver 1px"
            )
        )
    ),

    # main box
    material_card(
        style_more="margin: 10px; margin-top: 0px; ",
        shiny::fluidRow(

            # left column
            shiny::column(
                width=8,
                # visnet cell hierarchy plot
                visNetwork::visNetworkOutput(
                    "graph", width="100%", height="100vh-100px") %>%
                    shinycssloaders::withSpinner(
                        type=6, size=.5, proxy.height="250px"),
                style="border-right: solid silver 1px; "
            ),

            # right column
            shiny::column(
                width=4,
                shiny::fluidRow(
                    shiny::column(
                        width=12,
                        shiny::uiOutput('vn_slider_ui'),
                        style="text-align: center; margin: 10px; padding: 10px; border-bottom: solid silver 1px; "
                    ),
                    shiny::column(
                        width=12,
                        shiny::uiOutput('type_cpop_ui'),
                        style="text-align: center; margin: 10px"
                    ),
                    shiny::column(
                        width=12,
                        shiny::uiOutput('go_cpop_ui'),
                        style="text-align: center; padding-bottom: 15px; border-bottom: solid silver 1px; " # padding-top: 100px;
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
                        style="padding: 10px; margin: 10px; padding-bottom: 15px; text-align: center; border-bottom: solid silver 1px"
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
                )
            )
        )
    )
)

shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(disable=TRUE),
    shinydashboard::dashboardSidebar(disable=TRUE), body)


