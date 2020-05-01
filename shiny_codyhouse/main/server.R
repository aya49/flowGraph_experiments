# javascript function to control progress bar during network layouting
# from https://nz-stefan.shinyapps.io/cran-explorer/

sidel <- 350 # sidelength of plots
root_ <- "root"

server <- function(input, output, session) {
    values <- shiny::reactiveValues(root=root_)



    #### about ####
    output$about_ui <- renderUI({
        tagList(
            h1("flowGraph", style="align: center;"),
            p("flowGraph uses the SpecEnr feature to analyze differentially abundant cell populations across flow cytometry (FCM) samples of different phenotypic origins. SpecEnr negates the influence of overlapping cell populations on each other hence making independent analysis and finding of true differentially abundant cell populations possible via current statistical tests. See the flowGraph ", a(href="https://github.com/aya49/flowGraph", "package"), " and ", a(href="https://www.biorxiv.org/content/10.1101/837765v1", "paper"), " for more details.")
        )
    })



    #### input fg ####
    fg_obj <- shiny::reactive({
        fg
        # shiny::req(input$fg)
        # if (grepl("[.]rds",input$fg$datapath,ignore.case=TRUE))
        #     fg <- readRDS(input$fg$datapath)
        # if (grepl("[.]rdata",input$fg$datapath,ignore.case=TRUE))
        #     fg <- get(load(input$fg$datapath))
        # fg
    })

    # fg description
    output$fg_descript_ui <- shiny::renderUI({
        shiny::req(fg_obj())
        fg <- fg_obj()

        shiny::HTML(paste0(
            "<h4>flowGraph file contains</h4>",
            "<ul>",
            "<li>marker/gate(s): ",paste0(fg@markers,collapse=", "),"</li>",
            "<li>",length(fg@feat$node), " node and ",
            length(fg@feat$edge)," edge feature(s)</li>",
            "<li>",nrow(fg@meta)," samples x ",
            nrow(fg@graph$v)," cell population nodes and ",
            nrow(fg@graph$e)," edges.","</li>",
            "</ul>"
        ))
    })



    #### data table ####
    ip_p_thres <- reactive({ input$p_thres })
    ip_filter_adjust0 <- reactive({ input$adjust0 })
    ip_filter_es <- reactive({
        switch(input$effect_size, negligible=1, small=2, medium=3, large=4)
    })

    # number of significant cell pops in each summary statistic
    dt_signo <- reactive({
        shiny::req(fg_obj())
        fg <- fg_obj()
        filter_adjust0 <- ip_filter_adjust0()
        filter_es <- ip_filter_es()
        p_thres <- ip_p_thres()

        sdt <- fg@summary_desc$node
        signo <- purrr::map_int(seq_len(nrow(sdt)), function(x) {
            pp <- fg_get_summary(fg, index=x, filter_adjust0=filter_adjust0,
                                 filter_es=filter_es)
            sum(pp$values<p_thres)
        })
    })

    # summary statistics with >0 sig cell pops
    dt_indices <- reactive({
        req(dt_signo())
        signo <- dt_signo()

        which(signo>0)
    })
    dt_obj <- reactive({
        shiny::req(dt_indices())
        fg <- fg_obj()
        indices <- dt_indices()

        fg@summary_desc$node[indices,]
    })
    output$table_summary <- DT::renderDataTable({
        shiny::req(dt_obj())
        sdt <- dt_obj()

        DT::datatable(
            sdt, selection="single", style="bootstrap4", rownames=FALSE,
            options=list(pageLength=5, columnDefs=list(list(
                searchable=FALSE, columns.orderable=FALSE))))
    }, server=TRUE)



    #### GO: summary statistic selected ####
    # index of selected summary statistic,
    # dt_index() used to regulate outputs that depend on summary statistic
    dt_ind <- shiny::reactive({ input$table_summary_rows_selected })
    shiny::observeEvent(input$go, {
        if (is.null(dt_ind()))
            shinyWidgets::sendSweetAlert(
                session=session,
                type="warning",
                title="no summary statistic selected",
                text="please select a row from the summary statistic table",
                btn_labels="ok"
            )
        shiny::req(dt_ind())
        print("GO")

        values$go <- TRUE

        values$fg <- fg <- fg_obj()
        values$p_thres <- ip_p_thres()
        values$filter_adjust0 <- filter_adjust0 <- ip_filter_adjust0()
        values$filter_es <- filter_es <- ip_filter_es()

        values$index <- index <- dt_indices()[tail(dt_ind(),1)]
        values$pp <- pp <- fg_get_summary(
            fg, index=index,
            filter_adjust0=filter_adjust0, filter_es=filter_es)

        values$max_nodes_sig <- dt_signo()[index]

        # if feature is not SpecEnr, this is NONE,
        # else it is the actual and expected feature names
        nl <- "NONE"
        feat_ <- unlist(fg@summary_desc$node$feat[index])
        if (grepl("SpecEnr", feat_)) {
            if (feat_=="SpecEnr") {
                nl <- c("prop","expect_prop")
            } else {
                nl1 <- gsub("SpecEnr_", "", feat_)
                nl <- c(nl1, paste0("expect_", nl1))
            }
        }
        values$node_labels <- nl

        root_ <- input$root_name
        if (!is.null(root_))
            if (root_!="")
                values$root <- root_

        cpopi <- which.min(pp$values)
        values$cpop_ <- values$cpop <- cpop <- names(pp$values)[cpopi]
        if (values$cpop=="") values$cpop_ <- values$root
        values$cpop_min_ <- values$cpop_
        values$cpop_l <- fg@graph$v$phenolayer[cpopi]
    })
    show_boxes <- shiny::reactive({
        shiny::req(values$node_labels)
        values$node_labels[1]!="NONE"
    })
    # output$sum_stat_ui <- shiny::renderUI({
    #     shiny::req(values$cpop)
    #     fg <- values$fg
    #     index <- values$index
    #
    #     sm <- unlist(fg@summary_desc$node[index,])
    #     shiny::tagList(
    #         shiny::h3(paste0("summary statistics for: ", sm["feat"]))
    #     )
    # })




    #### qq plot ####
    output$qq_opts_ui <- shiny::renderUI({
        shiny::req(values$go)
        print("qq options ui")
        # shiny::req(values$cpop) # my way of not using conditionalPanel
        tags$div(
            shiny::checkboxInput(inputId="logged", value=FALSE, label="logged")
        )
    })
    output$plot_qq <- shiny::renderPlot({
        shiny::req(values$index)
        fg <- values$fg
        index <- values$index
        filter_adjust0 <- values$filter_adjust0
        filter_es <- values$filter_es
        logged <- input$logged; if (length(logged)==0) logged <- FALSE
        print(paste("qq: ", index, filter_adjust0, filter_es, logged))

        plot(fg_plot_qq(fg, index=index, main="", logged=logged,
                        filter_adjust0=filter_adjust0, filter_es=filter_es)) #shiny_plot=TRUE,
    })



    #### p vs difference plot ####
    output$plot_pVSdifference <- shiny::renderPlot({
        shiny::req(values$index)
        fg <- values$fg
        index <- values$index
        filter_adjust0 <- values$filter_adjust0
        filter_es <- values$filter_es
        print(paste("pvd: ", index, filter_adjust0, filter_es))

        plot(fg_plot_pVSdiff(
            fg, index=index, main="", label_max=0,
            filter_adjust0=filter_adjust0, filter_es=filter_es)) #shiny_plot=TRUE, nodes_max=Inf)
    })



    #### cell hierarchy plot ####
    output$vn_slider_ui <- shiny::renderUI({
        shiny::req(values$max_nodes_sig)
        max_n <- values$max_nodes_sig
        print(paste0("max nodes: ",max_n))
        shiny::sliderInput(
            inputId="max_nodes", label="max # of phenotypes to show",
            min=1, max=values$max_nodes_sig, step=5,
            value=min(20, max_n)
        )
    })
    ip_max_nodes <- shiny::reactive({ input$max_nodes })
    dt_graph <- shiny::reactive({
        shiny::req(ip_max_nodes())
        max_nodes_plot <- ip_max_nodes()
        fg <- values$fg
        index <- values$index
        p_thres <- values$p_thres
        filter_adjust0 <- values$filter_adjust0
        filter_es <- values$filter_es
        node_labels <- values$node_labels
        root <- values$root
        print(paste("gr: ", index, max_nodes_plot, p_thres, root, paste0(node_labels, collapse="/"), filter_adjust0, filter_es))

        gr <- fg_plot(fg, index=index, p_thres=p_thres,
                      label_max=max_nodes_plot, node_labels=node_labels,
                      filter_adjust0=filter_adjust0, filter_es=filter_es)
        gr$v$v_ind <- gr$v$label_ind
        gr$v$phenotype[gr$v$phenotype==""] <- gr$e$from[gr$e$from==""] <- root

        print(paste("gr again: ", sum(gr$v$label_ind), max_nodes_plot))
        gr
    })
    output$graph <- visNetwork::renderVisNetwork({
        shiny::req(dt_graph())
        gr <- dt_graph()

        plot_gr(gr, shiny_plot=TRUE)
    })



    ### cell pop search and go ###
    # textbox
    output$type_cpop_ui <- shiny::renderUI({
        shiny::req(values$pp)
        # gr <- dt_graph()
        pp <- values$pp
        cpop_m <- values$cpop_min_

        # head(gr$v)
        # gr$v$displayed <- ifelse(gr$v$v_ind, "shown", "hidden")
        # gr$v$label_shiny <- paste0("p=", signif(pp$values,3), " (effect-size: ",as.character(pp$cohensd_size),") ", signif(pp$m1,3), " vs ", signif(pp$m2,3))
        # shinysky::textInput.typeahead(
        #     id="type_cpop", placeholder="enter a cell population",
        #     local=gr$v, valueKey = "phenotype",
        #     tokens=c(),
        #     template = HTML("<p class='repo-language'>{{displayed}}</p> <p class='repo-name'>{{phenotype}}</p> <p class='repo-description'>{{label_shiny}}</p>")
        # )
        # shiny::tagAppendAttributes(
            shiny::textInput(
                inputId="type_cpop", label="enter a cell opulation",
                value=cpop_m
            )
            # , `data-proxy-click` = "go_cpop"
        # )
    })
    output$go_cpop_ui <- shiny::renderUI({
        shiny::req(values$pp)
        shiny::actionButton(inputId='go_cpop', label="go")
    })
    # fill textbox with node chosen in visnetwork
    shiny::observeEvent(input$graph_selected, {
        # shiny::req(input$graph_selected)
        shiny::updateTextInput(
            session, inputId="type_cpop", value=input$graph_selected)
    })



    #### GO: cell pop selected selected_cpop() ####
    shiny::observeEvent(input$go_cpop, {
        shiny::req(dt_graph())
        print("GO CPOP")
        gr <- dt_graph()
        root <- values$root
        text_cpop <- text_cpop_ <- input$type_cpop # init = NULL; empty = length()==1 ""

        if (!is.null(text_cpop_)) {
            if (text_cpop==root) text_cpop <- ""
            text_cpopi <- which(gr$v$phenotype==text_cpop)
            if (text_cpopi) {
                values$cpop <- text_cpop
                values$cpop_ <- text_cpop_
                values$cpop_l <- gr$v$phenolayer[text_cpopi]

                # select node in visnetwork with text in textbox
                if (gr$v$v_ind[text_cpopi]) {
                    shiny::isolate({
                        visNetwork::visNetworkProxy("plot_visnet") %>%
                            visNetwork::visSelectNodes(id=text_cpop)
                    })
                }
            } else {
                shinyWidgets::sendSweetAlert(
                    session=session,
                    type="warning",
                    title="cell population not found!",
                    text="try again :)",
                    html=TRUE, btn_labels="close", closeOnClickOutside=TRUE,
                    showCloseButton=TRUE
                )
            }
        }
    })

    # cell pop description
    output$cpop_info_ui <- shiny::renderUI({
        shiny::req(values$cpop)
        fg <- values$fg
        index <- values$index
        cpop <- values$cpop
        cpop_ <- values$cpop_
        node_labels <- values$node_labels
        pp <- values$pp

        cpopi <- which(names(pp$values)==cpop)
        sm <- unlist(fg@summary_desc$node[index,])
        m0a <- signif(pp$m1[cpopi],3)
        m0b <- signif(pp$m2[cpopi],3)

        print(paste("cpop INFO: ", index, cpop, cpop_, m0a, m0b, paste0(node_labels, collapse="/"), cpopi, paste0(sm, collapse="/")))

        fm <- ""
        if (show_boxes()) {
            m1a <- signif(mean(as.matrix(fg@feat$node[[node_labels[1]]])[pp$id1,cpopi]),3)
            m1b <- signif(mean(as.matrix(fg@feat$node[[node_labels[1]]])[pp$id2,cpopi]),3)
            m2a <- signif(mean(as.matrix(fg@feat$node[[node_labels[2]]])[pp$id1,cpopi]),3)
            m2b <- signif(mean(as.matrix(fg@feat$node[[node_labels[2]]])[pp$id2,cpopi]),3)
            fm <- paste0(" (actual/expect: ",m1a,"/",m2a," vs ",m1b,"/",m2b,")")
            print(fm)
        }

        shiny::HTML(paste0(
            "<h4>",sm["feat"],": ",cpop_,"</h4>",
            "<ul>",
            "<li>p-value: ",signif(pp$values[cpopi],3),"</li>",
            "<li>stat test: ",sm["test_name"],"</li>",
            "<li>",sm["class"],": ", m0a," vs ",m0b,"</li>",
            "<li>feat means: ",m0a," vs ",m0b,fm,"</li>",
        "</ul>"
        ))

    })



    #### plot boxplot ####
    output$box_opts_ui <- shiny::renderUI({
        shiny::req(values$index)
        print("boxplot options ui")
        shiny::req(values$go) # my way of not using conditionalPanel

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
            )
    })
    ip_paired <- shiny::reactive({ input$paired })
    ip_dotplot <- shiny::reactive({ input$dotplot })
    ip_outlier <- shiny::reactive({ input$outlier })
    # output$show_multibox <- shiny::reactive({ show_boxes() })
    # shiny::outputOptions(output, "show_multibox", suspendWhenHidden=FALSE)

    output$plot_box <- shiny::renderPlot({
        shiny::req(values$cpop)
        fg <- values$fg
        index <- values$index
        cpop <- values$cpop
        paired <- ip_paired()
        dotplot <- ip_dotplot()
        outlier <- ip_outlier()
        print(paste("box: ", cpop, paired, dotplot, outlier))
        head(show(fg))

        fg_plot_box(fg, type="node", index=index, node_edge=cpop, main="",
                         paired=paired, outlier=outlier, dotplot=dotplot)
    })

    output$plot_box1 <- shiny::renderPlot({
        shiny::req(values$cpop, show_boxes())
        fg <- values$fg
        index <- values$index
        cpop <- values$cpop
        paired <- ip_paired()
        dotplot <- ip_dotplot()
        outlier <- ip_outlier()
        node_labels <- values$node_labels
        print(paste("box1: ", cpop, paired, dotplot, outlier, node_labels[1]))

        fg_plot_box(
            fg, type="node", index=index, node_edge=cpop,
            main=node_labels[1],
            paired=paired, outlier=outlier, dotplot=dotplot,
            feature=node_labels[1])
    })

    output$plot_box2 <- shiny::renderPlot({
        shiny::req(values$cpop, show_boxes())
        fg <- values$fg
        index <- values$index
        cpop <- values$cpop
        paired <- ip_paired()
        dotplot <- ip_dotplot()
        outlier <- ip_outlier()
        node_labels <- values$node_labels
        print(paste("box2: ", cpop, paired, dotplot, outlier, node_labels[2]))

        fg_plot_box(
            fg, type="node", index=index, node_edge=cpop,
            main=node_labels[2],
            paired=paired, outlier=outlier, dotplot=dotplot,
            feature=node_labels[2])
    })

    output$box_plots_ui <- shiny::renderUI({
        shiny::div(
            shiny::plotOutput(outputId="plot_box"),
            shiny::conditionalPanel(
                condition=show_boxes(),
                shiny::plotOutput(outputId="plot_box1"),
                shiny::plotOutput(outputId="plot_box2")
            )
        )
    })
    output$box_click_ui <- shiny::renderUI({
        shiny::req(show_boxes())
        shiny::div(
            shiny::actionButton(inputId="click_box", label="show all boxplots")
        )
    })
    # shiny::observeEvent(input$click_box, {
    #     shiny::req(values$cpop)
    #     shinyWidgets::sendSweetAlert(
    #         session=session,
    #         title=paste0("SpecEnr + actual/expected boxplots: ", values$cpop_),
    #         text=tags$div(
    #             shiny::plotOutput("plot_box"),
    #             shiny::plotOutput("plot_box1"),
    #             shiny::plotOutput("plot_box2")
    #             #     shiny::plotOutput(
    #             #     outputId="plot_box",
    #             #     height=paste0(sidel,"px"), width=paste0(sidel,"px")),
    #             # shiny::plotOutput(
    #             #     outputId="plot_box1",
    #             #     height=paste0(sidel,"px"), width=paste0(sidel,"px")),
    #             # shiny::plotOutput(
    #             #     outputId="plot_box2",
    #             #     height=paste0(sidel,"px"), width=paste0(sidel,"px"))
    #         ),
    #         html=TRUE, btn_labels="close", width=paste0(sidel+40,"px")
    #     )
    # })




    #### plot scatterplots ####
    shinyFiles::shinyDirChoose(input, "fcs_dir", roots=start_dir, session=session) # fcs
    shinyFiles::shinyDirChoose(input, "ft_dir", roots=start_dir, session=session) # flowType
    output$fcs_dir_out <- shiny::renderPrint({
        if (is.integer(input$fcs_dir)) {
            cat("directory: FCS files (.fcs)")
        } else {
            parseDirPath(start_dir, input$fcs_dir)
        }
    })
    output$ft_dir_out <- renderPrint({
        if (is.integer(input$ft_dir)) {
            cat("directory: flowType files (.Rdata)")
        } else {
            parseDirPath(start_dir, input$ft_dir)
        }
    })
    ip_fcs_paths <- shiny::reactive({
        # shiny::req(input$fcs_dir)
        # fcs_dir <- shinyFiles::parseDirPath(start_dir, input$fcs_dir)
        # list.files(fcs_dir, full.names=T, pattern=".fcs", recursive=TRUE)
        list.files(fcs_dir, full.names=T, pattern=".fcs", recursive=TRUE)
    })
    ip_ft_paths <- shiny::reactive({
        # shiny::req(input$ft_dir)
        # ft_dir <- shinyFiles::parseDirPath(start_dir, input$ft_dir)
        # list.files(ft_dir, full.names=T, pattern=".fcs", recursive=TRUE)
        list.files(ft_dir, full.names=T, pattern=".fcs", recursive=TRUE)
    })

    gs_reactive <- shiny::reactive({
        shiny::req(values$cpop, ip_fcs_paths(), ip_ft_paths())
        index <- values$index
        cpop <- values$cpop
        pp <- values$pp
        fcs_paths <- ip_fcs_paths()
        ft_paths <- ip_ft_paths()
        print(paste("gs: ", cpop, fcs_paths[1], ft_paths[1]))

        # load fcs.Rdata or matrix m to choose median sample
        m <- fg@feat$node[[fg@summary_desc$node$feat[[index]]]]

        # choose median samples from each class label
        label_ind1 <- fg@meta$id[pp$id1]
        f_id1 <- label_ind1[which.min(abs(m[pp$id1,cpop]-median(m[pp$id1,cpop])))]
        label_ind2 <- fg@meta$id[pp$id2]
        f_id2 <- label_ind2[which.min(abs(m[pp$id2,cpop]-median(m[pp$id2,cpop])))]

        # load flowtype and fcs files
        # CHECK THIS, multiple matches
        fcs1 <- flowCore::read.FCS(fcs_paths[grep(f_id1, fcs_paths)[1]])
        fcs2 <- flowCore::read.FCS(fcs_paths[grep(f_id2, fcs_paths)[1]])

        ft1 <- get(load(ft_paths[grep(f_id1, ft_paths)[1]]))
        ft2 <- get(load(ft_paths[grep(f_id2, ft_paths)[1]]))

        # find best path
        igr <- fg@graph$e[,c("from","to")]
        igr <- igraph::graph_from_data_frame(igr)
        e_weights <- log(pp$values)[fg@graph$e$to]-min(log(pp$values))+1
        gpp <- igraph::V(igr)$name[igraph::shortest_paths(
            igr, from="", to=cpop, output="vpath",
            weights=e_weights)$vpath[[1]][-1]]
        gating_path_phenocode <-
            fg@graph$v$phenocode[match(gpp, fg@graph$v$phenotype)]
        # rchy <- RchyOptimyx::RchyOptimyx(
        #     pheno.codes=fg@graph$v$phenocode,
        #     phenotypeScores=-log(fg@summary$node[[index]]$values),
        #     startPhenotype=fg@graph$v$phenocode[fg@graph$v$phenotype==cpop],
        #     pathCount=1, trimPaths=FALSE)
        # gating_path_phenocode <- rchy@nodes[1,-1]

        # make gating set
        gates <- list(ft1@Thresholds, ft2@Thresholds)
        gs <- CytoML::GatingSet(as(list(fcs1, fcs2),"flowSet"))
        sampleNames(gs) <- names(gates) <- c(f_id1, f_id2)

        # leave only one threshold for each marker,,, because,,,
        # i'm not sure why there are two thresholds lol
        for (mi in seq_len(length(gates[[1]])))
            if (length(gates[[1]][[mi]])==2)
                gates[[1]][[mi]] <- gates[[1]][[mi]][1]

        scat_chans <- unlist(sapply (c("SSC-A","FSC-A"), function(x)
            grep(x, colnames(fcs1@exprs), value=T)))

        make_gs(gs, gates, markers=fg@markers,
                gating_path_phenocode, scat_chans[1])
    })

    # (Note: The hashing algorithm used in renderCachedPlot is xxHash-64.)
    # d <- data.frame(x = rnorm(400), y = rnorm(400))
    # system.time(digest::digest(d, "xxhash64"))
    output$plot_dens1 <- shiny::renderCachedPlot({
    # output$plot_dens1 <- shiny::renderPlot({
        shiny::req(gs_reactive())
        gs <- gs_reactive()
        flowWorkspace::plotGate(gs[[1]], xbin=70, gpar=list(ncol=1))
    }, cacheKeyExpr={ list("plot_dens1") })
    # })
    output$plot_dens2 <- shiny::renderCachedPlot({
    # output$plot_dens2 <- shiny::renderPlot({
        shiny::req(gs_reactive())
        gs <- gs_reactive()
        flowWorkspace::plotGate(gs[[2]], xbin=70, gpar=list(ncol=1))
    # })
    }, cacheKeyExpr={ list("plot_dens2") })

    output$click_sp_ui <- shiny::renderUI({
        shiny::req(values$go, ip_fcs_paths(), ip_ft_paths())
        shiny::actionButton(inputId="click_sp", label="show scatterplots")
    })
    shiny::observeEvent(input$click_sp, {
        shiny::req(gs_reactive())
        clayer <- values$cpop_l
        cpop_ <- values$cpop_

        shinyWidgets::sendSweetAlert(
            session=session,
            title=paste0("scatterplots: ", cpop_),
            text=tags$div(
                shiny::splitLayout(
                    cellWidths=rep(paste0(sidel,"px"),2),
                    shiny::plotOutput(
                        outputId="plot_dens1",
                        height=paste0(sidel*clayer,"px"), width=paste0(sidel,"px")),
                    shiny::plotOutput(
                        outputId="plot_dens2",
                        height=paste0(sidel*clayer,"px"), width=paste0(sidel,"px"))
                )
            ), html=TRUE, showCloseButton=TRUE,
            btn_labels="close", width=paste0(sidel*2+40,"px")
        )
    })

    # session$onSessionEnded(stopApp) # stop session when window closed
}
