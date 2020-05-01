# Shiny APP for flowGraph visulaization
# Written by: Mehrnoush Malek (MM), Alice Yue (AY)
# Revision: March 2020, AY


server <- function(input, output, session) {

    # callModule(module=serverChangeTheme, id="moduleChangeTheme")

    #### description ####
    output$about <- renderUI({
        tagList("flowGraph uses the SpecEnr feature to analyze differentially abundant cell populations across flow ctyometry (FCM) samples of different phenotypic origins. SpecEnr negates the influence of overlapping cell populations on each other hence making independent analysis and finding of true differentially abundant cell populations possible via current statitical tests. See the flowGraph ", a(href="https://github.com/aya49/flowGraph", "package"), " and ", a(href="https://www.biorxiv.org/content/10.1101/837765v1", "paper"), " for more details.")
    })

    #### FCS & flowType folder input ####
    shinyFiles::shinyDirChoose(input, "fcs_dir", roots=start_dir, session=session) # fcs
    # shinyFiles::shinyDirChoose(input, "fg_dir", roots=start_dir, session=session) # flowGraph
    shinyFiles::shinyDirChoose(input, "ft_dir", roots=start_dir, session=session) # flowType

    output$fcs_dir_out <- shiny::renderPrint({
        if (is.integer(input$fcs_dir)) {
            cat("< select a directory: FCS files (.fcs)")
        } else {
            parseDirPath(start_dir, input$fcs_dir)
        }
    })
    output$ft_dir_out <- renderPrint({
        if (is.integer(input$ft_dir)) {
            cat("< select a directory: flowType files (.Rdata)")
        } else {
            parseDirPath(start_dir, input$ft_dir)
        }
    })


    #### summary statistics table ####
    fg_obj <- reactive({
        req(input$fg)
        if (grepl("[.]rds",input$fg$datapath,ignore.case=TRUE))
            fg <- readRDS(input$fg$datapath)
        if (grepl("[.]rdata",input$fg$datapath,ignore.case=TRUE))
            fg <- get(load(input$fg$datapath))
        return(fg)
    })

    tindex <- shiny::reactive({ input$table_summary_rows_selected })
    output$table_summary <- DT::renderDataTable({
        req(fg_obj(), input$p_thres, input$effect_size, input$adjust0)

        filter_es <- switch(input$effect_size,
                            negligible=1, small=2, medium=3, large=4)
filter_adjust0 <- input$adjust0

        fg <<- fg_obj()
        print(show(fg))

        p_thres <- input$p_thres
        print(p_thres)

        shiny::withProgress(message='loading', value=0, {
            shiny::incProgress(.3, detail="generating summary statistics table")

            # get table
            sdt <- fg@summary_desc$node
            print(nrow(sdt))
            indices <<- which(purrr::map_int(seq_len(nrow(sdt)), function(x) {
                pp <- fg_get_summary(fg, index=x, filter_es=filter_es,
                                     filter_adjust0=filter_adjust0)
                sum(pp$values<p_thres)
            })>0)

            sdt <- sdt[indices,]
            print(nrow(sdt))
            DT::datatable(sdt, selection="single",
                          style="bootstrap4",
                          options=list(pageLength=5, columns.orderable=FALSE),
                          rownames=FALSE)
            # https://stackoverflow.com/questions/28274584/get-selected-row-from-datatable-in-shiny-app
            #         callback="function(table) {
            #       table.on('click.dt', 'tr', function() {
            #         table.$('tr.selected').removeClass('selected');
            #         $(this).toggleClass('selected');
            #     });
            # }"
            #         datatable(sdt, callback=JS(callback))
        })
    }, server=TRUE)

    # define index to global and print the selected indices
    output$table_row <- shiny::renderPrint({
        shiny::req(tindex())

        ri <- tindex()
        index <<- indices[tail(ri,1)]
        # if (length(ri))
        cat('row selected:', paste0(ri, collapse=", "),
            " (", paste0(index, collapse=", "), ")", sep="")
    })


    # # load example flowtype file
    # ft <- get(load(ft_paths[1]))
    #
    # # phenotypes
    # phenotype <- sapply(
    #     ft@PhenoCodes, flowTypeFilter::decodePhenotype,
    #     marker.names=fg@markers,
    #     partitions.per.marker=ft@PartitionsPerMarker,
    #     purrr::map_int(ft@Thresholds, length)==1)
    # print(paste("sample phenotypes:", head(phenotype)))


    #### cell hierarchy plot ####

    # selected_cpop <- reactive({ input$plot_hierarchy_selected })
    # selected_cpop <- reactive({ input$plot_pVSdifference_selected })
    selected_cpop <- reactive({ input$plot_visnet_selected })

    # output$plot_hierarchy <- renderGirafe({
    #     req(input$table_summary_rows_selected, input$p_thres, input$max_nodes)
    #     p_thres <- input$p_thres
    #     max_nodes <- input$max_nodes
    #     show_bgedges <- input$show_bgedges
    #
    #     withProgress(message='loading', value=0, {
    #         incProgress(.5, detail="generating cell hierarchy of significant phenotypes")
    #
    #         ri <- input$table_summary_rows_selected
    #
    #         # row index of selected row
    #         index <<- indices[tail(ri,1)]
    #         print(index)
    #
    #         # # load fcs.Rdata or matrix m
    #         # m <<- fg@feat$node[[fg@summary_desc$node$feat[[index]]]]
    #         # print(dim(m))
    #
    #         # plot cell hierarchy
    #         gr <- fg_plot(fg, index=index, p_thres=p_thres, show_bgedges=show_bgedges, label_max=max_nodes)#, node_labels=c("cells_ul_blood", "expect_cells_ul_blood"))
    #         gr$v$v_ind <- gr$v$label_ind
    #         gr$e$e_ind <- gr$e$from%in%gr$v$phenotype[gr$v$v_ind] &
    #             gr$e$to%in%gr$v$phenotype[gr$v$v_ind]
    #         gg_hier <- plot_gr(gr, shiny_plot=TRUE, bgedges=show_bgedges,
    #                            visNet_plot=FALSE)
    #
    #         ggiraph::girafe(
    #             code=print(gg_hier), width_svg=9, height_svg=5,
    #             options=list(
    #                 opts_selection(
    #                     type="single", css="fill:#FF3333;stroke:black;"),
    #                 opts_hover(css="fill:#FF3333;stroke:black;cursor:pointer;")
    #             ))
    #     })
    # })

    output$plot_visnet <- visNetwork::renderVisNetwork({
        req(tindex(), input$p_thres, input$max_nodes, input$adjust0)
        p_thres <- input$p_thres
        max_nodes <- input$max_nodes
        filter_es <- switch(input$effect_size,
                            negligible=1, small=2, medium=3, large=4)
        filter_adjust0 <- input$adjust0
        print(p_thres)
        print(max_nodes)
        print(filter_es)


        withProgress(message='loading', value=0, {
            incProgress(.5, detail="generating visNetwork cell hierarchy of significant phenotypes")

            node_labels <<- "NONE"
            feat_ <- unlist(fg@summary_desc$node$feat[index])
            if (grepl("SpecEnr", feat_)) {
                if (feat_=="SpecEnr") {
                    node_labels <<- c("prop","expect_prop")
                }
                else {
                    nl1 <- gsub("SpecEnr_", "", feat_)
                    node_labels <<- c(nl1, paste0("expect_", nl1))
                }
            }

            gr <- fg_plot(
                fg, index=index, p_thres=p_thres, label_max=max_nodes,
                node_labels=node_labels,
                filter_es=filter_es, filter_adjust0=filter_adjust0)
            gr$v$v_ind <- gr$v$label_ind

            if (!is.null(gr)) {
                plot_gr(gr, shiny_plot=TRUE)
            }
            else {
                showModal(modalDialog(
                    title="No significant cell populations found.",
                    "Change some settings?",
                    easyClose=TRUE
                ))
            }
        })
    })


    #### p vs difference plot ####
    output$plot_pVSdifference <- renderPlot({
        req(tindex())
        # req(tindex(), input$p_thres, input$max_nodes)
        # p_thres <- input$p_thres
        # max_nodes <- input$max_nodes

        withProgress(message='loading', value=0, {
            incProgress(.5, detail="generating cell hierarchy of significant phenotypes")

            # plot
            fg_plot_pVSdiff(fg, index=index, label_max=3) #shiny_plot=TRUE, nodes_max=Inf)

            # ggiraph::girafe(
            #     code=print(gp), width_svg=9, height_svg=5,
            #     options=list(
            #         opts_selection(type="single", css="fill:#FF3333;stroke:black;"),
            #         opts_hover(css="fill:#FF3333;stroke:black;cursor:pointer;")
            #     ))
        })
    })


    #### qq plot ####
    output$plot_qq <- renderPlot({
        req(tindex(), input$p_thres, input$max_nodes)
        p_thres <- input$p_thres
        max_nodes <- input$max_nodes

        withProgress(message='loading', value=0, {
            incProgress(.5, detail="generating cell hierarchy of significant phenotypes")

            # plot
            fg_plot_qq(fg, index=index) #shiny_plot=TRUE,
        })
    })


    #### plot boxplot ####
    # selected_label <- reactive({ input$plot_box_selected })
    output$plot_box <- renderPlot({
        req(selected_cpop())
        # cpop <- input$plot_hierarchy_selected
        cpop <- selected_cpop()
        paired <- input$paired
        dotplot <- input$dotplot
        outlier <- input$outlier
        print(paired)
        print(dotplot)
        print(outlier)

        fg_plot_box(
            fg, type="node", index=index, node_edge=cpop,
            paired=paired, outlier=outlier, dotplot=dotplot)
    })

    output$plot_box1 <- renderPlot({
        req(selected_cpop())
        # cpop <- input$plot_hierarchy_selected
        cpop <- selected_cpop()
        paired <- input$paired
        dotplot <- input$dotplot
        outlier <- input$outlier
        print(paired)
        print(dotplot)
        print(outlier)
        print(node_labels[1])
        if (node_labels[1]!="NONE")
            fg_plot_box(
                fg, type="node", index=index, node_edge=cpop,
                paired=paired, outlier=outlier, dotplot=dotplot,
                feature=node_labels[1])
    })

    output$plot_box2 <- renderPlot({
        req(selected_cpop())
        # cpop <- input$plot_hierarchy_selected
        cpop <- selected_cpop()
        paired <- input$paired
        dotplot <- input$dotplot
        outlier <- input$outlier
        print(paired)
        print(dotplot)
        print(outlier)
        print(node_labels[2])
        if (node_labels[1]!="NONE")
            fg_plot_box(
                fg, type="node", index=index, node_edge=cpop,
                paired=paired, outlier=outlier, dotplot=dotplot,
                feature=node_labels[2])
    })


    #### plot scatterplots ####
    gs_reactive <- reactive({
        req(selected_cpop(), input$ft_dir, input$fcs_dir)

        cpop <- selected_cpop()
        print(cpop)
        # label <- selected_label()
        # label <- unlist(fg@summary_desc$node$label2[index])
        # print(label)
        ft_dir <-  parseDirPath(start_dir, input$ft_dir)
        print(ft_dir)
        fcs_dir <- parseDirPath(start_dir, input$fcs_dir)
        print(fcs_dir)

        shiny::withProgress(message='loading', value=0, {
            shiny::incProgress(.2, detail="loading median samples")

            fcs_paths <- list.files(fcs_dir, full.names=T, pattern=".fcs", recursive=TRUE)
            ft_paths <- list.files(ft_dir, full.names=T, pattern=".Rdata", recursive=TRUE)

            # load fcs.Rdata or matrix m to choose median sample
            m <- fg@feat$node[[fg@summary_desc$node$feat[[index]]]]
            print(dim(m))

            pp <- fg_get_summary(fg, index=index)
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

            shiny::incProgress(.4, detail="finding best gating path")

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

            shiny::incProgress(.5, detail="creating gating set")

            # make gating set
            gates <- list(ft1@Thresholds, ft2@Thresholds)
            gs <- CytoML::GatingSet(as(list(fcs1, fcs2),"flowSet"))
            sampleNames(gs) <- names(gates) <- c(f_id1, f_id2)

            # leave only one threshold for each marker,,, because,,, i'm not sure why there are two lol
            for (mi in seq_len(length(gates[[1]])))
                if (length(gates[[1]][[mi]])==2)
                    gates[[1]][[mi]] <- gates[[1]][[mi]][1]

            scat_chans <- unlist(sapply (c("SSC-A","FSC-A"), function(x)
                grep(x, colnames(fcs1@exprs), value=T)))[1]

            markers <- fg@markers
            gs <- make_gs(gs, gates, markers=markers,
                          gating_path_phenocode, scat_chans[1])

            shiny::incProgress(.5, detail="plotting density plot")

            return(gs)
        })
    })
    output$plot_dens1 <- renderPlot({
        req(gs_reactive())
        gs <- gs_reactive()
        flowWorkspace::plotGate(gs[[1]], xbin=70, gpar=list(ncol=2))
    })
    output$plot_dens2 <- renderPlot({
        req(gs_reactive())
        gs <- gs_reactive()
        flowWorkspace::plotGate(gs[[2]], xbin=70, gpar=list(ncol=2))
    })

}

