# toggling outliers for box plot
# meaningful labels



# Shiny APP for flowGraph visulaization
# Written by: Mehrnoush Malek (MM), Alice Yue (AY)
# Revision: March 2020, AY

# runApp('/mnt/f/Brinkman group/current/Alice/flowGraph/flowType3/shiny')
# install.packages("/mnt/f/Brinkman group/current/Alice/flowtype_metric/src/flowTypeFilter/flowTypeFilterC_1.0_R_x86_64-pc-linux-gnu.tar.gz",repos=NULL, type="source")
setwd("/mnt/f/Brinkman group/current/Alice/flowtype_metric/src/shiny")
source("_func.R")

start_dir <<- c(data="/mnt/f/FCS data/Tobias Kollmann")
# rslt_dir <<- c(data="/mnt/f/Brinkman group/current/Alice/flowtype_metric/result")
# fcs_dir <<- "/mnt/f/FCS data/Tobias Kollmann/flowType_data/Bcell_panel"
fg_dir <<- "/mnt/f/Brinkman group/current/Alice/flowtype_metric/result/epichipc_bcell/fg.Rdata"
# fg <<- get(load(fg_dir))
yes_char <- "yes"

# library(flowType)
# library(flowCore)
library(flowDensity)
library(flowTypeFilterC)
library(flowGraph)
library(flowWorkspace)
library(shinyFiles)
library(shiny)
library(ggiraph)
library(ggplot2)
library(DT)

server <- function(input, output, session) {
    shinyFiles::shinyDirChoose(input, "fcs_dir", roots=start_dir, session=session) # fcs
    # shinyFiles::shinyDirChoose(input, "fg_dir", roots=start_dir, session=session) # flowGraph
    shinyFiles::shinyDirChoose(input, "ft_dir", roots=start_dir, session=session) # flowType

    output$fcs_dir_out <- renderPrint({
        if (is.integer(input$fcs_dir)) {
            cat("no fcs directory has been selected")
        } else {
            parseDirPath(start_dir, input$fcs_dir)
        }
    })
    output$ft_dir_out <- renderPrint({
        if (is.integer(input$ft_dir)) {
            cat("no flowType directory has been selected")
        } else {
            parseDirPath(start_dir, input$ft_dir)
        }
    })


    # print table
    output$table_summary = DT::renderDataTable({
        # req(input$fg, input$p_thres)
        req(input$p_thres)
        print(input$p_tres)

        withProgress(message='loading', value=0, {
            incProgress(.3, detail="generating summary statistics table")

            # if (grepl("[.]rds",input$fg$datapath,ignore.case=TRUE))
            #     fg <<- readRDS(input$fg$datapath)
            # if (grepl("[.]rdata",input$fg$datapath,ignore.case=TRUE))
            #     fg <<- get(load(input$fg$datapath))
            # print(show(fg))

            # get table
            sdt <- fg@summary_desc$node
            print(nrow(sdt))
            indices <<- which(purrr::map_int(seq_len(nrow(sdt)), function(x)
                sum(fg@summary$node[[x]]$values<input$p_thres))>1)
            sdt <- sdt[indices,]
            print(nrow(sdt))
            datatable(sdt, selection="single") # https://stackoverflow.com/questions/28274584/get-selected-row-from-datatable-in-shiny-app
            #         callback = "function(table) {
            #       table.on('click.dt', 'tr', function() {
            #         table.$('tr.selected').removeClass('selected');
            #         $(this).toggleClass('selected');
            #     });
            # }"
            #         datatable(sdt, callback=JS(callback))
        })
    }, server=TRUE)

    # print the selected indices
    output$table_row = renderPrint({
        shiny::req(input$table_summary_rows_selected)

        ri <- input$table_summary_rows_selected
        index <<- indices[tail(ri,1)]
        if (length(ri)) {
            cat('row selected:')
            cat(paste0(ri, collapse=", "), " (", paste0(index, collapse=", "), ")", sep="")
        }
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



    selected_cpop <- reactive({ input$plot_hierarchy_selected })

    output$plot_hierarchy <- renderGirafe({
        req(input$table_summary_rows_selected, input$p_thres, input$max_nodes)
        p_thres <- input$p_thres
        max_nodes <- input$max_nodes

        withProgress(message='loading', value=0, {
            incProgress(.5, detail="generating cell hierarchy of significant phenotypes")

            ri <- input$table_summary_rows_selected

            # row index of selected row
            index <<- indices[tail(ri,1)]
            print(index)

            # # load fcs.Rdata or matrix m
            # m <<- fg@feat$node[[fg@summary_desc$node$feat[[index]]]]
            # print(dim(m))

            # plot cell hierarchy
            gr <- fg_plot(fg, index=index, p_thres=p_thres, show_bgedges=FALSE, label_max=max_nodes)#, node_labels=c("cells_ul_blood", "expect_cells_ul_blood"))
            gr$v$v_ind <- gr$v$label_ind
            gr$e$e_ind <- gr$e$from%in%gr$v$phenotype[gr$v$v_ind] &
                gr$e$to%in%gr$v$phenotype[gr$v$v_ind]
            gg_hier <- plot_gr(gr, shiny_plot=TRUE, bgedges=FALSE)

            ggiraph::girafe(
                code=print(gg_hier), width_svg=9, height_svg=5,
                options=list(
                    opts_selection(
                        type = "single", css = "fill:#FF3333;stroke:black;"),
                    opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;")
                ))
        })
    })

    # plot boxplot
    selected_label <- reactive({ input$plot_box_selected })

    output$plot_box <- renderGirafe({
        req(input$plot_hierarchy_selected, input$paired, input$dotplot, input$outlier)
        cpop <- input$plot_hierarchy_selected
        paired <- input$paired==yes_char
        dotplot <- input$dotplot==yes_char
        outlier <- input$outlier==yes_char

        gg_box <- fg_plot_box(
            fg, type="node", index=index, node_edge=cpop,
            paired=paired, outlier=outlier, dotplot=dotplot, shiny_plot=TRUE)
        print(index)
        # ggiraph::girafe(code=print(gg_box))
        ggiraph::girafe(
            code=print(gg_box), #width_svg=3, height_svg=3,
            options=list(
                opts_selection(
                    type = "single", css = "fill:#FF3333;stroke:black;"),
                opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;")
            ))

    })

    # plot scatterplots
    output$plot_dens <- renderPlot({
        req(input$plot_hierarchy_selected, input$plot_box_selected,
            input$ft_dir, input$fcs_dir)

        cpop <- input$plot_hierarchy_selected
        print(cpop)
        label <- input$plot_box_selected
        print(label)
        ft_dir <-  parseDirPath(start_dir, input$ft_dir)
        print(ft_dir)
        fcs_dir <- parseDirPath(start_dir, input$fcs_dir)
        print(fcs_dir)


        withProgress(message='loading', value=0, {
            incProgress(.2, detail="loading median samples")

            ft_paths <- list.files(ft_dir, full.names=T, pattern=".Rdata", recursive=TRUE)
            fcs_paths <- list.files(fcs_dir, full.names=T, pattern=".fcs", recursive=TRUE)

            # load fcs.Rdata or matrix m to choose median sample
            m <<- fg@feat$node[[fg@summary_desc$node$feat[[index]]]]
            print(dim(m))

            label_ind <- which(fg@meta[[unlist(fg@summary_desc$node$class[index])]]==label)
            f_ind <- label_ind[which.min(abs(m[label_ind,cpop]-median(m[label_ind,cpop])))]
            f_id <- fg@meta$id[f_ind]

            # load flowtype and fcs files
            ft <- get(load(ft_paths[grepl(f_id, ft_paths)]))
            fcs <- flowCore::read.FCS(fcs_paths[grepl(f_id, fcs_paths)])

            incProgress(.4, detail="finding best gating path")

            # find best path
            rchy <- RchyOptimyx::RchyOptimyx(
                pheno.codes=fg@graph$v$phenocode,
                phenotypeScores=-log(fg@summary$node[[index]]$values),
                startPhenotype=fg@graph$v$phenocode[fg@graph$v$phenotype==cpop],
                pathCount=1, trimPaths=FALSE)
            gating_path_phenocode <- rchy@nodes[1,-1]

            incProgress(.5, detail="creating gating set")

            # make gating set
            gates <- list(ft@Thresholds)
            gs <- CytoML::GatingSet(as(fcs,"flowSet"))
            sampleNames(gs) <- names(gates) <- f_id

            # leave only one threshold for each marker,,, because,,, i'm not sure lol
            for (mi in seq_len(length(gates[[1]])))
                if (length(gates[[1]][[mi]])==2)
                    gates[[1]][[mi]] <- gates[[1]][[mi]][1]

            scat_chans <- unlist(sapply (c("SSC-A","FSC-A"), function(x)
                grep(x, colnames(fcs@exprs), value=T)))[1]

            gs <- make_gs(gs, gates, markers=fg@markers,
                          gating_path_phenocode, scat_chans[1])

            incProgress(.5, detail="plotting density plot")

            flowWorkspace::plotGate(gs[[1]], xbin=70, gpar=list(ncol=4))
        })
    })
}

