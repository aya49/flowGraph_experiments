# Shiny APP for flowGraph visulaization
# Written by: Mehrnoush Malek (MM), Alice Yue (AY)
# Revision: March 2020, AY

# runApp('/mnt/f/Brinkman group/current/Alice/flowGraph/flowType3/shiny')
setwd("/mnt/f/Brinkman group/current/Alice/flowtype_metric/src/flowType3/shiny")
source('3i-Helper.R')
source('slidingRchy-Helper.R')

# data_dir <- "/mnt/f/FCS data/Tobias Kollmann"
ft_dir <<- "/mnt/f/FCS data/Tobias Kollmann/flowType_results/Bcell_panel"
fcs_dir <<- "/mnt/f/FCS data/Tobias Kollmann/flowType_data/Bcell_panel"
fg_dir <<- "/mnt/f/Brinkman group/current/Alice/flowtype_metric/result/epichipc_bcell/fg.Rdata"
p_thres <<- .01
# root name
root_name <<- "live" # renderText({ input$root_name })
max_nodes <<- 100
paired <<- FALSE

library(flowType)
library(flowCore)
library(flowDensity)
library(flowTypeFilter)
library(rdist)
library(shinyFiles)
library(shiny)
library(shinydashboard)
library(fields)
library(visNetwork)
library(ggiraph)
library(ggplot2)
library(flowWorkspace)
library(DT)


server <- function(input, output, session) {
    
    # # choose files
    # volumes <- c(data=data_dir) # directory from where users start choosing files
    # shinyFiles::shinyDirChoose(input, "ft_folder", roots=volumes, session=session) # flowType
    # shinyFiles::shinyDirChoose(input, "fcs_folder", roots=volumes, session=session) # fcs
    # shinyFiles::shinyDirChoose(input, "fg_folder", roots=volumes, session=session) # fcs
    
    # # start analysis = go
    # shiny::observeEvent(input$go, {
    # get directories
    # ft_dir <<- paste0(shinyFiles::parseDirPath(volumes, input$ft_folder), "")
    print(paste("flowType folder:", ft_dir))
    ft_paths <<- list.files(ft_dir, full.names=T, pattern=".Rdata", recursive=TRUE)
    
    # fcs_dir <<- paste0(shinyFiles::parseDirPath(volumes, input$fcs_folder), "")
    print(paste("fcs folder:", fcs_dir))
    fcs_paths <<- list.files(fcs_dir, full.names=T, pattern=".fcs", recursive=TRUE)
    
    # fg <<- reactive({
    #     fg_dir <<- paste0(shinyFiles::parseDirPath(volumes, input$fg_folder), "")
    #     print("flowGraph folder:", fg_dir)
    #     flowGraph::fg_load(fg_dir)
    # })
    
    fg <<- get(load(fg_dir))
    print(show(fg))
    
    # get table
    sdt <- fg@summary_desc$node
    print(nrow(sdt))
    indices <- which(purrr::map_int(seq_len(nrow(sdt)), function(x) 
        sum(fg@summary$node[[x]]$values<p_thres)) >1)
    sdt <- sdt[indices,]
    print(nrow(sdt))
    
    # print table
    output$table_summary = DT::renderDataTable(    
        DT::datatable(
            sdt,
            selection="none",
            # select only one row
            callback="function(table) {
      table.on('click.dt', 'tr', function() {
                table.$('tr.selected').removeClass('selected');
                $(this).toggleClass('selected');            
                Shiny.onInputChange('rows',
                table.rows('.selected').data()[0][0]);
});
                }"
        ), server=TRUE)
    
    # row index of selected row
    index <- indices[input$table_summary_rows_selected]
    print(index)
    
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
    
    # load fcs.Rdata or matrix m
    m <- fg@feat$node[[fg@summary_desc$node$feat[[index]]]]
    print(dim(m))
    
    # plot cell hierarchy
    gr <- fg_plot(fg, index=index, p_thres=p_thres, show_bgedges=FALSE, label_max=max_nodes)
    gr$v$v_ind <- gr$v$label_ind
    gr$e$e_ind <- gr$e$from%in%gr$v$phenotype[gr$v$v_ind] &
        gr$e$to%in%gr$v$phenotype[gr$v$v_ind]
    gg_hier <- plot_gr(gr, shiny_plot=TRUE)
    
    selected_cpop <- reactive({ input$plot_hierarchy_selected })
    
    output$plot_hierarchy <- renderGirafe({
        ggiraph::girafe(
            code=print(gg_hier), width_svg=9, height_svg=5,
            options=list(
                opts_selection(
                    type = "single", css = "fill:#FF3333;stroke:black;"),
                opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;")
            ))
    })
    
    # plot boxplot
    selected_label <- reactive({ input$plot_box_selected })
    
    output$plot_box <- renderGirafe({
        cpop <- selected_cpop()
        print(cpop)
        if (!base::is.null(cpop)) {
            gg_box <- fg_plot_box(fg, cpop=cpop, paired=paired, shiny_plot=TRUE)
            gp <- ggiraph::girafe(
                code=print(gg_box), width_svg=3, height_svg=3,
                options=list(
                    opts_selection(
                        type = "single", css = "fill:#FF3333;stroke:black;"),
                    opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;")
                ))
            return(gp)
        }
    })
    
    # plot scatterplots
    
    
    
    output$plot_dens <- renderPlot({
        cpop <- selected_cpop()
        label <- selected_label()
        print(cpop)
        print(label)
        
        if (!base::is.null(cpop) & !base::is.null(label)) {
            
            label_ind <- which(fg@meta[[fg@summary_desc$node$class[index]]]==label)
            f_ind <- label_ind[which.min(abs(m[label_ind,cpop]-median))]
            f_id <- fg@meta$id[f_ind]
            
            # load flowtype and fcs files
            ft <- get(load(ft_paths[grepl(f_id, ft_paths)]))
            fcs <- flowCore::read.FCS(fcs_paths[grepl(f_id, fcs_paths)])
            
            # find best path
            pp <- fg@summary$node[[index]]
            rchy <- RchyOptimyx::RchyOptimyx(
                pheno.codes=fg@graph$v$phenocode, phenotypeScores=-log(pp$values),
                startPhenotype=fg@graph$v$phenocode[fg@graph$v$phenotype==cpop],
                pathCount=1, trimPaths=FALSE)
            gating_path_phenocode <- rchy@nodes[1,-1]
            
            # make gating set
            gating_path_phenocode_ <- lapply(gating_path_phenocode, function(p) 
                unlist(strsplit(p,"")))
            
            gs <- CytoML::GatingSet(as(fcs,"flowSet"))
            sampleNames(gs) <- f_id
            
            gates <- ft@Thresholds
            scat.chans <- unlist(sapply (c("SSC-A","FSC-A"), function(x) 
                grep(x, colnames(fs), value=T)))[1]
            for (i1 in 1:length(gating_path_phenocode_)) {
                cpop_vec <- gating_path_phenocode_[[i1]]
                ind_m <- which(cpop_vec!="0")
                if (names(ft@Thresholds)[ind_m]!="gate") {
                    fdens <- flowDensity(
                        fcs, 
                        channels=c(names(ft@Thresholds)[ind_m], scat.chans[1]), 
                        position=c(cpop_vec[ind_m]!="1",NA),
                        gates=c(gates[[ind_m]],NA))
                    filter_id <- 
                        paste0(fg@markers[ind_m], 
                               ifelse(cpop_vec[ind_m]!="1", "+", "-"))
                } else {
                    if (cpop_vec[ind_m]=="1") {
                        fdens_fun <- flowDensity::notSubFrame
                    } else {
                        fdens_fun <- flowDensity::flowDensity
                    }
                    fdens <- fdens_fun(
                        fcs, 
                        channels=colnames(gates[[ind_m]]),
                        position=c(T,T),
                        filter=as.matrix(gates[[ind_m]]))
                    filter_id <- paste0(ifelse(cpop_vec[ind_m]=="1", "Not ", ""),
                                        fg@markers[ind_m])
                }
                
                # mark off this marker for future gates
                if (i1<length(gating_path_phenocode_))
                    for (k1 in (i1+1):length(gating_path_phenocode_))
                        gating_path_phenocode_[[k1]][ind_m] <- "0"
                
                # apply to gating set
                poly <- list(polygonGate(filterId=filter_id, .gate=fdens@filter))
                names(poly) <- sampleNames(gs)
                
                nodeID <- flowWorkspace::gs_pop_add(
                    gs, poly, 
                    parent=tail(flowWorkspace::gs_get_pop_paths(gs),1))
                recompute(gs)
            }
            print(class(gs))
            
            return(plotGate(gs[[1]], xbin=70, gpar=list(ncol=4)))
        }
    })
    
}

