fg_dir <- input$fg$datapath
fg <- get(load(fg_dir))
ifrom <<- match(fg@graph$e$from, fg@graph$v$phenotype)
ito <<- match(fg@graph$e$to, fg@graph$v$phenotype)
sigv <- c("all","insig","sig")

server <- function(input, output, session) {
    values <- shiny::reactiveValues(fg=fg)



    #### about ####
    output$about_ui <- renderUI({
        tagList(
            h2("comparing edge features", style="align: center;")
        )
    })



    #### input fg ####
    fg_obj <- shiny::reactive({
        shiny::req(input$fg)
        get(load(input$fg$datapath))
        ifrom <<- match(fg@graph$e$from, fg@graph$v$phenotype)
        ito <<- match(fg@graph$e$to, fg@graph$v$phenotype)
    })

    # fg description
    output$fg_descript_ui <- shiny::renderUI({
        shiny::req(values$fg)
        fg <- values$fg

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



    #### data table & options ####
    ip_filter_es <- reactive({
        switch(input$effect_size, negligible=1, small=2, medium=3, large=4)
    })

    output$ui_layer <- shiny::renderUI({
        fg <- values$fg

        shiny::sliderInput(
            inputId="layer",
            label="layers to display (cell population edge points to)",
            min=1, max=max(fg@graph$v$phenolayer), step=1, value=c(1,4)
        )
    })

    output$ui_count <- shiny::renderUI({
        max_c <- 10000
        if (!is.null(values$mc_med)) max_c <- max(values$mc_med)
        fg <- values$fg

        shiny::sliderInput(
            inputId="med_count_thres",
            label="min cell count median (max between sample phenotypes)",
            min=1, max=max_c, step=5, value=50)
    })

    # number of significant cell pops in each summary statistic
    dt_signo <- reactive({
        shiny::req(values$fg)
        fg <- values$fg
        filter_adjust0 <- input$adjust0
        filter_es <- ip_filter_es()
        p_thres <- input$p_t1

        cand_rows <- grep("SpecEnr",fg@summary_desc$node$feat)

        sdt <- fg@summary_desc$node
        signo <- purrr::map_int(cand_rows, function(x) {
            pp <- fg_get_summary(fg, index=x, filter_adjust0=filter_adjust0,
                                 filter_es=filter_es)
            sum(pp$values<p_thres)
        })
        signo_ <- rep(0,nrow(sdt))
        signo_[cand_rows] <- signo
        signo_
    })

    # summary statistics with >0 sig cell pops
    dt_indices <- shiny::reactive({
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
    output$dt_ui <- shiny::renderUI({
        fg <- fg_obj()

        shiny::conditionalPanel(
            condition=length(fg)>0,
            DT::dataTableOutput("table_summary")
        )
    })



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

        # cf_se <- unlist(values$fg@summary_desc$node$feat[dt_ind()])
        # values$compare_feats <- compare_feats <-
        #     c(cf_se, ifelse(cf_se=="SpecEnr","prop",gsub("SpecEnr_","",cf_se)))

        values$fg <- fg <- fg_obj()

        values$index <- index <- dt_indices()[tail(dt_ind(),1)]

        # # if feature is not SpecEnr, this is NONE,
        # # else it is the specenr, actual and expected feature names
        # nl <- "NONE"
        feat_ <- unlist(fg@summary_desc$node$feat[index])
        # if (grepl("SpecEnr", feat_)) {
            if (feat_=="SpecEnr") {
                nl <- c(feat_, "prop","expect_prop")
            } else {
                nl1 <- gsub("SpecEnr_", "", feat_)
                nl <- c(feat_, nl1, paste0("expect_", nl1))
            }
        # }
        values$node_labels <- nl

        # get summaries
        values$pp0 <- fg_get_summary(
            fg, type="edge", summary_meta=unlist(append(
                "SpecEnr", fg@summary_desc$node[index,-1])))

        values$pp1 <- pp1 <- fg_get_summary(fg, index=index)

        values$pp2 <- fg_get_summary(
            fg, summary_meta=unlist(append(
                nl[2], fg@summary_desc$node[index,-1])))

        # get default cell population
        cpopi <- which.min(values$pp1$values)
        values$cpop <- cpop <- names(values$pp1$values)[cpopi]
        values$cpop_l <- fg@graph$v$phenolayer[cpopi]

        # get counts
        mc_med1 <- apply(as.matrix(fg@feat$node$count[pp1$id1,]),2,median)
        mc_med2 <- apply(as.matrix(fg@feat$node$count[pp1$id2,]),2,median)
        values$mc_med <- ifelse(mc_med1>mc_med2, mc_med1, mc_med2)

        # make df
        df <- data.frame(
            phenotype=paste0(fg@graph$e$from,"_",fg@graph$e$to),
            phenolayer=fg@graph$v$phenolayer[ito],
            eprop_change=abs(
                colMeans(as.matrix(fg@feat$edge$prop[pp1$id1,])) -
                    colMeans(as.matrix(fg@feat$edge$prop[pp1$id2,]))),
            mean_count=colMeans(as.matrix(fg@feat$node$count))[ito],
            xpfrom=pp1$values[ifrom],
            xpto=pp1$values[ito],
            ypfrom=pp2$values[ifrom],
            ypto=pp2$values[ito]
        )
        df$label_long <- paste0(
            df$phenotype,"\n- SpecEnr: ", df$xpfrom, " > ", df$xpto,
            "\n- raw: ", df$ypfrom, " > ", df$ypto)

        df$x <- log(df$xpto,10) - log(df$xpfrom,10)
        df$y <- log(df$ypto,10) - log(df$ypfrom,10)
        df$z <- log(values$pp0$values,10)

        values$df <- df

    })








    #### edge indices ####
    ip_einds0 <- shiny::reactive({
        req(values$df)

        p_t1 <- input$p_t1 # .05
        p_t2 <- input$p_t2 # .05
        cd_t1 <- input$cd_t1 # .5

        l_ <- input$layer
        l1 <- l_[1] # 1
        l2 <- l_[2] # 4

        sg <- c(input$xsigfrom, input$xsigto, input$ysigfrom, input$ysigto)
        med_count_thres <- input$med_count_thres

        filter_adjust0 <- input$filter_adjust0
        filter_es <- ip_filter_es()

        pp0 <- values$pp0
        pp1 <- values$pp1
        pp2 <- values$pp2

        mc_med <- values$mc_med

        # whether specenr should be significant based on how different actual
        # and expected raw values are; sigcands are those that are significantly different
        if (is.null(fg@etc$actualVSexpect))
            fg@etc$actualVSexpect <- list(node=list(), edge=list())
        df_ps <- fg@etc$actualVSexpect$node[[nl[1]]]

        sig1 <- df_ps$tpv1_paired<p_t1 & df_ps$tpv1<p_t1 & abs(df_ps$cd1)>=cd_t1
        sig2 <- df_ps$tpv2_paired<p_t1 & df_ps$tpv2<p_t1 & abs(df_ps$cd2)>=cd_t1
        sigcands <- sig1 | sig2

        # edge version of above
        df_ps_ <- fg@etc$actualVSexpect$edge$SpecEnr

        sig1_ <- df_ps_$tpv1_paired<p_t1 & df_ps_$tpv1<p_t1 & abs(df_ps_$cd1)>=cd_t1
        sig2_ <- df_ps_$tpv2_paired<p_t1 & df_ps_$tpv2<p_t1 & abs(df_ps_$cd2)>=cd_t1
        sigcands_ <- sig1_ | sig2_

        # get layer and count indices
        good_count <- mc_med > med_count_thres
        good_count_ <- good_count[ito] & good_count[ifrom]

        # get p-value indices
        pp0$values[(!sigcands_ | as.numeric(pp0$effect_size)>=filter_es |
                        pp0$adjust>filter_adjust0) & pp0$values<p_t1] <- p_t1
        # pp0sig <- pp0$values<p_t1

        pp1$values[(!sigcands | as.numeric(pp1$effect_size)>=filter_es |
                        pp1$adjust>filter_adjust0) & pp1$values<p_t1] <- p_t1
        pp1sig <- pp1$values<p_t1

        pp2sig <- pp2$values<p_t2


        df <- values$df
        f_inds <- df$phenolayer >= l1 & df$phenolayer <= l2 & good_count_
        return(f_inds)
    })

    ip_einds1 <- reactive({
        req(ip_einds0())
        f_inds <- ip_einds0()
        sg <- values$sg

        if (sg[1]==1) f_inds <- !pp1sig[ifrom] & f_inds
        if (sg[1]==2) f_inds <- pp1sig[ifrom] & f_inds
        if (sg[2]==1) f_inds <- !pp1sig[ito] & f_inds
        if (sg[2]==2) f_inds <- pp1sig[ito] & f_inds
        return(f_inds)
    })

    ip_einds2 <- reactive({
        shiny::req(ip_einds0())

        f_inds <- ip_einds0()
        sg <- values$sg
        if (sg[3]==1) f_inds <- !pp2sig[ifrom] & f_inds
        if (sg[3]==2) f_inds <- pp2sig[ifrom] & f_inds
        if (sg[4]==1) f_inds <- !pp2sig[ito] & f_inds
        if (sg[4]==2) f_inds <- pp2sig[ito] & f_inds
        return(f_inds)
    })

    ip_einds12 <- reactive({
        shiny::req(ip_einds2())
        ip_einds1() & ip_einds2()
    })

    #### static plots ####
    output$plot_SpecEnrVSpropdiff <- shiny::renderPlot({
        shiny::req(ip_einds1())
        df <- values$df
        df_inds <- ip_einds1()

        qp <- ggplot2::ggplot(df[df_inds,], ggplot2::aes(
            x=x,y=eprop_change, colour=eprop_change,
            size=mean_count, alpha=.3,stroke=1)) +
            ggplot2::geom_abline(intercept=0, slope=1) +
            ggplot2::ggtitle(paste0(
                "difference between child and parent ",
                "-log10 ","p-values")) +
            ggplot2::labs(
                x="SpecEnr p-value difference on edges",
                y="mean child/parent prop diff",
                col="mean child/parent prop diff",
                size="count mean (to)") +
            ggplot2::geom_point(shape=1) +
            ggplot2::scale_colour_gradientn(
                colours=c('blue','cyan','yellow','red'))

        qp
    })

    output$plot_rawVSpropdiff <- shiny::renderPlot({
        shiny::req(ip_einds1())
        df <- values$df
        df_inds <- ip_einds1()

        qp <- ggplot2::ggplot(df[df_inds,], ggplot2::aes(
            x=y,y=eprop_change, colour=eprop_change,
            size=mean_count, alpha=.3,stroke=1)) +
            ggplot2::geom_abline(intercept=0, slope=1) +
            ggplot2::ggtitle(paste0(
                "difference between child and parent ",
                "-log10 ","p-values")) +
            ggplot2::labs(
                x="raw p-value difference on edges",
                y="mean child/parent prop diff",
                col="mean child/parent prop diff",
                size="count mean (to)") +
            ggplot2::geom_point(shape=1) +
            ggplot2::scale_colour_gradientn(
                colours=c('blue','cyan','yellow','red'))

        qp
    })

    output$plot_edgeSpecEnrVSraw <- shiny::renderPlot({
        shiny::req(ip_einds1())
        df <- values$df
        df_inds <- ip_einds1()

        qp <- ggplot2::ggplot(df[df_inds,], ggplot2::aes(
            x=z,y=y, colour=eprop_change,
            size=mean_count)) +
            ggplot2::geom_abline(intercept=0, slope=1) +
            ggplot2::ggtitle(paste0("(",sum(df_inds),") ",
                "difference between child and parent ",
                "-log10 ","p-values")) +
            ggplot2::labs(
                x="SpecEnr p-value on edges",
                y="cells_ul_blood/prop p-value difference on edges",
                col="mean child/parent prop diff",
                size="count mean (to)") +
            ggplot2::scale_colour_gradientn(
                colours=c('blue','cyan','yellow','red')) +
                ggplot2::geom_point(shape=1,stroke=1)

        qp
    })

    #### interactive plot ####
    output$plot_SpecEnrVSraw <- ggiraph::renderGirafe({
        shiny::req(ip_einds12())
        df <- values$df
        df_inds <- ip_einds12()

        qp <- ggplot2::ggplot(df[df_inds,]) +
            ggiraph::geom_point_interactive(ggplot2::aes(
                x=x,y=y, colour=eprop_change, size=mean_count,
                alpha=.3,stroke=1, tooltip=label_long, data_id=phenotype)) +
            ggplot2::geom_abline(intercept=0, slope=1) +
            ggplot2::ggtitle(paste0("(",sum(f_inds),")",
                "difference between child and parent ",
                "-log10 ","p-values")) +
            ggplot2::labs(
                x="SpecEnr p-value difference on edges",
                y="cells_ul_blood/prop p-value difference on edges",
                col="mean child/parent prop diff",
                size="count mean (to)") +
            ggplot2::scale_colour_gradientn(
                colours=c('blue','cyan','yellow','red'))

        ggiraph::girafe(
            code = print(qp),
            width_svg = 6, height_svg = 8,
            options = list(
                opts_selection(
                    type = "multiple", css = "fill:#FF3333;stroke:black;"),
                opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;")
            ))
    })

    # reset interactive plot
    observeEvent(input$reset, {
        session$sendCustomMessage(type='plot_set', message=character(0))
    })

    # table
    selected_cpop <- reactive({
        input$plot_SpecEnrVSraw_selected
    })
    output$datatab <- renderTable({
        req(selected_cpop())
        df <- values$df

        out <- df[df$phenotype %in% selected_cpop(), -"label_long"]
        if( nrow(out) < 1 ) return(NULL)
        row.names(out) <- NULL
        out
    })





    ### cell pop search and go ###
    # fill textbox with node chosen in interactive plot
    shiny::observeEvent(selected_cpop(), {
        # shiny::req(input$graph_selected)
        shiny::updateTextInput(
            session, inputId="type_cpop", value=tail(selected_cpop(),1))
    })



    #### GO: cell pop selected selected_cpop() ####
    shiny::observeEvent(input$go_cpop, {
        shiny::req(input$type_cpop)
        print("GO CPOP")

        df <- values$df
        text_cpop <- text_cpop_ <- input$type_cpop # init = NULL; empty = length()==1 ""

        if (!is.null(text_cpop_)) {
            text_cpopi <- which(df$phenotype==text_cpop)
            if (text_cpopi) {
                values$cpop <- text_cpop
                values$cpop_l <- df$phenolayer[text_cpopi]
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
        # fg <- values$fg
        index <- values$index
        cpop <- values$cpop
        node_labels <- values$node_labels
        pp <- values$pp1

        cpopi <- which(names(pp$values)==cpop)
        sm <- unlist(fg@summary_desc$node[index,])
        m0a <- signif(pp$m1[cpopi],3)
        m0b <- signif(pp$m2[cpopi],3)

        print(paste("cpop INFO: ", index, cpop, m0a, m0b, paste0(node_labels, collapse="/"), cpopi, paste0(sm, collapse="/")))

        fm <- ""
            m1a <- signif(mean(as.matrix(fg@feat$node[[node_labels[1]]])[pp$id1,cpopi]),3)
            m1b <- signif(mean(as.matrix(fg@feat$node[[node_labels[1]]])[pp$id2,cpopi]),3)
            m2a <- signif(mean(as.matrix(fg@feat$node[[node_labels[2]]])[pp$id1,cpopi]),3)
            m2b <- signif(mean(as.matrix(fg@feat$node[[node_labels[2]]])[pp$id2,cpopi]),3)
            fm <- paste0(" (actual/expect: ",m1a,"/",m2a," vs ",m1b,"/",m2b,")")
            print(fm)

        shiny::HTML(paste0(
            "<h4>",sm["feat"],": ",cpop,"</h4>",
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

        fg_plot_box(fg, type="node", index=index, node_edge=cpop,
                    paired=paired, outlier=outlier, dotplot=dotplot)
    })

    output$plot_box1 <- shiny::renderPlot({
        shiny::req(values$cpop)
        fg <- values$fg
        index <- values$index
        cpop <- values$cpop
        paired <- ip_paired()
        dotplot <- ip_dotplot()
        outlier <- ip_outlier()
        node_labels <- values$node_labels
        print(paste("box1: ", cpop, paired, dotplot, outlier, node_labels[2]))

        fg_plot_box(
            fg, type="node", index=index, node_edge=cpop,
            paired=paired, outlier=outlier, dotplot=dotplot,
            feature=node_labels[2])
    })

    output$plot_box2 <- shiny::renderPlot({
        shiny::req(values$cpop)
        if (!show_boxes()) return(NULL)
        # fg <- values$fg
        index <- values$index
        cpop <- values$cpop
        paired <- ip_paired()
        dotplot <- ip_dotplot()
        outlier <- ip_outlier()
        node_labels <- values$node_labels
        print(paste("box2: ", cpop, paired, dotplot, outlier, node_labels[3]))

        fg_plot_box(
            fg, type="node", index=index, node_edge=cpop,
            paired=paired, outlier=outlier, dotplot=dotplot,
            feature=node_labels[3])
    })

    output$box_plots_ui <- shiny::renderUI({
        shiny::req(values$cpop)
        car <- bsplus::bs_carousel(id="box_carosel", use_indicators=TRUE) %>%
            bsplus::bs_append(shiny::plotOutput(outputId="plot_box")) %>%
            bsplus::bs_append(shiny::plotOutput(outputId="plot_box1")) %>%
            bsplus::bs_append(shiny::plotOutput(outputId="plot_box2"))
        return(car)
    })

    # session$onSessionEnded(stopApp) # stop session when window closed
}
