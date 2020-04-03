# Shiny APP for flowGraph visulaization
# Written by: Mehrnoush Malek (MM), Alice Yue (AY)
# Revision: March 2020, AY

# runApp('/mnt/f/Brinkman group/current/Alice/flowGraph/flowType3/shiny')
setwd("/mnt/f/Brinkman group/current/Alice/flowtype_metric/src/flowType3/shiny")
source('3i-Helper.R')
source('slidingRchy-Helper.R')
data_dir <- "/mnt/f/FCS data/Tobias Kollmann"

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


server <- function(input, output, session) {
    
    # choose files
    volumes <- c(data=data_dir) # directory from where users start choosing files
    shinyFiles::shinyDirChoose(input, "result", roots=volumes, session=session) # flowType
    shinyFiles::shinyDirChoose(input, "files", roots=volumes, session=session) # fcs
    shinyFiles::shinyFileChoose(input, "arguments", roots=volumes) # arguments
    shinyFiles::shinyFileChoose(input, "fcs", roots=volumes) # one fcs file
    
    # start analysis = go
    shiny::observeEvent(input$go, {
        # get directories
        ft_dir <<- paste0(shinyFiles::parseDirPath(volumes, input$result), "")
        print(ft_dir)
        ft_paths <- list.files(ft_dir, full.names=T, pattern=".Rdata")
        fcs_dir <<- paste0(shinyFiles::parseDirPath(volumes, input$files), "")
        print(fcs_dir)
        fcs_paths <- list.files(fcs_dir, full.names=T, pattern=".Rdata")
        
        # load example flowtype file
        ft <- get(load(ft_paths[1]))
        
        # load arguments
        arguments <- get(load(paste0(parseFilePaths(roots=volumes, input$arguments),"")))
        print(names(arguments)) # "prp" "part" "max" "part.v" "filt.v" "names" "threshold"; names=markers
        print(head(ft@CellFreqs))
        
        # phenotypes
        phenotype <- sapply(ft@PhenoCodes, flowTypeFilter::decodePhenotype, 
                            marker.names=arguments$names, 
                            partitions.per.marker=ft@PartitionsPerMarker, 
                            purrr::map_int(ft@Thresholds, length)==1)
        print(head(phenotype))
        
        # load fcs.Rdata or matrix mc
        withProgress(message='Analyzing: ', value=0, {
            incProgress(.3, detail="Calculating proportions")
            if (input$calc=="NO") { # load cell counts 
                mc <- c()
                cnt <- 1
                n <- 10
                for (ftp in ft_paths[1:n]) {
                    print(cnt)
                    ft <- get(load(ftp))
                    mc <- rbind(mc, ft@CellFreqs)
                    cnt <- cnt+1
                }
                colnames(mc) <- as.vector(phenotype)
                rownames(mc) <- gsub(basename(ft_paths)[1:n],pattern="_FT.Rdata",replacement="")
            } else {
                print("here")
                mc <- get(load(paste0(volumes,"/Allcounts_flowType.RData")))
            }
            incProgress(.5, detail="Calculating p-values and plotting.")
        })
        
        # rchyoptimyx output$plotgraph (boxplot inside)
        output$plotgraph <- renderVisNetwork({
            rchy.graph <- make.sliding.graph(
                arguments=arguments, ft=ft, counts=mc,
                method1=input$method1, method2=input$method2, parent="CD45+",
                count.lim=input$slider1, p.level=input$slider2)
            # returns a list: nodes, edges
            # nodes <- data.frame(id=merged.res@nodes[1,], # phenocode
                                # label=nodes.lbl, # phenotype
                                # color=cols, # colour
                                # shape="ellipse",
                                # level=phenolayer) # phenolayer
            # #You can add font.color, label, and font.size to edges
            # edges <- data.frame(
            #     from=unlist(lapply(merged.res@edges[1,],function(x) # phenocode
            #         unlist(strsplit(x,split = "~"))[1])),
            #     to=unlist(lapply(merged.res@edges[1,],function(x) # phenocode
            #         unlist(strsplit(x,split = "~"))[2])),
            #     weight=edge.weights,value=edge.weights)

            nodes <- rchy.graph$nodes
            print(nodes)
            nodes.id <<-nodes$id
            nodes.lbl <<-nodes$label
            edges <- rchy.graph$edges
            edges.to <<- edges$to
            edges.from <<- edges$from
            visNetwork(nodes, edges) %>%
                visNodes(shape="ellipse", size=55, label=nodes$label) %>%
                visEdges(arrows="to")%>%
                visOptions(highlightNearest=TRUE, 
                           nodesIdSelection=TRUE) %>%
                visInteraction(hover=TRUE,multiselect=TRUE) %>%
                visLayout(hierarchical=TRUE) %>%
                visHierarchicalLayout(parentCentralization=TRUE) %>%
                visEvents(select="function(nodes) {
                  Shiny.onInputChange('current_node_selection', nodes.nodes);
                  ;}")
            
        })
        
        selected_sample <- reactive({
            input$plot_selected
        })
        
        myNode <- reactiveValues(selected='', pops='')
        observeEvent(input$current_node_selection, {
            myNode$selected <<- as.vector(input$current_node_selection)
            myNode$pops <<- as.vector(nodes.lbl)[match(unlist(myNode$selected),as.vector(nodes.id) )]
            print(paste0("populations selected: ",paste(myNode$pops,collapse=", ")))
            
            # prop
            mp <- mc/mc[,1]
            mp <- mp*100
            print(dim(mp))
            
            # labels
            lbls <- sample (x=0:1, size=nrow(mc), replace=TRUE)
            
            # data frame with a cell population proportion, labels, filenames as columns
            pop.prop <-data.frame(
                props=c(mp[which(lbls==0), myNode$pops], mp[which(lbls==1), myNode$pops]),
                groups=c(rep(times=length(which(lbls==0)),"Group1"),
                         rep(times=length(which(lbls==1)),"Group2")),
                files=c(rownames(mp)[which(lbls==0)],
                        rownames(mp)[which(lbls==1)]))
            
            group1i <- which(pop.prop$groups=="Group1")
            group2i <- which(pop.prop$groups=="Group2")
            
            # extract only samples / rows with min/median/max proportions for the cell population
            sample.props <- rbind(
                pop.prop[group1i[which.min(pop.prop$props[group1i])],], # row in group 1 with min proportion
                pop.prop[group1i[which.min(abs(pop.prop$props[group1i]-
                                                   median(pop.prop$props[group1i])))],], # row in group 1 with median proportion
                pop.prop[group1i[which.max(pop.prop$props[group1i])],], # row in group 1 with max proportion
                pop.prop[group2i[which.min(pop.prop$props[group2i])],],
                pop.prop[group2i[which.min(abs(pop.prop$props[group2i]-
                                                   median(pop.prop$props[group2i])))],],
                pop.prop[group2i[which.max(pop.prop$props[group2i])],])
            sample.props$files <- paste(sample.props$files,sample.props$groups,sep="_g_")
            
            # You need to create 
            # output$boxplot
            if (input$random=="YES" | is.null(input$fcs)) {
                print("random samples")
                
                output$boxplot <- renderggiraph({
                    x <- ggiraph(code=print(
                        plot_results(cell.proportion=pop.prop, outliers=sample.props) 
                    ))
                    x <- girafe_options(
                        x, opts_selection(type="multiple", css="fill:#FF3333;stroke:black;"),
                        opts_hover(css="fill:#FF3333;stroke:black;cursor:pointer;"))
                    x
                })
            }else{
                print("selected sample")
                
                fcs.file <- parseFilePaths(roots=volumes, input$fcs)
                print(parseFilePaths(roots=volumes,input$fcs)$datapath)
                print(fcs.file)
                
                fcs.samples <- as.character(parseFilePaths(roots=volumes,input$fcs)$name)
                samp.ind <- which(gsub("_FT.Rdata","",dirname(ft_paths))==unlist(strsplit(fcs.samples,split=".Rdata"))[1])
                
                ft <- get(load(ft_paths[samp.ind]))
                samp.count<- ft@CellFreqs
                names(samp.count) <- colnames(mc)
                
                sample.gr <- "Group1"
                sample.data <- as.data.frame(props=samp.count[myNode$pops]*100/samp.count[1],files=unlist(strsplit(fcs.samples,split=".Rdata"))[1],groups=sample.gr)
                output$boxplot <-renderggiraph({x<-ggiraph(code=print(plot_results(cell.proportion=pop.prop,showsample=T,samplestoplotdata=sample.data,
                                                                                   outliers=sample.props)))
                x <- girafe_options(x, opts_selection(
                    type="multiple", css="fill:#FF3333;stroke:black;"),
                    opts_hover(css="fill:#FF3333;stroke:black;cursor:pointer;"))
                x
                })
                
                
            }
            
            ## START PATH pop.path from 1st layer to leaf
            
            # st = index of an edge from node selected
            st <- which(edges.to==unlist(myNode$selected))[1]
            print(st)
            
            # node 1
            selected.node <- unlist(myNode$selected)
            print(selected.node)
            
            # node 2
            start.path <- c(selected.node, as.vector(edges.from[st]))
            
            load(paste0(dirname(ft_dir),"/Settingup_argument.RData")) #??? why load arguments again... don't know scope lol
            
            # node n -- i'll just get the lowest p value path when i do it lol
            nxt <- tail(start.path,1)
            pop.path <- start.path
            while (nxt!=as.vector(edges.from[1])) {
                prnt <- which(edges.to==nxt)[1]
                nxt<- as.vector(edges.from[prnt])
                pop.path <-c( pop.path ,as.vector(nodes.id[which(nodes.id==nxt)]))
            }
            parent.names <- rev(nodes.lbl[match(pop.path,nodes.id)])
            parent.names <- parent.names[-length(parent.names)]
            pop.path <- rev(pop.path)[-1]
            
            # load FCS file and its gates
            fcs.file <- paste0(parseFilePaths(roots=volumes,input$fcs),"")
            print(parseFilePaths(roots=volumes,input$fcs)$datapath)
            print(fcs.file)
            
            fcs.samples <- as.character(parseFilePaths(roots=volumes,input$fcs)$name)
            samp.ind <- which(names(arguments$threshold)==unlist(strsplit(fcs.samples,split=".Rdata"))[1])
            
            gates <-  arguments$threshold[[samp.ind]]
            frame <- get(load(parseFilePaths(roots=volumes,input$fcs)$datapath))

            output$plotdens <- renderPlot({
                click.input <- input$boxplot_selected
                print(click.input)
                
                if (!is.null(click.input)) {
                    click.gr <- unlist(lapply(click.input, function(cl) unlist(strsplit(cl,"_g_"))[2]))
                    click.id <- unlist(lapply(click.input, function(cl) unlist(strsplit(cl,"_g_"))[1]))
                    print(paste("click is:",length(click.input)))
                    print(click.id)
                    print(click.gr)
                    
                    ft <- get(load(ft_paths[2]))
                    load(paste0(dirname(ft_dir),"/Settingup_argument.RData")) #??? why load arguments again... don't know scope lol
                    
                    gs <- make.gs(ft=ft, ids=tail(click.id,2),arguments=arguments,pop.path=pop.path)
                    print(class(gs))
                    
                    # arrange plots in a square * 2
                    vec <- c(1:(ceiling(sqrt(length(pop.path)))*2))
                    vec[which(vec>=length(gs_get_pop_paths(gs)))] <- 0
                    mat <-matrix(vec,nrow =ceiling(sqrt(length(pop.path))) ,byrow=T)
                    vec2 <- vec + max(vec,na.rm=T)
                    vec2[which(vec2==max(vec))]<-0
                    mat2 <-cbind(mat,matrix(vec2,nrow =ceiling(sqrt(length(pop.path))) ,byrow=T))
                    
                    if (length(gs)==2) { # number of samples
                        print("OK")
                        print(mat2)
                        plot.gates(gs,popstoplot=NULL,upto=NULL,plot.layout=mat2,
                                   show.stats=TRUE,contour.pops=NULL,
                                   density.overlay.pops=NULL,
                                    flowCut.channels=c("Time","BV510-A"),
                                   legend.cex=1.5,cell.size=c(1,2.5,3),
                                   show.palette=TRUE,offset=.5)
                        title(paste(click.gr,sep=" ",collapse="      -        "), outer=TRUE, line=2)
                        
                        
                    }else{
                        print("Fine")
                        print(mat)
                        plot.gates (gs,popstoplot=NULL,upto=NULL,plot.layout=mat,
                                    show.stats=TRUE,contour.pops=NULL,
                                    density.overlay.pops=NULL,
                                    flowCut.channels=c("Time","BV510-A"),
                                    legend.cex=1.5,cell.size=c(1,2.5,3),
                                    show.palette=TRUE,offset=.5)
                        title(paste(click.gr,sep="  ",collapse="  -   "), outer=TRUE, line=2)
                        
                        
                    }
                    
                    
                    
                }
            })
            
            
        })
    })
    
}
