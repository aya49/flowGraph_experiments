#Shiny APP for flowType/Rchy visulaization and checking outliers
#Written by Mehrnoush Malek (MM)
#Revise: March 2020, MM

# runApp('/mnt/f/Brinkman group/current/Alice/flowGraph/flowType3/shiny')
setwd("/mnt/f/Brinkman group/current/Alice/flowtype_metric/src/flowType3/shiny")
source('3i-Helper.R')
source('slidingRchy-Helper.R')

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


server<- function(input, output,session) {
  volumes = c(wd="~/data/semisupervised/")
  shinyDirChoose(input, "result", roots = volumes, session = session)
  shinyDirChoose(input, "files", roots = volumes, session = session)
  shinyFileChoose(input, "fcs",roots = volumes)
  observeEvent(input$go, {
    dir.path <<- paste0(parseDirPath(volumes,input$result),"")
    print(dir.path)
    ft.path<- list.files(dir.path,full.names = T, pattern = ".Rdata")
    files.path <<- paste0(parseDirPath(volumes,input$files),"")
    fs.path<- list.files(files.path,full.names = T, pattern = ".Rdata")
    tmp <- load(ft.path[1])
    ft <- get(tmp)
    rm(tmp)
    load(paste0(dirname(dir.path),"/Settingup_argument.RData"))
    print(dirname(dir.path))
    print(names(arguments))
    print(head(ft@CellFreqs))
    pop.names<-sapply(ft@PhenoCodes, flowTypeFilter::decodePhenotype, marker.names =arguments$names , partitions.per.marker=arguments$part.v,arguments$filt.v)
    print(head(pop.names))
    withProgress(message = 'Analyzing:--->', value = 0, {
      incProgress(.3, detail = "Calculating proportions")
      if (input$calc=="NO")
      {
        ft.counts<-c()
        cnt<-1
        n <- 200
        for ( ftp in ft.path[1:n])
        {
          print(cnt)
          tmp <- load(ftp)
          ft <- get(tmp)
          rm(tmp)
          ft.counts <-rbind( ft.counts, ft@CellFreqs)
          cnt <- cnt+1
        }
        colnames( ft.counts)<- as.vector(pop.names)
        rownames(ft.counts) <- gsub(basename(ft.path)[1:n],pattern = "_FT.Rdata",replacement = "")
       
        
      }else{
        print("here")
        tmp <- load(paste0(dirname(dir.path),"/Allcounts_flowType.RData"))
        ft.counts <- get(tmp)
        rm(tmp)
      }
      incProgress(.5, detail = "Calculating p-values and plotting.")
    })

    output$plotgraph <- renderVisNetwork({
      rchy.graph<- make.sliding.graph (arguments = arguments,ft = ft, counts = ft.counts,
                                       method1=input$method1,method2 = input$method2,parent="CD45+",
                                       count.lim = input$slider1,p.level = input$slider2)
      
      nodes <- rchy.graph$nodes
      print(nodes)
      nodes.id <<-nodes$id
      nodes.lbl <<-nodes$label
      edges <- rchy.graph$edges
      edges.to <<- edges$to
      edges.from <<- edges$from
      visNetwork(nodes, edges) %>%
        visNodes(shape = "ellipse", size=55, label=nodes$label) %>%
        visEdges(arrows = "to")%>%
        visOptions(highlightNearest=TRUE, 
                   nodesIdSelection = TRUE) %>%
        visInteraction(hover = TRUE,multiselect = TRUE) %>%
        visLayout(hierarchical = TRUE) %>%
        visHierarchicalLayout(parentCentralization=TRUE) %>%
        visEvents(select="function(nodes) {
                  Shiny.onInputChange('current_node_selection', nodes.nodes);
                  ;}")
      
    })
    
    selected_sample <- reactive({
      input$plot_selected
    })

    myNode <- reactiveValues(selected = '',pops='')
    observeEvent(input$current_node_selection, {
      myNode$selected <<- as.vector(input$current_node_selection)
      myNode$pops <<-as.vector(nodes.lbl)[match(unlist(myNode$selected),as.vector(nodes.id) )]
      print(paste0("populations selected: ",paste(myNode$pops,collapse=" ,")))
      ft.props <- ft.counts/ft.counts[,1]
      ft.props <- ft.props*100
      lbls <-sample (x = 0:1,size =nrow(ft.counts) ,replace = T)
    print(dim(ft.props))
      pop.prop <-data.frame(props=c(ft.props[which(lbls==0), myNode$pops],
                                    ft.props[which(lbls==1), myNode$pops]), 
                            groups=c(rep(times=length(which(lbls==0)),"Group1"),rep(times=length(which(lbls==1)),"Group2")),
                            files=c(rownames(ft.props)[which(lbls==0)],rownames(ft.props)[which(lbls==1)]))
      sample.props <- rbind(pop.prop[which(pop.prop$groups=="Group1")[which.min(pop.prop$props[which(pop.prop$groups=="Group1")])],],
                            pop.prop[which(pop.prop$groups=="Group1")[which.min(abs(pop.prop$props[which(pop.prop$groups=="Group1")]-
                           median(pop.prop$props[which(pop.prop$groups=="Group1")])))],],
                         pop.prop[which(pop.prop$groups=="Group1")[which.max(pop.prop$props[which(pop.prop$groups=="Group1")])],],
                         pop.prop[which(pop.prop$groups=="Group2")[which.min(pop.prop$props[which(pop.prop$groups=="Group2")])],],
                       pop.prop[which(pop.prop$groups=="Group2")[which.min(abs(pop.prop$props[which(pop.prop$groups=="Group2")]-
                           median(pop.prop$props[which(pop.prop$groups=="Group2")])))],],
                         pop.prop[which(pop.prop$groups=="Group2")[which.max(pop.prop$props[which(pop.prop$groups=="Group2")])],])
      sample.props$files <- paste(sample.props$files,sample.props$groups,sep="_g_")
      #You need to create 
      if (input$random=="YES" | is.null(input$fcs))
      {
        print("Random samples")
      output$boxplot <-renderggiraph({x<-ggiraph(code=print(plot_results(cell.proportion = pop.prop,outliers = sample.props)))
      x <- girafe_options(x, opts_selection(
        type = "multiple", css = "fill:#FF3333;stroke:black;"),
        opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;"))
      x
      })
      }else{
        print("selected sample")
        fcs.file <-paste0(parseFilePaths(roots = volumes,input$fcs),"")
        print(parseFilePaths(roots = volumes,input$fcs)$datapath)
        print(fcs.file)
       
        fcs.samples <- as.character(parseFilePaths(roots = volumes,input$fcs)$name)
        samp.ind <-which(gsub(dirname(ft.path),pattern = "_FT.Rdata",replacement = "")==unlist(strsplit(fcs.samples,split = ".Rdata"))[1])

          tmp <- load(ft.path[samp.ind])
          ft <- get(tmp)
          samp.count<- ft@CellFreqs
          names(samp.count)<-colnames(ft.counts)
         
          sample.gr <- "Group1"
        sample.data <- as.data.frame(props=samp.count[myNode$pops]*100/samp.count[1],files=unlist(strsplit(fcs.samples,split = ".Rdata"))[1],groups=sample.gr)
        output$boxplot <-renderggiraph({x<-ggiraph(code=print(plot_results(cell.proportion = pop.prop,showsample = T,samplestoplotdata = sample.data,
                                                                        outliers = sample.props)))
        x <- girafe_options(x, opts_selection(
          type = "multiple", css = "fill:#FF3333;stroke:black;"),
          opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;"))
        x
        })
        
        
      }
      
      st <-which(edges.to==unlist(myNode$selected))[1]
      print(st)
      selected.node <- unlist(myNode$selected)
      print(selected.node)
      start.path <- c(selected.node,as.vector(edges.from[st]))
      load(paste0(dirname(dir.path),"/Settingup_argument.RData"))
      nxt<- tail(start.path,1)
      pop.path <- start.path
      while(nxt!=as.vector(edges.from[1]))
      {
        prnt <- which(edges.to==nxt)[1]
        nxt<- as.vector(edges.from[prnt])
        pop.path <-c( pop.path ,as.vector(nodes.id[which(nodes.id==nxt)]))
      }
      parent.names <- rev(nodes.lbl[match(pop.path,nodes.id)])
      parent.names <- parent.names[-length(parent.names)]
      pop.path <- rev(pop.path)[-1]
  
     
        fcs.file <-paste0(parseFilePaths(roots = volumes,input$fcs),"")
        print(parseFilePaths(roots = volumes,input$fcs)$datapath)
        print(fcs.file)

        fcs.samples <- as.character(parseFilePaths(roots = volumes,input$fcs)$name)
        samp.ind <- which(names(arguments$threshold)==unlist(strsplit(fcs.samples,split = ".Rdata"))[1])
        gates <-  arguments$threshold[[samp.ind]]
        temp <- load(parseFilePaths(roots = volumes,input$fcs)$datapath)
        frame <- get(temp)
        rm(temp)
        # 
      
        click.input <- selected_sample()
        click.gr <- unlist(lapply(click.input, function(cl) unlist(strsplit(cl,"_g_"))[2]))
        click.id <- unlist(lapply(click.input, function(cl) unlist(strsplit(cl,"_g_"))[1]))
        
       
        if (input$random=="YES"){
        output$plotdens <- renderPlot({
          click.input <- input$boxplot_selected
          print(click.input)
          
          if(!is.null(click.input))
          {
            click.gr <- unlist(lapply(click.input, function(cl) unlist(strsplit(cl,"_g_"))[2]))
            click.id <- unlist(lapply(click.input, function(cl) unlist(strsplit(cl,"_g_"))[1]))
            print(paste("click is:",length(click.input)))
            print(click.id)
            print(click.gr)
            tmp <- load(ft.path[2])
            ft <- get(tmp)
            load(paste0(dirname(dir.path),"/Settingup_argument.RData"))
            gs <- make.gs(ft = ft, ids = tail(click.id,2),arguments = arguments,pop.path = pop.path)
            print(class(gs))
            vec <- c(1:(ceiling(sqrt(length(pop.path)))*2))
            vec[which(vec>=length(gs_get_pop_paths(gs)))] <-0
            mat <-matrix(vec,nrow =ceiling(sqrt(length(pop.path))) ,byrow = T)
            vec2 <- vec + max(vec,na.rm = T)
            vec2[which(vec2==max(vec))]<-0
            mat2 <-cbind(mat,matrix(vec2,nrow =ceiling(sqrt(length(pop.path))) ,byrow = T))
            source("~/code/plot.gate.R")
            if (length(gs)==2)
            {
              print("OK")
              print(mat2)
              plot.gates (gs,popstoplot=NULL,upto=NULL,plot.layout=mat2,show.stats=TRUE,contour.pops=NULL,density.overlay.pops=NULL,
                                     flowCut.channels=c("Time","BV510-A"),legend.cex=1.5,cell.size=c(1,2.5,3),show.palette=TRUE,offset=.5)
              title(paste(click.gr,sep=" ",collapse = "      -        "), outer = TRUE, line = 2)

              
              }else{
              print("Fine")
                print(mat)
             plot.gates (gs,popstoplot=NULL,upto=NULL,plot.layout=mat,show.stats=TRUE,contour.pops=NULL,density.overlay.pops=NULL,
                                     flowCut.channels=c("Time","BV510-A"),legend.cex=1.5,cell.size=c(1,2.5,3),show.palette=TRUE,offset=.5)
             title(paste(click.gr,sep="  ",collapse = "  -   "), outer = TRUE, line = 2)

             
             }
            
            
            
          }
        })
        }
        
    })
  })
  
}
