install.packages('igraph')
install.packages('network') 
install.packages('sna')
install.packages('ndtv')
install.packages('visNetwork')
install.packages('GGally')

libr('igraph')
libr('network') 
libr('sna')
libr('ndtv')
libr('visNetwork')
libr('GGally')
libr('ggplot2')





# random graph
net = rgraph(10, mode = "graph", tprob = 0.5)
net = network(net, directed = FALSE)

# vertex names
network.vertex.names(net) = letters[1:10]

# vertex colour/size
ggnet2(net, node.size = 6, node.color = "black", edge.size = 1, edge.color = "grey")
ggnet2(net, size = 6, color = "black", edge.size = 1, edge.color = "grey")
ggnet2(net, size = 6, color = rep(c("tomato", "steelblue"), 5)) #pass vector of colours
# This attribute can be passed to ggnet2 to indicate that the nodes belong to a group. All the user has to do is to pass the name of the vertex attribute to the color argument, which will find it in the list of vertex attributes and use it to map the colors of the nodes
net %v% "phono" = ifelse(letters[1:10] %in% c("a", "e", "i"), "vowel", "consonant")
ggnet2(net, color = "phono")
net %v% "color" = ifelse(net %v% "phono" == "vowel", "steelblue", "tomato")
ggnet2(net, color = "color")
ggnet2(net, color = "phono", palette = c("vowel" = "steelblue", "consonant" = "tomato")) #pass palettes
ggnet2(net, color = ifelse(net %v% "phono" == "vowel", "steelblue", "tomato"))

# vertex size
ggnet2(net, size = "phono")
ggnet2(net, size = "phono", size.palette = c("vowel" = 10, "consonant" = 1))
ggnet2(net, size = sample(0:2, 10, replace = TRUE), max_size = 9) # define only max size
ggnet2(net, size = "degree") # define by indegree
ggnet2(net, size = "degree", size.cut = 3) # cuts size into quartiles

# remove any isolated nodes
x = ggnet2(net, size = "degree", size.min = 1)
# remove all nodes
x = ggnet2(net, size = "degree", size.max = 1)
ggnet2(net, size = sample(0:2, 10, replace = TRUE), size.zero = TRUE) # if want 0-sized nodes plotted

# vertex placement (default: Fruchterman-Reingold force-directed algorithm)
ggnet2(net, mode = "circle")
ggnet2(net, mode = "kamadakawai")
ggnet2(net, mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75)) #arguments of algorithm
ggnet2(net, mode = "target", layout.par = list(niter = 100))
ggnet2(net, color = "phono", palette = "Set2")


# legend
ggnet2(net, alpha = "phono", alpha.legend = "Phonetics")
ggnet2(net, shape = "phono", shape.legend = "Phonetics")
ggnet2(net, color = "phono", color.legend = "Phonetics")
ggnet2(net, size = "degree", size.legend = "Centrality")
# additional legends have to be discrete_scale controllers, even when the scale applies to the size of the nodes

# control the colors of the nodes
ggnet2(net, color = "phono") +
  scale_color_brewer("", palette = "Set1",
                     labels = c("consonant" = "C", "vowel" = "V"),
                     guide = guide_legend(override.aes = list(size = 6)))
# control the size of the nodes
ggnet2(net, size = "degree") +
  scale_size_discrete("", range = c(5, 10), breaks = seq(10, 2, -2))
# style by theme of plot
ggnet2(net, color = "phono", legend.size = 12, legend.position = "bottom") +
  theme(panel.background = element_rect(color = "grey"))

# node labels
ggnet2(net, label = TRUE)
ggnet2(net, label = "phono")
ggnet2(net, label = 1:10)
# If label is a vector of values that does not contain exactly as many elements as the number of nodes in the graph, ggnet2 will label the nodes that match one of these values
ggnet2(net, label = c("a", "e", "i"), color = "phono", label.color = "black")
ggnet2(net, size = 12, label = TRUE, color = "black", label.color = "white")
ggnet2(net, size = 12, label = TRUE, label.size = 5) # control label size
ggnet2(net, label = TRUE, label.alpha = 0.75)
ggnet2(net, color = "grey15", size = 12, label = TRUE, label.color = "color") +
  theme(panel.background = element_rect(fill = "grey15")) # accept vector of labels and can change bg colour

# vertex shape
ggnet2(net, color = "phono", shape = 15)
ggnet2(net, color = "phono", shape = "phono")
ggnet2(net, alpha = "phono", alpha.palette = c("vowel" = 0.2, "consonant" = 1)) # can take manual palettes of shape and alpha
ggnet2(net, shape = "phono", shape.palette = c("vowel" = 19, "consonant" = 15))
ggnet2(net, shape = sample(1:10))
ggnet2(net, alpha = "phono")




## Example (2): Bipartite network
# weighted adjacency matrix
bip = data.frame(event1 = c(1, 2, 1, 0),
                 event2 = c(0, 0, 3, 0),
                 event3 = c(1, 1, 0, 4),
                 row.names = letters[1:4])
# weighted bipartite network
bip = network(bip,
              matrix.type = "bipartite",
              ignore.eval = FALSE,
              names.eval = "weights")

ggnet2(bip, label = TRUE) # default
col = c("actor" = "grey", "event" = "gold") # set colors for each mode
ggnet2(bip, color = "mode", palette = col, label = TRUE) # detect and color the mode

# edge labels
ggnet2(bip, color = "mode", palette = col, label = TRUE, edge.label = "weights")
ggnet2(bip, shape = "mode", edge.label = "weights", edge.label.color = "darkred")
ggnet2(bip, shape = "mode", edge.label = "weights", edge.label.size = 6)
set.edge.attribute(bip, "color", ifelse(bip %e% "weights" > 1, "black", "grey75"))
ggnet2(bip, shape = "mode", edge.label = "weights", edge.label.color = "color")
ggnet2(bip, shape = "mode", edge.label = "weights", edge.label.fill = NA) # color of that background, which is draw as a circle with geom_point, can be styled with edge.label.fill
# edge size/ colour
ggnet2(bip, color = "mode", palette = col, edge.size = "weights")
set.edge.attribute(bip, "color", ifelse(bip %e% "weights" > 1, "black", "grey75"))
ggnet2(bip, color = "mode", palette = col, edge.size = "weights", edge.color = "color")
# edge type
set.edge.attribute(bip, "lty", ifelse(bip %e% "weights" > 1, 1, 2))
ggnet2(bip, color = "mode", palette = col, edge.size = "weights", edge.lty = "lty")
# edge arrows (default hides them, workaround is to make edges shorter)
ggnet2(network(rgraph(10, tprob = 0.25), directed = TRUE),
       arrow.size = 12, arrow.gap = 0.025)

# ggnet2 supports this functionality by allowing the edge.color argument to take the c("color", "grey") value. The first value will tell  ggnet2 to color edges between nodes of the same group with the color of that group. The second value is the color to use for edges that connect nodes belonging to different groups.
ggnet2(net, color = "phono", palette = "Set1", edge.color = c("color", "grey50"))













## plotting pies as nodes
# http://stackoverflow.com/questions/22587959/using-pie-charts-as-vertices-in-graph-plots-and-specify-a-color-for-each-factor
libr(igraph)

g <- graph.ring(10)

values <- lapply(1:10, function(x) c(sample(1:10,3),0,0))

# make some unique bits
values[[7]][5] = 6
values[[9]][4] = 3

# default for all
V(g)$pie.color=list(heat.colors(5))

# make one stand out
V(g)[6]$pie.color=list(topo.colors(5))

# set.seed() keeps the vertices in the same place each plot
set.seed(1492)

plot(g,
     vertex.shape="pie", 
     vertex.pie=values,
     vertex.size=seq(10, 30, length=10), 
     vertex.label=NA)

