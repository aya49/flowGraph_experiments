# these are plots generated for a report on 2020-04-27 to match specenr results with sebastiano's cells per ml of blood ones

result_dir <- "/mnt/f/Brinkman group/current/Alice/flowtype_metric/result/epichipc_bcell"
fg <- get(load(paste0(result_dir,"/fg.Rdata")))

show(fg)

# sebastiano's cell poulations
cpops <- c("CD38+CD138+IgM+",
           "CD38+CD138+IgD-IgM+CD10-",
           "CD38+CD138+IgD-IgM+CD27-CD10-",
           "CD38+CD138+IgD-IgM+CD10-CD38CD10-",
           "CD38-CD138-CD20+IgM-CD27-CD10-")

# top specenr filtered cell populations
cpops <- c("CD38+CD138+IgM+CD10+",
           "CD66+CD38+CD138+SSCACD19+",
           "CD66-CD38-CD138-IgM-CD27-CD10+",
           "CD38-CD138-CD20-IgD+CD27-",
           "CD66+CD38+IgM+CD10+",
           "CD38+CD138+CD20+IgM-CD27-CD10-")

# boxplots
for (cpop in cpops) {
    fg_plot_box(fg, index=3, node_edge=cpop,
                outlier=T, dotplot=F,
                path=paste0(result_dir,"/temp/", cpop, "_3_out.png"))
    fg_plot_box(fg, index=3, node_edge=cpop,
                outlier=F, dotplot=F,
                path=paste0(result_dir,"/temp/", cpop, "_3_noout.png"))
    fg_plot_box(fg, index=4, node_edge=cpop,
                outlier=T, dotplot=F,
                path=paste0(result_dir,"/temp/", cpop, "_4_out.png"))
    fg_plot_box(fg, index=4, node_edge=cpop,
                outlier=F, dotplot=F,
                path=paste0(result_dir,"/temp/", cpop, "_4_noout.png"))
}


# relationship between median0 and adjust0
fg_adjust0_ <- function(m, id1, id2, adjust0_lim=c(-.1, .1)) {
    m <- as.matrix(m)
    xl <- min(sum(id1), length(id1))
    yl <- min(sum(id2), length(id2))
    purrr::map_dbl(seq_len(ncol(m)), function(ci)
        min(sum(m[id1,ci]<adjust0_lim[2] & m[id1,ci]>adjust0_lim[1])/xl,
            sum(m[id2,ci]<adjust0_lim[2] & m[id2,ci]>adjust0_lim[1])/yl))
}

fg_med0_ <- function(m, id1, id2) {
    m <- as.matrix(m)
    purrr::map_dbl(seq_len(ncol(m)), function(ci)
        max(abs(c(median(m[id1,ci]), median(m[id2,ci])))))
}


pp <- fg_get_summary(fg, index=4)
med0 <- fg_med0_(fg@feat$node$SpecEnr_cells_ul_blood, pp$id1, pp$id2)
a0 <- fg_adjust0_(fg@feat$node$SpecEnr_cells_ul_blood, pp$id1, pp$id2)
names(med0) <- names(a0) <- fg@graph$v$phenotype

plot(sort(a0))
plot(sort(med0))

plot(a0[order(a0)], col=ifelse(med0[order(a0)]<.1,"red","blue"))


gr <- fg_plot(fg, index=4, filter_adjust0=.4, filter_es=2, p_thres=.6)
gr$v$v_ind <- F

vs <- fg@edge_list$parent[["CD38+CD138+IgM+CD10+"]]
vs <- append(vs, fg@edge_list$child[["CD38+CD138+IgM+CD10+"]])
vs <- append(vs, fg@edge_list$child[["CD38+CD138+IgM+CD10-"]])
vs <- append(vs, c("CD38+CD138+IgM+CD10+", "CD38+CD138+IgM+CD10-"))
vs <- unique(vs)
gr$sv$v_ind <- gr$v$phenotype%in%vs



