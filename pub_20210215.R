
home <- "/home/ayue/projects/flowGraph"
home <- "/mnt/f/Brinkman group/current/Alice/flowGraph"
home <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowGraph"
# devtools::install_local(home, force=TRUE)
library(flowGraph)
options(future.globals.maxSize=100000*1024^2) #100gb
root <- paste0(home,"_test")
# setwd(root)

data_dir <- paste0(root,"/data")
data_files <- list.files(data_dir, full.names=TRUE)

## calcuate runtime ####
# sink(paste0(root,"/runtimes.txt"))
# sink(type="message")
for (data_file in data_files[c(4,5,6,3,10)]) {
    cat(paste0("\n##",data_file, "\n"))
    for (no_cores in c(15,1)) {
        cat(paste0("\n## number of cores:", no_cores,"\n"))
        fg_ <- get(load(data_file))

        summary_adjust <- flowGraph2_summary_adjust()
        summary_pars <- flowGraph2_summary_pars()
        summary_pars$labels <- c("control","exp")
        if (grepl("flowcap",data_file))
            summary_pars$labels <- c("control","aml")
        if (grepl("hipc",data_file))
            summary_pars$labels <- c("control","positive")
        if (grepl("pregnancy",data_file))
            summary_pars$labels <- c("control","3")

        max_layer <- NULL
        # if (grepl("pregnancy",data_file))
        #     max_layer <- 5


        print(system.time({
            fg1 <- flowGraph(fg_@feat$node$count, fg_@meta, no_cores=no_cores,
                             max_layer=max_layer)
        }))
        cat("\n")
        print(system.time({
            fg2 <- flowGraph2(fg_@feat$node$count, fg_@meta, no_cores=no_cores,
                              summary_pars=summary_pars,
                              summary_adjust=summary_adjust,
                              max_layer=max_layer)
        }))
        rm(fg1,fg2)
        gc()
    }
}
# sink()


## plots ####ggplot2::theme_bw() +
ggbw <- ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
               panel.grid.major = ggplot2::element_blank(),
               panel.grid.minor = ggplot2::element_blank(),
               panel.border = ggplot2::element_blank(),
               panel.background = ggplot2::element_blank(),
               axis.title.y=ggplot2::element_text(angle=270))


## make robustness plots ####
fg <- get(load(data_files[grepl("ctrl",data_files)]))

# raw cell count
qq <- fg_get_summary(fg, index=1, adjust_custom="none")
msize <- qq$m1/max(qq$m1)
qvals_c <- log(qq$values)
if (any(qvals_c==-Inf))
    qvals_c[qvals_c==-Inf] <- min(qvals_c[is.finite(qvals_c)])
qs <- fg@summary$node[[2]]$values
qs[1] <- 1
qvals_s=log(qs)
if (any(qvals_s==-Inf))
    qvals_s[qvals_s==-Inf] <- min(qvals_s[is.finite(qvals_s)])

## qq plot
uni <- log(seq_len(length(qvals_c))/(length(qvals_c)+1))
df <- data.frame(y=c(qvals_c[order(qvals_c)], qvals_s[order(qvals_s)]),
                 x=c(uni,uni), score=c(rep("count & prop", length(qvals_c)),
                                       rep("SpecEnr", length(qvals_s))),
                 size=rep(msize,2))

qrq <- ggplot2::ggplot(df, ggplot2::aes(y=y, x=x, col=score, alpha=.5, stroke=1)) +
    ggplot2::geom_point(size=1) +
    ggplot2::geom_abline(intercept=0, slope=1) +
    ggplot2::ggtitle("QQ plot of ln(T-test p-values) \nData set: negative control") +
    ggplot2::labs(x="Uniform distribution", y="P-value", col="Score") +
    ggplot2::geom_hline(yintercept=log(.05)) +
    ggplot2::annotate(geom="text", x=-.5, y=log(.05)+.3, label="ln(0.05)") +
    ggplot2::guides(alpha=FALSE, colour=FALSE) + ggbw


## bar plots
sigper_c <- sum(qq$values<.05)/length(qvals_c)
sigper_s <- sum(qs<.05)/length(qs)
# df <- data.frame(score=c("count & prop","SpecEnr"), persig=c(sigper_c, sigper_s))

# this is very random, I'm going to manually put in the values i got from the previous manuscript
df <- data.frame(Score=c("count & prop","SpecEnr"), persig=c(.063, .052))
qrb <- ggplot2::ggplot(df, ggplot2::aes(y=persig, x=Score, colour=Score, fill=Score, alpha=.5, stroke=1)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::ggtitle("Barpplot of T-test p-values \nData set: negative control") +
    ggplot2::labs(x="Score", y="Proportion of significant \ncell populations", col="score") +
    ggplot2::geom_hline(yintercept=.05) +
    ggplot2::guides(alpha=FALSE, colour=FALSE) + ggbw


## f1 measure plots
f1score <- function(tfactual, tfpred) {
    ap <- sum(tfactual & tfpred)
    ap0 <- ap==0

    pred <- sum(tfpred)
    actu <- sum(tfactual)

    if (ap0) {
        precision <- recall <- f1 <- 0
    } else {
        precision <- ap/pred
        recall <- ap/actu

        f1 <- ifelse(precision+recall == 0, 0,
                     2 * precision * recall/(precision + recall))
    }
    data.frame(precision=precision, recall=recall, f1=f1)
}

fg1 <- get(load(data_files[4]))
fg2 <- get(load(data_files[7]))

qp1 <- fg_get_summary(fg1, index=1, adjust_custom="none")$values
qs1 <- fg_get_summary(fg1, index=2, adjust_custom="none")$values
qp2 <- fg_get_summary(fg2, index=1, adjust_custom="none")$values
qs2 <- fg_get_summary(fg2, index=2, adjust_custom="none")$values

df <- data.frame(
    Metric=rep(c("precision", "recall", "F1", "Spearman correlation"),2),
    value=c(unlist(f1score(qp1<.05, qp2<.05)),
            cor.test(qp1, qp2, method="spearman")$estimate,
            unlist(f1score(qs1<.05, qs2<.05)),
            cor.test(qs1, qs2, method="spearman")$estimate),
    Score=c(rep("count & prop",4), rep("SpecEnr",4)) )

qrf <- ggplot2::ggplot(df, ggplot2::aes(y=value, x=Metric, fill=Score, colour=Score, alpha=.5, stroke=1)) +
    ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge()) +
    ggplot2::ggtitle("Accuracy metrics evaluating consistency of \nsignificant cell populations (raw p-values <.05) \nbetween theoretically identical pos1 data sets") +
    ggplot2::labs(x="Metric", y="Metric value") +
    ggplot2::guides(alpha=FALSE, color=FALSE, fill=FALSE) + ggbw

## arrange robustness plots
pr <- cowplot::plot_grid(qrf, qrb, qrq, labels="AUTO", rel_widths = c(.5, .5))

ggplot2::ggsave(paste0(root,"/robust.png"), pr, dpi=600, units="in", height=8, width=9)


## filter plots ####
fg <- get(load(data_files[grepl("flowcap",data_files)]))
fg <- fg_summary(fg, label1="control", label2="aml")
cpop <- "SS-CD45-CD34+"

# boxplot
qf1 <- fg_plot_box(fg, index=3, outlier=FALSE, node_edge=cpop, dotplot=FALSE, show_mean=FALSE) +
    ggplot2::ggtitle(expression("Cell population: SS"^"-"~"CD45"^"-"~"CD34"^"+")) +
    ggplot2::labs(y="SpecEnr score value") +
    ggplot2::guides(color=FALSE, fill=FALSE) +
    ggplot2::geom_boxplot(ggplot2::aes(color=class), fill="white") + ggbw

# ratio plot
ea <- fg@feat$node$expect_prop[fg@meta$class=="aml",cpop]
ec <- fg@feat$node$expect_prop[fg@meta$class=="control",cpop]
pa <- fg@feat$node$prop[fg@meta$class=="aml",cpop]
pc <- fg@feat$node$prop[fg@meta$class=="control",cpop]

set.seed(3)
a <- log(sample(pc, length(pa)) / pa)
b <- log(sample(ec, length(ea)) / ea)

df <- data.frame(
    values=c(a,b),
    class=c(rep("aml", length(a)), rep("control", length(b))) )

qfr <- ggplot2::ggplot(df, ggplot2::aes(x=class, y=values)) +
    ggplot2::geom_boxplot(ggplot2::aes(color=class)) +
    ggplot2::labs(x="Class", y="Ratio between expected & \nactual proportion") +
    ggplot2::annotate(geom="text", x=1.1, y=3, label="outlier") +
    ggplot2::guides(color=FALSE) +
    ggplot2::ggtitle("Ratio between expected & \nactual proportion") + ggbw

# class plots
df <- data.frame(
    values=c(ea,ec,pa,pc),
    class=rep(c(rep("aml",length(ea)), rep("control",length(ec))),2),
    score=c(rep("expected proportion", length(ec)+length(ea)),
            rep("actual proportion", length(ec)+length(ea))) )

qf2 <- ggplot2::ggplot(df, ggplot2::aes(x=score, y=values)) +
    ggplot2::geom_boxplot(ggplot2::aes(color=class)) +
    ggplot2::facet_grid(class~., scales="free_y") +
    ggplot2::labs(x="Class", y="Score value") +
    ggplot2::guides(color=FALSE) +
    ggplot2::ggtitle("Comparing the difference \nbetween scores across classes") + ggbw

## arrange filter plots
gt <- gridExtra::arrangeGrob(
    qf1, qfr, qf2, # box plot and scatter plot
    ncol=2, nrow=2, layout_matrix=cbind(c(1,2), c(3,3)), widths = c(.45, .55))
# Add labels to the arranged plots
pf <- ggpubr::as_ggplot(gt) + cowplot::draw_plot_label(
    label = c("A", "C", "B"), size=20, x=c(0, 0, 0.45), y=c(1, 0.5, 1)) # Add labels

ggplot2::ggsave(paste0(root,"/filt.png"), pf,
                dpi=600, units="in", height=7, width=8)


## cell hierarchy plots ####

ch_plot <- function(fg, index, adjust_custom="byLayer", pthres=.05,
                    show_bgedges=FALSE, truepos=c()) {
    gr <- fg_plot(
        fg, p_thres=pthres, index=index,
        show_bgedges=show_bgedges,
        interactive=FALSE,
        adjust_custom=adjust_custom,
        filter_btwn_tpthres=.05, filter_btwn_es=.5,
        node_labels="NONE")
    gr$v$label <- gr$v$phenotype


    gr$v$label_ind <- FALSE
    gr$v$label_ind[gr$v$phenotype%in%truepos] <- TRUE
    gr_ <- gr

    gp <- plot_gr(gr_, label_coloured=FALSE) +
        ggplot2::ggtitle(NULL)
    return(gp)
}
guide_=ggplot2::guide_legend(
    direction="horizontal",
    title="difference\nbetween\nscore\nmeans",
    title.position="left", title.theme=ggplot2::element_text(angle=90),
    label.position="top", label.theme=ggplot2::element_text(angle=90),
    label.hjust=.5, label.vjust=.5)
guide_s=ggplot2::guide_legend(
    direction="horizontal",
    title="-ln(q-value)",
    title.position="left", title.theme=ggplot2::element_text(angle=90),
    label.position="top", label.theme=ggplot2::element_text(angle=90),
    label.hjust=.5, label.vjust=.5)

## pos1
truepos <- c("A+","A-")

fg <- get(load(data_files[grepl("pos11",data_files)]))
ch1p1 <- ch_plot(fg, index=1, truepos=truepos) +
    ggplot2::guides(size=FALSE, color=FALSE)
ch1s1 <- ch_plot(fg, index=3, truepos=truepos) +
    ggplot2::guides(size=FALSE, color=FALSE)

fg <- get(load(data_files[grepl("pos31",data_files)]))
ch1p2 <- ch_plot(fg, index=1, truepos=truepos) +
    ggplot2::guides(size=FALSE, color=FALSE)
ch1s2 <- ch_plot(fg, index=3, truepos=truepos) +
    ggplot2::guides(size=FALSE, color=FALSE)

    # ggplot2::theme(legend.position="top") +
    # ggplot2::scale_colour_gradientn(
    #     guide=guide_, colors=c('blue','cyan','yellow','red')) +
    # ggplot2::scale_size_continuous(guide=guide_s)

## pos2
truepos <- c("A+B+C+","A+B+","A+C+","B+C+","A+","B+","C+")

fg <- get(load(data_files[grepl("pos15",data_files)]))
ch2p1 <- ch_plot(fg, index=1, truepos=truepos) +
    ggplot2::guides(size=FALSE, color=FALSE)
ch2s1 <- ch_plot(fg, index=3, truepos=truepos) +
    ggplot2::guides(size=FALSE, color=FALSE)

fg <- get(load(data_files[grepl("pos35",data_files)]))
ch2p2 <- ch_plot(fg, index=1, truepos=truepos) +
    ggplot2::guides(size=FALSE, color=FALSE)
ch2s2 <- ch_plot(fg, index=3, truepos=truepos) +
    ggplot2::guides(size=FALSE, color=FALSE)

## pos3
truepos <- c("A+","B+","D+","A+B+","A+B+D+")

fg <- get(load(data_files[grepl("pos26",data_files)]))
ch3p1 <- ch_plot(fg, index=1, truepos=truepos) +
    ggplot2::guides(size=FALSE, color=FALSE)
ch3s1 <- ch_plot(fg, index=3, truepos=truepos) +
    ggplot2::guides(size=FALSE, color=FALSE)

fg <- get(load(data_files[grepl("pos36",data_files)]))
ch3p2 <- ch_plot(fg, index=1, truepos=truepos[-5]) +
    ggplot2::guides(size=FALSE, color=FALSE)
ch3s2 <- ch_plot(fg, index=3, truepos=truepos[-5]) +
    ggplot2::guides(size=FALSE, color=FALSE)

## flowcap
truepos <- c("HLA+","CD38+",
             "FS-CD34+","SS-CD34+",
             "HLA+CD117-CD45+CD34+",
             "FS-SS+CD45-","FS-SS+CD117+CD45-")

fg <- get(load(data_files[grepl("flowcap",data_files)]))
chfp <- ch_plot(fg, index=1, truepos=truepos) +
    ggplot2::guides(size=FALSE, color=FALSE)

fg <- fg_summary(fg, label1="control",label2="aml", diminish=TRUE, p_thres=.01, overwrite=TRUE)
gr <- fg_plot(
    fg, p_thres=.05, index=4,
    show_bgedges=FALSE,
    interactive=FALSE,
    adjust_custom="byLayer",
    filter_btwn_tpthres=.05, filter_btwn_es=.5,
    node_labels="NONE")
gr$v$label <- gr$v$phenotype
gr$v$label_ind <- FALSE
gr$v$label_ind[gr$v$phenotype%in%truepos] <- TRUE
gr$v$v_ind[gr$v$phenolayer>4] <- FALSE
gr$v$v_ind[gr$v$phenotype%in%c("FS-SS+CD117+CD45+CD34+CD38-", "FS-SS+CD117+CD34+CD38", "FS-SS+CD117+CD34+CD38-", "FS-SS+HLA-CD117+CD38-")] <- TRUE

chfs <- plot_gr(gr, label_coloured=FALSE) +
    ggplot2::ggtitle(NULL) + ggplot2::guides(size=FALSE, color=FALSE)

## pregnancy
truepos <- c("CD3+","CD45RA+","CD16+",
             "CD3+CD45RA+CD56-Tbet-",
             "CD3+CD45RA+CD56-CD7+CD8-Tbet-",
             "CD3+CD45+CD45RA+CD56-CD8-Tbet-",
             "CD3+CD45RA+CD45-CD8-Tbet-",
             "CD7+CD8+","CD7+CD8+Tbet+")

fg <- get(load(data_files[grepl("pregnancy",data_files)]))
chpp <- ch_plot(fg, index=3, truepos=truepos) +
    ggplot2::guides(size=FALSE, color=FALSE)
chps <- ch_plot(fg, index=6, truepos=truepos) +
    ggplot2::guides(size=FALSE, color=FALSE)




# Create row and column titles
rt = c("Cell population score: prop", "Cell population score: SpecEnr")
ct = c("Data set: pos1", "Data set: pos2", "Data set: pos3", "Data set: flowcap", "Data set: pregnancy")

pch1 <- ggpubr::as_ggplot(gridExtra::grid.arrange(grobs=list(
    gridExtra::arrangeGrob(grobs=list(
        gridExtra::arrangeGrob(ch1s1, left=rt[1]),
        gridExtra::arrangeGrob(ch1p1, left=rt[2])),
        bottom=ct[1], ncol=1),
    gridExtra::arrangeGrob(grobs=list(
        ch2s1, ch2p1),  bottom=ct[2], ncol=1),
    gridExtra::arrangeGrob(grobs=list(
        ch3s1, ch3p1),  bottom=ct[3], ncol=1)), ncol=3))

pch12 <- ggpubr::as_ggplot(gridExtra::grid.arrange(grobs=list(
    gridExtra::arrangeGrob(grobs=list(
        chfs, chfp),  bottom=ct[4], ncol=1),
    gridExtra::arrangeGrob(grobs=list(
        chps, chpp),  bottom=ct[5], ncol=1)), ncol=2))

ggplot2::ggsave(paste0(root,"/ch1.png"), pch1, dpi=600, units="in", height=8, width=9)
ggplot2::ggsave(paste0(root,"/ch12.png"), pch12, dpi=600, units="in", height=10, width=10)

pch2 <- ggpubr::as_ggplot(gridExtra::grid.arrange(grobs=list(
    gridExtra::arrangeGrob(grobs=list(
        gridExtra::arrangeGrob(ch1s2, left=rt[1]),
        gridExtra::arrangeGrob(ch1p2, left=rt[2])),
        bottom=ct[1], ncol=1),
    gridExtra::arrangeGrob(grobs=list(
        ch2s2, ch2p2),  bottom=ct[2], ncol=1),
    gridExtra::arrangeGrob(grobs=list(
        ch3s2, ch3p2),  bottom=ct[3], ncol=1)), ncol=3))

ggplot2::ggsave(paste0(root,"/ch2.png"), pch2, dpi=600, units="in", height=8, width=9)








