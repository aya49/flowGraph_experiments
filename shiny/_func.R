# given a gating set, gates (list of flowType@Thresholds), a vector of phenocodes on a gating strategy path, and a marker to start gating plots on,
# outputs the gating strategy with gates for plotting

make_gs <- function(gs, gates, markers=NULL, gating_path_phenocode, scat_chans) {
    fs <- flowWorkspace::gs_pop_get_data(gs)
    if (base::is.null(markers)) markers <- names(gates[[1]])
    gating_path_phenocode_ <- lapply(gating_path_phenocode, function(p)
        unlist(strsplit(p,"")))

    alt_gate <- scat_chans
    for (i1 in 1:length(gating_path_phenocode_)) {
        cpop_vec <- gating_path_phenocode_[[i1]]
        ind_m <- which(cpop_vec!="0")
        if (names(gates[[1]])[ind_m]!="gate") {
            gatenum <- as.numeric(cpop_vec[ind_m])
            fdens <- flowCore::fsApply(fs, function(fcs) {
                fcs_id <- flowCore::identifier(fcs)
                fden <- flowDensity::flowDensity(
                    fcs,
                    channels=c(names(gates)[ind_m], alt_gate),
                    position=c(gatenum>1,NA),
                    gates=c(gates[[fcs_id]][[ind_m]][max(1,gatenum-1)],NA))
                if (gatenum>1 & length(gates[[ind_m]])>(gatenum-1))
                    fden <- flowDensity::flowDensity(
                        fden@flow.frame,
                        channels=c(names(gates)[ind_m], alt_gate),
                        position=c(FALSE, NA),
                        gates=c(gates[[fcs_id]][[ind_m]][gatenum],NA))
                return(fden)
            })
            filter_id <-
                paste0(markers[ind_m],
                       ifelse(cpop_vec[ind_m]!="1",
                              rep("+",gatenum-1), "-"))
        } else {
            if (cpop_vec[ind_m]=="1") {
                fdens_fun <- flowDensity::notSubFrame
            } else {
                fdens_fun <- flowDensity::flowDensity
            }
            fdens <- flowCore::fsApply(fs, function(fcs) {
                fcs_id <- flowCore::identifier(fcs)
                fdens_fun(
                    fcs,
                    channels=colnames(gates[[fcs_id]][[ind_m]]),
                    position=c(TRUE, TRUE),
                    filter=as.matrix(gates[[fcs_id]][[ind_m]]))
            })
            filter_id <- paste0(
                ifelse(cpop_vec[ind_m]=="1", "Not ", ""),
                markers[ind_m])
        }

        # mark off this marker for future gates
        if (i1<length(gating_path_phenocode_))
            for (k1 in (i1+1):length(gating_path_phenocode_))
                gating_path_phenocode_[[k1]][ind_m] <- "0"

        alt_gate <- names(gates)[ind_m]

        # apply to gating set
        poly <- lapply(fdens, function(fden)
            flowCore::polygonGate(filterId=filter_id,.gate=fden@filter))
        names(poly) <- flowWorkspace::sampleNames(gs)

        nodeID <- flowWorkspace::gs_pop_add(
            gs, poly, parent=tail(flowWorkspace::gs_get_pop_paths(gs),1))
        flowWorkspace::recompute(gs)
    }
    return(gs)
}
