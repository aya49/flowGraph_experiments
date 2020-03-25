
## build flowGraph library
setwd("/home/joann/EPICHIPC/.Ruserdata/ayue/flowGraph")

## this set of data must be done on the server; 4 cores!
no_cores = 1
## import libraries and connection to AWS system ###################
# setwd("/home/joann/EPICHIPC/.Ruserdata/smontante")
source("/home/joann/EPICHIPC/.Ruserdata/smontante/Utils.R")
source("/home/joann/EPICHIPC/.Ruserdata/ayue/src/00_data_epichipc_functions.R")

#----- Call libraries -----
#*** load data
library(aws.s3)
library(aws.signature)
library(aws.ec2metadata)
#***

library(data.table)
library(parallel)
library(stringr)
library(readxl)
library(RchyOptimyx)
library(plyr)
library(qvalue)

#---- Call the lists of the Amazon S3 buckets ----
bucketlist() #*** buckets are a large folder, we are interested in the epichipc-main folder



# -------------- get flowtype results paths -----
#*** only save data in this folder; don't touch Rdata/integration; can't delete files on the server, must email sofia
#*** this function lists all the paths in a folder
paths_all_temp <- get_bucket_df(
    "epichipc-main", #list bucket name
    # prefix = "Flow_Cytometry/Rdata", #list the name of the specific folder
    # prefix = "Flow_Cytometry/Rdata/FlowAnalysis/FlowType_results", # all flowtype results
    prefix = "Flow_Cytometry/Rdata/FlowAnalysis/Alice_Yue_data", #list the name of the specific folder
    max = Inf #this allows an infinite number of objects (default is  1000)
)
paths_all_temp <- paths_all_temp$Key
# [1] "Flow_Cytometry/Rdata/FlowAnalysis/Alice_Yue_data/"
# [2] "Flow_Cytometry/Rdata/FlowAnalysis/Alice_Yue_data/flowType_results_response_Bcells.csv"
# [3] "Flow_Cytometry/Rdata/FlowAnalysis/Alice_Yue_data/rchyoptimyx_df_Bcells_full_data_cells_ul_blood.csv"

# PREPROSSED flowtype data for all samples (0s removed)
# class: vaccine_response
## there is a BCELL and a MYELOID panel
## only samples with above 2.50ml antibody concentration is kept
path_data <- paths_all_temp[grepl("flowType_results_response_Bcells", paths_all_temp)]
df_flowType_results_Bcells_final <- aws.s3::s3read_using(
    FUN = read.csv,
    colClasses="character",
    object = sprintf("s3://epichipc-main/%s", path_data)
)
# columns: X       samples_name phenotype_names cell_counts Pheno_codes vaccine_response randomization_group_vec sex_vec cells_ul_blood cell_frequencies


## convert flowType values into matrix
fg_dir <- "Flow_Cytometry/Rdata/FlowAnalysis/Alice_Yue_data/fg_bcell.Rdata"
if (!aws.s3::object_exists(object=sprintf("s3://epichipc-main/%s",path_data))) {
    sns <- unique(df_flowType_results_Bcells_final$samples_name)
    propANDul <- do.call(rbind, purrr::map(sns, function(sn) {
        inds <- df_flowType_results_Bcells_final$samples_name==sn
        as.numeric(
            append(df_flowType_results_Bcells_final$cell_frequencies[inds],
                   df_flowType_results_Bcells_final$cells_ul_blood[inds]))
    }))
    mprop <- propANDul[,seq_len(ncol(propANDul)/2)]
    if (max(mprop==100)) mprop <- mprop/100
    mul <- propANDul[,(1+ncol(propANDul)/2):ncol(propANDul)]

    # sample meta
    meta <- df_flowType_results_Bcells_final[!duplicated(df_flowType_results_Bcells_final$samples_name),]
    colnames(meta)[2] = "id"
    rownames(mul) <- rownames(mp) <- sns

    # phenotype meta
    phenocodes <- unique(df_flowType_results_Bcells_final$Pheno_codes)
    markers <- c("CD66", "CD38", "CD138", "CD20", "IgD", "IgM", "CD27", "CD10", "SSCACD19", "SSCACD34", "CD38CD10") # last 3 is labelled "Not____"
    phenotypes <- purrr::map(phenocodes, flowType::decodePhenotype, markers, 2) # unique(df_flowType_results_Bcells_final$phenotype_names)
    colnames(mprop) <- colnames(mul) <- phenotypes

    # create specenr features and add in cells/ml blood feature
    fg <- flowGraph(mprop, meta=meta, no_cores=no_cores, calculate_summary=FALSE)
    # fg <- fg_rm_feature(fg, type="node", feature="count") # cant
    fg <- fg_add_feature(fg, type="node", feature="cells_ul_blood",
                         m=mul, overwrite=FALSE)
    fg <- fg_add_feature(
        fg, type="node", feature="expect_cells_ul_blood", overwrite=FALSE,
        m=fg@feat$node$cells_ul_blood[,1]*fg@feat$node$expect_prop)
    SE_cb <- log(fg@feat$node$cells_ul_blood/fg@feat$node$expect_cells_ul_blood)
    e0 <- as.matrix(fg@feat$node$expect_cells_ul_blood)==0
    SE_cb[is.nan(as.matrix(SE_cb))] <- 0
    SE_cb[e0] <- log(as.matrix(fg@feat$node$cells_ul_blood)[e0])
    SE_cb[fg@feat$node$cells_ul_blood==0] <- 0
    fg <- fg_add_feature(fg, type="node", feature="SpecEnr_cells_ul_blood",
                         m=SE_cb, overwrite=TRUE)

    # save features
    aws.s3::s3save(fg, object=paste0("s3://epichipc-main/",fg_dir))
}
aws.s3::s3load(sprintf("s3://epichipc-main/%s",fg_dir)) # fg


# calculate summary statistics
fg <- fg_summary(fg, no_cores=no_cores, test_name="wilcox",
                 cls="vaccine_response", adjust_custom="none",
                 diminish=FALSE,
                 node_features=c("prop","SpecEnr", "cells_ul_blood",
                                 "SpecEnr_cells_ul_blood"),
                 overwrite=FALSE)
fg <- fg_summary(fg, no_cores=no_cores, test_name="wilcox_diminish",
                 cls="vaccine_response", adjust_custom="none",
                 node_features=c("prop","SpecEnr", "cells_ul_blood",
                                 "SpecEnr_cells_ul_blood"),
                 overwrite=FALSE)

fg <- fg_summary(fg, no_cores=no_cores, test_name="wilcox_byLayer",
                 cls="vaccine_response",
                 diminish=FALSE,
                 node_features=c("prop","SpecEnr", "cells_ul_blood",
                                 "SpecEnr_cells_ul_blood"),
                 overwrite=FALSE)
fg <- fg_summary(fg, no_cores=no_cores,
                 cls="vaccine_response",
                 node_features=c("prop","SpecEnr", "cells_ul_blood",
                                 "SpecEnr_cells_ul_blood"),
                 overwrite=FALSE)

fg <- fg_summary(fg, no_cores=no_cores, test_name="wilcox_BH",
                 cls="vaccine_response", adjust_custom="BH",
                 diminish=FALSE,
                 node_features=c("prop","SpecEnr", "cells_ul_blood",
                                 "SpecEnr_cells_ul_blood"),
                 overwrite=FALSE)
fg <- fg_summary(fg, no_cores=no_cores, test_name="wilcox_BH_diminish",
                 cls="vaccine_response", adjust_custom="BH",
                 node_features=c("prop","SpecEnr", "cells_ul_blood",
                                 "SpecEnr_cells_ul_blood"),
                 overwrite=FALSE)


st <- function(x) tryCatch(stats::shapiro.test(x)$p.value, error=function(e) 0)
test_custom <- function(x,y) {
    if (st(x)>0.05 & st(y)>0.05)
        return(tryCatch(stats::t.test(x,y)$p.value, error=function(e) 1))
    return(tryCatch(stats::wilcox.test(x,y)$p.value, error=function(e) 1))
}
fg <- fg_summary(fg, no_cores=no_cores, test_name="wilt",
                 cls="vaccine_response", adjust_custom="none",
                 diminish=FALSE, test_custom=test_custom,
                 node_features=c("prop","SpecEnr", "cells_ul_blood",
                                 "SpecEnr_cells_ul_blood"),
                 overwrite=FALSE)
fg <- fg_summary(fg, no_cores=no_cores, test_name="wilt_diminish",
                 cls="vaccine_response", adjust_custom="none",
                 test_custom=test_custom,
                 node_features=c("prop","SpecEnr", "cells_ul_blood",
                                 "SpecEnr_cells_ul_blood"),
                 overwrite=FALSE)

fg <- fg_summary(fg, no_cores=no_cores, test_name="wilt_byLayer",
                 cls="vaccine_response", test_custom=test_custom,
                 diminish=FALSE,
                 node_features=c("prop","SpecEnr", "cells_ul_blood",
                                 "SpecEnr_cells_ul_blood"),
                 overwrite=FALSE)
fg <- fg_summary(fg, no_cores=no_cores, test_name="wilt_byLayer_diminish",
                 cls="vaccine_response", test_custom=test_custom,
                 node_features=c("prop","SpecEnr", "cells_ul_blood",
                                 "SpecEnr_cells_ul_blood"),
                 overwrite=FALSE)

fg <- fg_summary(fg, no_cores=no_cores, test_name="wilt_BH",
                 cls="vaccine_response", test_custom=test_custom,
                 diminish=FALSE, adjust_custom="BH",
                 node_features=c("prop","SpecEnr", "cells_ul_blood",
                                 "SpecEnr_cells_ul_blood"),
                 overwrite=FALSE)
fg <- fg_summary(fg, no_cores=no_cores, test_name="wilt_BH_diminish",
                 cls="vaccine_response", test_custom=test_custom,
                 adjust_custom="BH",
                 node_features=c("prop","SpecEnr", "cells_ul_blood",
                                 "SpecEnr_cells_ul_blood"),
                 overwrite=FALSE)

mt <- function(x,y) tryCatch(stats::mood.test(x,y)$p.value, error=function(e) 1)
fg <- fg_summary(fg, no_cores=no_cores, test_name="mood",
                 cls="vaccine_response", test_custom=mt,
                 adjust_custom="none", diminish=FALSE,
                 node_features=c("prop","SpecEnr", "cells_ul_blood",
                                 "SpecEnr_cells_ul_blood"),
                 overwrite=FALSE)

qt <- function(x,y)
    min( sum(x<.1 & x>-.1)/length(x), sum(y<.1 & y>-.1)/length(y))
fg <- fg_summary(fg, no_cores=no_cores, test_name="test0",
                 cls="vaccine_response", test_custom=qt,
                 adjust_custom="none", diminish=FALSE,
                 node_features=c("prop","SpecEnr", "cells_ul_blood",
                                 "SpecEnr_cells_ul_blood"),
                 overwrite=TRUE)



aws.s3::s3save(fg, FUN=save, object=sprintf("s3://epichipc-main/%s",fg_dir))

