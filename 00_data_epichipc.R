## this set of data must be done on the server
## import libraries and connection to AWS system ###################
setwd("/home/joann/EPICHIPC/.Ruserdata/smontante")
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


#------------- Bcells
# rchyoptimyx_df_Bcells_full_data_freqs <- generate_rchyoptimyx_df(
#     df_flowType_results_Bcells_final,
#     type_sex = "All",type_group = "All",type_counts = "freqs")

path_rchy <- paths_all_temp[grepl("rchyoptimyx_df_Bcells_full_data_cells_ul_blood", paths_all_temp)]
df_rchy <- aws.s3::s3read_using(
    FUN = read.csv,
    colClasses="character",
    object = sprintf("s3://epichipc-main/%s", path_rchy)
) # unadjusted p-values; rchyoptimyx input; cells_ul_blood


# with standard p-value
rchyoptimyx_results(
    df_rchy,
    panel = "Bcells",
    n_path_counts = 2,
    trim_level = 0,
    n_starting_phenotypes = 2,
    select_specific_starting_phenotypes = "None"
)







############ get boxplot for the poster/paper ########################

### myeloid panel #####
get_boxplots_phenotypes(flowType_results_response_myeloid,"CD64+_gd+_CD123-_",manual_selection=F)


# analyzing p-values before adjustment
inds_1<-which(rchyoptimyx_df_myeloid_full_data_cells_ul_blood$p_value==1)

inds_005<-which(rchyoptimyx_df_myeloid_full_data_cells_ul_blood$p_value<0.05)

inds_01_02<-which((rchyoptimyx_df_myeloid_full_data_cells_ul_blood$p_value>0.05) & (rchyoptimyx_df_myeloid_full_data_cells_ul_blood$p_value<0.2))

hist(rchyoptimyx_df_myeloid_full_data_cells_ul_blood$p_value,nclass = 20)
hist(rchyoptimyx_df_myeloid_full_data_freqs$p_value,nclass=20)

# analyzing p-values after adjustment

rchyoptimyx_df_myeloid_full_data_cells_ul_blood_test_after_adjust<-rchyoptimyx_df_myeloid_full_data_cells_ul_blood
rchyoptimyx_df_myeloid_full_data_freqs_after_adjust<-rchyoptimyx_df_myeloid_full_data_freqs


rchyoptimyx_df_myeloid_full_data_cells_ul_blood_test_after_adjust$p_value<-(qvalue(rchyoptimyx_df_myeloid_full_data_cells_ul_blood_test_after_adjust$p_value,fdr.level = 0.25,pi0.method="bootstrap"))$qvalues
rchyoptimyx_df_myeloid_full_data_cells_ul_blood_test_after_adjust$p_value<-p.adjust(rchyoptimyx_df_myeloid_full_data_cells_ul_blood_test_after_adjust$p_value,method = "BH")

rchyoptimyx_df_myeloid_full_data_freqs_after_adjust$p_value<-(qvalue(rchyoptimyx_df_myeloid_full_data_freqs_after_adjust$p_value,fdr.level = 0.25,pi0.method="bootstrap"))$qvalues
rchyoptimyx_df_myeloid_full_data_freqs_after_adjust$p_value<-p.adjust(rchyoptimyx_df_myeloid_full_data_freqs_after_adjust$p_value,method = "BH")


hist(rchyoptimyx_df_myeloid_full_data_cells_ul_blood_test_after_adjust$p_value,nclass = 20)
hist(rchyoptimyx_df_myeloid_full_data_freqs_after_adjust$p_value,nclass = 20)


inds<-which(rchyoptimyx_df_myeloid_full_data_cells_ul_blood_test_after_adjust$p_value<0.6)
tail(rchyoptimyx_df_myeloid_full_data_cells_ul_blood_test_after_adjust[inds,],200)

length(negative_group_values) # 93
length(positive_group_values) # 35

### Bcell panel #####
get_boxplots_phenotypes(flowType_results_response_Bcells,"CD38+_CD138+_IgM+_",manual_selection=F)


# analyzing p-values before adjustment
inds_1<-which(rchyoptimyx_df_Bcells_full_data_cells_ul_blood$p_value==1)

inds_005<-which(rchyoptimyx_df_Bcells_full_data_cells_ul_blood$p_value<0.05)

inds_01_02<-which((rchyoptimyx_df_Bcells_full_data_cells_ul_blood$p_value>0.05) & (rchyoptimyx_df_Bcells_full_data_cells_ul_blood$p_value<0.2))

hist(rchyoptimyx_df_Bcells_full_data_freqs$p_value,nclass = 20)

hist(rchyoptimyx_df_Bcells_full_data_cells_ul_blood$p_value,nclass = 20)

# analyzing p-values after adjustment

rchyoptimyx_df_Bcells_full_data_cells_ul_blood_test_after_adjust<-rchyoptimyx_df_Bcells_full_data_cells_ul_blood
rchyoptimyx_df_Bcells_full_data_freqs_test_after_adjust<-rchyoptimyx_df_Bcells_full_data_freqs

rchyoptimyx_df_Bcells_full_data_cells_ul_blood_test_after_adjust$p_value<-(qvalue(rchyoptimyx_df_Bcells_full_data_cells_ul_blood_test_after_adjust$p_value,fdr.level = 0.25,pi0.method="bootstrap"))$qvalues
rchyoptimyx_df_Bcells_full_data_cells_ul_blood_test_after_adjust$p_value<-p.adjust(rchyoptimyx_df_Bcells_full_data_cells_ul_blood_test_after_adjust$p_value,method = "BH")

rchyoptimyx_df_Bcells_full_data_freqs_test_after_adjust$p_value<-(qvalue(rchyoptimyx_df_Bcells_full_data_freqs_test_after_adjust$p_value,fdr.level = 0.25,pi0.method="bootstrap"))$qvalues
rchyoptimyx_df_Bcells_full_data_freqs_test_after_adjust$p_value<-p.adjust(rchyoptimyx_df_Bcells_full_data_freqs_test_after_adjust$p_value,method = "BH")


hist(rchyoptimyx_df_Bcells_full_data_freqs_test_after_adjust$p_value,nclass = 20)
hist(rchyoptimyx_df_Bcells_full_data_cells_ul_blood_test_after_adjust$p_value,nclass = 20)

inds<-which(rchyoptimyx_df_Bcells_full_data_cells_ul_blood_test_after_adjust$p_value<0.05)
rchyoptimyx_df_Bcells_full_data_cells_ul_blood_test_after_adjust[inds,]

length(negative_group_values) # 121
length(positive_group_values) # 48




####################################################### Info about the data ##############################################


#------- number of samples

# Bcells panel
samples_to_esclude<-get_info(flowType_results_response_Bcells,"NonGran_",n_samples_to_remove = 20)
SubjectID_to_esclude<-extract_ID(samples_to_esclude,type_ID = "SubjectID")


inds<-which(flowType_results_response_Bcells$vaccine_response=="negative")
samples_names_negative<-flowType_results_response_Bcells$samples_name[inds]
uniquie_samples_negative<-unique(samples_names_negative)
length(uniquie_samples_negative) # 120
inds<-which(flowType_results_response_Bcells$vaccine_response=="positive")
samples_names_positive<-flowType_results_response_Bcells$samples_name[inds]
uniquie_samples_positive<-unique(samples_names_positive)
length(uniquie_samples_positive)# 47

# myeloid panel
samples_to_esclude<-get_info(flowType_results_response_myeloid,"NonGran_",n_samples_to_remove = 20)
SubjectID_to_esclude<-extract_ID(samples_to_esclude,type_ID = "SubjectID")

inds<-which(flowType_results_response_myeloid$vaccine_response=="negative")
samples_names_negative<-flowType_results_response_myeloid$samples_name[inds]
negative_cells_ul_blood_myeloid<-flowType_results_response_myeloid$cells_ul_blood[inds]

uniquie_samples_negative<-unique(samples_names_negative)
length(uniquie_samples_negative) # 122

inds<-which(flowType_results_response_myeloid$vaccine_response=="positive")
samples_names_positive<-flowType_results_response_myeloid$samples_name[inds]
positive_cells_ul_blood_myeloid<-flowType_results_response_myeloid$cells_ul_blood[inds]

uniquie_samples_positive<-unique(samples_names_positive)
length(uniquie_samples_positive) # 47


#-------- number of subjects

all_subjects_ID<-extract_ID(unique(flowType_results_response_Bcells$samples_name),type_ID="SubjectID") # 84 IDs (2 samples) + 5 IDs (1 sample)=167 samples
unique_subjectsIDs<-unique(all_subjects_ID)
length(unique_subjectsIDs) # 86 subjects ID


all_subjects_ID<-extract_ID(unique(flowType_results_response_myeloid$samples_name),type_ID="SubjectID") # 81 IDs (2 samples) + 5 IDs (1 sample)= 167 samples
unique_subjectsIDs<-unique(all_subjects_ID)
length(unique_subjectsIDs) # 86 subjects ID

#-------- flow cytometry data info
df_myeloid_freqs<-import_all_csv_files_same_type(list_path_all_myeloid_results_folder,type="freqs",panel = "myeloid",ML_like = "None")
df_bcells_freqs<-import_all_csv_files_same_type(list_path_all_Bcells_results_folder,type="freqs",panel = "Bcells",ML_like = "None")

samples_myeloid<-unique(df_myeloid_freqs$samples_name)
samples_bcells<-unique(df_bcells_freqs$samples_name)
samples_myeloid_df<-as.data.frame(samples_myeloid)
samples_bcells_df<-as.data.frame(samples_bcells)

path_output<-"Flow_Cytometry/Rdata/FlowAnalysis/samples_passed_QC/samples_myeloid_df.csv"
aws.s3::s3write_using(samples_myeloid_df,FUN = write.csv,object = sprintf("s3://epichipc-main/%s",path_output))

path_output<-"Flow_Cytometry/Rdata/FlowAnalysis/samples_passed_QC/samples_bcells_df.csv"
aws.s3::s3write_using(samples_bcells_df,FUN = write.csv,object = sprintf("s3://epichipc-main/%s",path_output))

length(samples_myeloid)
length(samples_bcells)
