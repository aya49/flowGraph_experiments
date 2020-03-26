#--------------------- load libraries ---------------------------------
library (knitr)
library(aws.s3)
library(aws.signature)
library(aws.ec2metadata)
library(data.table)
library(ggplot2)
library(dbplyr)
# ------------------------------ Access content s3 ---------------------
access_content_s3<-function(){
    ACCESS_KEY <- "AKIA2SPT4T3BF5QWESDH"
    SECRET_ACCESS_KEY <- "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"
    my_bucket <- "epichipc-main"
    Sys.setenv("AWS_ACCESS_KEY_ID" = ACCESS_KEY,
               "AWS_SECRET_ACCESS_KEY" = SECRET_ACCESS_KEY, "AWS_DEFAULT_REGION" = "us-east-1", "AWS_SESSION_TOKEN"="")
    aws.signature::locate_credentials()
}



#---------------------------------- import functions -------------------------------
import_files_antibodytiter<-function(bucket='epichipc-main'){
    # import data 
    antibodytiter_data_batch_1 <<- s3read_using(FUN = read.csv, object = "s3://epichipc-main/Antibody_Titers/Rdata/DMC_BCH_EPIC-HIPC_aHBsAb_Batch1.csv", header = TRUE, row.names = 1, na.strings = c("", "NA","QNS","#N/A"),check.names=F)
    antibodytiter_data_batch_2 <<- s3read_using(FUN = read.csv, object = "s3://epichipc-main/Antibody_Titers/Rdata/BCH_EPIC-HIPC_aHBsAb measurements_Batch 2_14JUN2019.csv", header = TRUE, row.names = NULL, na.strings = c("", "NA","QNS","#N/A"),sep=";",check.names=F)
}

import_files_results_auto_analysis<-function(bucket='epichipc-main'){
    # import data 
    final_counts_batch_1_plate <<- s3read_using(FUN = read.csv, object = "s3://epichipc-main/Flow_Cytometry/Automated_analysis_results/Results_Gambia_myeloid_final_new/1st_batch/Stats/final_counts_file.csv", header = TRUE, row.names = NULL, na.strings = c("", "NA","QNS","#N/A"),check.names=F)
    final_counts_batch_2_plate_1 <<- s3read_using(FUN = read.csv, object = "s3://epichipc-main/Flow_Cytometry/Automated_analysis_results/Results_Gambia_myeloid_final_new/2nd_batch/plate_1_results/Stats/final_counts_file.csv", header = TRUE, row.names = NULL, na.strings = c("", "NA","QNS","#N/A"),check.names=F)
    final_counts_batch_2_plate_2 <<- s3read_using(FUN = read.csv, object = "s3://epichipc-main/Flow_Cytometry/Automated_analysis_results/Results_Gambia_myeloid_final_new/2nd_batch/plate_2_results/Stats/final_counts_file.csv", header = TRUE, row.names = NULL, na.strings = c("", "NA","QNS","#N/A"),check.names=F)
    final_counts_batch_2_plate_3 <<- s3read_using(FUN = read.csv, object = "s3://epichipc-main/Flow_Cytometry/Automated_analysis_results/Results_Gambia_myeloid_final_new/2nd_batch/plate_3_results/Stats/final_counts_file.csv", header = TRUE, row.names = NULL, na.strings = c("", "NA","QNS","#N/A"),check.names=F)
    final_counts_batch_2_plate_4 <<- s3read_using(FUN = read.csv, object = "s3://epichipc-main/Flow_Cytometry/Automated_analysis_results/Results_Gambia_myeloid_final_new/2nd_batch/plate_4_results/Stats/final_counts_file.csv", header = TRUE, row.names = NULL, na.strings = c("", "NA","QNS","#N/A"),check.names=F)
    final_counts_batch_2_plate_6 <<- s3read_using(FUN = read.csv, object = "s3://epichipc-main/Flow_Cytometry/Automated_analysis_results/Results_Gambia_myeloid_final_new/2nd_batch/plate_6_results/Stats/final_counts_file.csv", header = TRUE, row.names = NULL, na.strings = c("", "NA","QNS","#N/A"),check.names=F)
    
    final_counts_batch_3_plate_1 <<- s3read_using(FUN = read.csv, object = "s3://epichipc-main/Flow_Cytometry/Automated_analysis_results/Results_Gambia_myeloid_final_new/3rd_batch/Plate_1_samples/Stats/final_counts_file.csv", header = TRUE, row.names = NULL, na.strings = c("", "NA","QNS","#N/A"),check.names=F)
    final_counts_batch_3_plate_2 <<- s3read_using(FUN = read.csv, object = "s3://epichipc-main/Flow_Cytometry/Automated_analysis_results/Results_Gambia_myeloid_final_new/3rd_batch/Plate_2_samples/Stats/final_counts_file.csv", header = TRUE, row.names = NULL, na.strings = c("", "NA","QNS","#N/A"),check.names=F)
    final_counts_batch_3_plate_3 <<- s3read_using(FUN = read.csv, object = "s3://epichipc-main/Flow_Cytometry/Automated_analysis_results/Results_Gambia_myeloid_final_new/3rd_batch/Plate_3_samples/Stats/final_counts_file.csv", header = TRUE, row.names = NULL, na.strings = c("", "NA","QNS","#N/A"),check.names=F)
    final_counts_batch_3_plate_4 <<- s3read_using(FUN = read.csv, object = "s3://epichipc-main/Flow_Cytometry/Automated_analysis_results/Results_Gambia_myeloid_final_new/3rd_batch/Plate_4_samples/Stats/final_counts_file.csv", header = TRUE, row.names = NULL, na.strings = c("", "NA","QNS","#N/A"),check.names=F)
    final_counts_batch_3_plate_5_p1 <<- s3read_using(FUN = read.csv, object = "s3://epichipc-main/Flow_Cytometry/Automated_analysis_results/Results_Gambia_myeloid_final_new/3rd_batch/Plate_5_samples/FCS files samples plate 1/Stats/final_counts_file.csv", header = TRUE, row.names = NULL, na.strings = c("", "NA","QNS","#N/A"),check.names=F)
    final_counts_batch_3_plate_5_p2 <<- s3read_using(FUN = read.csv, object = "s3://epichipc-main/Flow_Cytometry/Automated_analysis_results/Results_Gambia_myeloid_final_new/3rd_batch/Plate_5_samples/FCS files samples plate 2/Stats/final_counts_file.csv", header = TRUE, row.names = NULL, na.strings = c("", "NA","QNS","#N/A"),check.names=F)
    
}
#get_bucket(bucket = my_bucket, check_region=FALSE) # access to a specific folder of the s3 system 


#---------------------------- get results counts as cell concentration --------------
get_results_counts_as_cell_concentration<-function(counts_df){
    inds<-grep("Beads",counts_df$Population)
    beads_counts_all_samples<-counts_df[inds,]
    samples_name<-unique(counts_df$name)
    all_modified_df<-lapply(1:length(samples_name),function(i){
        sample_name_x<-samples_name[i]
        inds<-grep(sample_name_x,counts_df$name)
        counts_df_sample_x<-counts_df[inds,]
        inds<-grep("Margin|Time|Singlets|Beads",counts_df_sample_x$Population)
        counts_df_sample_x<-counts_df_sample_x[-inds,]
        inds<-grep(sample_name_x,beads_counts_all_samples$name)
        beads_count_sample_x<-beads_counts_all_samples[inds,]$Count
        counts_df_sample_x$Count<-counts_df_sample_x$Count/beads_count_sample_x
        if(any(is.infinite(counts_df_sample_x$Count))){
            stop("counts_df_sample_x$Count/beads_count_sample_x has given an inf value")
        }
        counts_df_sample_x$Count<-counts_df_sample_x$Count*480000
        counts_df_sample_x$Count<-counts_df_sample_x$Count/200
        counts_df_sample_x$Count<- counts_df_sample_x$Count/112.5
        counts_df_sample_x$Count<-round(counts_df_sample_x$Count,2)
        ind<-grep("ParentCount",colnames(counts_df_sample_x))
        counts_df_sample_x<-counts_df_sample_x[,-ind]
        return(counts_df_sample_x)
    })
    counts_df_ul_blood<-as.data.frame(rbindlist(all_modified_df))
    
}
# ------------- get single final count df file (all batches) ---------------------------------
get_single_final_counts_df<-function(){
    final_counts_batch_1_plate_ul_blood<-get_results_counts_as_cell_concentration(counts_df=final_counts_batch_1_plate)
    final_counts_batch_2_plate_1_ul_blood<-get_results_counts_as_cell_concentration(counts_df=final_counts_batch_2_plate_1)
    final_counts_batch_2_plate_2_ul_blood<-get_results_counts_as_cell_concentration(counts_df=final_counts_batch_2_plate_2)
    final_counts_batch_2_plate_3_ul_blood<-get_results_counts_as_cell_concentration(counts_df=final_counts_batch_2_plate_3)
    final_counts_batch_2_plate_4_ul_blood<-get_results_counts_as_cell_concentration(counts_df=final_counts_batch_2_plate_4) 
    final_counts_batch_2_plate_6_ul_blood<-get_results_counts_as_cell_concentration(counts_df=final_counts_batch_2_plate_6) 
    final_counts_batch_3_plate_1_ul_blood<-get_results_counts_as_cell_concentration(counts_df=final_counts_batch_3_plate_1) 
    final_counts_batch_3_plate_2_ul_blood<-get_results_counts_as_cell_concentration(counts_df=final_counts_batch_3_plate_2)
    final_counts_batch_3_plate_3_ul_blood<-get_results_counts_as_cell_concentration(counts_df=final_counts_batch_3_plate_3) 
    final_counts_batch_3_plate_4_ul_blood<-get_results_counts_as_cell_concentration(counts_df=final_counts_batch_3_plate_4) 
    final_counts_batch_3_plate_5_p1_ul_blood<-get_results_counts_as_cell_concentration(counts_df=final_counts_batch_3_plate_5_p1)
    final_counts_batch_3_plate_5_p2_ul_blood<-get_results_counts_as_cell_concentration(counts_df=final_counts_batch_3_plate_5_p2) 
    Sample_info<-rep("Batch_1-plate",nrow(final_counts_batch_1_plate_ul_blood))
    final_counts_batch_1_plate_ul_blood<-cbind(final_counts_batch_1_plate_ul_blood,Sample_info)
    Sample_info<-rep("Batch_2-plate_1",nrow(final_counts_batch_2_plate_1_ul_blood))
    final_counts_batch_2_plate_1_ul_blood<-cbind(final_counts_batch_2_plate_1_ul_blood,Sample_info)
    Sample_info<-rep("Batch_2-plate_2",nrow(final_counts_batch_2_plate_2_ul_blood))
    final_counts_batch_2_plate_2_ul_blood<-cbind(final_counts_batch_2_plate_2_ul_blood,Sample_info)
    Sample_info<-rep("Batch_2-plate_3",nrow(final_counts_batch_2_plate_3_ul_blood))
    final_counts_batch_2_plate_3_ul_blood<-cbind(final_counts_batch_2_plate_3_ul_blood,Sample_info)
    Sample_info<-rep("Batch_2-plate_4",nrow(final_counts_batch_2_plate_4_ul_blood))
    final_counts_batch_2_plate_4_ul_blood<-cbind(final_counts_batch_2_plate_4_ul_blood,Sample_info)
    Sample_info<-rep("Batch_2-plate_6",nrow(final_counts_batch_2_plate_6_ul_blood))
    final_counts_batch_2_plate_6_ul_blood<-cbind(final_counts_batch_2_plate_6_ul_blood,Sample_info)
    Sample_info<-rep("Batch_3-plate_1",nrow(final_counts_batch_3_plate_1_ul_blood))
    final_counts_batch_3_plate_1_ul_blood<-cbind(final_counts_batch_3_plate_1_ul_blood,Sample_info)
    Sample_info<-rep("Batch_3-plate_2",nrow(final_counts_batch_3_plate_2_ul_blood))
    final_counts_batch_3_plate_2_ul_blood<-cbind(final_counts_batch_3_plate_2_ul_blood,Sample_info)
    Sample_info<-rep("Batch_3-plate_3",nrow(final_counts_batch_3_plate_3_ul_blood))
    final_counts_batch_3_plate_3_ul_blood<-cbind(final_counts_batch_3_plate_3_ul_blood,Sample_info)
    Sample_info<-rep("Batch_3-plate_4",nrow(final_counts_batch_3_plate_4_ul_blood))
    final_counts_batch_3_plate_4_ul_blood<-cbind(final_counts_batch_3_plate_4_ul_blood,Sample_info)
    Sample_info<-rep("Batch_3-plate_5",nrow(final_counts_batch_3_plate_5_p1_ul_blood))
    final_counts_batch_3_plate_5_p1_ul_blood<-cbind(final_counts_batch_3_plate_5_p1_ul_blood,Sample_info)
    Sample_info<-rep("Batch_3-plate_5",nrow(final_counts_batch_3_plate_5_p2_ul_blood))
    final_counts_batch_3_plate_5_p2_ul_blood<-cbind(final_counts_batch_3_plate_5_p2_ul_blood,Sample_info)
    df_list<-list(final_counts_batch_1_plate_ul_blood,final_counts_batch_2_plate_1_ul_blood,final_counts_batch_2_plate_2_ul_blood,
                  final_counts_batch_2_plate_3_ul_blood,final_counts_batch_2_plate_4_ul_blood,final_counts_batch_2_plate_6_ul_blood,
                  final_counts_batch_3_plate_1_ul_blood,final_counts_batch_3_plate_2_ul_blood,final_counts_batch_3_plate_3_ul_blood,
                  final_counts_batch_3_plate_4_ul_blood,final_counts_batch_3_plate_5_p1_ul_blood,final_counts_batch_3_plate_5_p2_ul_blood)
    final_df_counts_all_batches_ul_blood<-as.data.frame(rbindlist(df_list))
    return(final_df_counts_all_batches_ul_blood)
}

#------------------------------------------------ fix duplicated files -------------------------
fix_duplicated_csv<-function(df,path_s3_output,panel,type,write_to_s3=F){
    df[, ] <- lapply(df[, ], as.character)
    if(type=="counts"){
        unique_names<-unique(df$name)
    }else{
        unique_names<-unique(df$samples_name)
    }
    list_fixed_df_all_samples<-list()
    i<-0
    if(type=="counts"){
        for(name in unique_names){
            i<-i+1
            inds<-grep(name,df$name)
            if(panel=="myeloid"){
                check_number<-32
            }else if(panel=="Bcells"){
                check_number<-25
            }
            df_current_sample<-df[inds,]
            print(nrow(df_current_sample))
            if(nrow(df_current_sample)!=check_number){
                print("Warning:duplicated samples found!")
                print(name)
                print(nrow(df_current_sample))
                df_current_sample_fixed<-df_current_sample[(check_number+1):nrow(df_current_sample),]
                print("duplicated samples fixed!")
                print(nrow(df_current_sample_fixed))
            }else{
                print("no duplication")
                df_current_sample_fixed<-df_current_sample
                
            }
            list_fixed_df_all_samples[[i]]<-df_current_sample_fixed
        }
        df_fixed<-as.data.frame(rbindlist(list_fixed_df_all_samples))
    }else if(type=="freqs"){
        for(name in unique_names){
            i<-i+1
            inds<-grep(name,df$samples_name)
            df_current_sample<-df[inds,]
            if(nrow(df_current_sample)!=1){
                print("Warning:duplicated samples found!")
                print(name)
                print(nrow(df_current_sample))
                df_current_sample_fixed<-df_current_sample[2,]
                print("duplicated samples fixed!")
                print(nrow(df_current_sample_fixed))
            }else{
                print("no duplication")
                df_current_sample_fixed<-df_current_sample
                
            }
            list_fixed_df_all_samples[[i]]<-df_current_sample_fixed
        }
        df_fixed<-as.data.frame(rbindlist(list_fixed_df_all_samples))
    }
    if(write_to_s3==T){
        aws.s3::s3write_using(df_fixed,
                              FUN = write.csv,
                              object = path_s3_output
        )
    }
    
}

fix_csv_file<-function(path_csv_file_to_fix,type_file,read_only,name_sample,write_to_s3,path_samples_to_add="None",path_new_file="None"){
    df <- aws.s3::s3read_using(
        FUN = read.csv, 
        object =  path_csv_file_to_fix
    )
    df[, ] <- lapply(df[, ], as.character)
    if(read_only==T){
        print(nrow(df))
        return(df)
    }else{
        if(type_file=="Description"){
            ind<-grep(name_sample,df$Samples)
            # df$Status[ind]<-"Flagged"
            # df$Description[ind]<-"poor quality"
            # df$Criteria[ind]<-"Second Criteria"
            print(nrow(df))
            df$Samples[ind]<-"FCT_G002B_4L67.fcs"
            print(nrow(df))
            print(df)
            if(write_to_s3==T){
                aws.s3::s3write_using(df,FUN = write.csv,object = path_csv_file_to_fix)
            }
        }
        if(type_file=="counts"){
            inds<-grep(name_sample,df$name)
            print(nrow(df))
            #df<-df[-inds,]
            df$name[inds]<-"FCT_G002B_4L67.fcs"
            print(nrow(df))
            print(df)
            if(write_to_s3==T){
                aws.s3::s3write_using(df,FUN = write.csv,object = path_csv_file_to_fix)
            }
        }
        if(type_file=="freqs"){
            inds<-grep(name_sample,df$samples_name)
            print(nrow(df))
            #df<-df[-inds,]
            df$samples_name[inds]<-"FCT_G002B_4L67.fcs"
            print(nrow(df))
            print(df)
            if(write_to_s3==T){
                aws.s3::s3write_using(df,FUN = write.csv,object = path_csv_file_to_fix)
            }
        }
        if(type_file=="add_samples"){
            df_new_samples <- aws.s3::s3read_using(
                FUN = read.csv, 
                object =  path_samples_to_add
            )
            print(nrow(df))
            df_fixed<-rbind(df,df_new_samples)
            print(nrow(df_fixed))
            if(write_to_s3==T){
                aws.s3::s3write_using(df_fixed,FUN = write.csv,object = path_csv_file_to_fix)
            }
        }
        if(type_file=="replace_whole_file"){
            new_df <- aws.s3::s3read_using(
                FUN = read.csv, 
                object =  path_new_file
            )
            print(nrow(new_df))
            if(write_to_s3==T){
                aws.s3::s3write_using(new_df,FUN = write.csv,object = path_csv_file_to_fix)
            }
        }
    }
}
#------------- check_number_samples----------------

check_number_samples<-function(list_paths_csv_to_check,info_file,real_unique=F){
    vec_n_samples<-c()
    for(file in list_paths_csv_to_check){
        df <- aws.s3::s3read_using(
            FUN = read.csv, 
            object =  sprintf("s3://epichipc-main/%s",file)
        )
        df[, ] <- lapply(df[, ], as.characdf_current_pop_Not_V1ter)
        if(real_unique==T){
            unique_samples_names<-unique(df$samples_name)
        }else{
            unique_samples_names<-df$samples_name
        }
        print(length(unique_samples_names))
        i<-0
        ind_info<-grep(info_file,file)
        print(length(ind_info))
        
        if(length(ind_info)>0){
            #print(unique_samples_names)
            print(df$samples_name)
            
        }
        for(name in unique_samples_names){
            i<-i+1
            inds<-grep(name,df$samples_name)
            if(length(inds)>1){
                #print("threre is a duplication")
                df_selected_current_name<-df[inds[1],] # we take only the first one
                if(length(ind_info)>0){
                    print(df_selected_current_name[inds,])
                }
            }else{
                #print("no_duplication")
                df_selected_current_name<-df[inds,]
            }
        }
        vec_n_samples<-append(vec_n_samples,length(unique_samples_names))
    }
    total_n_samples<-sum(vec_n_samples)
    print(total_n_samples)
}

#----------------------------------- get info df ---------------
# function to extract various info from the df
get_info<-function(df,type_info,batch,plate){
    if(batch!="all"){
        inds_batch<-grep(batch,df$Batch)
        df<-df[inds_batch,]
    }
    if(plate!="all"){
        inds_plate<-grep(plate,df$Plate)
        df<-df[inds_plate,]
    }
    print(paste0("Number samples selected df:",nrow(df)))
    if(type_info=="V1_vs_V2"){
        inds_V1<-grep("^V1$",df$vec_visits)
        inds_V2<-grep("^V2$",df$vec_visits)
        df_V1<-df[inds_V1,]
        df_V2<-df[inds_V2,]
        print("number of samples based on visits")
        print(paste0("V1:",nrow(df_V1)))
        print(paste0("Not_V1:",nrow(df_V2)))
        print("ID based on visits")
        print("---- V1 IDs:")
        print(df_V1[,1])
        print("---- Not_V1 IDs:")
        print(df_V2[,1])
    }
    
}

#----------------------------------- add visits info to df ---------------------------
# function to add visit columns to the gated pop df
add_visits_columns<-function(df){
    info_samples_IDs <- aws.s3::s3read_using(
        FUN = read.csv, 
        object =  "s3://epichipc-main/Flow_Cytometry/Rdata/FlowAnalysis/Samples_info_plates_all_batches_complete.csv" 
    )
    info_samples_IDs[, ] <- lapply(info_samples_IDs[, ], as.character)
    # we find the visitID of V1
    inds<-grep("Visit 1|V1",info_samples_IDs$Visit.Num)
    visitID_V1<-info_samples_IDs$Visit.ID[inds]
    # we find the visitID of V2
    inds<-grep("Visit 2|V2",info_samples_IDs$Visit.Num)
    visitID_V2<-info_samples_IDs$Visit.ID[inds]
    # add visit colum to gated pop df
    vec_visits<-rep("a",nrow(df))
    df<-cbind(df,vec_visits)
    df[, ] <- lapply(df[, ], as.character)
    string_v1<-paste0(visitID_V1,collapse = "|")
    string_v2<-paste0(visitID_V2,collapse = "|")
    inds<-grep(string_v1,df[,1])
    df$vec_visits[inds]<-"V1"
    inds<-grep(string_v2,df[,1])
    df$vec_visits[inds]<-"V2"
    inds<-grep("a",df$vec_visits)
    df$vec_visits[inds]<-"None"
    return(df)
}

#----------------------------------- add visits info to df  flowType---------------------------

#---------------------- extract subject ID or visit ID -----------

# function to extract ID from full samples name

extract_ID<-function(vec_samples_names,type_ID){
    stringsplitted<-strsplit(vec_samples_names,"_")
    vec_IDs<-c()
    for(s in stringsplitted){
        if(type_ID=="SubjectID"){
            subjectid<-s[2]
            vec_IDs<-append(vec_IDs,subjectid)
        }else if(type_ID=="visitID"){
            visitid<-s[3]
            vec_IDs<-append(vec_IDs,visitid)
        }
    }
    return(vec_IDs)
    
}

#-------------------- get info flowtype_results df ------------------

# function to get several info about the df
get_info<-function(flowtype_results_df,info_phenotype,find_samples_to_remove=T,n_samples_to_remove){
    if(info_phenotype!="None"){
        inds<-which(flowtype_results_df$phenotype_names==info_phenotype)
        flowtype_results_df_selected_pheno<-flowtype_results_df[inds,]
    }
    if(find_samples_to_remove==T){
        sorted_cells_ul_blood<-sort(flowtype_results_df_selected_pheno$cells_ul_blood)
        print(sorted_cells_ul_blood)
        string<-paste0(sorted_cells_ul_blood[1:n_samples_to_remove],collapse = "|")
        inds<-grep(string,flowtype_results_df_selected_pheno$cells_ul_blood)
        samples_to_esclude<-flowtype_results_df_selected_pheno$samples_name[inds]
        return(samples_to_esclude)
    }
    
}

