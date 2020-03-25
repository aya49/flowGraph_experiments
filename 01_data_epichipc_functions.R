############################################## import libraries and connection to AWS system ###################
setwd("/home/joann/EPICHIPC/.Ruserdata/smontante")
source("/home/joann/EPICHIPC/.Ruserdata/smontante/Utils.R")
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
#---- Calls the lists of the Amazon S3 buckets ----
bucketlist() #*** buckets are a large folder, we are interested in the epichipc-main folder


############################################ definition of import functions ############################################

#------ function to import counts and freq files Bcells separately -----------------------


#*** function to import the csv files of counts or frequencies
import_csv_file<-function(list_path_all_results_folder,batch,plate,type){

    inds_counts_files<-grep("final_counts",list_path_all_results_folder)
    vec_all_counts_files_results<-list_path_all_results_folder[inds_counts_files]
    inds_freqs_files<-grep("final_freqs",list_path_all_results_folder)
    vec_all_freqs_files_results<-list_path_all_results_folder[inds_freqs_files]
    if(type=="counts"){
        ind<-grep(batch,vec_all_counts_files_results)
        if(plate!="None"){
            vec_all_counts_files_results<-vec_all_counts_files_results[ind]
            ind<-grep(plate,vec_all_counts_files_results)
        }
        df <- aws.s3::s3read_using(
            FUN = read.csv,
            object = sprintf("s3://epichipc-main/%s",vec_all_counts_files_results[ind])
        )
    }else{
        ind<-grep(batch,vec_all_freqs_files_results)
        if(plate!="None"){
            vec_all_freqs_files_results<-vec_all_freqs_files_results[ind]
            ind<-grep(plate,vec_all_freqs_files_results)
        }
        df <- aws.s3::s3read_using(
            FUN = read.csv,
            object = sprintf("s3://epichipc-main/%s",vec_all_freqs_files_results[ind])
        )
    }
    return(df)
}

#------ function to import counts and freq files Bcells in a unique dataframe -----------------------
import_all_csv_files_same_type<-function(list_path_all_results_folder,type,panel,ML_like="None"){

    inds_counts_files<-grep("final_counts",list_path_all_results_folder)
    vec_all_counts_files_results<-list_path_all_results_folder[inds_counts_files]
    inds_freqs_files<-grep("final_freqs",list_path_all_results_folder)
    vec_all_freqs_files_results<-list_path_all_results_folder[inds_freqs_files]
    if(type=="counts"){
        list_all_df<-list()
        c<-0
        for(f in vec_all_counts_files_results){
            c<-c+1
            stringsplitted<-strsplit(f,"/")[[1]]
            ind_batch<-grep("Batch|batch",stringsplitted)
            batch<-stringsplitted[ind_batch]
            ind_plate<-grep("Plate|plate",stringsplitted)
            if(length(ind_plate)!=0){
                plate<-stringsplitted[ind_plate]
                if(length(plate)>1){ # case of 3rd batch plate 5 (because there two subplate)
                    plate<-plate[1]
                }
            }else{
                plate<-"None"
            }
            df <- aws.s3::s3read_using(
                FUN = read.csv,
                object = sprintf("s3://epichipc-main/%s",f)
            )
            # in Batch 2 plate 5 Bcell panel,we select only 1 sample FCT_G216J_4L42.fcs
            if(panel=="Bcells"){
                ind_temp<-grep("Batch_2/Plate_5",f)
                if(length(ind_temp)>0){
                    ind<-grep("FCT_G216J_4L42.fcs",df$name)
                    df<-df[ind,]
                }
            }

            # we remove the X column
            ind_X<-grep("^X$",colnames(df))
            if(length(ind_X)>0){ # one column contains the row names
                df<-df[,-ind_X]
            }

            # add column of cell/ul blood
            inds<-which(df$Population=="Beads")
            inds_1<-which(df$Population=="Time") # first population of each sample
            if(panel=="myeloid"){
                n_populations<-31 # 32 populations in each sample included time,so +31, or 31 populations included Singlets, so +30
                if(length(inds_1)==0){
                    inds_1<-which(df$Population=="Singlets") # first population of each sample is Singlets (Time removed for this plate)
                    n_populations<-30
                }
            }else if(panel=="Bcells"){
                n_populations<-24 # 25 populations in each sample included time,so +24, or 24 populations included Singlets, so +23
                if(length(inds_1)==0){
                    inds_1<-which(df$Population=="Singlets") # first population of each sample is Singlets (Time removed for this plate)
                    n_populations<-23
                }
            }else{
                stop("error input: incorrect panel value")
            }
            final_column_cells_ulbloods<-c()
            s<-0
            for(indstime in inds_1){
                s<-s+1 # counter of the sample
                counts_all_pops_sample_x<-df$Count[indstime:(indstime+n_populations)]

                ind_bead_current_sample<-inds[s]
                counts_beads_sample_x<-df$Count[ind_bead_current_sample]
                temp_counts<-counts_all_pops_sample_x/counts_beads_sample_x
                C<-480000
                D<-200
                dilution<-112.5
                temp_counts<-temp_counts*C
                temp_counts<-temp_counts/D
                final_counts<-temp_counts/dilution
                final_column_cells_ulbloods<-append(final_column_cells_ulbloods,final_counts)
            }
            final_column_cells_ulbloods<-format(final_column_cells_ulbloods,scientific = F)
            final_column_cells_ulbloods<-as.numeric(final_column_cells_ulbloods)
            final_column_cells_ulbloods<-round(final_column_cells_ulbloods,2)

            df<-cbind(df,final_column_cells_ulbloods)

            ind<-grep("final_column_cells_ulbloods",colnames(df))
            colnames(df)[ind]<-"cells_ul_blood"
            # add plate and batch
            vec_batch<-rep(batch,nrow(df))
            df<-cbind(df,vec_batch)

            vec_plate<-rep(plate,nrow(df))
            df<-cbind(df,vec_plate)
            inds<-grep("vec_batch|vec_plate",colnames(df))
            colnames(df)[inds]<-c("Batch","Plate")

            list_all_df[[c]]<-df
        }
        df<-as.data.frame(rbindlist(list_all_df))

        # we can also represent the counts in a Machine learning like structure (row=samples,cols=features).
        if(ML_like!="None"){
            df$name<-as.character(df$name)
            unique_samples<-unique(df$name)
            list_df_reshaped<-list()
            for (i in 1:length(unique_samples)){
                current_sample_name<-unique_samples[i]
                inds_current_sample<-grep(current_sample_name,df$name)
                # print(inds_current_sample)
                # print(length(inds_current_sample))
                df_current_sample<-df[inds_current_sample,]
                check_n_pops<-25
                if(panel=="myeloid"){
                    check_n_pops<-32
                }
                if(nrow(df_current_sample)>check_n_pops){
                    df_current_sample$name[(check_n_pops+1):nrow(df_current_sample)]<-sprintf("%s_dup",current_sample_name)
                }
                if(ML_like=="counts"){
                    df_current_sample<-df_current_sample[,c("name","Population","Count")]
                    df_current_sample_reshaped<-reshape(df_current_sample,idvar = "name",timevar = "Population",direction = "wide")
                }else if(ML_like=="ul_blood"){
                    df_current_sample<-df_current_sample[,c("name","Population","cells_ul_blood")]
                    df_current_sample_reshaped<-reshape(df_current_sample,idvar = "name",timevar = "Population",direction = "wide")
                }
                ind<-grep("dup",df_current_sample_reshaped$name)
                df_current_sample_reshaped$name[ind]<-current_sample_name
                list_df_reshaped[[i]]<-df_current_sample_reshaped
            }
            # list_df_reshaped<-mclapply(1:length(unique_samples),function(i){
            #   current_sample_name<-unique_samples[i]
            #   inds_current_sample<-grep(current_sample_name,df$name)
            #   df_current_sample<-df[inds_current_sample,]
            #   if(ML_like=="counts"){
            #     df_current_sample<-df_current_sample[,c("name","Population","Count")]
            #     df_current_sample_reshaped<-reshape(df_current_sample,idvar = "name",timevar = "Population",direction = "wide")
            #   }else if(ML_like=="ul_blood"){
            #     df_current_sample<-df_current_sample[,c("name","Population","cells_ul_blood")]
            #     df_current_sample_reshaped<-reshape(df_current_sample,idvar = "name",timevar = "Population",direction = "wide")
            #   }
            #   return(df_current_sample_reshaped)
            # },mc.cores = 40)
            df<-as.data.frame(rbindlist(list_df_reshaped,fill = T))
        }
    }else{
        list_all_df<-list()
        c<-0
        for(f in vec_all_freqs_files_results){
            c<-c+1
            stringsplitted<-strsplit(f,"/")[[1]]
            ind_batch<-grep("Batch|batch",stringsplitted)
            batch<-stringsplitted[ind_batch]
            ind_plate<-grep("Plate|plate",stringsplitted)
            if(length(ind_plate)!=0){
                plate<-stringsplitted[ind_plate]
                if(length(plate)>1){ # case of 3rd batch plate 5 (because there two subplate)
                    plate<-plate[1]
                }
            }else{
                plate<-"None"
            }
            df <- aws.s3::s3read_using(
                FUN = read.csv,
                object = sprintf("s3://epichipc-main/%s",f)
            )
            # in Batch 2 plate 5,we select only 1 sample FCT_G216J_4L42.fcs
            if(panel=="Bcells"){
                ind_temp<-grep("Batch_2/Plate_5",f)
                if(length(ind_temp)>0){
                    ind<-grep("FCT_G216J_4L42.fcs",df$samples_name)
                    df<-df[ind,]
                }
            }
            # we remove the X column
            ind_X<-grep("^X$",colnames(df))
            if(length(ind_X)>0){
                df<-df[,-ind_X]
            }
            # round the frequencies
            df[,2:ncol(df)]<-round(df[,2:ncol(df)],2)
            # add plate and batch info
            vec_batch<-rep(batch,nrow(df))
            df<-cbind(df,vec_batch)
            vec_plate<-rep(plate,nrow(df))
            df<-cbind(df,vec_plate)
            inds<-grep("vec_batch|vec_plate",colnames(df))
            colnames(df)[inds]<-c("Batch","Plate")
            list_all_df[[c]]<-df
        }
        df<-as.data.frame(rbindlist(list_all_df,fill = T))
        #---- rename colnames freqs df---
        if(panel=="Bcells"){
            columns_name_Bcells<-c("All cells","Atypical B cells","Blasts","CD10+","CD10-","CD19+ B cells",
                                   "CD19+CD20+","CD19+CD20-","Time","Granulocytes","IGD-IGM- B cells","IGD+IGM- B cells","IGD-IGM+ B cells",
                                   "IGD+IGM+ B cells","IGM","Immature transition B cells","Live cells","Naive B cells","Non granulocytes","Plasmablasts",
                                   "Singlets","Switched memory B cells-","Unswitched memory B cells","plasma cells","root")
            inds<-grep("samples_name|Batch|Plate",colnames(df))
            colnames(df)[-inds]<-columns_name_Bcells
        }else if(panel=="myeloid"){
            columns_name_myeloid<-c("All cells","B cells","Basophils","CD11b+CD16+ Mature Neutrophils","CD11b+CD16- Granulocytes","CD11b-CD16+ Immature Neutrophils 2",
                                    "CD11-CD16- Immature Neutrophils 1","CD3+ T cells","CD3-","CD45+CD66- (Non Granulocytes)","CD56 Hi NK","CD56- dim CD16- NK","CD56- dim CD16+ NK",
                                    "CD56+CD16+ NKT cells","CD56-CD16+ NK","CD56-CD16- cells","CD64+","CD64-","Classical Monocytes","Time","Granulocytes",
                                    "HLADR+CD14+ Monocytes","HLADR+ CD14-","HLADR-CD14-","Live cells","Non classical Monocytes","Singlets","gd T cells","gd- T cells","mDC",
                                    "pDC","root")
            inds<-grep("samples_name|Batch|Plate",colnames(df))
            colnames(df)[-inds]<-columns_name_myeloid

        }


    }
    df[, ] <- lapply(df[, ], as.character)
    return(df)
}

#----- function import vaccine response file -----------------------
#*** import data to create the class: negative/positive
import_vaccine_response<-function(list_path_antibody_files){
    ind<-grep("AB_FullExport.csv",list_path_antibody_files)
    path_ab_full_export<-list_path_antibody_files[ind][1]
    antibody_full_export <- aws.s3::s3read_using(
        FUN = read.csv,
        object = sprintf("s3://epichipc-main/%s",path_ab_full_export)
    )
    colnames(antibody_full_export)<-c("Subject_ID","Visit_ID","VisitNum","Randomization_group","Sex","Conclusion_Result_mIU_ml","ConcentrationActual","ConcentrationImputed")
    full_antibodytiters_df<-antibody_full_export
    full_antibodytiters_df[, ] <- lapply(full_antibodytiters_df[, ], as.character)
    #------- find subject IDs that have concentration <LOQ (at least one of V1 or V3)
    inds<-grep("V4|MAT",full_antibodytiters_df$VisitNum)
    full_antibodytiters_df_temp<-full_antibodytiters_df[-inds,]
    inds<-which(full_antibodytiters_df_temp$ConcentrationImputed>2.50)
    all_subjcetID_inf_LOQ<-full_antibodytiters_df_temp$Subject_ID[-inds]
    all_subjcetID_inf_LOQ<-unique(all_subjcetID_inf_LOQ) # V1 or V3 or both have concentration < 2.50
    string_all_subjcetID_inf_LOQ<-paste0(all_subjcetID_inf_LOQ,collapse = "|")
    #----------- add differences between concentration after vaccine and the baseline (V1)
    unique_sample_names<-unique(full_antibodytiters_df$Subject_ID)
    # add V3-V1, V4-V1
    list_modified_df<-list()
    i<-0
    for(sample_name in unique_sample_names){
        i<-i+1
        inds_current_sample<-grep(sample_name,full_antibodytiters_df$Subject_ID)
        antibodytiters_df_current_sample<-full_antibodytiters_df[inds_current_sample,]
        antibodytiters_df_current_sample<-antibodytiters_df_current_sample[,c("Subject_ID","Visit_ID","VisitNum","Randomization_group","Sex","ConcentrationActual","ConcentrationImputed")]
        ind<-which(antibodytiters_df_current_sample$ConcentrationImputed=="")
        if(length(ind)>0){
            antibodytiters_df_current_sample<-antibodytiters_df_current_sample[-ind,]
        }
        ind<-which(antibodytiters_df_current_sample$ConcentrationImputed=="QNS")
        if(length(ind)>0){
            antibodytiters_df_current_sample<-antibodytiters_df_current_sample[-ind,]
        }
        inds_v1<-grep("V1",antibodytiters_df_current_sample$VisitNum)
        inds_v3<-grep("V3",antibodytiters_df_current_sample$VisitNum)
        inds_v4<-grep("V4",antibodytiters_df_current_sample$VisitNum)
        V3_current_sample<-as.numeric(antibodytiters_df_current_sample$ConcentrationImputed[inds_v3])
        V4_current_sample<-as.numeric(antibodytiters_df_current_sample$ConcentrationImputed[inds_v4])
        # we esclude incomplete cases ( where one of the V is absent)
        if((length(inds_v3)!=0)&(length(inds_v4)!=0)&(length(inds_v1)!=0)){
            baseline_current_sample<-as.numeric(antibodytiters_df_current_sample$ConcentrationImputed[inds_v1])
            diff_V3_V1<-V3_current_sample-baseline_current_sample
            diff_V4_V1<-V4_current_sample-baseline_current_sample
            diff_V4_V3<-V4_current_sample-V3_current_sample

            # print(V4_current_sample)
            # print(V3_current_sample)
            # print(baseline_current_sample)

            antibodytiters_df_current_sample<-rbind(antibodytiters_df_current_sample,antibodytiters_df_current_sample[1,])
            row.names(antibodytiters_df_current_sample)<-NULL
            antibodytiters_df_current_sample[nrow(antibodytiters_df_current_sample),"Subject_ID"]<-sample_name
            antibodytiters_df_current_sample[nrow(antibodytiters_df_current_sample),"VisitNum"]<-"V3-V1"
            antibodytiters_df_current_sample[nrow(antibodytiters_df_current_sample),"Visit_ID"]<-"None"
            antibodytiters_df_current_sample[nrow(antibodytiters_df_current_sample),"ConcentrationImputed"]<-diff_V3_V1
            antibodytiters_df_current_sample<-rbind(antibodytiters_df_current_sample,antibodytiters_df_current_sample[1,])
            antibodytiters_df_current_sample[nrow(antibodytiters_df_current_sample),"Subject_ID"]<-sample_name
            antibodytiters_df_current_sample[nrow(antibodytiters_df_current_sample),"VisitNum"]<-"V4-V1"
            antibodytiters_df_current_sample[nrow(antibodytiters_df_current_sample),"Visit_ID"]<-"None"
            antibodytiters_df_current_sample[nrow(antibodytiters_df_current_sample),"ConcentrationImputed"]<-diff_V4_V1
            antibodytiters_df_current_sample<-rbind(antibodytiters_df_current_sample,antibodytiters_df_current_sample[1,])
            antibodytiters_df_current_sample[nrow(antibodytiters_df_current_sample),"Subject_ID"]<-sample_name
            antibodytiters_df_current_sample[nrow(antibodytiters_df_current_sample),"VisitNum"]<-"V4-V3"
            antibodytiters_df_current_sample[nrow(antibodytiters_df_current_sample),"Visit_ID"]<-"None"
            antibodytiters_df_current_sample[nrow(antibodytiters_df_current_sample),"ConcentrationImputed"]<-diff_V4_V3
        }
        list_modified_df[[i]]<-antibodytiters_df_current_sample
    }
    full_antibodytiters_df<-as.data.frame(rbindlist(list_modified_df))


    #----- generate df reporting subjectID with positive response
    inds<-grep("V3-V1",full_antibodytiters_df$VisitNum)
    full_antibodytiters_df_differences<-full_antibodytiters_df[inds,]
    group_antibody_response<-rep("a",nrow(full_antibodytiters_df_differences))
    inds_pos<-which(full_antibodytiters_df_differences$ConcentrationImputed>=0)
    group_antibody_response[inds_pos]<-"positive"
    group_antibody_response[-inds_pos]<-"negative"
    full_antibodytiters_df_differences<-cbind(full_antibodytiters_df_differences,group_antibody_response)
    subjectIDs<-full_antibodytiters_df_differences$Subject_ID
    vaccine_response<-as.character(full_antibodytiters_df_differences$group_antibody_response)
    randomization_groups<-as.character(full_antibodytiters_df_differences$Randomization_group)
    Sex<-as.character(full_antibodytiters_df_differences$Sex)
    df_vaccine_response<-cbind(subjectIDs,vaccine_response,randomization_groups,Sex)
    df_vaccine_response<-as.data.frame(df_vaccine_response)
    #----- Add LOQ info
    inds<-grep(string_all_subjcetID_inf_LOQ,df_vaccine_response$subjectIDs)
    LOQ_info<-rep("a",nrow(df_vaccine_response))
    LOQ_info[inds]<-"inf_LOQ"
    LOQ_info[-inds]<-"great_LOQ"
    df_vaccine_response<-cbind(df_vaccine_response,LOQ_info)
    return(list(full_antibodytiters_df=full_antibodytiters_df,df_vaccine_response=df_vaccine_response))
}

#------------------ function to remove duplicates and generate coupled IDs df -------------

generate_final_df<-function(df,type_data,select_only_coupled_ID,ML_like=T){
    print("################# fixing the df of gated pops ####################")
    if(type_data=="counts"){
        if(ML_like==F){
            unique_populations<-unique(df$Population)
            list_all_df_fixed<-list() # list all df without duplicates,with only coupled ID, for each pop
            p<-0
            for(pop in unique_populations){
                p<-p+1
                print(sprintf("-------------- current pop analyzed: %s ---------------",pop))
                pop_regex<-sprintf("^%s$",pop)
                # deal with special characters
                inds<-grep(pop_regex,df$Population)
                if(length(inds)==0){ # because probably there are special char in the name of the pop
                    inds<-which(df$Population==pop)
                }
                df_current_pop<-df[inds,]
                df_current_pop<-df_current_pop[,c("name","Population","Count","cells_ul_blood","Batch","Plate","vec_visits")]
                print(paste0("n. samples current pop original:",nrow(df_current_pop)))
                ##################  we remove the duplicates from the current pop###############
                unique_names<-unique(df_current_pop$name)
                list_all_df_names<-list()
                i<-0
                d<-0 # counter of the duplicates
                for(name in unique_names){
                    #print(name)
                    i<-i+1
                    inds<-grep(name,df_current_pop$name)
                    #print(inds)
                    if(length(inds)>1){
                        d<-d+1
                        print("threre is a duplication")
                        print(df_current_pop[inds,c("name","Batch","Plate")],row.names=F)
                        df_current_pop_selected_current_name<-df_current_pop[inds[2],] # we take only the second one (outside plate 5)

                    }else{
                        #print("no_duplication")
                        df_current_pop_selected_current_name<-df_current_pop[inds,]
                        # if((df_current_pop_selected_current_name$Batch=="Batch_2") && (df_current_pop_selected_current_name$Plate=="Plate_5")){
                        #   print("no_duplication")
                        #   print(df_current_pop_selected_current_name)
                        # }
                    }
                    list_all_df_names[[i]]<-df_current_pop_selected_current_name

                }

                df_current_pop<-as.data.frame(rbindlist(list_all_df_names))
                df_current_pop$cells_ul_blood<-as.numeric(df_current_pop$cells_ul_blood)
                df_current_pop$Count<-as.numeric(df_current_pop$Count)
                print(paste0("n. samples duplicated:",d))
                print(paste0("n. samples current pop no duplicates ID:",nrow(df_current_pop)))
                print("exclusion of internal duplicates ID (same batch and plate) and duplicates ID across different plate and batches")
                if(select_only_coupled_ID==T){
                    ############################### select only coupled IDs #########
                    # coupled ID = same subjectID,different visitID
                    inds_v1<-grep("^V1$",df_current_pop$vec_visits)
                    samples_names_V1_current_pop<-df_current_pop$name[inds_v1]
                    stringsplitted<-strsplit(samples_names_V1_current_pop,"_")
                    list_all_df<-list()
                    i<-0
                    for(s in stringsplitted){
                        i<-i+1
                        subjectID<-s[2]
                        inds<-grep(subjectID,df_current_pop$name)
                        if(length(inds)==2){
                            df_current_pop_coupled_ID<-df_current_pop[inds,]
                            list_all_df[[i]]<-df_current_pop_coupled_ID
                        }

                    }
                    df_current_pop<-as.data.frame(rbindlist(list_all_df))
                    print(paste0("n. samples current pop only coupled ID:",nrow(df_current_pop)))

                }
                list_all_df_fixed[[p]]<-df_current_pop
            }
            final_df<-as.data.frame(rbindlist(list_all_df_fixed))
        }else if(ML_like==T){
            print(paste0("n. samples current pop original:",nrow(df)))
            ##################  we remove the duplicates ###############
            unique_names<-unique(df$name)
            list_all_df_names<-list()
            i<-0
            d<-0 # counter of the duplicates
            for(name in unique_names){
                i<-i+1
                inds<-grep(name,df$name)
                if(length(inds)>1){
                    d<-d+1
                    #print("threre is a duplication")
                    print(df[inds,c("name")])
                    df_selected_current_name<-df[inds[2],] # we take only the second one (outside plate 5)
                }else{
                    #print("no_duplication")
                    df_selected_current_name<-df[inds,]
                }
                list_all_df_names[[i]]<-df_selected_current_name
            }
            df<-as.data.frame(rbindlist(list_all_df_names))
            print(paste0("n. samples duplicated:",d))
            print(paste0("n. samples current pop no duplicates ID:",nrow(df)))
            if(select_only_coupled_ID==T){
                ############################### select only coupled IDs #########
                inds_v1<-grep("^V1$",df$vec_visits)
                samples_names_V1<-df$samples_name[inds_v1]
                stringsplitted<-strsplit(samples_names_V1,"_")
                list_all_df<-list()
                i<-0
                for(s in stringsplitted){
                    i<-i+1
                    subjectID<-s[2]
                    inds<-grep(subjectID,df$samples_name)
                    if(length(inds)==2){
                        df_coupled_ID<-df[inds,]
                        list_all_df[[i]]<-df_coupled_ID
                    }

                }
                df<-as.data.frame(rbindlist(list_all_df))
                print(paste0("n. samples current pop only coupled ID:",nrow(df)))
            }
            final_df<-df
        }
    }else{
        print(paste0("n. samples current pop original:",nrow(df)))
        ##################  we remove the duplicates ###############
        unique_names<-unique(df$samples_name)
        list_all_df_names<-list()
        i<-0
        d<-0 # counter of the duplicates
        for(name in unique_names){
            i<-i+1
            inds<-grep(name,df$samples_name)
            if(length(inds)>1){
                d<-d+1
                print("there is a duplication")
                if(ML_like==T){
                    print(df[inds,c("samples_name")])
                }else{
                    print(df[inds,c("samples_name","Batch","Plate")])
                }
                df_selected_current_name<-df[inds[2],] # we take only the second one (outside plate 5)
            }else{
                #print("no_duplication")
                df_selected_current_name<-df[inds,]
            }
            list_all_df_names[[i]]<-df_selected_current_name
        }
        df<-as.data.frame(rbindlist(list_all_df_names))
        print(paste0("n. samples duplicated:",d))
        print(paste0("n. samples current pop no duplicates ID:",nrow(df)))
        if(select_only_coupled_ID==T){
            ############################### select only coupled IDs #########
            inds_v1<-grep("^V1$",df$vec_visits)
            samples_names_V1<-df$samples_name[inds_v1]
            stringsplitted<-strsplit(samples_names_V1,"_")
            list_all_df<-list()
            i<-0
            for(s in stringsplitted){
                i<-i+1
                subjectID<-s[2]
                inds<-grep(subjectID,df$samples_name)
                if(length(inds)==2){
                    df_coupled_ID<-df[inds,]
                    list_all_df[[i]]<-df_coupled_ID
                }

            }
            df<-as.data.frame(rbindlist(list_all_df))
            print(paste0("n. samples current pop only coupled ID:",nrow(df)))
        }
        final_df<-df
    }
    # ############# remove samples with no info about the visit Number #####
    # inds<-grep("None",final_df$vec_visits)
    # if(length(inds)>0){
    #   final_df<-final_df[-inds,]
    # }
    return(final_df)
}
#######################################################################################################################################
##################################### functions about the correlation analysis between the gated pop and antibody titer ###############
#######################################################################################################################################

# function to calculate the correlation scores with the antibody response for each population
counts_comparison_gated_pop_vs_antibody_titer<-function(df,full_antibodytiters_df,type_visit_antibody,info_pop,type_visit_gated_pop){
    if(type_visit_antibody=="V1"){
        stop("Error input type_visit_antibody: comparison againt V1 does not make sense,because V1 is the baseline,you are comparing the baseline agaist the baseline")
    }
    ###################  Selection antibody concentration #######################
    # We want to compare V1 against the selected type of antibody concentration (V3,V4,V3-V1,V4-V1)
    type_visit_antibody_regex<-sprintf("^%s$",type_visit_antibody)
    inds<-grep(type_visit_antibody_regex,full_antibodytiters_df$Visit)
    antibodytitier_df_selected<-full_antibodytiters_df[inds,]
    print(sprintf("n. samples antibody titer of %s:%s",type_visit_antibody,nrow(antibodytitier_df_selected)))
    inds_None<-grep("None",full_antibodytiters_df$Visit_ID)
    full_antibodytiters_df_no_differences<-full_antibodytiters_df[-inds_None,]
    print(paste0("n. samples antibody titer total before selection:",nrow(full_antibodytiters_df_no_differences)))
    # V1 is the baseline against which we compare one of the possible antibody concentration and V2 of the gated pop
    inds<-grep("^V1$",full_antibodytiters_df$Visit)
    antibodytitier_df_V1<-full_antibodytiters_df[inds,]
    print(paste0("n. samples antibody titer V1:",nrow(antibodytitier_df_V1)))
    ###################  Selection current gated population counts #######################
    # we select the V2 counts
    inds_V2<-grep("^V2$",df$vec_visits)
    df_V2<-df[inds_V2,]
    inds_v1<-grep("^V1$",df$vec_visits)
    df_v1<-df[inds_v1,]

    print(paste0("n. samples df v1 selected:",nrow(df_v1)))
    print(paste0("n. samples df v2 selected:",nrow(df_V2)))

    inds_pop<-grep(info_pop,df_V2$Population)
    if(length(inds_pop)==0){
        inds_pop<-which(df_V2$Population==info_pop)
    }
    df_current_pop_V2<-df_V2[inds_pop,]
    print(paste0("n. samples df selected pop V2:",nrow(df_current_pop_V2)))
    inds_pop<-grep(info_pop,df_v1$Population)
    if(length(inds_pop)==0){
        inds_pop<-which(df_v1$Population==info_pop)
    }
    df_current_pop_v1<-df_v1[inds_pop,]
    print(paste0("n. samples df selected pop v1:",nrow(df_current_pop_v1)))
    ###################### select same number of samples in both df #################

    SubjectIDs_antibody_visit_selected<-antibodytitier_df_selected$Subject_ID
    string<-paste0(SubjectIDs_antibody_visit_selected,collapse = "|")
    inds<-grep(string,df_current_pop_V2[,1])
    df_current_pop_V2<-df_current_pop_V2[inds,]
    print(paste0("n. samples df selected pop v2 fixed:",nrow(df_current_pop_V2)))
    SubjectIDs_antibody_V1<-antibodytitier_df_V1$Subject_ID
    string<-paste0(SubjectIDs_antibody_V1,collapse = "|")
    inds<-grep(string,df_current_pop_v1[,1])
    df_current_pop_v1<-df_current_pop_v1[inds,]
    print(paste0("n. samples df selected pop v1 fixed:",nrow(df_current_pop_v1)))

    SubjectIDs_gated_pop<-extract_ID(df_current_pop_v1$name,"SubjectID")
    string<-paste0(SubjectIDs_gated_pop,collapse = "|")
    inds<-grep(string,antibodytitier_df_V1$Subject_ID)
    antibodytitier_df_selected_v1<-antibodytitier_df_V1[inds,]
    print(paste0("n. samples df antibody titer selected fixed v1:",nrow(antibodytitier_df_selected_v1)))

    # the subject ID of  the antibody titer selected(V3,V4,V3-V1 etc...) must be the same of the V2 of the gated pop
    SubjectIDs_gated_pop<-extract_ID(df_current_pop_V2$name,"SubjectID")
    string<-paste0(SubjectIDs_gated_pop,collapse = "|")
    inds<-grep(string,antibodytitier_df_selected$Subject_ID)
    antibodytitier_df_selected<-antibodytitier_df_selected[inds,]
    print(paste0("n. samples df antibody titer selected fixed not v1:",nrow(antibodytitier_df_selected))) # V3 or V4
    ########################## calculate correlation score (linear relationship)############################
    cells_ul_blood_v1<-as.numeric(df_current_pop_v1$cells_ul_blood)
    true_counts_v1<-as.numeric(df_current_pop_v1$Count)
    Concentration_antibody_v1<-as.numeric(antibodytitier_df_selected_v1$ConcentrationActual)
    correlation_score_ul_blood_v1_gated_pop_v1_Ab<-cor(cells_ul_blood_v1,Concentration_antibody_v1)
    vec_v1<-rep("V1/V1",length(cells_ul_blood_v1))
    cells_ul_blood_V2<-as.numeric(df_current_pop_V2$cells_ul_blood)
    true_counts_V2<-as.numeric(df_current_pop_V2$Count)
    Concentration_antibody_selected<-as.numeric(antibodytitier_df_selected$ConcentrationActual)
    correlation_score_ul_blood_v2_gated_pop_vs_v_selected_Ab<-cor(cells_ul_blood_V2,Concentration_antibody_selected)
    vec_v2<-rep(sprintf("V2/%s",type_visit_antibody),length(cells_ul_blood_V2))
    df_scatter_v1<-as.data.frame(cbind(cells_ul_blood_v1,Concentration_antibody_v1,vec_v1))
    df_scatter_Ab_selected_V2<-as.data.frame(cbind(cells_ul_blood_V2,Concentration_antibody_selected,vec_v2))
    colnames(df_scatter_v1)<-c("cells_ul_blood","Concentration_antibody_selected","Visit_cells/visit_Ab")
    colnames(df_scatter_Ab_selected_V2)<-c("cells_ul_blood","Concentration_antibody_selected","Visit_cells/visit_Ab")
    df_scatter_final<-rbind(df_scatter_v1,df_scatter_Ab_selected_V2)
    df_scatter_final$cells_ul_blood<-as.numeric(as.character(df_scatter_final$cells_ul_blood))
    df_scatter_final$Concentration_antibody_selected<-as.numeric(as.character(df_scatter_final$Concentration_antibody_selected))
    df_scatter_final<<-df_scatter_final
    print(correlation_score_ul_blood_v1_gated_pop_v1_Ab)
    print(correlation_score_ul_blood_v2_gated_pop_vs_v_selected_Ab)


    ######################### calculate percentage protected vs not protected ################
    df_scatter_final_selected_v<-df_scatter_final[df_scatter_final$`Visit_cells/visit_Ab`!="V1/V1",]
    tot_samples<-nrow(df_scatter_final_selected_v)
    inds<-which(df_scatter_final_selected_v$Concentration_antibody_selected>10)# 10 is the CoP.
    n_samples_protected<-length(inds)
    percentage_protection<-(n_samples_protected/tot_samples)*100
    percentage_protection<-round(percentage_protection,2)
    print(sprintf("percentage protection %s:%s",type_visit_antibody,percentage_protection))
    ######################### calculate percentage out of range, undetectable ################
    inds<-which(df_scatter_final_selected_v$Concentration_antibody_selected<2.50)
    n_samples_out_of_range<-length(inds)
    percentage_out_of_range<-(n_samples_out_of_range/tot_samples)*100
    percentage_out_of_range <-round(percentage_out_of_range,2)
    inds<-which(df_scatter_final_selected_v$Concentration_antibody_selected==0)
    n_samples_undetectable<-length(inds)
    percentage_undetectable<-(n_samples_undetectable/tot_samples)*100
    percentage_undetectable<-round(percentage_undetectable,2)
    print(sprintf("percentage undetactable %s:%s",type_visit_antibody,percentage_undetectable))
    print(sprintf("percentage out of range %s:%s",type_visit_antibody,percentage_out_of_range))
    # ###################### make a scatterplot ########################
    p<-ggplot(df_scatter_final,aes(x=cells_ul_blood,y=Concentration_antibody_selected,color=`Visit_cells/visit_Ab`)) + geom_point() + ggtitle(sprintf("%s",info_pop))
    show(p)
}
# function to compare the counts of the populations at V1 vs V2 and see how the antibody response changes based on the counts of the population
counts_comparison_different_visits_gated_pop<-function(df,info_pop,type_data="counts",type_visual="single_pop",type_test="paired"){
    if(type_data=="counts"){
        if(type_visual=="single_pop"){
            print(sprintf("################# visualize pop: %s ####################",info_pop))
            inds_pop_selected<-grep(info_pop,df$Population)
            final_df_pop_selected<-df[inds_pop_selected,]
            #---------- boxplot between the v1 and V2 ------------------
            # with outliers
            # p <- ggplot(final_df_pop_selected, aes(x=vec_visits, y=cells_ul_blood)) + geom_boxplot() + ggtitle(sprintf("%s",info_pop))
            # without outliers
            sts <- boxplot.stats(final_df_pop_selected$cells_ul_blood)$stats # a vector of length 5, containing the extreme of the lower whisker, the lower ‘hinge’, the median, the upper ‘hinge’ and the extreme of the upper whisker.
            p <- ggplot(final_df_pop_selected, aes(x=vec_visits, y=cells_ul_blood)) + geom_boxplot(outlier.shape = NA) + ggtitle(sprintf("%s",info_pop))
            p <- p + coord_cartesian(ylim = c(sts[1],sts[5]*1.05))
            #---------- Mann withnee test between the v1 and not v1 ------------------
            inds_V1<-grep("^V1$",final_df_pop_selected$vec_visits)
            inds_V2<-grep("^V2$",final_df_pop_selected$vec_visits)
            df_current_pop_V1<-final_df_pop_selected[inds_V1,]
            df_current_pop_V2<-final_df_pop_selected[inds_V2,]
            score<-wilcox.test(df_current_pop_V1$cells_ul_blood,df_current_pop_V2$cells_ul_blood)
            show(p)
            print(score)
        }
        if(type_visual=="all_pops"){
            #-------- remove pops we don't want to plot ---------
            pops_to_remove<-c("All cells","Time","Singlets","Beads","Live cells")
            string<-paste0(pops_to_remove,collapse = "|")
            inds_pops_to_remove<-grep(string,df$Population)
            df<-df[-inds_pops_to_remove,]
            #----------- set correct class for each variable -------
            df$Population<-factor(df$Population)
            df$vec_visits<-factor(df$vec_visits)
            df$cells_ul_blood<-as.numeric(df$cells_ul_blood)
            #---------- Mann withnee test between the v1 and V2 for each pop ------------------
            final_scores_all_pops<-c()
            unique_populations<-unique(df$Population)
            for(pop in unique_populations){
                print(pop)
                pop_regex<-sprintf("^%s$",pop)
                inds_current_pop<-grep(pop_regex,df$Population)
                if(length(inds_current_pop)==0){
                    inds_current_pop<-which(df$Population==pop)
                }
                df_current_pop<-df[inds_current_pop,]
                df_current_pop$vec_visits<-as.character(df_current_pop$vec_visits)
                inds_V1<-grep("^V1$",df_current_pop$vec_visits)
                inds_V2<-grep("^V2$",df_current_pop$vec_visits)
                df_current_pop_V1<-df_current_pop[inds_V1,]
                df_current_pop_V2<-df_current_pop[inds_V2,]
                if(type_test!="paired"){
                    # not paired wilcox test
                    score<-wilcox.test(df_current_pop_V1$cells_ul_blood,df_current_pop_V2$cells_ul_blood,paired = F)
                    final_scores_all_pops<-append(final_scores_all_pops,sprintf("%s:%s",pop,format(round(score$p.value, 5), nsmall = 5)))

                }else{
                    # paired wilcox test
                    # Note: in this case the two df must have equal size
                    stringsplitted<-strsplit(df_current_pop_V2$name,"_")
                    sujectIDs_V2<-c()
                    for(s in stringsplitted){
                        sujectID<-s[2]
                        sujectIDs_V2<-append(sujectIDs_V2,sujectID)
                    }
                    string<-paste0(sujectIDs_V2,collapse = "|")
                    inds<-grep(string,df_current_pop_V1$name)
                    df_current_pop_V1<-df_current_pop_V1[inds,]
                    score<-wilcox.test(df_current_pop_V1$cells_ul_blood,df_current_pop_V2$cells_ul_blood,paired = T)
                    final_scores_all_pops<-append(final_scores_all_pops,sprintf("%s:%s",pop,format(round(score$p.value, 5), nsmall = 5)))

                }

            }
            print(final_scores_all_pops)
            #---------- make the paired boxplots of all pops ------------------
            # with outliers
            # p<-ggplot(df,aes(x=Population,y=cells_ul_blood,fill=vec_visits)) + geom_boxplot()  + theme(axis.text.x = element_text(face="bold",size=7,angle = 90)) + scale_fill_manual(values = c("blue","red"))
            # show(p)
            # without outliers
            upper_limit=100
            print(unique(df$vec_visits))
            p<-ggplot(df,aes(x=Population,y=cells_ul_blood,fill=vec_visits)) + geom_boxplot(outlier.shape = NA) +
                coord_cartesian(ylim = c(0,upper_limit)) + theme(axis.text.x = element_text(face="bold",size=7,angle = 90)) +
                scale_fill_manual(labels=c("pre-vaccine","post-vaccine"),values = c("blue","red")) + ggtitle("Gated populations: V1 vs V2")

            show(p)
            #annotate(geom="text", x=dflabels$vec_positions_x, y=dflabels$vec_positions_y, label=dflabels$vec_scores,color="red")
        }



    }else if(type_data=="freqs"){
        if(type_visual=="all_pops"){
            col_names<-colnames(final_df)
            inds<-grep("samples_name|vec_visits|Batch|Plate",col_names)
            col_names<-col_names[-inds]
            melted_df<-reshape2::melt(final_df,id.vars=c("samples_name","vec_visits"),measure.vars=col_names)
            melted_df$value<-as.numeric(melted_df$value)
            #---- remove pops -----
            pops_to_remove<-c("All cells","Time","Singlets","Beads","Live cells")
            string<-paste0(pops_to_remove,collapse = "|")
            inds_pops_to_remove<-grep(string,melted_df$variable)
            melted_df<-melted_df[-inds_pops_to_remove,]
            #---------- make the paired boxplots of all pops ------------------
            # with outliers
            # p<-ggplot(df,aes(x=Population,y=cells_ul_blood,fill=vec_visits)) + geom_boxplot()  + theme(axis.text.x = element_text(face="bold",size=7,angle = 90)) + scale_fill_manual(values = c("blue","red"))
            # show(p)
            # without outliers
            p<-ggplot(melted_df,aes(x=variable,y=value,fill=vec_visits)) + geom_boxplot(outlier.shape = NA) +
                theme(axis.text.x = element_text(face="bold",size=7,angle = 90)) + scale_fill_manual(labels=c("pre-vaccine","post-vaccine"),values = c("blue","red")) +
                ggtitle("Gated populations: V1 vs V2")
            show(p)
            #---------- Mann withnee test between the v1 and V2 for each pop ------------------
            final_scores_all_pops<-c()
            unique_populations<-unique(melted_df$variable)
            for(pop in unique_populations){
                pop_regex<-sprintf("^%s$",pop)
                inds_current_pop<-grep(pop_regex,melted_df$variable)
                if(length(inds_current_pop)==0){
                    inds_current_pop<-which(melted_df$variable==pop)
                }
                df_current_pop<-melted_df[inds_current_pop,]
                df_current_pop$vec_visits<-as.character(df_current_pop$vec_visits)
                inds_V1<-grep("^V1$",df_current_pop$vec_visits)
                inds_V2<-grep("^V2$",df_current_pop$vec_visits)
                df_current_pop_V1<-df_current_pop[inds_V1,]
                df_current_pop_V2<-df_current_pop[inds_V2,]
                if(type_test!="paired"){
                    # not paired wilcox test
                    score<-wilcox.test(df_current_pop_V1$value,df_current_pop_V2$value,paired = F)
                    final_scores_all_pops<-append(final_scores_all_pops,sprintf("%s:%s",pop,format(round(score$p.value, 5), nsmall = 5)))
                }
            }
            print(final_scores_all_pops)
        }
    }
}

#-------------------------- function to show the boxplot of the counts for all pops---------------------------------
boxplots_all_counts<-function(df,info_pop,threshold_checking){
    #-------- remove pops we don't want to plot ---------
    pops_to_remove<-c("CD64-","Time","Singlets","Beads")
    string<-paste0(pops_to_remove,collapse = "|")
    inds_pops_to_remove<-grep(string,df$Population)
    df<-df[-inds_pops_to_remove,]
    #----------- set correct class for each variable -------
    df$Population<-factor(df$Population)
    df$cells_ul_blood<-as.numeric(df$cells_ul_blood)
    df$Count<-as.numeric(df$Count)
    #----- make boxplot ------
    upper_limit=100
    # with outlier
    # myeloid panel
    # p<<-ggplot(df,aes(x=Population,y=Count)) + geom_boxplot()  + theme(axis.text.x = element_text(face="bold",size=7,angle = 90)) +
    # ggtitle("boxplot pops myeloid panel:1379 samples analyzed,average of 49 outliers (3.5%) ")
    # b cell panel
    p<<-ggplot(df,aes(x=Population,y=Count)) + geom_boxplot()  + theme(axis.text.x = element_text(face="bold",size=7,angle = 90)) +
        ggtitle("boxplot pops B cell panel:1475 samples analyzed, average of 61 outliers (4.1%) ")
    # without outliers
    # p<-ggplot(df,aes(x=Population,y=cells_ul_blood)) + geom_boxplot(outlier.shape = NA) +
    #   coord_cartesian(ylim = c(0,upper_limit)) + theme(axis.text.x = element_text(face="bold",size=7,angle = 90))
    show(p)

    #--- check outliers
    inds<-which(df$Population==info_pop)
    df<-df[inds,]
    inds<-which(df$Count>threshold_checking)
    name_samples_outliers<-df$name[inds]
    batches<-df$Batch[inds]
    Plates<-df$Plate[inds]
    counts<-df$Count[inds]
    print(name_samples_outliers)
    print(batches)
    print(Plates)
    print(counts)
    #---- get other info---
    # get average number of outliers
    out<-ggplot_build(p)$data[[1]]
    list_outliers<-as.data.frame(out)$outliers
    n_pop<-length(list_outliers)
    vec_n_outliers<-c()
    for (i in 1:n_pop){
        n_outlier_current_pop<-length(list_outliers[[i]])
        vec_n_outliers<-append(vec_n_outliers,n_outlier_current_pop)
    }
    mean_n_outliers<-mean(vec_n_outliers)
    print(mean_n_outliers)

}


#----------------------------- function to generate the flowType df with frequencies ------

generate_flowtype_results_frequencies<-function(df_flowType_results){
    unique_samples_name<-unique(df_flowType_results$samples_name)
    df_flowType_results$cell_counts<-as.numeric(df_flowType_results$cell_counts)
    print("starting frequencies generation")
    list_all_df<-mclapply(1:length(unique_samples_name),function(i){
        current_sample_id<-unique_samples_name[i]
        print(paste0("current_sample:",current_sample_id))
        inds<-grep(current_sample_id,df_flowType_results$samples_name)
        df_flowType_results_current_sample<-df_flowType_results[inds,]
        cell_frequencies<-df_flowType_results_current_sample$cell_counts/df_flowType_results_current_sample$cell_counts[1]
        cell_frequencies<-cell_frequencies*100
        df_flowType_results_current_sample<-cbind(df_flowType_results_current_sample,cell_frequencies)
        return(df_flowType_results_current_sample)
    },mc.cores = 4)
    print("Combining all df....")
    df_flowType_results<-as.data.frame(rbindlist(list_all_df))
    df_flowType_results$cell_frequencies<-round(df_flowType_results$cell_frequencies,3)
    return(df_flowType_results)
}

#--------------------------------- functions to generate the df ready to be used by RchyOptimyx ---------------
# type indicate the type of analysis conducted:
#1) full data or >LOD
# 2) all randomization groups togheter (all vaccine togheter) or specific randomization groups
# 3) Sex: Male or Female patients
combine_response_flowtype_results<-function(df_flowType_results,df_vaccine_response,type_data="great_LOQ"){
    df_vaccine_response[, ] <- lapply(df_vaccine_response[, ], as.character)
    df_flowType_results[, ] <- lapply(df_flowType_results[, ], as.character)
    All_subjectIDs<-as.character(df_vaccine_response$subjectIDs)
    if(type_data=="great_LOQ"){
        inds<-which(df_vaccine_response$LOQ_info=="great_LOQ")
        All_subjectIDs<-df_vaccine_response$subjectIDs[inds]
    }
    print(sprintf("all selected subject ID in vaccine response:%s",paste0(All_subjectIDs,collapse = ",")))
    print("-------- assign response to each subject ID in flowType results ----------")
    list_modifed_df<<-mclapply(1:length(All_subjectIDs),function(i){
        subject_id<-All_subjectIDs[i]
        print(paste0("current_subject_ID_under_analysis:",subject_id))
        inds_current_id<-grep(subject_id,df_flowType_results$samples_name)
        if(length(inds_current_id)>0){
            print("the subject id has been found in the flowType results")
            df_flowType_results_current_id<-df_flowType_results[inds_current_id,]
            inds_current_id_response<-grep(subject_id,df_vaccine_response$subjectIDs)

            response<-df_vaccine_response$vaccine_response[inds_current_id_response] # vaccine response for the current subject ID
            vaccine_response<-rep(response,nrow(df_flowType_results_current_id))

            randomization_group<-df_vaccine_response$randomization_groups[inds_current_id_response]
            randomization_group_vec<-rep(randomization_group,nrow(df_flowType_results_current_id))

            sex<-df_vaccine_response$Sex[inds_current_id_response]
            sex_vec<-rep(sex,nrow(df_flowType_results_current_id))

            df_flowType_results_current_id<-cbind(df_flowType_results_current_id,vaccine_response,randomization_group_vec,sex_vec)
            return(df_flowType_results_current_id)
        }
    },mc.cores = 4)
    flowType_results_response<-as.data.frame(rbindlist(list_modifed_df))
    return(flowType_results_response)
}

##*** THIS IS THE FUNCTION; takes in the df
generate_rchyoptimyx_df<-function(flowType_results_response,type_data="All",type_group="All",type_sex="All",type_counts="counts"){
    if(type_group!="All"){
        inds<-grep(type_group,flowType_results_response$randomization_group_vec)
        flowType_results_response<-flowType_results_response[inds,]
    }
    if(type_sex!="All"){
        inds<-grep(sex,flowType_results_response$sex_vec)
        flowType_results_response<-flowType_results_response[inds,]
    }
    flowType_results_response$cell_counts<-as.numeric(flowType_results_response$cell_counts)
    flowType_results_response$cell_frequencies<-as.numeric(flowType_results_response$cell_frequencies)
    flowType_results_response$cells_ul_blood<-as.numeric(flowType_results_response$cells_ul_blood)

    # I need to perform a wilcox.test for each population
    print("-------- wilcox.test testing for each phenotype ----------")
    unique_phenotypes<-unique(flowType_results_response$phenotype_names)
    list_results<<-mclapply(1:length(unique_phenotypes),function(i){
        pheno<-unique_phenotypes[i]
        print(paste0("phenotype under analysis:",pheno))
        print(paste0("n_phenotype:",i))
        inds_current_pheno<-which(flowType_results_response$phenotype_names==pheno)
        df_current_pheno<-flowType_results_response[inds_current_pheno,]
        df_current_pheno$vaccine_response<-factor(df_current_pheno$vaccine_response)
        if(type_counts=="counts"){
            test_current_pheno<-wilcox.test(df_current_pheno$cell_counts~df_current_pheno$vaccine_response)
        }else if(type_counts=="freqs"){
            # test for normality
            inds_neg<-which(df_current_pheno$vaccine_response=="negative")
            negative_freqs<-df_current_pheno$cell_frequencies[inds_neg]
            inds_pos<-which(df_current_pheno$vaccine_response=="positive")
            positive_freqs<-df_current_pheno$cell_frequencies[inds_pos]
            p.value_shapiro_negative<-tryCatch({(shapiro.test(negative_freqs))$p.value},error=function(e){0})
            p.value_shapiro_positive<-tryCatch({(shapiro.test(positive_freqs))$p.value},error=function(e){0})
            print(p.value_shapiro_negative)
            print(p.value_shapiro_positive)
            # execution of statistical test
            if((p.value_shapiro_negative>0.05) && (p.value_shapiro_positive>0.05)){
                print("performing t test")
                test_current_pheno<-t.test(df_current_pheno$cell_frequencies~df_current_pheno$vaccine_response)
            }else{
                print("performing wilcoxon test")
                test_current_pheno<-wilcox.test(df_current_pheno$cell_frequencies~df_current_pheno$vaccine_response)
            }
        }else if(type_counts=="cells_ul_blood"){
            # test for normality
            inds_neg<-which(df_current_pheno$vaccine_response=="negative")
            negative_cells_ul_blood<-df_current_pheno$cells_ul_blood[inds_neg]
            inds_pos<-which(df_current_pheno$vaccine_response=="positive")
            positive_cells_ul_blood<-df_current_pheno$cells_ul_blood[inds_pos]
            p.value_shapiro_negative<-tryCatch({(shapiro.test(negative_cells_ul_blood))$p.value},error=function(e){0})
            p.value_shapiro_positive<-tryCatch({(shapiro.test(positive_cells_ul_blood))$p.value},error=function(e){0})
            print(p.value_shapiro_negative)
            print(p.value_shapiro_positive)
            # execution of statistical test
            if((p.value_shapiro_negative>0.05) && (p.value_shapiro_positive>0.05)){
                print("performing t test")
                test_current_pheno<-t.test(df_current_pheno$cells_ul_blood~df_current_pheno$vaccine_response)
            }else{
                print("performing wilcoxon test")
                test_current_pheno<-wilcox.test(df_current_pheno$cells_ul_blood~df_current_pheno$vaccine_response)
            }
        }
        p.value<-test_current_pheno$p.value
        current_pheno_code<-unique(df_current_pheno$Pheno_codes)
        if(is.na(p.value)==T){
            p.value<-"None"
            print(p.value)
        }else{
            print(p.value)
        }
        if(is.null(pheno)==T){
            stop("Null found")
        }
        if(is.null(p.value)==T){
            stop("Null found")
        }
        if(is.null(current_pheno_code)==T){
            stop("Null found")
        }
        return(list(phenotypes=pheno,p_value=p.value,pheno_code=current_pheno_code))
    },mc.cores = 4)
    if(length(list_results)!=length(unique(flowType_results_response$phenotype_names))){
        stop("length(list_results) has less phenotypes than the expected ones")
    }
    final_df<-ldply(list_results, data.frame)
    if(nrow(final_df)!=length(unique(flowType_results_response$phenotype_names))){
        stop("nrow(final_df) has less phenotypes than the expected ones")
    }

    # fix p_value of phenotypes with frequencies == 0
    inds<-which(final_df$p_value=="None")
    if(length(inds)>0){
        #final_df$p_value[inds]<-runif(1, min=0.2, max=1)
        final_df<-final_df[-inds,]
    }
    # define correct class for each column
    final_df$p_value<-as.numeric(final_df$p_value)
    final_df$phenotypes<-as.character(final_df$phenotypes)
    final_df$pheno_code<-as.character(final_df$pheno_code)
    # check presence of root live (it must be present)
    ind<-grep("root_live",final_df$phenotypes)
    #add root live if absent
    if(length(ind)==0){
        pheno_code_root_live<-paste0(rep(0,nchar(as.character(final_df$pheno_code[1]))),collapse = "")
        root_live_row<-c("root_live",0.8,pheno_code_root_live)
        names(root_live_row)<-c("phenotypes","p_value","pheno_code")
        final_df<-rbind(root_live_row,final_df)
    }
    final_df$p_value<-as.numeric(final_df$p_value)
    final_df$phenotypes<-as.character(final_df$phenotypes)
    final_df$pheno_code<-as.character(final_df$pheno_code)
    return(final_df)
}
#-------------------------- function to generate the RchyOptimyx results --------------------------------
#*** flowtype result!
rchyoptimyx_results<-function(rchyoptimyx_df,panel,n_path_counts=2,trim_level=0,n_starting_phenotypes=1,
                              select_specific_starting_phenotypes="None",adjust_p=T){
    if(panel=="Bcells"){
        pheno.codes<-rchyoptimyx_df$pheno_code
        phenotypeScores<-rchyoptimyx_df$p_value
        pheno.names<-rchyoptimyx_df$phenotypes
        # adjust pvalue (Bonferroni)
        if(adjust_p==T){
            phenotypeScores<-p.adjust(phenotypeScores,method = "BH")
            #phenotypeScores<-(qvalue(phenotypeScores,fdr.level = 0.25,pi0.method="bootstrap"))$qvalues
        }
        # define the correct markers vector
        MarkerNames<-c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W","CD19","CD138","CD10","CD14","HLADR","CD20","IgD","CD27","CD38","CD66","CD34","viability dye","IgM","Time")
        filter_markers_bcells<-"SSCACD19"
        filter_markers_blast<-"SSCACD34"
        filter_markers_trans<-"CD38CD10"
        PropMarkers<-c(16,15,8,12,13,19,14,9)
        marker.names<-MarkerNames[PropMarkers]
        marker.names<-c(marker.names,filter_markers_bcells,filter_markers_blast,filter_markers_trans)
        print("All marker names:")
        print(marker.names)
        # we select the vector of starting phenotypes with pvalue < 0.05 (thus the significant phenotypes)
        inds<-which(phenotypeScores<0.05)
        selected_pvals<<-phenotypeScores[inds]
        selected_pheno<<-as.character(pheno.codes[inds])
        selected_pheno_names<<-as.character(pheno.names[inds])

        # print(selected_pheno_names)
        # print(selected_pvals)
        # We select the starting_phenotypes that we consider
        if(select_specific_starting_phenotypes=="None"){
            # we select the n last starting phenotypes
            selected_pheno<-selected_pheno[(length(selected_pheno)-(n_starting_phenotypes-1)):length(selected_pheno)] # we select the last x phenotypes
            selected_pheno_names<-selected_pheno_names[(length(selected_pheno_names)-(n_starting_phenotypes-1)):length(selected_pheno_names)]

        }else{
            # we select specific starting phenotypes
            inds<-grep(select_specific_starting_phenotypes,selected_pheno_names)
            selected_pheno_names<-selected_pheno_names[inds]
            print(selected_pheno_names)
            selected_pheno<-selected_pheno[inds]
        }
        print(sprintf("starting_phenotypes:%s",paste0(selected_pheno_names,collapse = ",")))
        print(paste0("n_starting_phenotypes_selected:",length(selected_pheno)))
        res<-RchyOptimyx(pheno.codes, -log10(phenotypeScores), startPhenotype=selected_pheno[1],factorial(n_path_counts),trimPaths = F,trim.level=trim_level )
        if(n_starting_phenotypes>1){
            for (i in 2:length(selected_pheno)){
                temp<-RchyOptimyx(pheno.codes, -log10(phenotypeScores), startPhenotype=selected_pheno[i],factorial(n_path_counts),trimPaths = F,trim.level=trim_level)
                res=merge(res,temp)
            }
            plot(res, phenotypeScores=-log10(phenotypeScores), phenotypeCodes=pheno.codes, marker.names=marker.names, ylab='-log10(Pvalue)')
        }else{
            plot(res, phenotypeScores=-log10(phenotypeScores), phenotypeCodes=pheno.codes, marker.names=marker.names, ylab='-log10(Pvalue)')
        }
    }else if(panel=="myeloid"){
        pheno.codes<-rchyoptimyx_df$pheno_code
        phenotypeScores<-rchyoptimyx_df$p_value
        pheno.names<-rchyoptimyx_df$phenotypes
        if(adjust_p==T){
            #phenotypeScores<-p.adjust(phenotypeScores,method = "BH")
            phenotypeScores<-(qvalue(phenotypeScores,fdr.level = 0.25,pi0.method="bootstrap"))$qvalues
            # in the setting of exploratory discovery where one is interested in finding candidate hypothesis to be further validated as a results of future research

        }
        # define the correct markers vector
        MarkerNames<-c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W","CD11c","CD64","CD14","HLADR","CD56","CD11b","CD3","CD123","CD66","gdTCR","Viability dye","CD16","CD45","Time")
        filter_markers_basophils<-"basophils"
        filter_markers_NonGran<-"NonGran"
        filter_markers_DRneg<-"HLDR-CD14-"
        filter_markers_DRpos<-"HLDR+CD14-"
        filter_markers_granulocytes<-"granulocytes"
        filter_markers_monocytes<-"monocytes"
        PropMarkers<-c(8,12,18,13,16,11,14,7)
        marker.names<-MarkerNames[PropMarkers]
        marker.names<-c(marker.names,filter_markers_basophils,filter_markers_NonGran,filter_markers_DRneg,filter_markers_DRpos,filter_markers_granulocytes,filter_markers_monocytes)
        print("All marker names:")
        print(marker.names)
        # we select the vector of starting phenotypes with pvalue < 0.05 (thus the significant phenotypes)
        inds<-which(phenotypeScores<0.30)
        selected_pvals<-phenotypeScores[inds]
        selected_pheno<-as.character(pheno.codes[inds])
        selected_pheno_names<-as.character(pheno.names[inds])
        print(selected_pheno_names)
        print(selected_pvals)
        # We select the starting_phenotypes that we consider
        if(select_specific_starting_phenotypes=="None"){
            # we select the n last starting phenotypes
            selected_pheno<-selected_pheno[(length(selected_pheno)-(n_starting_phenotypes-1)):length(selected_pheno)] # we select the last x phenotypes
            selected_pheno_names<-selected_pheno_names[(length(selected_pheno_names)-(n_starting_phenotypes-1)):length(selected_pheno_names)]
        }else{
            # we select specific starting phenotypes
            inds<-grep(select_specific_starting_phenotypes,selected_pheno_names)
            selected_pheno_names<-selected_pheno_names[inds]
            print(selected_pheno_names)
            selected_pheno<-selected_pheno[inds]
        }
        print(sprintf("starting_phenotypes:%s",paste0(selected_pheno_names,collapse = ",")))
        print(paste0("n_starting_phenotypes_selected:",length(selected_pheno)))
        res<-RchyOptimyx(pheno.codes, -log10(phenotypeScores), startPhenotype=selected_pheno[1],factorial(n_path_counts),trimPaths = F,trim.level=trim_level )
        if(n_starting_phenotypes>1){
            for (i in 2:length(selected_pheno)){
                temp<-RchyOptimyx(pheno.codes, -log10(phenotypeScores), startPhenotype=selected_pheno[i],factorial(n_path_counts),trimPaths = F,trim.level=trim_level)
                res=merge(res,temp)
            }
            plot(res, phenotypeScores=-log10(phenotypeScores), phenotypeCodes=pheno.codes, marker.names=marker.names, ylab='-log10(Pvalue)')
        }else{
            plot(res, phenotypeScores=-log10(phenotypeScores), phenotypeCodes=pheno.codes, marker.names=marker.names, ylab='-log10(Pvalue)')
        }
    }


}

#------------------ function to remove pops -----------

remove_pops<-function(pops_to_remove,df){
    inds<-grep(pops_to_remove,colnames(df))
    if(length(inds)>0){
        df<-df[,-inds]
    }else{
        print("pops not found")
    }
    return(df)
}

#######################################################################################################################################
##################################### functions about the  Analysis of the antibody titers data #######################################
#######################################################################################################################################

get_info_antibody_data<-function(df_vaccine_response){
    n_subjects_full_data<-nrow(df_vaccine_response)
    print(paste0("n_subjects_full_data:",n_subjects_full_data))
    inds<-grep("great_LOQ",df_vaccine_response$LOQ_info)
    pristine_data_vaccine_response<-df_vaccine_response[inds,]
    n_subjects_pristine_data<-nrow(pristine_data_vaccine_response)
    print(paste0("n_subjects_pristine_data:",n_subjects_pristine_data))

}
#######################################################################################################################################
##################################### functions about the data for the poster #######################################
#######################################################################################################################################

get_boxplots_phenotypes<-function(flowType_results_response,phenotype="None",manual_selection=T,type_data="freqs"){
    inds<-which(flowType_results_response$phenotype_names==phenotype)
    flowType_results_response_phenotype<<-flowType_results_response[inds,]
    nrow(flowType_results_response_phenotype)
    negative_group_inds<-which(flowType_results_response_phenotype$vaccine_response=="negative")
    if(type_data=="freqs"){
        negative_group_values<<-flowType_results_response_phenotype$cell_frequencies[negative_group_inds]
    }else{
        negative_group_values<<-flowType_results_response_phenotype$cells_ul_blood[negative_group_inds]
    }
    print(length(negative_group_values))
    median_negative_group_values<-median(negative_group_values)
    positive_group_inds<-which(flowType_results_response_phenotype$vaccine_response=="positive")
    if(type_data=="freqs"){
        positive_group_values<<-flowType_results_response_phenotype$cell_frequencies[positive_group_inds]
    }else{
        positive_group_values<<-flowType_results_response_phenotype$cells_ul_blood[positive_group_inds]

    }
    print(length(positive_group_values))
    median_positive_group_values<-median(positive_group_values)
    print(paste0("median_positive_group_values:",median_positive_group_values))
    print(paste0("median_negative_group_values:",median_negative_group_values))
    #---- find points near median ----
    sts_negative_group_values <- boxplot.stats(negative_group_values)$stats # a vector of length 5, containing the extreme of the lower whisker, the lower ‘hinge’, the median, the upper ‘hinge’ and the extreme of the upper whisker.
    lower_hinge<-sts_negative_group_values[2]
    inds<-which(negative_group_values<median_negative_group_values)
    negative_group_values_temp<-negative_group_values[inds]
    inds<-which(negative_group_values_temp>lower_hinge)
    negative_group_values_temp<-negative_group_values_temp[inds]
    point_near_median_negative_group<-max(negative_group_values_temp)
    if(manual_selection==T){
        # point_near_median_negative_group<-3.794 # myeloid panel
        point_near_median_negative_group<-0.433 # bcell panel
    }
    if(type_data=="freqs"){
        ind_point_near_median<-which(flowType_results_response_phenotype$cell_frequencies==point_near_median_negative_group)

    }else{
        ind_point_near_median<-which(flowType_results_response_phenotype$cells_ul_blood==point_near_median_negative_group)

    }
    flowType_results_response_phenotype_temp_negative_group<<-flowType_results_response_phenotype[ind_point_near_median,]


    sts_positive_group_values <- boxplot.stats(positive_group_values)$stats # a vector of length 5, containing the extreme of the lower whisker, the lower ‘hinge’, the median, the upper ‘hinge’ and the extreme of the upper whisker.
    lower_hinge<-sts_positive_group_values[2]
    inds<-which(positive_group_values<median_positive_group_values)
    positive_group_values_temp<-positive_group_values[inds]
    inds<-which(positive_group_values_temp>lower_hinge)
    positive_group_values_temp<-positive_group_values_temp[inds]
    point_near_median_positive_group<-max(positive_group_values_temp)
    if(manual_selection==T){
        #point_near_median_positive_group<-1.737 # myeloid panel
        point_near_median_positive_group<-0.109 # bcell panel

    }
    if(type_data=="freqs"){
        ind_point_near_median<-which(flowType_results_response_phenotype$cell_frequencies==point_near_median_positive_group)

    }else{
        ind_point_near_median<-which(flowType_results_response_phenotype$cells_ul_blood==point_near_median_positive_group)

    }
    flowType_results_response_phenotype_temp_positive_group<<-flowType_results_response_phenotype[ind_point_near_median,]
    #----- make plot----------------
    if(type_data=="freqs"){
        sts <<- boxplot.stats(flowType_results_response_phenotype$cell_frequencies)$stats # a vector of length 5, containing the extreme of the lower whisker, the lower ‘hinge’, the median, the upper ‘hinge’ and the extreme of the upper whisker.
        p<-ggplot(flowType_results_response_phenotype,aes(x=vaccine_response,y=cell_frequencies)) + geom_boxplot(outlier.shape = NA)
        p <- p + coord_cartesian(ylim = c(sts[1],sts[5]*1.10))
        p<- p + geom_point(aes(x=2,y=point_near_median_positive_group,color="red"))
        p<-p+labs(x="vaccine response group",y="proportion")
        p<-p+geom_point(aes(x=1,y=point_near_median_negative_group,color="red")) +
            theme(legend.position = "none",axis.text = element_text(size=15),axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16))
        p<- p+ scale_x_discrete(labels=c("V3-V1<0", "V3-V1>0"))

    }else{
        sts <<- boxplot.stats(flowType_results_response_phenotype$cells_ul_blood)$stats # a vector of length 5, containing the extreme of the lower whisker, the lower ‘hinge’, the median, the upper ‘hinge’ and the extreme of the upper whisker.
        p<-ggplot(flowType_results_response_phenotype,aes(x=vaccine_response,y=cells_ul_blood)) + geom_boxplot(outlier.shape = NA)
        p <- p + coord_cartesian(ylim = c(sts[1],sts[5]*1.10))
        p<- p + geom_point(aes(x=2,y=point_near_median_positive_group,color="red"))
        p<-p+labs(x="vaccine response group",y="cells/ul blood")
        p<-p+geom_point(aes(x=1,y=point_near_median_negative_group,color="red")) +
            theme(legend.position = "none",axis.text = element_text(size=15),axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16))
        p<- p+ scale_x_discrete(labels=c("V3-V1<0", "V3-V1>0"))

    }
    show(p)
    print(paste0("point_near_median_positive_group:",point_near_median_positive_group))
    print(paste0("point_near_median_negative_group:",point_near_median_negative_group))
    sample_name_positive<-flowType_results_response_phenotype_temp_positive_group$samples_name
    print(paste0("sample_positive_group:",sample_name_positive))
    sample_name_negative<-flowType_results_response_phenotype_temp_negative_group$samples_name
    print(paste0("sample_negative_group:",sample_name_negative))

}



#------------ function to generate the cells/ul blood -----
generate_flowtype_results_cells_ul_blood<-function(df,df_flowType_results){
    unique_samples_name<-unique(df_flowType_results$samples_name)
    list_all_df<<-mclapply(1:length(unique_samples_name),function(i){
        current_sample<-unique_samples_name[i]
        print(current_sample)
        inds<-grep(current_sample,df$name)
        if(length(inds)>0){
            df_current_sample<-df[inds,]
            ind<-which(df_current_sample$Population=="Beads")
            beads_value<-as.numeric(df_current_sample$Count[ind])
            inds<-grep(current_sample,df_flowType_results$samples_name)
            df_flowType_results_current_sample<-df_flowType_results[inds,]
            df_flowType_results_current_sample$cell_counts<-as.numeric(df_flowType_results_current_sample$cell_count)
            tryCatch({cells_ul_blood<-df_flowType_results_current_sample$cell_counts/beads_value},
                     warning=function(w){
                         print(df_flowType_results_current_sample)
                         print(beads_value)
                         print(df_current_sample)
                         print(i)
                         stop()

                     })
            cells_ul_blood<-cells_ul_blood*480000
            cells_ul_blood<-cells_ul_blood/200
            cells_ul_blood<-cells_ul_blood/112.5
            df_flowType_results_current_sample<-cbind(df_flowType_results_current_sample,cells_ul_blood)

            return(df_flowType_results_current_sample)
        }
    },mc.cores = 4)
    df_flowType_results<-as.data.frame(rbindlist(list_all_df))
    df_flowType_results$cells_ul_blood<-round(df_flowType_results$cells_ul_blood,3)

    return(df_flowType_results)
}


#------------ function to filter immunophenotypes -------------------
filter_immunophenotypes<-function(flowType_results_response){
    unique_phenotypes<-unique(flowType_results_response$phenotype_names)
    list_results_temp<<-mclapply(1:length(unique_phenotypes),function(i){
        pheno<-unique_phenotypes[i]
        print(paste0("phenotype under analysis:",pheno))
        print(paste0("n_phenotype:",i))
        inds_current_pheno<-which(flowType_results_response$phenotype_names==pheno)
        df_current_pheno<-flowType_results_response[inds_current_pheno,]
        if(length(unique(df_current_pheno$cells_ul_blood))==1){
            return("None")
        }else{
            return(df_current_pheno)
        }

    },mc.cores = 4)
    inds<-which(list_results_temp=="None")
    if(length(inds)>0){
        list_results_temp<-list_results_temp[-inds]
    }
    final_df<-ldply(list_results_temp, data.frame)
    return(final_df)
}
