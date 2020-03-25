#### import code flowTypeFilter and libraries ######
# library(devtools)
# source('./flowTypeFilter_code/R/flowTypeFilterEdited_copy.R')
# source('./flowTypeFilter_code/R/HelperFunc.R')
install.packages("/mnt/f/Brinkman group/current/Alice/flowtype_metric/flowTypeFilter_code.tar.gz",repos=NULL, type="source")
#load_all('/home/rstudio/Code_Bcells_Myeloid_cells/flowTypeFilter_code/package/flowTypeFilterC')
library(flowCore)
library(stringr)
library(data.table)
library(flowTypeFilterC)
library(parallel)
############################# functions to use flowType #########

import_fs<-function(path_fcs_files){
  paths<-list.files(path_fcs_files,recursive = T,full.names = T)
  names_samples<-list.files(path_fcs_files,recursive = T,full.names = F)
  frames<-mclapply(paths, read.FCS,mc.cores = 8)
  fs<-as(frames, "flowSet")
  sampleNames(fs)<-names_samples
  return(fs)
}


import_thresholds<-function(path_df_thresholds){
  paths<-list.files(path_df_thresholds,recursive = T,full.names = T)
  df<-read.csv(paths,header = T)
  colnames(df)[1]<-"samples_name"
  return(df)
}

import_filters<-function(path_filter_objects,pop){
  paths<-list.files(path_filter_objects,recursive = T,full.names = T)
  ind<-grep(pop,paths)
  list_filter<-readRDS(paths[ind])
  return(list_filter)
}

get_indices_marker<-function(df_f,marker_name,channel_name="None",mode=1){
  ind<-which(df_f$desc==marker_name)
  if((length(ind)==0)||(length(ind)>1)){
    ind<-which(df_f$name==channel_name)
  }
  if(mode==2){
    ind<-which(df_f$desc==marker_name)
    if(length(ind)==0){
      ind<-which(df_f$name==marker_name)
    }
  }
  return(as.integer(ind))
}
get_markers_filter_df<-function(df_f,df_filters){
  markers<-colnames(df_filters)
  ind_marker_1<-get_indices_marker(df_f,markers[1],mode = 2)
  ind_marker_2<-get_indices_marker(df_f,markers[2],mode = 2)
  marker_1<-as.character(df_f$desc[ind_marker_1])
  marker_2<-as.character(df_f$desc[ind_marker_2])
  if(is.na(marker_1)==T){
    marker_1<-as.character(df_f$name[ind_marker_1])
  }
  if(is.na(marker_2)==T){
    marker_2<-as.character(df_f$name[ind_marker_2])
  }
  marker_1<-str_remove(marker_1,"-")
  marker_2<-str_remove(marker_2,"-")
  filter_markers<-paste0(marker_2,marker_1)
  return(filter_markers)
}
get_df_from_results<-function(list_results){
  list_df<-list()
  samples_names<-names(list_results)
  for(i in 1:length(list_results)){
    name_current_sample<-samples_names[i]
    df_i<-as.data.frame(list_results[[i]][c(1,2,4)]) # for each sample, we save counts,pheno_names and pheno_codes
    samples_name<-rep(name_current_sample,nrow(df_i))
    df_i<-cbind(samples_name,df_i)
    list_df[[i]]=df_i
  }
  df_flowType_results<-as.data.frame(rbindlist(list_df))
  return(df_flowType_results)
}
execute_flowType_filter<-function(fs,path_df_thresholds,path_filter_objects,path_output,panel){
  if(panel=="Bcells"){
    # import df_trehsholds current plate
    df_thresholds<<-import_thresholds(path_df_thresholds=path_df_thresholds)
    # import filters objects current plate
    list_filters_bcells<-import_filters(path_filter_objects=path_filter_objects,"Bcells")
    list_filters_blast<-import_filters(path_filter_objects=path_filter_objects,"blast")
    list_filters_trans<-import_filters(path_filter_objects=path_filter_objects,"trans")
    
    #---------- get flowType results--------------
    list_results<<-mclapply(1:length(fs),function(i){
      f<-fs[[i]]
      name_current_sample<-identifier(f)
      print(name_current_sample)
      #---- get info current flowFrame
      df_f<-pData(parameters(f))
      MarkerNames<-as.vector(df_f$desc)
      inds<-which(is.na(MarkerNames)==T)
      MarkerNames[inds]<-df_f$name[inds]
      #----- get filters df
      ind<-which(names(list_filters_bcells)==name_current_sample)
      df_filters_bcells<-list_filters_bcells[[ind]]
      df_filters_blast<-list_filters_blast[[ind]]
      df_filters_trans<-list_filters_trans[[ind]]
      # print(df_filters_bcells)
      # print(df_filters_blast)
      # print(df_filters_trans)
      #---- get thresholds of interest
      ind<-which(df_thresholds$samples_name==name_current_sample)
      thr_non_gran_CD66<-as.numeric(strsplit(as.character(df_thresholds$Non.granulocytes[ind]),",")[[1]][1])
      thr_pc_CD38<-as.numeric(strsplit(as.character(df_thresholds$plasma.cells[ind]),",")[[1]][1])
      thr_pc_CD138<-as.numeric(strsplit(as.character(df_thresholds$plasma.cells[ind]),",")[[1]][2])
      thr_CD20<-as.numeric(strsplit(as.character(df_thresholds$CD19.CD20..1[ind]),",")[[1]][2])
      thr_IgD<-as.numeric(strsplit(as.character(df_thresholds$IGD.IGM..B.cells[ind]),",")[[1]][1])
      thr_IgM<-as.numeric(strsplit(as.character(df_thresholds$IGD.IGM..B.cells[ind]),",")[[1]][2])
      thr_CD27<-as.numeric(strsplit(as.character(df_thresholds$Plasmablasts[ind]),",")[[1]][2])
      thr_CD10<-as.numeric(strsplit(as.character(df_thresholds$CD10..1[ind]),",")[[1]][1])
      thr_CD27_2<-as.numeric(strsplit(as.character(df_thresholds$Naive.B.cells[ind]),",")[[1]][1])
      thr_IgD_2<-as.numeric(strsplit(as.character(df_thresholds$IGD.IGM..B.cells[ind]),",")[[1]][2])
      
      ind_CD66<-get_indices_marker(df_f,"CD66")
      ind_CD38<-get_indices_marker(df_f,"CD38")
      ind_CD138<-get_indices_marker(df_f,"CD138")
      ind_CD20<-get_indices_marker(df_f,"CD20")
      ind_IgD<-get_indices_marker(df_f,"IgD")
      ind_IgM<-get_indices_marker(df_f,"IgM")
      ind_CD27<-get_indices_marker(df_f,"CD27")
      ind_CD10<-get_indices_marker(df_f,"CD10")
      
      # print(thr_non_gran_CD66)
      # print(thr_pc_CD138)
      # print(thr_CD20)
      # print(thr_IgD)
      # print(thr_IgM)
      # print(thr_CD27)
      # print(thr_CD10)
      # print(thr_CD27_2)
      # print(thr_IgD_2)
      
      # print(ind_CD66)
      # print(ind_CD38)
      # print(ind_CD20)
      # print(ind_CD138)
      # print(ind_IgD)
      # print(ind_IgM)
      # print(ind_CD27)
      # print(ind_CD10)
      #---set arguments value for flowType filter
      list_thresholds<-list(thr_non_gran_CD66,thr_pc_CD38,thr_pc_CD138,thr_CD20,c(thr_IgD,thr_IgD_2),thr_IgM,c(thr_CD27,thr_CD27_2),thr_CD10,df_filters_bcells,df_filters_blast,df_filters_trans)
      PropMarkers<<-c(ind_CD66,ind_CD38,ind_CD138,ind_CD20,ind_IgD,ind_IgM,ind_CD27,ind_CD10) # only indices markers with simple thresholds
      filter_markers_bcells<-get_markers_filter_df(df_f,df_filters_bcells)
      filter_markers_blast<-get_markers_filter_df(df_f,df_filters_blast)
      filter_markers_trans<-get_markers_filter_df(df_f,df_filters_trans)
      PartitionsPerMarker_vec<-rep(2,11)
      # execute flowType filter
      print("execution of flowTypefilter")
      output_flowtype_filter<-flowTypeFilterC::flowTypeFilter(Frame=f,PropMarkers=PropMarkers,Methods='Filters',Thresholds=list_thresholds,MarkerNames=MarkerNames,verbose = F,cores = 10,MaxMarkersPerPop=6,PartitionsPerMarker = PartitionsPerMarker_vec) # the deepest pop has 6 markers
      print("decodification of the phenotypes codes")
      phenotype.names<<-unlist(lapply(output_flowtype_filter@PhenoCodes,function(x){return(flowTypeFilterC::decodePhenotype(x,c("CD66","CD38","CD138","CD20","IgD","IgM","CD27","CD10",filter_markers_bcells,filter_markers_blast,filter_markers_trans),output_flowtype_filter@PartitionsPerMarker,c(rep(FALSE,8),T,T,T)))}))
      phenotype.names[1]<-"root_live"
      cell_counts<-output_flowtype_filter@CellFreqs
      Pheno_codes<-output_flowtype_filter@PhenoCodes
      return(list(phenotype_names=phenotype.names,cell_counts=cell_counts,sample_name=name_current_sample,Pheno_codes=Pheno_codes))
    },mc.cores = 10)
  }else if(panel=="myeloid"){
    # import df_trehsholds current plate
    df_thresholds<<-import_thresholds(path_df_thresholds=path_df_thresholds)
    # import filters objects current plate
    list_filters_basophils<-import_filters(path_filter_objects=path_filter_objects,"basophils")
    list_filters_CD66negCD45pos<-import_filters(path_filter_objects=path_filter_objects,"CD66negCD45pos")
    list_filters_DRn_14n<-import_filters(path_filter_objects=path_filter_objects,"DRn_14n")
    list_filters_DRp_14n<-import_filters(path_filter_objects=path_filter_objects,"DRp_14n")
    list_filters_granulocytes<-import_filters(path_filter_objects=path_filter_objects,"granulocytes")
    list_filters_monocytes<-import_filters(path_filter_objects=path_filter_objects,"monocytes")
    #---------- get flowType results--------------
    list_results<<-mclapply(1:length(fs),function(i){
      f<-fs[[i]]
      name_current_sample<-identifier(f)
      print(name_current_sample)
      #---- get info current flowFrame
      df_f<-pData(parameters(f))
      MarkerNames<-as.vector(df_f$desc)
      inds<-which(is.na(MarkerNames)==T)
      MarkerNames[inds]<-df_f$name[inds]
      #----- get filters df
      ind<-which(names(list_filters_basophils)==name_current_sample)
      df_filters_basophils<-list_filters_basophils[[ind]]
      df_filters_CD66negCD45pos<-list_filters_CD66negCD45pos[[ind]]
      df_filters_DRn_14n<-list_filters_DRn_14n[[ind]]
      df_filters_DRp_14n<-list_filters_DRp_14n[[ind]]
      df_filters_granulocytes<-list_filters_granulocytes[[ind]]
      df_filters_monocytes<-list_filters_monocytes[[ind]]
      
      # print(df_filters_basophils)
      # print(df_filters_CD66negCD45pos)
      # print(df_filters_DRn_14n)
      # print(df_filters_DRp_14n)
      # print(df_filters_granulocytes)
      # print(df_filters_monocytes)
      
      # #---- get thresholds of interest
      ind<-which(df_thresholds$samples_name==name_current_sample)
      thr_CD64<-as.numeric(df_thresholds$CD64.[ind])
      thr_mature_neutrophils_CD11B<-as.numeric(strsplit(as.character(df_thresholds$CD11b.CD16..Mature.Neutrophils[ind]),",")[[1]][1])
      thr_mature_neutrophils_CD16<-as.numeric(strsplit(as.character(df_thresholds$CD11b.CD16..Mature.Neutrophils[ind]),",")[[1]][2])
      thr_CD3<-as.numeric(df_thresholds$CD3.[ind])
      thr_gd<-as.numeric(strsplit(as.character(df_thresholds$gd.T.cells[ind]),",")[[1]][2])
      thr_CD16_NKT<-as.numeric(strsplit(as.character(df_thresholds$CD56.CD16..NKT.cells[ind]),",")[[1]][1])
      thr_CD56NKT<-as.numeric(strsplit(as.character(df_thresholds$CD56.CD16..NKT.cells[ind]),",")[[1]][2])
      thr_CD16_NK<-as.numeric(strsplit(as.character(df_thresholds$CD56.CD16..cells[ind]),",")[[1]][1])
      thr_CD56_NK<-as.numeric(strsplit(as.character(df_thresholds$CD56.CD16..cells[ind]),",")[[1]][2])
      thr_CD123<-as.numeric(strsplit(as.character(df_thresholds$B.cells[ind]),",")[[1]][1])
      thr_CD11c<-as.numeric(strsplit(as.character(df_thresholds$B.cells[ind]),",")[[1]][2])
      thr_CD16_monocytes<-as.numeric(strsplit(as.character(df_thresholds$Classical.Monocytes[ind]),",")[[1]][2])
      
    
      ind_CD64<-get_indices_marker(df_f,"CD64",channel_name = "Alexa Fluor 700-A")
      ind_CD11B<-get_indices_marker(df_f,"CD11b",channel_name = "BV786-A")
      ind_CD16<-get_indices_marker(df_f,"CD16",channel_name = "FITC-A")
      if(length(ind_CD16)>1){
        ind_CD16<-ind_CD16[2]
      }
      ind_CD3<-get_indices_marker(df_f,"CD3",channel_name = "PE-CF594-A")
      ind_gd<-get_indices_marker(df_f,"gd",channel_name = "PE-A")
      if(length(ind_gd)==0){
        ind_gd<-get_indices_marker(df_f,"gd TCR")
      }
      ind_CD56<-get_indices_marker(df_f,"CD56",channel_name = "BV650-A")
      ind_CD123<-get_indices_marker(df_f,"CD123",channel_name = "PE-Cy7-A")
      ind_CD11c<-get_indices_marker(df_f,"CD11c",channel_name = "APC-A")

      # print(thr_CD64)
      # print(thr_mature_neutrophils_CD11B)
      # print(thr_mature_neutrophils_CD16)
      # print(thr_CD3)
      # print(thr_gd)
      # print(thr_CD16_NKT)
      # print(thr_CD56NKT)
      # print(thr_CD16_NK)
      # print(thr_CD56_NK)
      # print(thr_CD123)
      # print(thr_CD11c)
      # print(thr_CD16_monocytes)
      
      print(ind_CD64)
      print(ind_CD11B)
      print(ind_CD16)
      print(ind_CD3)
      print(ind_gd)
      print(ind_CD56)
      print(ind_CD123)
      print(ind_CD11c)
      #---set arguments value for flowType filter
      list_thresholds<-list(thr_CD64,thr_mature_neutrophils_CD11B,c(thr_mature_neutrophils_CD16,thr_CD16_NKT,thr_CD16_NK,thr_CD16_monocytes),thr_CD3,thr_gd,c(thr_CD56NKT,thr_CD56_NK),thr_CD123,thr_CD11c,
                            df_filters_basophils,df_filters_CD66negCD45pos,df_filters_DRn_14n,df_filters_DRp_14n,df_filters_granulocytes,df_filters_monocytes)
      PropMarkers<-c(ind_CD64,ind_CD11B,ind_CD16,ind_CD3,ind_gd,ind_CD56,ind_CD123,ind_CD11c) # only indices markers with simple thresholds
      PartitionsPerMarker_vec<-rep(2,14)

      # execute flowType filter
      print("execution of flowTypefilter")
      output_flowtype_filter<-flowTypeFilterC::flowTypeFilter(Frame=f,PropMarkers=PropMarkers,Methods='Filters',Thresholds=list_thresholds,MarkerNames=MarkerNames,verbose = T,cores = 10,MaxMarkersPerPop=5,PartitionsPerMarker = PartitionsPerMarker_vec) # the deepest pop has 5 markers
      # the output_flowtype_filter variable should contain everything. The counts,the phenocodes,the thresholds and filters.
      print("decodification of the phenotypes codes")
      phenotype.names<-unlist(lapply(output_flowtype_filter@PhenoCodes,function(x){return(flowTypeFilterC::decodePhenotype(x,c("CD64","CD11B","CD16","CD3","gd","CD56","CD123","CD11c","basophils","NonGran","HLDR-CD14-","HLDR+CD14-","granulocytes","monocytes"),output_flowtype_filter@PartitionsPerMarker,c(rep(FALSE,8),T,T,T,T,T,T)))}))
      phenotype.names[1]<-"root_live"
      cell_counts<-output_flowtype_filter@CellFreqs
      Pheno_codes<-output_flowtype_filter@PhenoCodes
      return(list(phenotype_names=phenotype.names,cell_counts=cell_counts,sample_name=name_current_sample,Pheno_codes=Pheno_codes))
    },mc.cores = 1)
    
  }
  vec_samples_name<<-sapply(1:length(list_results),function(i){
    sample_name_i<-list_results[[i]]$sample_name
    return(sample_name_i)
  })
  names(list_results)<-vec_samples_name
  #---- get final df ---------
  df_flowType_results<-get_df_from_results(list_results)
  df_flowType_results[, ] <- lapply(df_flowType_results[, ], as.character)
  write.csv(df_flowType_results,path_output,row.names = F)
  return(df_flowType_results)
  }

########### function to combine all the flowtype results ###########
combine_all_flowtype_dfs<-function(paths_all_flowType_results,type_data="great_LOQ",type_visit="V2",n_cores=12){
  list_all_df<-list()
  # These are the subjects IDs with the best quality regarding the antibody concentration (reported on the server)
  subject_IDs_great_LOQ<-"G014K|G016D|G018C|G026E|G028G|G033F|G053B|G063J|G067C|G068A|G073K|G078D|G089J|G093G|G094J|G102K|G113F|G118B|G124C|G134F|G144E|G150K|G152D|G154G|G159F|G174A|G175E|G181C|G183G|G200J|G208B|G216J|G233G|G241A|G260G|G264A|G290B|G306H|G308A|G311F|G312J|G346D|G366B|G388G|G404K|G411J|G425G|G433D|G453H|G458E|G468K|G483J|G492K|G499D|G501C|G502A|G511G|G521E|G530H|G532K|G533C|G537J|G542C|G556D|G557E|G567B|G572E|G575G|G576C|G582J|G587H|G616F|G617K|G625J|G634G|G650H|G655F|G659D|G667A|G685K|G688E|G691G|G700D|G703E|G704C|G713J"
  info_df<-read.csv(file = "/home/rstudio/data/Samples_info_plates_all_batches_complete.csv")
  inds<-grep("V2|Visit 2",info_df$Visit.Num)
  visitID_v2<-as.character(info_df$Visit.ID[inds])
  string_visit_v2<-paste0(visitID_v2,collapse = "|")
  list_all_df<<-mclapply(1:length(paths_all_flowType_results),function(i){
    current_df<-read.csv(paths_all_flowType_results[i],colClasses = "character")
    if(type_data=="great_LOQ"){
      inds<-grep(subject_IDs_great_LOQ,current_df$samples_name)
      if(length(inds)!=0){
        current_df<-current_df[inds,]
      }else{
        return(i)
      }
    }
    if(type_visit=="V2"){
      inds<-grep(string_visit_v2,current_df$samples_name)
      if(length(inds)!=0){
        current_df<-current_df[inds,]
      }else{
        return(i)
      }
    }
    ind<-grep("Batch_2/Plate_5/",paths_all_flowType_results[i])
    if(length(ind)!=0){
      ind<-grep("FCT_G216J_4L42.fcs",current_df[,1])
      current_df<-current_df[ind,]
    }
    return(current_df)
  },mc.cores = n_cores )
 vec_classes<-sapply(list_all_df, class)
 inds<-which(vec_classes!="data.frame")
 list_all_df<-list_all_df[-inds]
 df<-as.data.frame(rbindlist(list_all_df))
 return(df)
}
#########################################################################################################################
################# execute functions #####################################################################################
#########################################################################################################################
fs_Bcells<-import_fs("/home/rstudio/data/flowType_data/Bcell_panel/Batch_4/Plate_3/fcs_files")
fs_myeloid<-import_fs("/home/rstudio/data/flowType_data/Myeloid_panel/Batch_1/fcs_files")
length(fs_Bcells)
length(fs_myeloid)
# I change fs_Bcells or fs_myeloid everytime I execute the function on a different set of fcs file (fcs file of each plate)
# the fcs files are stored on the bioinformatic server
########################## Bcell panel ######################################################
#############################################################################################

# Batch 1
start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_1/batch_1/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_1/batch_1/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_1/batch_1_replaced/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_1/batch_1_replaced/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 2 plate 1
start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_2/Plate_1/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_2/Plate_1/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 2 plate 2
start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_2/Plate_2/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_2/Plate_2/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")

end_time <- Sys.time()
running_time<-end_time - start_time
# batch 2 plate 3
start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_2/Plate_3/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_2/Plate_3/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time

# batch 2 plate 4
start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_2/Plate_4/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_2/Plate_4/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time


# batch 2 plate 5
start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_2/Plate_5/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_2/Plate_5/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time

# batch 2 plate 6
start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_2/Plate_6/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_2/Plate_6/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time

# batch 3 plate 1
start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_3/Plate_1/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_3/Plate_1/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time

# batch 3 plate 2

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_3/Plate_2/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_3/Plate_2/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time

# batch 3 plate 3

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_3/Plate_3/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_3/Plate_3/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 3 plate 4

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_3/Plate_4/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_3/Plate_4/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time

# batch 3 plate 5

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_3/Plate_5/plate_5/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_3/Plate_5/plate_5/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_3/Plate_5/plate_6/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_3/Plate_5/plate_6/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time

# batch 4 plate 1

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_4/Plate_1/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_4/Plate_1/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 4 plate 2

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_4/Plate_2/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_4/Plate_2/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 4 plate 3

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_4/Plate_3/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_4/Plate_3/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 4 plate 4
start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_4/Plate_4/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_4/Plate_4/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 4 plate 5

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_Bcells,path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_4/Plate_5/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Bcell_panel/Batch_4/Plate_5/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "Bcells")
end_time <- Sys.time()
running_time<-end_time - start_time

######################### Myeloid panel ##############################################################################
######################################################################################################################

# batch 1
start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_1/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_1/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time

# batch 2 plate 1
start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_1/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_1/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time

# batch 2 plate 2

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_2/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_2/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 2 plate 3

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_3/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_3/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time

# batch 2 plate 4

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_4/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_4/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time


# batch 2 plate 6

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_6/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_6/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 3 plate 1

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_3/Plate_1/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_3/Plate_1/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 3 plate 2

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_3/Plate_2/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_3/Plate_2/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 3 plate 3

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_3/Plate_3/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_3/Plate_3/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 3 plate 4

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_3/Plate_4/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_3/Plate_4/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 3 plate 5

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_3/Plate_5/Plate_1/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_3/Plate_5/Plate_1/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_3/Plate_5/Plate_2/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_3/Plate_5/Plate_2/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 4 plate 1

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_4/Plate_1/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_4/Plate_1/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 4 plate 2

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_4/Plate_2/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_4/Plate_2/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 4 plate 3

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_4/Plate_3/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_4/Plate_3/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time
# batch 4 plate 4

start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_4/Plate_4/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_4/Plate_4/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time

# batch 4 plate 5
start_time <- Sys.time()
df_flowType_results<-execute_flowType_filter(fs_myeloid,path_df_thresholds ="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_4/Plate_5/df_thresholds" ,path_filter_objects = "/home/rstudio/data/flowType_data/Myeloid_panel/Batch_4/Plate_5/filter_objects",
                                             path_output ="/home/rstudio/results/flowType_results/df_flowType_results.csv", panel = "myeloid")
end_time <- Sys.time()
running_time<-end_time - start_time




########## combine all dfs #########
paths_all_flowType_results_Bcells<-list.files("/home/rstudio/data/FlowType_results_all_batches_Bcells/",recursive = T,full.names = T)
paths_all_flowType_results_myeloid<-list.files("/home/rstudio/data/FlowType_results_all_batches_myeloid/",recursive = T,full.names = T)
# NOTE: I combine the results of each plate in one df,but I performed a filtering of the samples considering the samples with an LOQ greater than 2.5
# The LOQ is the limit of quantification related to the quality of measurement of the antibody concentration that I need to create the groups in the server.
# So,in summary, not all samples I have been used,only 86  subjects out of 711.
final_flowtyperesults_df_bcells_great_LOQ_v2<-combine_all_flowtype_dfs(paths_all_flowType_results_Bcells,type_data = "great_LOQ",type_visit="V2",n_cores=17)
final_flowtyperesults_df_bcells_great_LOQ_all_visits<-combine_all_flowtype_dfs(paths_all_flowType_results_Bcells,type_data="great_LOQ",type_visit="All",n_cores=17)


final_flowtyperesults_df_myeloid_great_LOQ_v2<-combine_all_flowtype_dfs(paths_all_flowType_results_myeloid,type_data = "great_LOQ",type_visit = "V2",n_cores = 17) 
final_flowtyperesults_df_myeloid_great_LOQ_all_visits<-combine_all_flowtype_dfs(paths_all_flowType_results_myeloid,type_data = "great_LOQ",type_visit = "All",n_cores = 17) 

#------------- check df
unique(final_flowtyperesults_df_myeloid_great_LOQ$samples_name)
unique(final_flowtyperesults_df_myeloid_great_LOQ_all_visits$samples_name)


unique(final_flowtyperesults_df_bcells_great_LOQ$samples_name)



#------------ write the selected flow type df
write.csv(final_flowtyperesults_df_bcells_great_LOQ_v2,"/home/rstudio/results/flowType_results/final_flowtyperesults_df_bcells_LOQ_v2.csv",row.names = F)
write.csv(final_flowtyperesults_df_bcells_great_LOQ_all_visits,"/home/rstudio/results/flowType_results/final_flowtyperesults_df_bcells_great_LOQ_all_visits.csv",row.names = F)

write.csv(final_flowtyperesults_df_myeloid_great_LOQ_v2,"/home/rstudio/results/flowType_results/final_flowtyperesults_df_myeloid_great_LOQ_v2.csv",row.names = F)
write.csv(final_flowtyperesults_df_myeloid_great_LOQ_all_visits,"/home/rstudio/results/flowType_results/final_flowtyperesults_df_myeloid_great_LOQ_all_visits.csv",row.names = F)


# once we have the proportions of each pop for each sample(using flowType),we find the samples that are positive or negative to the vaccine response (using the info in the AWS).
# Take all the pop x of samples with positive response,the pop x of samples with negative,we test if the mean is different.
# and so on for the other populations. # in the end we will have a correlation score for each phenotype. We will use Rhychptomix for the final plot.
