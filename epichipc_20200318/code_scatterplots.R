########################## import libraries ###########################
library("flowCore")
library("flowDensity")
library("ggcyto")
library(flowCut)
library(flowWorkspace)
library('MASS')
library(sp) # to make SpatialPolygon object
library(rgeos) # to use the function gIntersect,gDifference ecc..,for polygon comparison
library('flowPeaks') # Note : you need the library libgsl-dev
library(stringr)
#library(cytoUtils)
library(CytoML) # to use GatingSetFlowjo
library(data.table) # to use rbindlist
library(reshape2) # to use facet_wrap
library(plyr) # to use revalue
library(readxl)
library(parallel) # to enable parallization
library(RColorBrewer)

#########################################################################################################
########################################### Definitions of  functions ###################################
#########################################################################################################

#------------------------ function import fcs files -----
# fcs file that contains cleaned live cells (flowType input)

import_fcs_file<-function(path_fcs_file){
  f<-read.FCS(path_fcs_file)
  return(f)
}


#-------  function to execute gating of flowtype phenotypes --------------
# f= flowFrame
gating_flowtype_pheno<-function(f,channel_1="None",channel_2="None",phenotype_pheno="None",
                                gate_channel_1="None",gate_channel_2="None",name_gated_pop="None",
                                position_ch_1="pos",position_ch_2="pos",filter_pop=F,filter_object="none"){
  if(filter_pop==F){
    # we take the gates from the df_thresholds file
    gate_channel_1 <- gate_channel_1
    gate_channel_2 <- gate_channel_2
    if(is.na(gate_channel_1)==T){
      position_1<-NA
    }else{
      if(position_ch_1=="pos"){
        position_1<-T
      }else{
        position_1<-F
      }
    }
    if(is.na(gate_channel_2)==T){
      position_2<-NA
    }else{
      if(position_ch_2=="pos"){
        position_2<-T
      }else{
        position_2<-F
      }
    }
    
    # we gate the population
    print(channel_1)
    print(channel_2)
    print(position_1)
    print(position_2)
    print(gate_channel_1)
    print(gate_channel_2)
    cell_pop <- flowDensity(f, channels = c(channel_1, channel_2), position = c(position_1,position_2), gates = c(gate_channel_1, gate_channel_2))
  }else if(filter_pop==T){
    cell_pop <- flowDensity(f,channels = c(channel_1, channel_2), position = c(F,F),filter=filter_object)
    
  }

  plotDens(f, channels = c(channel_1, channel_2),  main = phenotype_pheno, cex.lab = 1, cex.axis = 1, cex.main=2)
  lines(cell_pop@filter,lwd=2)
  text(mean(cell_pop@filter[,1]), mean(cell_pop@filter[,2]), labels = name_gated_pop, cex = 2)
  f_gated<-getflowFrame(cell_pop)
  return(f_gated)
}








#########################################################################################################
########################################### execute functions ###########################################
#########################################################################################################

########### myeloid  panel################

#-------- positive group ----------------

f<-import_fcs_file(path_fcs_file="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_4/Plate_5/fcs_files/FCT_G511G_4C72.fcs")
f@parameters@data$desc[c(15,19,12,7)]<-c("CD66","CD45","CD11b","CD11c")

# first gate
#  it is a filter:we load the filter object from the flow type folder
filter_list<-readRDS("/home/rstudio/data/flowType_data/Myeloid_panel/Batch_4/Plate_5/filter_objects/filter_CD66negCD45pos.rds")
filter_df<-filter_list[[56]]
channel_1_gate_value<-max(filter_df[,1])
channel_2_gate_value<-max(filter_df[,2])

f_non_gran<-gating_flowtype_pheno(f,channel_1="BV711-A",channel_2="V450-A",phenotype_pheno="All cells",
                                   gate_channel_1=channel_1_gate_value,gate_channel_2 = channel_2_gate_value,
                                   position_ch_1="neg",position_ch_2="neg",name_gated_pop="Non Gran",
                                   filter_pop = T,filter_object = filter_df)
# second gate

f_CD11c_pos<-gating_flowtype_pheno(f_non_gran,channel_1="BV786-A",channel_2="APC-A",phenotype_pheno="NonGran",
                                   gate_channel_1=NA,gate_channel_2 = 2.39,position_ch_1="none",position_ch_2="pos",
                                   name_gated_pop="CD11c+")

# third gate

f_CD11b_pos<-gating_flowtype_pheno(f_CD11c_pos,channel_1="BV786-A",channel_2="APC-A",phenotype_pheno="NonGranCD11c+",
                                   gate_channel_1=2.36,gate_channel_2 = NA,position_ch_1="pos",position_ch_2="none",
                                   name_gated_pop="CD11b+")
#-------- negative group ------------------
f<-import_fcs_file(path_fcs_file="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_4/Plate_1/fcs_files/FCT_G556D_4U37.fcs")

# first gate
#  it is a filter:we load the filter object from the flow type folder
filter_list<-readRDS("/home/rstudio/data/flowType_data/Myeloid_panel/Batch_4/Plate_1/filter_objects/filter_CD66negCD45pos.rds")
filter_df<-filter_list[[55]]
channel_1_gate_value<-max(filter_df[,1])
channel_2_gate_value<-max(filter_df[,2])

f_non_gran<-gating_flowtype_pheno(f,channel_1="BV711-A",channel_2="V450-A",phenotype_pheno="All cells",
                                  gate_channel_1=channel_1_gate_value,gate_channel_2 = channel_2_gate_value,
                                  position_ch_1="neg",position_ch_2="neg",name_gated_pop="Non Gran",
                                  filter_pop = T,filter_object = filter_df)
# second gate

f_CD11c_pos<-gating_flowtype_pheno(f_non_gran,channel_1="BV786-A",channel_2="APC-A",phenotype_pheno="NonGran",
                                   gate_channel_1=NA,gate_channel_2 = 2.00,position_ch_1="none",position_ch_2="pos",
                                   name_gated_pop="CD11c+")

# third gate

f_CD11b_pos<-gating_flowtype_pheno(f_CD11c_pos,channel_1="BV786-A",channel_2="APC-A",phenotype_pheno="NonGranCD11c+",
                                   gate_channel_1=2.75,gate_channel_2 = NA,position_ch_1="pos",position_ch_2="none",
                                   name_gated_pop="CD11b+")



########### Bcell panel ################
f<-import_fcs_file(path_fcs_file="/home/rstudio/data/flowType_data/Bcell_panel/Batch_4/Plate_4/fcs_files/FCT_G688E_5U39.fcs")

#---------- positive sample ---------

# first gate

f_IgMpos<-gating_flowtype_pheno(f,channel_1="FITC-A",channel_2="PE-Cy5-A",phenotype_pheno="All cells",
                                 gate_channel_1=2.03,gate_channel_2 = NA,position_ch_1="pos",position_ch_2="none",
                                 name_gated_pop="IgM+")

# second gate

f_CD38pos<-gating_flowtype_pheno(f_IgMpos,channel_1="PE-Cy5-A",channel_2="FITC-A",phenotype_pheno="IgM+",
                                 gate_channel_1=2.13,gate_channel_2 = NA,position_ch_1="pos",position_ch_2="none",
                                 name_gated_pop="CD38+")

# third gate 
f_138pos<-gating_flowtype_pheno(f_CD38pos,channel_1="V450-A",channel_2="FITC-A",phenotype_pheno="IgM+CD38+",
                                 gate_channel_1=1.78,gate_channel_2 = NA,position_ch_1="pos",position_ch_2="none",
                                 name_gated_pop="CD138+")

#---------- negative sample ----------
f<-import_fcs_file(path_fcs_file="/home/rstudio/data/flowType_data/Bcell_panel/Batch_3/Plate_5/plate_5/fcs_files/FCT_G483J_7F23.fcs")


# first gate

f_IgMpos<-gating_flowtype_pheno(f,channel_1="FITC-A",channel_2="PE-Cy5-A",phenotype_pheno="All cells",
                                gate_channel_1=2.16,gate_channel_2 = NA,position_ch_1="pos",position_ch_2="none",
                                name_gated_pop="IgM+")

# second gate

f_CD38pos<-gating_flowtype_pheno(f_IgMpos,channel_1="PE-Cy5-A",channel_2="FITC-A",phenotype_pheno="IgM+",
                                 gate_channel_1=2.12,gate_channel_2 = NA,position_ch_1="pos",position_ch_2="none",
                                 name_gated_pop="CD38+")

# third gate 
f_138pos<-gating_flowtype_pheno(f_CD38pos,channel_1="V450-A",channel_2="FITC-A",phenotype_pheno="IgM+CD38+",
                                gate_channel_1=0.62,gate_channel_2 = NA,position_ch_1="pos",position_ch_2="none",
                                name_gated_pop="CD138+")
