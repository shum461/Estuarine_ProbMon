



library(tidyverse)
library(lubridate)
library(fuzzyjoin)
library(stringr)
library(readxl)
library(janitor)
library(rio)

source("Salinity_regime_funs.R")
source("file_load.R")
options(scipen = 999)

#----------------------------------------

# look for unit issues, generally metals are in ppm and everything else in ppb
# 1 ug/g = 1 mg/kg = 1 ppm
# 1 ug/kg = 1 ppb 
# For ESB convert ppb back to ppm
#----------------------------------------

# ------Load lookup tables------
#uu <- "\U00B5"
#Encoding(df$Units) = "latin1" 


#________________Read in excel files, combine into one df________________________
# mostly used for multiple sed chem EDD files per year


# Tox test data is stores on multiple excel sheets. Read in all tab & sheets, combine row wise

Stations_tox=pin_get("EstProbMon_Stations", board = "rsconnect")


CEDS_Stations=pin_get("EstProbMon_DCLS", board = "rsconnect") %>%
filter(Year %in% 2016:2019)%>%
#dplyr::select(Fdt_Sta_Id,CBP_NAME,Ana_Sam_Mrs_Container_Id_Desc,Year)%>%
filter(!Ana_Sam_Mrs_Container_Id_Desc %in% "EB")


#====================== Toxicity testing data=======================


# several Station IDs don't match well
# "Sample description" is Ches Bay IDs but there are several naming errors with each year
# In some cases sites with duplicate tox tests (i.e. -S1 and -S2) do not have matching duplicate water chem samples (only "R" )
# PERCENT_SURVIVAL col is just for the "A" replicate. Use "MEAN_PERCENT_SURVIVAL" for overall test
# 2-JMS087.11	VA16-033A, VA06-0083A


Hyalella=c("TN-16-312","TN-17-238","TN-17-296","TN-18-593","TN-19-495","TN-19-523")


Tox_files_path="C:/Users/vvt97279/Documents/RStudio_Test/Toxicity_data"


readfilefun=function(folderpath){

file_list=list.files(folderpath,full.names=TRUE)
 
df=map_dfr(1:length(file_list),~rio::import_list(file_list[.x],rbind=TRUE)%>%
        mutate('Standard Deviation'=as.numeric('Standard Deviation'))%>%
        mutate_if(is.character,list(~na_if(., "N/A")))) %>%
        janitor::clean_names(.,"screaming_snake")%>%
        mutate(Year=lubridate::year(START_DATE)) %>% 
        #distinct(TEST_NUMBER,SAMPLE_DESCRIPTION,Year,.keep_all = T)
       mutate(Organism=ifelse(TEST_NUMBER %in% Hyalella,"Hyalella","Leptocheirus"))


return(df)

}

#====================================  


#====================================
  New_Files=readfilefun(Tox_files_path) %>%
  mutate(Ana_Sam_Mrs_Container_Id_Desc=
         case_when(
           str_detect(SAMPLE_DESCRIPTION,"-S1")~"S1",
           str_detect(SAMPLE_DESCRIPTION,"-S2") ~"S2",
           TRUE ~ "R")) %>%
  mutate(CBP_NAME=dupe_finder_fun(SAMPLE_DESCRIPTION)) %>%
  filter(!str_detect(CBP_NAME,c("Control|PRESS")))%>%
  drop_na(MEAN_SURVIVAL_PERCENT)%>%
mutate(CBP_NAME=
         case_when(
  CBP_NAME=="VA09-0005A"~"VA19-0005A",
  CBP_NAME=="VA90-0013A"~ "VA19-0013A",
  CBP_NAME=="VA09-0013A"~ "VA19-0013A",
  CBP_NAME=="VA19-004B"~ "VA19-0004B",
  CBP_NAME=="VA19-002A"~ "VA19-0002A",
  CBP_NAME=="VA19-0041"~ "VA19-0041A",
  CBP_NAME=="VA19-003A"~ "VA19-0003A",
  CBP_NAME=="VA18-033A"~ "VA18-0033A",
TRUE ~ CBP_NAME
))
                    
full_join(.,CEDS_Stations,by = c("Year","CBP_NAME")) 

# missing tox ? VA16-033A

# R's with 2 tox samples
c("VA18-0038A","VA18-0041A")


New_Files[which(!New_Files$CBP_NAME %in% CEDS_Stations$CBP_NAME),]


Non_match=New_Files %>%
group_by(Organism)%>%
anti_join(.,CEDS_Stations) %>% 
  dplyr::select(CBP_NAME,Year,Ana_Sam_Mrs_Container_Id_Desc) %>%
  filter(!str_detect(CBP_NAME,c("Control|PRESS")))


# purrr instead of rio
readfile_purr_fun=function(folderpath){

file_list=list.files(folderpath,full.names=TRUE)%>%
as_tibble() 

mutate(file_list,sheet_name = map(value,readxl::excel_sheets))%>% 
unnest(cols = c(sheet_name)) %>%
mutate(Tox_data=map2(value, sheet_name,~ readxl::read_excel(.x, .y))) %>% 
     
         
         
mutate()
         
#map_df(Tox_data,dplyr::)

}





file_list=list.files(Tox_files_path,full.names=TRUE)

Sheets=file_list%>%
as_tibble() %>%
mutate(sheet_name = map(value,readxl::excel_sheets))%>% 
unnest(cols = c(sheet_name))%>%
group_by(value) %>%
mutate(Tox_data=map2(value, sheet_name,~ readxl::read_excel(.x, .y))) %>%
mutate(Tox_data = map(Tox_data, ~  
                        mutate(.x,'Standard Deviation'=as.numeric('Standard Deviation'))%>%
                        mutate_if(is.character,list(~na_if(.x, "N/A")))%>%
                        mutate_at(Tox_data,7,as.numeric)))
                        

Files=rio::import_list(Tox_files_list,rbind=TRUE) %>% 
  rename(files="_file")%>% 
  mutate_if(is.character,list(~na_if(., "N/A")))%>%
  janitor::clean_names(.,"screaming_snake")%>%
  mutate(Year=lubridate::year(START_DATE)) %>%

%>%


  map_dfr(value
  ,~rio::import_list(value[.x],rbind=TRUE))

#==============

Tox_files_list=list.files(Tox_files_path,full.names=TRUE) %>%
set_names(nm=c("2017","2018","2019"))

yyy=rio::import_list(Tox_files_list,rbind=TRUE) %>% 
mutate(Year=lubridate::year(`Start Date`)) %>% 
rename(files="_file") 

%>%

sheet_name = map(Tox_files_list,~readxl::excel_sheets(.x))





#===================================



non_detects=c(NA,"ND","nd")

# duplicates cause issues when joining by ID
# Remove -S1 and -S2 so that salinities can be the same for dupes. 
# NCCA revisits will have differnt salinities but will also have different dates 
dupe_finder_fun=function(x){  
  str_replace(x,c("-S1|-S2"),"")
}



#=================================================================================


# Column names can be different than CAS, REsult, Units e.g. 'Real_results' would work
# creates Class column and converts ug/l to mg/l for metals
mutate_metals_fun=function(df,CAS,Result,Analyte,Units){

CAS=enquo(CAS) 
Result=enquo(Result)
Units=enquo(Units)
Analyte=enquo(Analyte)

df %>%
mutate(Class=case_when(
  !!Analyte %in% metals_CAS$NAME~ "METAL",
  !!CAS %in% PCB_Congener_CAS$CAS  & str_detect(!!Analyte,"Surr")==F | !!Analyte %in% "2,3,4,4',5-PeCB1,2"~ "PCB",
  !!Analyte %in% Pesticides_CAS$Name | !!Analyte %in% "Hexachlorocyclohexane" ~ "PESTICIDE",
  !!Analyte %in% PAHs_CAS$PAH ~ "PAH",
  
  TRUE ~ NA_character_),
Result_SMH=case_when(Class=="METAL" & !!Units %in% c("ug/Kg-dry","µg/Kg-dry") ~ !!Result/1000,
                     TRUE ~ !!Result ),
Units_SMH=case_when(Class=="METAL" ~ "mg/Kg-dry",
                    !!Units %in% "ug/Kg-dry" ~ "µg/Kg-dry",
           TRUE ~ !!Units ))

}




#===========add bottom salinity and salinity regime to metals data===================

Salinities=Toxicity_Salinitiy_fun(CEDS_Sites_Salinities,years=2017:2019,output="df")


Unique_Analytes=data.frame("Analyte"=unique(Sed_Chem_NEW$Analyte))
#### Working area#######

PAH_adj_values=PAH_values[1:15,] %>% 
left_join(PAHs_CAS)

Pesticides_CAS_adj=Pesticides_CAS %>%
drop_na(CAS)

PAH_Type=PAHs_CAS %>% 
dplyr::select(TYPE,PAH)


                       
# Non matching analytes
#PAH_adj_values[which(!PAH_adj_values$PAH %in% Unique_Analytes$analytes),]
#metals_CAS[which(!metals_CAS$NAME %in% Unique_Analytes$analytes),]
#PCBs_CAS[which(!PCBs_CAS$PCB_NAME%in% Unique_Analytes$analytes),]
#Pesticides_CAS[which(!Pesticides_CAS$Name %in% Unique_Analytes$analytes),]
#PAHs_CAS[which(!PAHs_CAS$PAH %in% Unique_Analytes$analytes),]


#=====================================================================================
# join chemistry data to salinities by siteID
# Some names don't match
#Sed_Chem_2017_18_19 from file_load
#(df,CAS,Result,Analyte,Units)

Full= function(x){

Sed_Chem_Data = x %>% 
mutate_if(is.character, ~gsub('[^ -~]', '', .))%>%
mutate_if(is.character,~gsub("<U\\+00B4>", "", .))%>%
mutate_metals_fun(CAS,Result,Analyte,Units) %>%
left_join(PAH_Type,by=c("Analyte"="PAH")) 


Totals=Sed_Chem_Data %>%
group_by(ClientSampID,Year)%>%
summarise(Total_LMW_PAHs = sum(Result_SMH[which(Class=="PAH" & TYPE=="LMW")]),
          Total_HMW_PAHs=sum(Result_SMH[which(Class=="PAH" & TYPE=="HMW")]),
          Total_PAHs=Total_LMW_PAHs+Total_HMW_PAHs,
          Total_PCBs=sum(Result_SMH[which(Class=="PCB")]),
          Total_Chlorodane= sum(Result_SMH[which(CAS %in% Total_Chlorodanes)]),
          Total_DDT= sum(Result_SMH[which(CAS %in% Total_DDTs)]))%>%
          gather(Analyte,Result_SMH, Total_LMW_PAHs:Total_DDT) %>%
          mutate(Class="Totals",Units_SMH="µg/Kg-dry")

Modified_Sed_Chem_Data =
bind_rows(Sed_Chem_Data,Totals) %>%
left_join(ERM_PEC_PEL,by="Analyte")%>%
mutate(ERM_Q=Result_SMH/ERM,PEC_Q=Result_SMH/PEC,
       ERM_PEC_min_Q=Result_SMH/ERM_PEC_min,
       PEL_Q=Result_SMH/PEL,ERL_Q=Result_SMH/ERL)%>%
  group_by(ClientSampID,Year)

return(Modified_Sed_Chem_Data)
}          
#===========

try=Full(Sed_Chem_2017_18_19)

#===========================




#distinct(Year,ClientSampID)%>%
#group_by(Year)%>%
#tally()

 
Station_data_2011_2020= read_csv("App_lookup_Tables/Station_data_2011_2020.csv") %>%
janitor::clean_names(.,"screaming_snake") 



Sediment_Chem_18_19= try %>%
filter(Year %in% 2018:2019)

Stations_18_19=Station_data_2011_2020 %>% 
filter(YEAR %in% 2018:2019) %>%
select(1:4)

"2A" "3B" "5A" "2B"



#============= Hyland and global Quotients=====================
# means >0.1 indicate 75% chance of observing toxicity           
#Hylands=c("Nickel","Dieldrin","4,4-DDD","4,4-DDT","Total_PAHs","Total_LMW_PAHs","Total_HMW_PAHs")  
Hyland_values = read_csv("App_lookup_Tables/Hyland_values.csv")
Global_values = read_csv("App_lookup_Tables/Global_values.csv")

All_Calculated=try %>% 
#filter(str_detect(ClientSampID,"PRESS"))%>%
#filter(ClientSampID %in% "PRESS-10") %>%
group_by(ClientSampID,Year) %>% 
 summarise(ERM_Q_Globalmean=mean(ERM_Q[which(Analyte %in% Global_values$Global)],na.rm=TRUE),
           PEC_Q_Globalmean=mean(PEC_Q[which(Analyte %in% Global_values$Global)],na.rm=TRUE),
           Min_Global=mean(ERM_PEC_min_Q[which(Analyte %in% Global_values$Global)],na.rm=TRUE),
           PEL_Q_Globalmean=mean(PEL_Q[which(Analyte %in% Global_values$Global)],na.rm=TRUE),
ERM_Hyland=mean(ERM_Q[which(Analyte %in% Hyland_values$Hyland)],na.rm=TRUE),
PEC_Hyland=mean(PEC_Q[which(Analyte %in% Hyland_values$Hyland)],na.rm=TRUE),
ERM_PEC_min_Hyland=mean(ERM_PEC_min_Q[which(Analyte %in% Hyland_values$Hyland)],na.rm=TRUE),
PEL_Hyland=mean(PEL_Q[which(Analyte %in% Hyland_values$Hyland)],na.rm=TRUE)                
     )



#=======================================================
Press_Calc=All_Calculated %>%
filter(str_detect(ClientSampID,"PRESS"))%>%
#filter(ClientSampID %in% "PRESS-10") %>%
dplyr::select(ClientSampID,Year,PEC_Hyland,ERM_Hyland)

Press_full=try%>%
filter(str_detect(ClientSampID,"PRESS"))%>%
  

#============================================================

try %>%
filter(Analyte %in% "Manganese" & Result_SMH >0) %>%
ggplot() + 
  geom_hline(yintercept = 1100,size=1.25,color="red") +
  geom_point(aes(x=ClientSampID,y=Result_SMH),color="blue",size=3)+
  geom_point(data=Press_full %>% filter(Analyte %in% "Manganese"  & Result_SMH >0),
             aes(x=ClientSampID,y=Result_SMH),color="green",size=3)+
  facet_wrap(~Year)+
  #ggtitle(paste("Count of Exceedances=",Title))+
  theme_bw()



pRESS_TRY %>% 
  filter(Year %in% 2019) %>% 
  dplyr::select(PEC_Hyland,PEC_Q_Globalmean)

        
PRESS10=try %>%
filter(ClientSampID %in% "PRESS-10", Year %in% 2019) %>%
  dplyr::select(Analyte,PEC_Q)
##############################################################################


plot_fun=function(df,x){

x=enquo(x)
  
Title= df %>% 
  filter(str_detect(ClientSampID,"PRESS") & Class=="METAL") %>% 
  tally(!!x > 1)
  
df %>% 
  drop_na(!!x)%>%
  filter(str_detect(ClientSampID,"PRESS") & Class=="METAL") %>% 
  ggplot() + 
  geom_hline(yintercept = 1,size=1.25) +
  geom_point(aes(x=ClientSampID,y=!!x,color=Analyte),size=3)+
  facet_wrap(~Year)+
  ggtitle(paste("Count of Exceedances=",Title))+
  theme_bw() 

}


try_cols=c("ERM_Q","PEC_Q","PEL_Q","ERL_Q")

 
plot_fun(try,ERL_Q)


