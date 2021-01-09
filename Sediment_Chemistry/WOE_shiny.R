



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

Tox_df=readfilefun(Tox_files_path)%>%
  drop_na(MEAN_SURVIVAL_PERCENT)%>%
  filter(!str_detect(SAMPLE_DESCRIPTION,"Control"))%>%
  mutate(Ana_Sam_Mrs_Container_Id_Desc=case_when(
    str_detect(SAMPLE_DESCRIPTION,"-S1")~"S1",
    str_detect(SAMPLE_DESCRIPTION,"-S2") ~"S2",
    TRUE~"R"))%>%
  #add_count(SAMPLE_DESCRIPTION, wt = length(unique(DATE_COLLECT)))%>% 
  #group_by(SITE_ID)%>%
  mutate(Ana_Sam_Mrs_Container_Id_Desc=
           case_when(
             SAMPLE_DESCRIPTION %in% "NCA20_VA-10001"~
               case_when(
                 TEST_NUMBER == "TN-20-513" ~"V1",
                 TRUE ~"V2"),
             SAMPLE_DESCRIPTION %in% "NCA20_VA-10002"~
               case_when(
                 LABORATORY_ACCESSION_NUMBER == "AT0-752" ~"V2",
                 TRUE ~"V1"),
             TRUE~Ana_Sam_Mrs_Container_Id_Desc
           ))%>% 
  mutate(CBP_NAME=dupe_finder_fun(SAMPLE_DESCRIPTION))%>%
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
             CBP_NAME=="PRESS-Alt4"~ "PRESS-4",
             CBP_NAME=="PRESS-ALt2"~ "PRESS-Alt2",
             TRUE ~ CBP_NAME
           ))

CEDS_Stations=pin_get("EstProbMon_DCLS", board = "rsconnect") %>%
  filter(Year %in% 2016:2020)%>%
  filter(!Ana_Sam_Mrs_Container_Id_Desc %in% "EB")


# readfilefun is only using files with "Toxicity" in file name
# list.files(getwd(),pattern="(Toxicity).*\\.(xlsx|xls)$", full.names=TRUE) 

Tox_2016_2020=Tox_df %>%
  left_join(.,stations_fix2)%>%
  mutate(SMH_BIO_SIG=ifelse(CONTROL_CORRECTED_SURVIVAL>80,"N","Y"))%>%
  Don_TOX_fun() %>%
  NCCA_TOX_fun()%>%
  Cal_TOX_fun()

pin(Tox_2016_2020,"EstProbMon_ToxTests_2016_20",board = "rsconnect")

not=tibble(STATIONS=not_in_fun((Tox_df2$STATIONS,stations_fix2$CBP_NAME))
           matches=stringdist_left_join(not,stations_fix2,by=c("STATIONS"="CBP_NAME"),max_dist = 1)
           

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
  mutate_if(is.character, ~gsub('[^ -~]', '', .))%>%
  mutate_if(is.character,~gsub("<U\\+00B4>", "", .))%>%
  mutate_at(vars(contains('Date')|contains('date')),~lubridate::mdy(.))%>%
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
left_join(PAHs_CAS,by=c("Analyte"="PAH","CAS")) 


Totals=Sed_Chem_Data %>%
group_by(Fdt_Sta_Id,Ana_Sam_Mrs_Container_Id_Desc,Year)%>%
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

#==================== sed chem cleaning =======================================================================================


Sed_Chem_2017_18_19 = readr::read_csv("Sediment_Chemistry/Sed_Chem_2017_18_19.csv") %>%
#mutate_metals_fun(CAS,Result,Analyte,Units)%>%
mutate(Ana_Sam_Mrs_Container_Id_Desc=
           case_when(
           str_detect(ClientSampID,"-S1")~"S1",
            str_detect(ClientSampID,"-S2") ~"S2",
           TRUE ~"R"
           ))   %>%
  left_join(stations_fix,by=c("ClientSampID_SMH"='CBP_NAME'))%>%
  #janitor::clean_names(.,"screaming_snake")%>% 
  mutate(DATE_COLLECTED=mdy(DateCollected),YEAR=lubridate::year(DATE_COLLECTED))%>%
 dplyr::select(SITE_ID=ClientSampID,Fdt_Sta_Id,CBP_NAME=ClientSampID_SMH,Ana_Sam_Mrs_Container_Id_Desc,DATE_COLLECT=DATE_COLLECTED,
               PARAMETER=Analyte,RESULT=Result,UNIT=Units,CAS_NO=CAS,MDL,QA_CODES=Qualifier,YEAR)
 

#============ 2016 2015 sed chem============================

VA_2015_SEDCHEM_ENHANCEMENT = read_excel("Sediment_Chemistry/VA_2015_SEDCHEM_ENHANCEMENT_RAW (Draft data).xlsx")%>%
mutate(Ana_Sam_Mrs_Container_Id_Desc=
         case_when(
           str_detect(SITE_ID,"-S1")~"S1",
           str_detect(SITE_ID,"-S2") ~"S2",
           TRUE ~"R"
         ),CBP_NAME=dupe_finder_fun(SITE_ID))%>%
  left_join(a)%>%
  rename(Fdt_Sta_Id=STATION_ID)%>%
  mutate(DATE_COLLECT=lubridate::mdy(DATE_COLLECT))%>%
  mutate(RESULT=as.numeric(RESULT),MDL=as.numeric(MDL))%>%
  dplyr::select(SITE_ID,Fdt_Sta_Id,CBP_NAME,Ana_Sam_Mrs_Container_Id_Desc,DATE_COLLECT,PARAMETER,RESULT,UNIT,CAS_NO,MDL,QA_CODES=QC_CODE)%>%
  mutate(YEAR=2015)


Copy_of_NCCA15 = read_excel("Sediment_Chemistry/Copy of NCCA15_VA_SedChem_PHYSIS_data(DRAFT).10.25.2016.hjs.xlsx")%>%
left_join(a,by=c("SITE_ID"="value")) %>%
rename(Fdt_Sta_Id=STATION_ID)%>%
add_count(SITE_ID, wt = length(unique(DATE_COLLECT)))%>% 
group_by(SITE_ID)%>%
mutate(Ana_Sam_Mrs_Container_Id_Desc=
case_when(
  n==1 ~"R",
 DATE_COLLECT == first(DATE_COLLECT) ~"V2",
 DATE_COLLECT == last(DATE_COLLECT) ~"V1"
 )) %>%
mutate(DATE_COLLECT=as.Date(DATE_COLLECT))%>%
dplyr::select(SITE_ID,Fdt_Sta_Id,CBP_NAME,Ana_Sam_Mrs_Container_Id_Desc,DATE_COLLECT,PARAMETER,RESULT,UNIT,CAS_NO,MDL,QA_CODES)%>%
mutate(YEAR=2015)
  

sed_chem_2015=bind_rows(VA_2015_SEDCHEM_ENHANCEMENT,Copy_of_NCCA15)

#===================2016===================================================================================================
SedChem2016_Condolidated = read_excel("Sediment_Chemistry/SedChem2016 - Condolidated.xls") %>%
mutate(CBP_NAME=dupe_finder_fun(ClientSampID)) %>%
left_join(Station_data_2011_2020[,1:2])%>%
rename(Fdt_Sta_Id=STATION_ID)%>%
mutate(RESULT=as.numeric(str_replace_all(R_Rslt,c("'|,"),"")),MDL=as.numeric(str_replace_all(R_MDL,c("'|,"),"")),DATE_COLLECT=as.Date(DateCollected),YEAR=2016,Ana_Sam_Mrs_Container_Id_Desc=
           case_when(
             str_detect(ClientSampID,"-S1")~"S1",
             str_detect(ClientSampID,"-S2") ~"S2",
             TRUE ~"R"
           ))%>%
  
dplyr::select(SITE_ID=ClientSampID,Fdt_Sta_Id,CBP_NAME,Ana_Sam_Mrs_Container_Id_Desc,DATE_COLLECT,PARAMETER=R_Analyte,RESULT,
              UNIT=R_Units,CAS_NO=CAS,MDL,QA_CODES=R_Qual,YEAR)




Sediment_Chem_2015_2019=bind_rows(SedChem2016_Condolidated,sed_chem_2015,Sed_Chem_2017_18_19)
write.csv(Sediment_Chem_2015_2019,"Sediment_Chemistry/Sediment_Chem_2015_2019.csv")



try=Full(Sed_Chem_2017_18_19)

#===========================

dupe_finder_fun=function(x){  
  str_replace(x,c("-S1|-S2"),"")
}



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

#==========================================================================

#====================== cleaning BUg data ==================================

# Need to do some string matching of ODU's "Station names" to Dons CBP_Names (i.e with letters added)
stations_fix=CEDS_Stations[ ,c(1:2)]

str_right =function(string, n) {
  substr(string, nchar(string) - (n - 1), nchar(string))
}
# ===================2019 bugs===============================

#------------------- EMAP and MAIA----------------------------------

NCA_MAIA_EMAP2019 =  read_excel("~/RStudio_Test/EMAP_MAIA/NCA_MAIA_EMAP2019.xlsx", 
                                sheet = "NCA_MAIA_EMAP2019")%>%
  rename(Fdt_Sta_Id=`DEQ Stion ID` ) %>%
  group_by(Fdt_Sta_Id)%>%
  add_tally()%>%
  mutate(duplicates=make.unique(STATION,sep="-S"),
         Ana_Sam_Mrs_Container_Id_Desc=
           case_when(
             n==2 & str_detect(duplicates,"-S1")~"S1",
             n==2 & str_detect(duplicates,"-S1")==FALSE ~"S2",
             n==1 ~"R"
           ))   %>%
  left_join(stations_fix) %>%
  dplyr::select(Fdt_Sta_Id,sort(colnames(.)),-STATION,-n,-duplicates)


#------------------ AMBI and MAMBI -----------------------------------

AMBI_ProbMon_2019= read_excel("~/RStudio_Test/AMBI/AMBI_ProbMon_2019.xlsx") %>%
  rename(STATIONS=`Stations( & Rep)`)%>%
  filter(!str_detect(STATIONS,"PR"))%>%
  separate(STATIONS, c("A", "B","C"),"-") %>% 
  mutate(STATION=paste0(A,"-",B)) %>%
  stringdist_left_join(stations_fix,by=c("STATION"="CBP_NAME"), max_dist = 1)%>%
  distinct(Fdt_Sta_Id,C,AMBI,.keep_all = T)%>%
  group_by(Fdt_Sta_Id) %>%
  add_tally()%>%
  mutate(duplicates=make.unique(Fdt_Sta_Id,sep="-S"),
         Ana_Sam_Mrs_Container_Id_Desc=
           case_when(
             n==2 & str_detect(duplicates,"-S1")~"S1",
             n==2 & str_detect(duplicates,"-S1")==FALSE ~"S2",
             n==1 ~"R"
           )) %>%
  dplyr::select(Fdt_Sta_Id,sort(colnames(.)),-STATION,-A,-B,-C,-n,-duplicates)%>%
  left_join(NCA_MAIA_EMAP2019,by=c("Fdt_Sta_Id","Ana_Sam_Mrs_Container_Id_Desc"))

#--------------------------------------------------------------------------


#Individuals   =  or TOTAL_ABUND
#Taxa          = N_SP_REP_EPI or TOTAL_SPECIES
#Shannon H'   
#Gleason D   
#Pielou Eveness   
#Tubificidae   
#Spionidae   

# N_SP_REP_EPI	and TOTAL_ABUND_EPI match Dons numbers

# ches bay bibi
Ches_bay_2019_SMH=read_excel("~/RStudio_Test/VARIB/VARIBI5NCA_2019.xlsx") %>%
  rename(Fdt_Sta_Id=`DEQ Station ID`)%>%
  filter(!str_detect(STATION,"PR"))%>%
  left_join(stations_fix)%>%
  distinct(Fdt_Sta_Id,ABUN_MSQ,.keep_all = T)%>%
  group_by(Fdt_Sta_Id)%>%
  add_tally() %>%
  mutate(duplicates=make.unique(Fdt_Sta_Id,sep="-S"),
         Ana_Sam_Mrs_Container_Id_Desc=
           case_when(
             n==2 & str_detect(duplicates,"-S1")~"S1",
             n==2 & str_detect(duplicates,"-S1")==FALSE ~"S2",
             n==1 ~"R"
           ))%>%
  dplyr::select(Fdt_Sta_Id,sort(colnames(.)),-STATION,-n,-duplicates)


#---------------------------------------------------------------------------------------------------------------

Raw_Bugs_2019=Ches_bay_2019_SMH %>%
  left_join(AMBI_ProbMon_2019,by=c("Fdt_Sta_Id","Ana_Sam_Mrs_Container_Id_Desc"))%>%
  mutate(Year=2019) 


Clean_bugs_2019=Raw_Bugs_2019%>%
  bug_scores_fun(Basin=BASIN,CB=B_IBI,MAIA=MAIA_INDEX,EMAP=EMAP_INDEX,MAMBI=`M-AMBI`) %>%
  dplyr::select(CBP_NAME,Year,Ana_Sam_Mrs_Container_Id_Desc,BASIN,B_IBI,MAIA_INDEX,EMAP_INDEX,MAMBI=`M-AMBI`,Matrix_score)
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------



# ===================  2018 bugs  ===============================


NCA_MAIA_EMAP2018 = read_excel("~/RStudio_Test/EMAP_MAIA/NCA_MAIA_EMAP18new.xlsx", 
                               sheet = "NCA_MAIA_EMAP18new")%>%
  stringdist_left_join(stations_fix,by=c("STATION"="CBP_NAME"), max_dist = 1)%>%
  distinct(STATION,DATE,REP_NUM,.keep_all = T)%>%
  group_by(Fdt_Sta_Id)%>%
  add_tally()%>%
  mutate(duplicates=make.unique(Fdt_Sta_Id,sep="-S"),
         Ana_Sam_Mrs_Container_Id_Desc=
           case_when(
             n==2 & str_detect(duplicates,"-S1")~"S1",
             n==2 & str_detect(duplicates,"-S1")==FALSE ~"S2",
             n==1 ~"R"
           )) %>%
  dplyr::select(Fdt_Sta_Id,sort(colnames(.)),-STATION,-n,-duplicates)


# AMBI and MAMBI 
AMBI_ProbMon_2018= read_excel("~/RStudio_Test/AMBI/5. NCA2018_M_AMBI.xlsx")%>%
  stringdist_left_join(stations_fix,by=c("Stations"="CBP_NAME"), max_dist = 1) %>%
  distinct(Stations,AMBI,.keep_all = T)%>%
  group_by(Fdt_Sta_Id) %>%
  add_tally()%>%
  mutate(duplicates=make.unique(Fdt_Sta_Id,sep="-S"),
         Ana_Sam_Mrs_Container_Id_Desc=
           case_when(
             n==2 & str_detect(duplicates,"-S1")~"S1",
             n==2 & str_detect(duplicates,"-S1")==FALSE ~"S2",
             n==1 ~"R"
           )) %>%
  dplyr::select(Fdt_Sta_Id,sort(colnames(.)),-Stations,-n,-duplicates)%>%
  left_join(NCA_MAIA_EMAP2018,by=c("Fdt_Sta_Id","Ana_Sam_Mrs_Container_Id_Desc"))


# Ches bay Bibi
Ches_bay_2018_SMH = read_excel("~/RStudio_Test/VARIB/VARIBI5NCA_2018.xlsx") %>%
  stringdist_left_join(eighteen,by=c("STATION"="CBP_NAME"), max_dist = 1) %>%
  distinct(CBP_NAME,DATE,REP_NUM,B_IBI,.keep_all = T)%>%
  group_by(Fdt_Sta_Id)%>%
  add_tally() %>%
  mutate(duplicates=make.unique(Fdt_Sta_Id,sep="-S"),
         Ana_Sam_Mrs_Container_Id_Desc=
           case_when(
             n==2 & str_detect(duplicates,"-S1")~"S1",
             n==2 & str_detect(duplicates,"-S1")==FALSE ~"S2",
             n==1 ~"R"
           )) %>%
  dplyr::select(Fdt_Sta_Id,sort(colnames(.)),-STATION,-n,-duplicates)

#-------------------------------------------------------------------------------------------------------
Raw_Bugs_2018=Ches_bay_2018_SMH %>%
left_join(AMBI_ProbMon_2018,by=c("Fdt_Sta_Id","Ana_Sam_Mrs_Container_Id_Desc"))%>%
select(-DATE.y,-REP_NUM.y)%>%
rename(DATE=DATE.x,REP_NUM=REP_NUM.x)%>%
mutate(Year=2018)


Clean_bugs_2018=Raw_Bugs_2018 %>%
bug_scores_fun(Basin=BASIN,CB=B_IBI,MAIA=MAIA_INDEX,EMAP=EMAP_INDEX,MAMBI=`M-AMBI`) %>%
dplyr::select(CBP_NAME,Year,Ana_Sam_Mrs_Container_Id_Desc,BASIN,B_IBI,MAIA_INDEX,EMAP_INDEX,MAMBI=`M-AMBI`,Matrix_score)



# ============== 2017 bugs=========================================

stations_fix2 =read_csv("stations_fix.csv")
# ===========2017 EMAP===============
NCCA_POT_EMAP2017 <- read_excel("EMAP_MAIA/NCCA&POT_EMAP2017.xls") %>%
mutate(EMAP_INDEX=EMAP_VP_INDEX2,EMAP_STATUS=NA)%>%
dplyr::select(STATION,EMAP_INDEX,EMAP_STATUS)

#=========== 2017 MAIA=================
NCCA_POT_MAIA2017 <- read_excel("EMAP_MAIA/NCCA&POT_MAIA2017.xlsx", 
                                sheet = "NCA_MAIA_EMAP17")%>%
left_join(NCCA_POT_EMAP2017)%>%
rename(Fdt_Sta_Id='DEQ           Station Id.')%>%
mutate(STATION=press_names(STATION))%>%
stringdist_left_join(stations_fix2,by="Fdt_Sta_Id",max_dist = 1)%>%
distinct(CBP_NAME,TOTAL_DEN_REP,.keep_all = T)%>%
group_by(CBP_NAME)%>%
add_tally()%>%
  mutate(duplicates=make.unique(CBP_NAME,sep="-S"),
         Ana_Sam_Mrs_Container_Id_Desc=
           case_when(
             n==2 & str_detect(duplicates,"-S1")~"S1",
             n==2 & str_detect(duplicates,"-S1")==FALSE ~"S2",
             n==1 ~"R"
           ))%>%
  dplyr::select(Fdt_Sta_Id=Fdt_Sta_Id.x,sort(colnames(.)),-STATION,-n,-duplicates)

#==================== AMBI==================
NCA_2017b_M_AMBI = read_excel("AMBI/NCA_2017b_M-AMBI.xls")%>%
slice(3:n())%>%
mutate(Stations=press_names(Stations))%>%
  mutate(Stations=case_when(
    Stations %in% "VA17-0003"~paste0(Stations,"C"),
    Stations %in% c("VA17-0014","VA17-0006","VA17-0023","VA17-0031","VA17-0039","VA17-0043","VA17-0048")~paste0(Stations,"B"),
    str_detect(Stations,"PRE") ~Stations,
    TRUE~paste0(Stations,"A")))%>%
left_join(stations_fix2,by=c("Stations"="CBP_NAME")) %>%
rename(CBP_NAME=Stations)%>%
distinct(Stations,AMBI,.keep_all = T)%>%
group_by(CBP_NAME) %>%
add_tally()%>%
  mutate(duplicates=make.unique(CBP_NAME,sep="-S"),
         Ana_Sam_Mrs_Container_Id_Desc=
           case_when(
             n==2 & str_detect(duplicates,"-S1")~"S1",
             n==2 & str_detect(duplicates,"-S1")==FALSE ~"S2",
             n==1 ~"R"
           )) %>%
  dplyr::select(Fdt_Sta_Id,sort(colnames(.)),-Stations,-n,-duplicates)
 

  


press_names=function(x){
  case_when(
    str_detect(x,"VA17-PR01")~"PRESS-1",
    str_detect(x,"VA17-PR02")~"PRESS-Alt2",
    str_detect(x,"VA17-PR03")~"PRESS-Alt3",
    str_detect(x,"VA17-PR04")~"PRESS-4",
    str_detect(x,"VA17-PR05")~"PRESS-5",
    str_detect(x,"VA17-PR06|VA17-PR10") ~"PRESS-10",
    str_detect(x,"VA17-PR07|VA17-PR11") ~"PRESS-11",
    str_detect(x,"VA17-PR08|VA17-PR12") ~"PRESS-12",
    TRUE~x)
}


# Ches bay BIBI
VARIB_2017_f = read_excel("~/RStudio_Test/VARIB/VABIBINCCA_&POT_2017.xls") %>%
  stringdist_left_join(stations_fix %>% filter(str_detect(CBP_NAME,"VA17")),by=c("STATION"="CBP_NAME"), max_dist = 1) %>%
  filter(!str_detect(STATION,"PR"))%>%
  distinct(CBP_NAME,DATE,REP_NUM,B_IBI,.keep_all = T)%>%
  group_by(Fdt_Sta_Id)%>%
  add_tally() %>%
  #dplyr::select(`DEQ Station ID`,n) %>%
  mutate(duplicates=make.unique(Fdt_Sta_Id,sep="-S"),
         Ana_Sam_Mrs_Container_Id_Desc=
           case_when(
             n==2 & str_detect(duplicates,"-S1")~"S1",
             n==2 & str_detect(duplicates,"-S1")==FALSE ~"S2",
             n==1 ~"R"
           )) %>%
  mutate(Year="2017",CB_weight=ifelse(CBPSEG=="CBAYS",0,2))%>%
  dplyr::select(CBP_NAME,Fdt_Sta_Id,Year,CB_IBI=B_IBI,CB_weight,Ana_Sam_Mrs_Container_Id_Desc,SHANNON,HABITAT,CBPSEG,BASIN,WATER_BODY  
                ,STRATUM,-n,-duplicates) 




#========================== 2016 bugs=====================================

sixteen=
  stations_fix %>% 
  filter(str_detect(CBP_NAME,"VA16")) %>% 
  mutate(letters=str_right(CBP_NAME,1))

#CB bibi
VARIB_2016_f = read_excel("~/RStudio_Test/VARIB/VARIBI5NCA_2016ModifiedMediomastus.XLS", 
                          skip = 1)%>%
  mutate(STATION=str_remove(STATION, "0"))%>%
  stringdist_left_join(sixteen,by=c("STATION"="CBP_NAME"), max_dist = 1) %>%
  distinct(STATION,DATE,REP_NUM,B_IBI,.keep_all = T)%>%
  group_by(Fdt_Sta_Id)%>%
  add_tally() %>%
  #dplyr::select(`DEQ Station ID`,n) %>%
  mutate(duplicates=make.unique(Fdt_Sta_Id,sep="-S"),
         Ana_Sam_Mrs_Container_Id_Desc=
           case_when(
             n==2 & str_detect(duplicates,"-S1")~"S1",
             n==2 & str_detect(duplicates,"-S1")==FALSE ~"S2",
             n==1 ~"R"
           )) %>%
  mutate(Year="2016",CB_weight=ifelse(CBPSEG=="CSTBY",0,2))%>%
  dplyr::select(CBP_NAME,Fdt_Sta_Id,Year,CB_IBI=B_IBI,CB_weight,Ana_Sam_Mrs_Container_Id_Desc,SHANNON,HABITAT,CBPSEG,BASIN,WATER_BODY  
                ,STRATUM,-n,-duplicates) 


#=======================================================================================
# VARIB_2015 =read_excel("~/RStudio_Test/VARIB/VARIBI5NCA_2015.XLS")



CB_IBI_2016_2019=rbind(VARIB_2016_f,VARIB_2017_f,VARIB_2018_f,VARIB_2019_f)
#pin(CB_IBI_2016_2019,"EstProbMon_Bugs_2016_2019",board = "rsconnect")

write.csv(CB_IBI_2016_2019,"Bugs/Bugs_2016_2019.csv")

#pin(Sediment_Chem_2015_2019,"EstProbMon_Sed_Chem_2015_2019",board = "rsconnect")
#=======================================================================================





# check names any(!names(EDD_2017) %in% c(names(EDD_2018),names(EDD_2019))


#EDD_2019=read_csv("Chemistry_data/EDD_2019_combined_original.csv")%>%
#  mutate_if(is.character, ~gsub('[^ -~]', '', .))%>%
#  mutate(Result=as.numeric(ifelse(Result %in% non_detects,0,Result)))%>%
#  mutate(ClientSampID_SMH=dupe_finder_fun(ClientSampID))


#EDD_2018=read_csv("Chemistry_data/EDD_2018_combined_original.csv")%>%
#  mutate_if(is.character, ~gsub('[^ -~]', '', .))%>%
#  mutate(Result=as.numeric(ifelse(Result %in% non_detects,0,Result)))%>%
#  mutate(ClientSampID_SMH=dupe_finder_fun(ClientSampID))


#EDD_2017=read_csv("Chemistry_data/EDD_2017_combined_original.csv")%>%
#  mutate_if(is.character, ~gsub('[^ -~]', '', .))%>%
#  mutate(Result=as.numeric(ifelse(Result %in% non_detects,0,Result))) %>%
#  mutate(ClientSampID_SMH=dupe_finder_fun(ClientSampID))

#EDD_2016 <- read_csv("Chemistry_data/EDD_2016.csv")%>%
#  mutate_if(is.character, ~gsub('[^ -~]', '', .))%>%
#  rename_all(.,~gsub('R_', '', .))  %>%
#  rename("Result"="Rslt","Qualifier"="Qual","AnalysisDate"="AnalDate") %>%
#  mutate(Result=ifelse(Result=="ND" | Result=="nd" ,0,Result))

#Master_171819 = rbind(EDD_2019,EDD_2018,EDD_2017)%>%
#distinct(Year,SampID,Analyte,.keep_all=TRUE)

#write.csv(Master_171819,"Sed_Chem_2017_18_19.csv")



