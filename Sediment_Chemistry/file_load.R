
#_____________ read lookup tables___________________

# Metals, pesticides, PCBs and PAHs - ERM,PEC,PEL values in "values"
# check units are in mg/kg (ppm) for metals and ug/kg (ppb) for PAH & PCB

ERM_PEC_PEL = read_csv("App_lookup_Tables/ERM_PEC_PEL.csv")%>%
mutate(ERM_PEC_min=ifelse(ERM>0 & ERM<PEC,ERM,PEC))

#~~~~~~~~~~~~~  DEQ Station IDs (CEDS) with NCCA or Bay Program IDs  ~~~
Station_data_2011_2020= read_csv("App_lookup_Tables/Station_data_2011_2020.csv") %>%
  janitor::clean_names(.,"screaming_snake") 

#~~~~~~~~~~~~~~ metals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
metals_values =read_csv("App_lookup_Tables/metals_values.csv") 
metals_CAS = read_csv("App_lookup_Tables/metals_CAS.csv")
metals_full=left_join(metals_CAS,metals_values,by=c("NAME","CAS"))
#~~~~~~~~~~~~~~~ pesticides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pesticide_values =read_csv("App_lookup_Tables/pesticide_values.csv") %>%
  mutate_if(is.character, ~gsub('[^ -~]', '', .))

Pesticides_CAS = read_csv("App_lookup_Tables/Pesticides_CAS.csv")

#~~~~~~~~~~~~~~~~ PCBs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PCBs_CAS = read_csv("App_lookup_Tables/PCBs_CAS.csv")%>%
  mutate_if(is.character, ~gsub('[^ -~]', '', .))

PCB_values=tibble(PCB="Total_PCBs",ERM=180,PEC=676,PEL=189)

PCB_Congener_CAS = read_csv("App_lookup_Tables/PCB_Congener_CAS.csv")
#~~~~~~~~~~~~~~~~ PAHs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PAH_values = read_csv("App_lookup_Tables/PAH_values.csv") %>% 
dplyr::select(1:4)

PAHs_CAS= read_csv("App_lookup_Tables/PAHs_CAS.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Total_DDTs=c("72-54-8","72-55-9","50-29-3","53-19-0","3424-82-6","789-02-6")
Total_Chlorodanes=c("5103-71-9","5103-74-2","26880-48-8")
Total_PCBs=PCBs_CAS$CAS


#============Sediment Chemistry Original Files ====================================

Sed_Chem_2017_18_19 = read_csv("Sed_Chem_2017_18_19.csv")



#====================== cleaning BUg data ==================================

# Need to do some string matching of ODU's "Station names" to Dons CBP_Names (i.e with letters added)
stations_fix=CEDS_Stations[ ,c(1:2)]

str_right =function(string, n) {
  substr(string, nchar(string) - (n - 1), nchar(string))
}
# ===================2019 bugs===============================

# EMAP and MAIA

NCA_MAIA_EMAP2019 =  read_excel("~/RStudio_Test/EMAP_MAIA/NCA_MAIA_EMAP2019.xlsx", 
                                                    sheet = "NCA_MAIA_EMAP2019")%>%
dplyr::select(STATION,`DEQ Stion ID`,STATUS,MAIA_INDEX,MAIA_Status=STATUS,EMAP_INDEX,
              EMAP_STATUS,TOTAL_IND,TOTAL_SPECIES,PIELOUS_EVENNESS_REP,TUBIFICIDAE,SPIONIDAE)%>%
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
  dplyr::select(everything(),-STATION,-n,-duplicates)

# AMBI and MAMBI 
AMBI_ProbMon_2019= read_excel("~/RStudio_Test/AMBI/AMBI_ProbMon_2019.xlsx") %>%
  rename(STATIONS=`Stations( & Rep)`)%>%
  filter(!str_detect(STATIONS,"PR"))%>%
  separate(STATIONS, c("A", "B","C"),"-") %>% 
  mutate(STATION=paste0(A,"-",B)) %>%
 stringdist_left_join(stations_fix,by=c("STATION"="CBP_NAME"), max_dist = 1)%>%
distinct(STATION,C,AMBI,.keep_all = T)%>%
  group_by(STATION) %>%
  add_tally()%>%
  mutate(duplicates=make.unique(STATION,sep="-S"),
         Ana_Sam_Mrs_Container_Id_Desc=
           case_when(
             n==2 & str_detect(duplicates,"-S1")~"S1",
             n==2 & str_detect(duplicates,"-S1")==FALSE ~"S2",
             n==1 ~"R"
           )) %>%
  dplyr::select(Fdt_Sta_Id,CBP_NAME,Ana_Sam_Mrs_Container_Id_Desc,everything(),-A,-B,-C,-n,-duplicates)%>%
  left_join(NCA_MAIA_EMAP2019,by=c("Fdt_Sta_Id","Ana_Sam_Mrs_Container_Id_Desc"))


           
#Individuals   =  or TOTAL_ABUND
#Taxa          = N_SP_REP_EPI or TOTAL_SPECIES
#Shannon H'   
#Gleason D   
#Pielou Eveness   
#Tubificidae   
#Spionidae   

# N_SP_REP_EPI	and TOTAL_ABUND_EPI match Dons numbers

# ches bay bibi
Bugs_2019_SMH=read_excel("~/RStudio_Test/VARIB/VARIBI5NCA_2019.xlsx") %>%
group_by(`DEQ Station ID`)%>%
add_tally() %>%
mutate(duplicates=make.unique(`DEQ Station ID`,sep="-S"),
       Ana_Sam_Mrs_Container_Id_Desc=
case_when(
  n==2 & str_detect(duplicates,"-S1")~"S1",
  n==2 & str_detect(duplicates,"-S1")==FALSE ~"S2",
  n==1 ~"R"
)) %>%
mutate(Year="2019")%>%
dplyr::select(Fdt_Sta_Id=`DEQ Station ID`,Year,B_IBI,Ana_Sam_Mrs_Container_Id_Desc,SHANNON,HABITAT,CBPSEG,BASIN,WATER_BODY  
              ,STRATUM,-n,-duplicates) %>%
left_join(AMBI_ProbMon_2019,by=c("Fdt_Sta_Id","Ana_Sam_Mrs_Container_Id_Desc")) %>%
bug_scores_fun(Basin=BASIN,CB=B_IBI,MAIA=MAIA_INDEX,EMAP=EMAP_INDEX,MAMBI=`M-AMBI`)

x,CB,MAIA,EMAP,MAMBI
# EMAP & MAIA 





# ===================  2018 bugs  ===============================

# Ches bay Bibi
VARIB_2018_f = read_excel("~/RStudio_Test/VARIB/VARIBI5NCA_2018.xlsx") %>%
stringdist_left_join(eighteen,by=c("STATION"="CBP_NAME"), max_dist = 1) %>%
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
  mutate(Year="2018",CB_weight=ifelse(CBPSEG=="CSTBY",0,2))%>%
  dplyr::select(CBP_NAME,Fdt_Sta_Id,Year,CB_IBI=B_IBI,CB_weight,Ana_Sam_Mrs_Container_Id_Desc,SHANNON,HABITAT,CBPSEG,BASIN,WATER_BODY  
                ,STRATUM,-n,-duplicates) 


# ============== 2017 bugs=========================================
  
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

