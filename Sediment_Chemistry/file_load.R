
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

