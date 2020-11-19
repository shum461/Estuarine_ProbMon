



#~~~~~~~~~~~~~  DEQ Station IDs (CEDS) with NCCA or Bay Program IDs  ~~~
Station_data_2011_2020= read_csv("Station_data_2011_2020.csv") %>%
janitor::clean_names(.,"screaming_snake") 

#---------------------------------------------------------------------------------


#_____________Salinity Regime Functions______________________

## Set salinity regime based on bottom salinity 
#See llanso_and_dauer._2002._Chesapeake Bay B-IBIfor cut offs ##

Sal_Regime_fun=function(y){
  case_when(
    is.numeric(y)==FALSE ~ NA_character_,
    y <= 0.5 ~ "TF",
    y > 0.5 & y <= 5 ~ "OH",
    y >5 & y <= 12 ~ "LM",
    y > 12 & y <= 18 ~ "HM",
    y > 18 & y <= 30 ~ "PO",
    y > 30 ~ "EH",
    TRUE ~ NA_character_)
  
}

#_________________________
## WQS salinity regimes ##
Sal_WQS_fun=function(y){
  case_when(
    is.numeric(y)==FALSE ~ NA_character_,
    y <= 0.5 ~ "Tidal Freshwater",
    y >0.5 & y <=5 ~ "Transitional Zone",
    y >5 & y <=30 ~ "Estuarine Waters",
    y>30 ~ "Oceanic Waters",
    TRUE ~ NA_character_)
}

#_________________________
# Ches Bay protocol only for high meso & polyhaline #
Silt_clay_fun=function(salinity,grainsize){
  case_when(
    is.numeric(grainsize)==FALSE ~ NA_character_,
    salinity <= 12 ~ NA_character_,
    grainsize >40 ~ "mud",
    grainsize<= 40 ~"sand")
}

#______________________________________-

# Pull percent sand from CEDS PART analysis
Percent_Sand_fun=function(percent_sand){
  
case_when(
  is.numeric(percent_sand)==FALSE ~ NA_character_,
  percent_sand<= 60 ~"Mud",
  percent_sand > 60 ~ "Sand",
  TRUE ~ NA_character_)
}

# Pull TOC from CEDS PART analysis
TOC_fun=function(percent_TOC){
  
  case_when(
    is.numeric(percent_TOC)==FALSE ~ NA_character_,
    percent_TOC<= 20 ~"Good",
    percent_TOC >20 & percent_TOC <= 60 ~ "Fair",
    percent_TOC >60 ~"Poor",
    TRUE ~ NA_character_)
}

# bottom DO NCCA criteria
NCCA_DO_fun=function(x){
  
  case_when(
    is.numeric(x)==FALSE ~ NA_character_, 
    x <= 2 ~ "Poor",
    x >2 & x <= 5 ~"Fair",
    x > 5 ~"Good"
  )}

NCCA_Chla_fun=function(x){
  
  case_when(
    is.numeric(x)==FALSE ~ NA_character_, 
    x < 5 ~ "Good",
    x >= 5 & x <= 20 ~"Fair",
    x > 20 ~"Poor"
  )}



#------------------------------ DCLS Analytes -------------------------------------------------------------------------

# Total Dissolved Nitrogen plus particulate nitrogen is used to calculate total nitrogen (TN=TDN+PN).
# Parameter short names (Pg_Parm_Short_Name) are blank for some rows and shouldn't be used to filter data. Use Pg_Storet_Code instead

# NITROGEN TOTAL, FIELD FILTERED, DISSOLVED,WTR MG/L =  49571
# NITROGEN PARTICULATE, FIELD FILT., SUSP., WTR MG/L =  49570
# PHOSPHOROUS TOTAL, FIELD FILTRED, DISSLVD,WTR MG/L =  49572
# PHOSPHOROUS PARTICULATE, FIELD FILT.,SUSP,WTR MG/L =  49567
# PERCENT SAND IN SEDIMENT ON A DRY WEIGHT BASIS =      82007
# CARBON, ORGANIC, IN BED MATERIAL (GM/KG AS C)  =      00687
# CHLOROPHYLL-A UG/L SPECTROPHOTOMETRIC ACID. METH.=    32211
 
# Need to make more concise & clean up with purrr

# param_names=c("Total_Nitrogen","Total_Phosphorus","Percent_Sand","TOC")  
# Parm_Codes=list(c("49571","49570"), c("49572","49567"),"8207","00687")
  
DCLS_Analytes_fun=function(x){
  
 
TN= x %>%
    group_by(Ana_Sam_Mrs_Container_Id_Desc,Ana_Sam_Fdt_Id,Fdt_Sta_Id,Fdt_Date_Time)%>%
    filter(Pg_Storet_Code %in% c("49571","49570"))%>%
    summarise(Total_Nitrogen=sum(Ana_Value))
  
TP= x %>%
    group_by(Ana_Sam_Mrs_Container_Id_Desc,Ana_Sam_Fdt_Id,Fdt_Sta_Id,Fdt_Date_Time)%>%
    filter(Pg_Storet_Code %in% c("49572","49567"))%>%
    summarise(Total_Phosphorus=sum(Ana_Value))
  
Sand= x %>%
    group_by(Ana_Sam_Mrs_Container_Id_Desc,Fdt_Sta_Id,Fdt_Date_Time)%>%
    filter(Pg_Storet_Code %in% "82007")%>%
    summarise(Percent_Sand=max(Ana_Value))%>%
    mutate(Subtrate_Type=Percent_Sand_fun(Percent_Sand))
  
TOC= x %>%
    group_by(Ana_Sam_Mrs_Container_Id_Desc,Fdt_Sta_Id,Fdt_Date_Time)%>%
    filter(Pg_Storet_Code %in% "00687")%>%
    summarise(TOC=max(Ana_Value))%>%
    mutate(TOC_Quality=TOC_fun(TOC))
  
Chla= x %>%
  group_by(Ana_Sam_Mrs_Container_Id_Desc,Fdt_Sta_Id,Fdt_Date_Time)%>%
  filter(Pg_Storet_Code %in% "32211")%>%
  summarise(Chlorophyll=max(Ana_Value)) %>%
  mutate(Chlorophyll_Quality=NCCA_Chla_fun(Chlorophyll))
  
  #data.frame(plyr::join_all(list(TN,TP,Sand), type='left'))
  Nutrients=left_join(TN,TP)
  Sand_Toc=left_join(Sand,TOC)
  
  Total=left_join(Nutrients,Sand_Toc,by=c("Ana_Sam_Mrs_Container_Id_Desc","Fdt_Sta_Id","Fdt_Date_Time")) %>%
  left_join(Chla)
           
  return(Total)          
}

# ==============CEDS all C2 ever---- use for Tox testing determinations===============






#____________________Just for toxicity salinities____________________
# x= latest C2 file from CEDS  
Toxicity_Salinitiy_fun=function(x,years,output="csv"){
  
  tox_sals=x %>%
   janitor::clean_names(.,"screaming_snake") %>%
    drop_na(SALINITY) %>%
    mutate(DATE_TIME=mdy(DATE_TIME),YEAR=year(DATE_TIME)) %>%
    filter(YEAR %in% years)%>%
    group_by(STATION_ID,DATE_TIME) %>% 
    arrange(desc(DEPTH)) %>% 
    slice(1) %>%
    left_join(Station_data_2011_2020,by=c("STATION_ID","YEAR")) %>%
    #filter(YEAR.y %in% years) %>%
    mutate(SALINITY_REGIME=Sal_Regime_fun(SALINITY),
           STANDARDS_SALINITY=Sal_WQS_fun(SALINITY)) %>%
    dplyr::select(YEAR,STATION_ID,CBP_NAME,DATE_TIME,SALINITY,SALINITY_REGIME,STANDARDS_SALINITY) %>%
    arrange(desc(DATE_TIME)) 
  
  df= data.frame(tox_sals)
  
  
  df.csv=write.csv(tox_sals,paste0("Toxicity_Testing_Salinity/Toxicity_Test_Salinity_Regimes_",Sys.Date(),".csv"))
  
  ifelse(output=="csv",return(df.csv),return(df))
  
}




