



#~~~~~~~~~~~~~  DEQ Station IDs (CEDS) with NCCA or Bay Program IDs  ~~~
Station_data_2011_2020= read_csv("Station_data_2011_2020.csv") %>%
janitor::clean_names(.,"screaming_snake") 

#---------------------------------------------------------------------------------

dupe_finder_fun=function(x){  
  str_replace(x,c("-S1|-S2"),"")
}

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



# =================== Reading file fun ===============================================


# Multiple excel sheets with multiple tabs. Gross


#====================== Toxicity testing data=======================


# several Station IDs don't match well
# "Sample description" is Ches Bay IDs but there are several naming errors with each year
# In some cases sites with duplicate tox tests (i.e. -S1 and -S2) do not have matching duplicate water chem samples (only "R" )
# PERCENT_SURVIVAL col is just for the "A" replicate. Use "MEAN_PERCENT_SURVIVAL" for overall test
# 2-JMS087.11	VA16-033A, VA06-0083A


Hyalella=c("TN-16-312","TN-17-238","TN-17-296","TN-18-593","TN-19-495","TN-19-523")


Tox_files_path="C:/Users/vvt97279/Documents/RStudio_Test/CEDS_data"


readfilefun=function(folderpath){
  
  file_list= list.files(folderpath,pattern="(Toxicity).*\\.(xlsx|xls)$", full.names=TRUE) 
  
  df=map_dfr(1:length(file_list),~rio::import_list(file_list[.x],rbind=TRUE)%>%
               mutate('Standard Deviation'=as.numeric('Standard Deviation'))%>%
               mutate_if(is.character,list(~na_if(., "N/A")))) %>%
    janitor::clean_names(.,"screaming_snake")%>%
    mutate(Year=lubridate::year(START_DATE)) %>% 
    #distinct(TEST_NUMBER,SAMPLE_DESCRIPTION,Year,.keep_all = T)
    mutate(Organism=ifelse(TEST_NUMBER %in% Hyalella,"Hyalella","Leptocheirus"))
  
  
  return(df)
  
}

#=================================================== 

#--------------Don Toxicity Scoring------------------ 

# -	No data available		
# 0	No significant toxicity		
# 1	Slight significant toxicity (corrected survivorship ≥ 75%),              1 spp		
# 2	Moderate 1 spp (50< corrected survivorship < 75) or slight 2 or more spp mortality                                    		
# 3	Severe 1 ssp (corrected survivorship ≤ 50) or moderate 2 or more spp toxicity		

Don_TOX_fun=function(x){
  
  df= x %>%
    mutate(Don_TOX_Desc=
             case_when(
               SIGNFICANTLY_DIFFERENT_Y_OR_N=="N" & BIOLOGICALLY_SIGNIFICANT_80_PERCENT_Y_OR_N=="N" ~ 0,
               SIGNFICANTLY_DIFFERENT_Y_OR_N | BIOLOGICALLY_SIGNIFICANT_80_PERCENT_Y_OR_N %in% "Y" & CONTROL_CORRECTED_SURVIVAL > 75 ~ 1,
               SIGNFICANTLY_DIFFERENT_Y_OR_N | BIOLOGICALLY_SIGNIFICANT_80_PERCENT_Y_OR_N %in% "Y" & CONTROL_CORRECTED_SURVIVAL <=75 & CONTROL_CORRECTED_SURVIVAL >50 ~ 2,
               SIGNFICANTLY_DIFFERENT_Y_OR_N | BIOLOGICALLY_SIGNIFICANT_80_PERCENT_Y_OR_N %in% "Y" & CONTROL_CORRECTED_SURVIVAL <50 ~ 3,
               TRUE ~ NA_real_))
  }

               #------------NCCA Toxicity Scoring-------------------

#  (0)  	 If Ecological and Statistical significance both equal N, "Good"!			
#  (1) 	   If only one equals N, "Fair".			
#  (2) 	   If both equal Y, and control-corrected survivorship ≥ 50%, "Poor"!			
#  (3)  	 If both equal Y, and control-corrected survivorship < 50%, "Very Poor"!			
  


#SIGNFICANTLY_DIFFERENT_Y_OR_N 
#BIOLOGICALLY_SIGNIFICANT_80_PERCENT_Y_OR_N
#CONTROL_CORRECTED_SURVIVAL


NCCA_TOX_fun=function(x){
  
df= x %>%
mutate(NCCA_TOX_Desc=
case_when(
  SIGNFICANTLY_DIFFERENT_Y_OR_N=="N" & BIOLOGICALLY_SIGNIFICANT_80_PERCENT_Y_OR_N=="N" ~ "Good",
  SIGNFICANTLY_DIFFERENT_Y_OR_N=="Y" & BIOLOGICALLY_SIGNIFICANT_80_PERCENT_Y_OR_N=="N" ~ "Fair",
  SIGNFICANTLY_DIFFERENT_Y_OR_N=="N" & BIOLOGICALLY_SIGNIFICANT_80_PERCENT_Y_OR_N=="Y" ~ "Fair",
  SIGNFICANTLY_DIFFERENT_Y_OR_N=="Y" & BIOLOGICALLY_SIGNIFICANT_80_PERCENT_Y_OR_N=="Y" & CONTROL_CORRECTED_SURVIVAL >= 50 ~ "Poor",
  SIGNFICANTLY_DIFFERENT_Y_OR_N=="Y" & BIOLOGICALLY_SIGNIFICANT_80_PERCENT_Y_OR_N=="Y" & CONTROL_CORRECTED_SURVIVAL < 50 ~ "Very Poor",
  TRUE ~ NA_character_
  ),NCCA_TOX_Score=
  as.numeric(case_when(
    NCCA_TOX_Desc=="Good"~ 0,
    NCCA_TOX_Desc=="Fair"~ 1,
    NCCA_TOX_Desc=="Poor"~ 2,
    NCCA_TOX_Desc=="Very Poor"~ 3,
    TRUE ~ NA_real_)))

return(df)

}


#------------- California Toxicity categories -------------------------

# control normalized =data from station/control data * 100
# http://ftp.sccwrp.org/pub/download/DOCUMENTS/TechnicalReports/777_CASQO_TechnicalManual.pdf
# https://oehha.ca.gov/media/downloads/ecotoxicology/general-info/marinetox3.pdf
# Compiled from U.S. EPA 1994b, ASTM 2000e
# http://ftp.sccwrp.org/pub/download/TOOLS/SQO/MLOE_Assessment_categories_table.pdf
# Welch’s t-test
#Nontoxic,
#Low, 
#Moderate, 
#High
Cal_TOX_fun=function(x){
  
df= x %>%
mutate(Cal_TOX_Desc=case_when(
  MEAN_SURVIVAL_PERCENT >80 ~ "Nontoxic",
  CONTROL_CORRECTED_SURVIVAL >75 & SIGNFICANTLY_DIFFERENT_Y_OR_N =="N" ~"Nontoxic",
  CONTROL_CORRECTED_SURVIVAL >75 & SIGNFICANTLY_DIFFERENT_Y_OR_N =="Y" ~"Low",
  CONTROL_CORRECTED_SURVIVAL <75 & CONTROL_CORRECTED_SURVIVAL >50 & SIGNFICANTLY_DIFFERENT_Y_OR_N =="Y" ~"Moderate",
  CONTROL_CORRECTED_SURVIVAL <50 ~"High"))

return(df)

}
  
 

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




