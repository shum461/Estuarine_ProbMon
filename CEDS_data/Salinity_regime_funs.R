


#~~~~~~~~~~~~~  DEQ Station IDs (CEDS) with NCCA or Bay Program IDs  ~~~
Station_data_2011_2020= read_csv("Station_data_2011_2020.csv") %>%
janitor::clean_names(.,"screaming_snake") 

#---------------------------------------------------------------------------------

dupe_finder_fun=function(x){  
  str_replace(x,c("-S1|-S2"),"")
}

#=======================================

not_in_fun=function(x,y){
  
  x[which(!x %in% y)]
}
#=======================================

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

# dissolved inorganic N
DIN_fun=function(x){
  
case_when(
is.numeric(x)==FALSE ~ NA_character_, 
x <=0.1 ~"Good",
x >0.1 & x <0.5 ~"Fair",
x >=0.5  ~"Poor",
          TRUE ~ NA_character_)
}

# orthophosphate
DIP_fun=function(x){
  
  case_when(
    is.numeric(x)==FALSE ~ NA_character_, 
    x <=0.01 ~"Good",
    x >0.01 & x <0.5 ~"Fair",
    x >=0.05  ~"Poor",
    TRUE ~ NA_character_)
}


Entero_fun=function(x){
case_when(
    is.numeric(x)==FALSE ~ NA_character_, 
    x <=10 ~"Good",
    x >10 & x <104 ~"Fair",
    x >=104  ~"Poor",
    TRUE ~ NA_character_)
}

Ecoli_fun=function(x){
  case_when(
    is.numeric(x)==FALSE ~ NA_character_, 
    x <=24 ~"Good",
    x >24 & x <235 ~"Fair",
    x >=235  ~"Poor",
    TRUE ~ NA_character_)
}


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

# NITRATE NITROGEN, DISSOLVED (MG/L AS N)= 00618
# NITRITE NITROGEN, DISSOLVED (MG/L AS N)= 00613
# NITROGEN, AMMONIA, TOTAL DISSOLVED (MG/L AS N)= 00608
# PHOSPHORUS, DISSOLVED ORTHOPHOSPHATE (MG/L AS P)=00671


# ENTEROCOCCI- ME-MF N0/100ML = 31649
# E. COLI - MTEC-MF N0/100ML  = 31648

# Need to make more concise & clean up with purrr

# param_names=c("Total_Nitrogen","Total_Phosphorus","Percent_Sand","TOC")  
# Parm_Codes=list(c("49571","49570"), c("49572","49567"),"8207","00687")
  
DCLS_Analytes_fun=function(x){
  

Summarized_Analytes= x %>% 
  filter(Ana_Sam_Mrs_Container_Id_Desc %in% c("S1","S2","R")) %>% # get rid of blanks
    group_by(Ana_Sam_Mrs_Container_Id_Desc,Fdt_Sta_Id,Fdt_Date_Time) %>% # one value for each replicate, station id and date_time,Unit
       summarise(Total_Nitrogen=sum(Ana_Value[Pg_Storet_Code %in% c("49571","49570")]),
              Total_Phosphorus=sum(Ana_Value[Pg_Storet_Code %in% c("49572","49567")]),
              Dis_Inorg_N=sum(Ana_Value[Pg_Storet_Code %in% c("00618","00613","00608")]),
              Dis_Inorg_P=sum(Ana_Value[Pg_Storet_Code %in% "00671"]),
              TOC=sum(Ana_Value[Pg_Storet_Code %in% "00687"]),
              Chlorophyll=sum(Ana_Value[Pg_Storet_Code %in% "32211"]),
              Enterococci=sum(Ana_Value[Pg_Storet_Code %in% "31649"]),
              Ecoli=sum(Ana_Value[Pg_Storet_Code %in% "31648"]),                    
              Percent_Sand=sum(Ana_Value[Pg_Storet_Code %in% "82007"])) %>%

  mutate(TOC_Quality=TOC_fun(TOC),
          DIN_Quality=DIN_fun(Dis_Inorg_N),
          DIP_Quality=DIP_fun(Dis_Inorg_P),
          Chla_Quality=NCCA_Chla_fun(Chlorophyll),
          Entero_Quality=Entero_fun(Enterococci),
          Ecoli_Quality=Ecoli_fun(Ecoli),
          Subtrate_Type=Percent_Sand_fun(Percent_Sand)) 
  
  
   
return(Summarized_Analytes) 
  
             
}



# =================== Reading file fun ===============================================


# Multiple excel sheets with multiple tabs. Gross


#====================== Toxicity testing data=======================


# several Station IDs don't match well
# "Sample description" is Ches Bay IDs but there are several naming errors with each year
# In some cases sites with duplicate tox tests (i.e. -S1 and -S2) do not have matching duplicate water chem samples (only "R" )
# PERCENT_SURVIVAL col is just for the "A" replicate. Use "MEAN_PERCENT_SURVIVAL" for overall test
# 2-JMS087.11	VA16-033A, VA06-0083A


Hyalella=c("TN-16-312","TN-17-238","TN-17-296","TN-18-593","TN-19-495",
           "TN-19-523","TN-20-466","TN-20-553")

Tox_files_path="C:/Users/vvt97279/Documents/RStudio_Test/ToxTesting"

#======================================================================
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

#==================================================
              
      
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
               SIGNFICANTLY_DIFFERENT_Y_OR_N=="Y" | BIOLOGICALLY_SIGNIFICANT_80_PERCENT_Y_OR_N =="Y" & CONTROL_CORRECTED_SURVIVAL > 75 ~ 1,
               SIGNFICANTLY_DIFFERENT_Y_OR_N=="Y" | BIOLOGICALLY_SIGNIFICANT_80_PERCENT_Y_OR_N =="Y" & CONTROL_CORRECTED_SURVIVAL <=75 & CONTROL_CORRECTED_SURVIVAL >50 ~ 2,
               SIGNFICANTLY_DIFFERENT_Y_OR_N=="Y" | BIOLOGICALLY_SIGNIFICANT_80_PERCENT_Y_OR_N =="Y" & CONTROL_CORRECTED_SURVIVAL <50 ~ 3,
               TRUE ~ NA_real_))
  }
>>>>>>> 0d25fbb73e93243260059170dcc0adfde969b4fa

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
  
#=====================================================================================
=======
  
  
  #SIGNFICANTLY_DIFFERENT_Y_OR_N 
  #BIOLOGICALLY_SIGNIFICANT_80_PERCENT_Y_OR_N
  #CONTROL_CORRECTED_SURVIVAL
  >>>>>>> 0d25fbb73e93243260059170dcc0adfde969b4fa


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
#==========================================


#============================= B-IBI ===================================================


#  Don's bug scores and weighted average
#  notice that EMAP and MAMBI weights are always 1. 
#  Ches bay and MAIA weights change based on basin

bug_scores_fun=function(x,Basin,CB,MAIA,EMAP,MAMBI){

  Basin=enquo(Basin)
  CB=enquo(CB)  
  MAIA=enquo(MAIA)
  EMAP=enquo(EMAP)
  MAMBI=enquo(MAMBI)
  
  x %>% 
  mutate(
  CB_WOE=
  case_when(
   is.numeric(!!CB)==FALSE ~ NA_real_, 
    !!CB >= 2.95 ~0,
    !!CB < 2.95 & !!CB >= 2.65 ~1,
    !!CB < 2.65 & !!CB >= 2.05 ~2,
    !!CB < 2.05 ~3
  ),
  MAIA_WOE=case_when(
    is.numeric(!!MAIA)==FALSE ~ NA_real_, 
    !!MAIA >= 2.95 ~0,
    !!MAIA < 2.95 & !!MAIA >= 2.65 ~1,
    !!MAIA < 2.65 & !!MAIA >= 2.05 ~2,
    !!MAIA < 2.05 ~3
  ),
  EMAP_WOE=case_when( 
    is.numeric(!!EMAP)==FALSE ~ NA_real_, 
    !!EMAP > 0.1 ~ 0,
    !!EMAP < 0.1 & !!EMAP >= -0.1 ~ 1,
    !!EMAP < -0.1 & !!EMAP >= -0.996 ~ 2,
    !!EMAP< -0.996 ~3
  ),
MAMBI_WOE=
case_when( 
is.numeric(!!MAMBI)==FALSE ~ NA_real_, 
!!MAMBI < 0.20 ~3,		
!!MAMBI >=0.20 & !!MAMBI < 0.39 ~ 2,	
!!MAMBI >=0.39 & !!MAMBI < 0.53	~1,
!!MAMBI >=0.53 ~0
),
CB_Weight=ifelse(str_detect(!!Basin,"VA Coastal Bays"),0,2),
MAIA_Weight=ifelse(str_detect(!!Basin,"VA Coastal Bays"),2,1),
EMAP_Weight=1,
MAMBI_Weight=1,
Matrix_score=sum((CB_WOE*CB_Weight),(MAIA_WOE*MAIA_Weight),
                 (EMAP_WOE*EMAP_Weight),(MAMBI_WOE*MAMBI_Weight))/ sum(CB_Weight,MAIA_Weight,EMAP_Weight,MAMBI_Weight)
)

}


#================================================================================================



#================  Sediment Chemistry Funs ===========================================================

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

#=============================================================================================
#============================Sediment Chemistry=================================================================

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

#=============================================================================================



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




