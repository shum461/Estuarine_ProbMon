


source("Sediment_Chemistry/file_load.R")
options(scipen = 999)

#PRESS_Sites <- read_csv("PRESS_Sites.csv")
#pin(stations_summary,"EstProbMon_Stations_Summary",board="rsconnect")


#===================== To DO ===============================================
# include fish tissue data from NCCA sites
# Inquire about sampling all parameters at a true reference site
# look into adjusting total DDT based on TOC

# Chemicals without ERM, PEC values should be compared to 95th percentile values

# Elevated chemistry was defined as 6 or more chemicals exceeding ERM guidelines, a
# mean ERMQ above 0.5, or one or more chemicals at concentrations high enough to likely be
# associated with biological effects


#================= Adjust metals based on PART size=========================
#  decreasing grain size and increasing metal concentrations
#  An example of normalizing a bulk sediment concentration for a metal to the fine fraction for
#  a sample with 84 mg/kg of lead and 60% fines (40% silt + 20% clay) is 84 mg Pb/kg ÷ 0.60 kg
#  fines /kg sediment = 140 mg lead / kg of fines.
#

#================ Adjust PAH, PCB, Pesticides based on TOC percent==========
# expressed on an assumed dry weight normalized basis at 1% organic carbon
# To convert the study site TPAH concentration to a dry weight concentration normalized to 1%,
# divide the 7,300 ug/kg value by 5 (5% TOC content) = 1,460 ug TPAH/kg at 1% TOC. On the
# common basis of 1% TOC, the study site TPAH concentration is less than the TEC
# concentration (1,460 ug/kg study site vs. 1,610 ug/kg TEC)


#=====================load pins from rsconnect===============================
Sed_Chem=pin_get("EstProbMon_Sed_Chem_2015_2019",board="rsconnect")
DCLS_params=pin_get("EstProbMon_DCLS_Params",board="rsconnect")

stations_summary=pin_get("EstProbMon_Stations_Summary",board="rsconnect")

#====== Clean names, join chemicals to lookup table screening values =========
Sed_Chem_all=Sed_Chem %>% 
mutate(Ana_Sam_Mrs_Container_Id_Desc=gsub("V","S",Ana_Sam_Mrs_Container_Id_Desc))%>%
mutate(RESULT=ifelse(RESULT<MDL,0,RESULT))%>%
mutate(CBP_NAME=case_when(
    str_detect(CBP_NAME,"VA15-0023A")~ "VA15-0023D",
    str_detect(CBP_NAME,"VA15-0028A")~ "VA15-0028B",
    str_detect(CBP_NAME,"VA16-034A") ~"VA16-034D",
    str_detect(CBP_NAME,"VA17-0021")~"VA17-0021A",
     str_detect(CBP_NAME,"VA15-0023A")~"VA15-0023C",
  str_detect(CBP_NAME,"VA15-0028A")~"VA15-0028B",
  str_detect(CBP_NAME,"VA16-034A")~"VA16-034D",
TRUE~ CBP_NAME))%>%
mutate_metals_fun(CAS=CAS_NO,Result=RESULT,Analyte=PARAMETER,UNIT)%>%
left_join(PAHs_CAS,by=c("PARAMETER"="PAH"))%>%
Full(.)%>%
left_join(DCLS_params %>% select(CBP_NAME,Salinity_regime),by=c("CBP_NAME"))%>%
mutate(Salinity_regime=
         case_when(
  str_detect(CBP_NAME,"PRESS")~"TF",
  str_detect(CBP_NAME,"VA15-0023D")~"HM",
  TRUE~Salinity_regime))%>%
mutate_at(vars(ERM_Q,ERM_PEC_min_Q,PEC_Q,PEL_Q,ERL_Q),list(Over=over_the_limit))#%>%
#group_by(CBP_NAME,Ana_Sam_Mrs_Container_Id_Desc,YEAR)%>%
#nest()

d=Sed_Chem_all %>% mutate(Global_SMH=map(data, ~.x %>% 
                                summarise_at(vars(ERM_Q,ERM_PEC_min_Q,PEC_Q,PEL_Q,ERL_Q),max)
                                   list(~mean(.,na.rm=T), ~max(.,na.rm=T)))})

                          
                          summarise_all(funs(mean, median                          
                          
                          %>%
mutate(LRM=exp(B0+B1*log10(Result_SMH))%%(1+exp(B0+B1*log10(Result_SMH))),
       LRM_Prob=ifelse(LRM==0,0,0.11+(0.33*LRM)+(0.4*LRM^2)))

#gapminder_nested %>% 
#  mutate(avg_lifeExp = map_dbl(data, ~{mean(.x$lifeExp)}))

Sed_Chem_LRM=Sed_Chem_all %>%
  filter(CBP_NAME=="VA19-0030A",PARAMETER=="Arsenic")%>%
  select(LRM,LRM_Prob)
  mutate(LRM=exp(B0+B1)*log10(Result_SMH))
  
  filter(YEAR %in% 2019) %>%
  group_by(CBP_NAME,Ana_Sam_Mrs_Container_Id_Desc,YEAR) %>%
  nest()%>%
  mutate(LRM_max=map_dbl(data,~.x %>% pluck("LRM") %>% max(.,na.rm=T)),
         LRM_Prob_max=map_dbl(data,~.x %>% pluck("LRM_Prob") %>% max(.,na.rm=T)))



LRM_fun=function(x,Result){

exp(x$B0+x$B1)*log10(x$Result)%%(1+exp(x$B0+x$B1)*log10(x$Result))

}




Result_SMH==0 ~ NA
LRM_25<MDL~NA
case_when(
is.na(B0)==T | is.na(B1)==T) ~NA,

#=============Evaluating quotients======================================

# Certain chemicals ar dropped during Hylands averaging of ERMs hence the lookup tables
# Dons Global averages dropped slightly less chems
# SMH Global is all chems

Global_param_drops=c("Dieldrin","Heptachlor epoxide","gamma-BHC")
Hyland_values = read_csv("App_lookup_Tables/Hyland_values.csv")
Global_values = read_csv("App_lookup_Tables/Global_values.csv")

over_the_limit=function(x){
case_when(x>=1~"Exceedance",TRUE~NA_character_)}

#--------------------------------------------------------------
# 3 different ways of collecting mean quotients
# Hyland is the one used in the sed chem WOE score 


 quos(PARAMETER %in% Hyland_values$Hyland)))
       
Sed_Chem_all %>%
pmap(list(data,Type=c("SMH","Don"),
  filters=c(quos(!PARAMETER %in% Global_param_drops),
  quos(PARAMETER %in% Global_param_drops))),
    ~{data %>% filter(.,!!!filters) %>%
    mutate(.,Quotient_Type=Type)})




Global_SMH=Sed_Chem_all %>%  
group_by(CBP_NAME,Ana_Sam_Mrs_Container_Id_Desc,YEAR)%>%
summarise_at(vars(ERM_Q,ERM_PEC_min_Q,PEC_Q,PEL_Q,ERL_Q), list(~mean(.,na.rm=T), ~max(.,na.rm=T)))%>%
mutate(Type="SMH")

Global_Don=Sed_Chem_all %>%  
group_by(CBP_NAME,Ana_Sam_Mrs_Container_Id_Desc,YEAR)%>%
filter(!PARAMETER %in% Global_param_drops)%>%
summarise_at(vars(ERM_Q,ERM_PEC_min_Q,PEC_Q,PEL_Q,ERL_Q), list(~mean(.,na.rm=T),~max(.,na.rm=T)))%>%
mutate(Type="Don")
  
Hyland=Sed_Chem_all %>%  
group_by(CBP_NAME,Ana_Sam_Mrs_Container_Id_Desc,YEAR)%>%
filter(PARAMETER %in% Hyland_values$Hyland)%>%
summarise_at(vars(ERM_Q,ERM_PEC_min_Q,PEC_Q,PEL_Q,ERL_Q), list(~mean(.,na.rm=T),~max(.,na.rm=T)))%>%
mutate(Type="Hyland")
#----------------------------------------------------------------------


DCLS_params_f=DCLS_params  %>% filter(Year %in% 2015:2020)%>%
select(CBP_NAME,Year,Salinity_regime)#%>%
  #mutate(Ana_Sam_Mrs_Container_Id_Desc=gsub("V","S",Ana_Sam_Mrs_Container_Id_Desc))

#-----------------------------------------------------------------------

Quotients=bind_rows(Hyland,Global_Don,Global_SMH)%>%
ungroup()%>%
mutate(Year=as.numeric(YEAR))%>%
left_join(DCLS_params_f,by=c("CBP_NAME","Year"))%>%
mutate(Salinity_regime=ifelse(str_detect(CBP_NAME,"PRESS"),"TF",Salinity_regime))%>%
Hyland_Scores_fun(ERM_Q_mean) %>%
Hyland_Scores_fun(PEC_Q_mean) %>%
Hyland_Scores_fun(ERM_PEC_min_Q_mean)%>%
WOE_Sed_Chem_fun(Salinity_regime)
  

Hyland_means=Quotients %>% 
  filter(Type %in% "Hyland")  
 
#========================================================

# disgusting -but group_map and other attempts just aren't working

#cols_to_map=c("ERM_Q","ERM_PEC_min_Q","PEC_Q","PEL_Q","ERL_Q")
#Sed_Chem_all %>% 
#group_by(CBP_NAME,Ana_Sam_Mrs_Container_Id_Desc,YEAR)%>%
#nest()%>%
#mutate(Tops=map(data[cols_to_map],~mean(.,na.rm=T)))

Sed_Chem_all %>%
group_by(CBP_NAME,Ana_Sam_Mrs_Container_Id_Desc,YEAR) %>%
nest()%>%

 
function(df,col) {
df %>%
  mutate(Top5=map(data, ~.x %>% top_n(n = 5, wt = .data$ERM_Q)%>%
  select(.data$PARAMETER,.data$ERM_Q))) 

}

ERM5=Sed_Chem_all %>% 
Top_5_fun(.data$ERM_Q)

PEC5=Sed_Chem_all %>% 
Top_5_fun(PEC_Q)

ERM_PEC_min5=Sed_Chem_all %>% 
Top_5_fun(ERM_PEC_min_Q)

PEL5=Sed_Chem_all %>% 
Top_5_fun(PEL_Q)

ERL5=Sed_Chem_all %>% 
Top_5_fun(ERL_Q)


# combine into one df  
Top5s=reduce(list(ERM5,PEC5,ERM_PEC_min5,PEL5,ERL5),left_join)%>%
mutate(Year=as.numeric(YEAR))%>%
select(-YEAR)%>%
left_join(DCLS_params_f,by=c("CBP_NAME","Year"))

# --------------shiny integration---------- pull out cols by site
# Top5s %>% 
# filter(CBP_NAME=="VA16-004A") %>% 
# unnest()

# pin(Top5s,"EstProbMon_Sed_Chem_Top5",board="rsconnect")

#===============================================================================================
#=====================Functions=================================================================

# Which quotient to use depends on the bottom salinity of the site
# Tidal fresh uses PEC, Transitional i.e. oligohaline uses min benchmark b/w ERM and PEC,
# regular estuarine/marine use ERM benchmarks

WOE_Sed_Chem_fun=function(x,sal){

  sal=enquo(sal)

x %>%
  mutate(WOE_Score=case_when(
  !!sal=="TF"~ PEC_Q_mean_Hyland_score,
  !!sal=="OH"~ ERM_PEC_min_Q_mean_Hyland_score,
  TRUE~ERM_Q_mean_Hyland_score),
WOE_Desc=case_when(
    !!sal=="TF"~ PEC_Q_mean_Hyland_desc,
    !!sal=="OH"~ ERM_PEC_min_Q_mean_Hyland_desc,
    TRUE~ERM_Q_mean_Hyland_desc)
  
)}

#============= Top 5 quotients and their parameters==============

# pull top 5 quotients of selected benchmark quotients (e.g. ERM_Q PEC_Q etc.)
# will rename 
Top_5_fun= function(df,benchmark){
  
  benchmark=ensym(benchmark)
  
  df %>%
    group_by(CBP_NAME,Ana_Sam_Mrs_Container_Id_Desc,YEAR) %>% 
    top_n(n = 5, wt = !!benchmark) %>% 
    select(!!benchmark,PARAMETER) %>% 
    nest(data=c(!!benchmark,PARAMETER))%>%
    mutate(data2=map(data,~.x %>% rename(!!str_c(as_label(ensym(benchmark)),"_Top_5_Params"):=PARAMETER)))%>%
    rename(!!str_c(as_label(ensym(benchmark)),"_Top_5"):=data2)%>%
    select(-data)
  
}
 
# ------------Test that all PAH names in lookup table join to dataset by name --------------------

#  Sed_Chem %>% 
#  mutate_metals_fun(CAS=CAS_NO,Result=RESULT,Analyte=PARAMETER,Units=UNIT)%>% 
#  filter(Class %in% "PAH")%>%
#  anti_join(PAHs_CAS,by=c("PARAMETER"="PAH"))
#-------------------------------------------------------------------------------------------------

# ==================== Sediment Chemistry =========================

# Adds class of METAL,PAH,PCB based on lookup table CAS number or parameter Names
# Adds ppm,ppb or ppt units descriptions 
# converts (ug/) ppb to (mg/kg) ppm for metals and mg/kg to g/kg for TOC
mutate_metals_fun=function(df,CAS,Result,Analyte,Units){
  
  CAS=enquo(CAS) 
  Result=enquo(Result)
  Units=enquo(Units)
  Analyte=enquo(Analyte)
  
  df %>%
    mutate_if(is.character, ~gsub('[^ -~]', '', .))%>%
    mutate_if(is.character,~gsub("<U\\+00B4>", "", .))%>%
    #mutate_at(vars(contains('Date')|contains('date')),~lubridate::mdy(.))%>%
    mutate(Class=case_when(
      !!Analyte %in% metals_CAS$NAME~ "METAL",
      str_detect(!!Analyte,"DD")==T ~ "PESTICIDE",
      !!CAS %in% PCB_Congener_CAS$CAS  & str_detect(!!Analyte,"Surr")==F | !!Analyte %in% "2,3,4,4',5-PeCB1,2"~ "PCB",
      !!Analyte %in% Pesticides_CAS$Name | !!Analyte %in% "Hexachlorocyclohexane" ~ "PESTICIDE",
      !!Analyte %in% PAHs_CAS$PAH ~ "PAH",
      TRUE ~ NA_character_))%>%
     mutate(
       Result_SMH=case_when(Class=="METAL" & !!Units %in% c("ug/Kg-dry","µg/Kg-dry") ~ !!Result/1000,
                           .$PARAMETER=="Organic Carbon, Total" & !!Units %in% "mg/Kg-dry" ~ !!Result/1000,
                           TRUE ~ !!Result),
      Units_SMH=case_when(Class %in% "METAL" ~ "mg/Kg-dry",
                          #!!Units %in% "ug/Kg-dry" ~ "µg/Kg-dry",
                          PARAMETER=="Organic Carbon, Total" & !!Units %in% "mg/Kg-dry"~"g/Kg-dry",
                          TRUE ~ !!Units),
      Units_desc=case_when(
        Units_SMH %in% c("ug/Kg-dry","µg/Kg-dry","ng/dry g")~"ppb",
        Units_SMH %in% "mg/Kg-dry"~"ppm",
        Units_SMH %in% "g/Kg-dry"~"ppt")
        
  )
}

#====================================================================
Full= function(x){
  
  #Sed_Chem_Data = x %>% 
   # mutate_if(is.character, ~gsub('[^ -~]', '', .))%>%
   # mutate_if(is.character,~gsub("<U\\+00B4>", "", .))%>%
   # mutate_metals_fun(CAS,Result,Analyte,Units) %>%
   # left_join(PAHs_CAS,by=c("Analyte"="PAH","CAS")) 
  
Totals=x %>%
    group_by(CBP_NAME,Ana_Sam_Mrs_Container_Id_Desc,YEAR)%>%
    summarise(Total_LMW_PAHs = sum(Result_SMH[which(Class=="PAH" & TYPE=="LMW")]),
              Total_HMW_PAHs=sum(Result_SMH[which(Class=="PAH" & TYPE=="HMW")]),
              Total_PAHs=Total_LMW_PAHs+Total_HMW_PAHs,
              Total_PCBs=sum(Result_SMH[which(Class=="PCB")]),
              Total_Chlorodane= sum(Result_SMH[which(CAS %in% Total_Chlorodanes)]),
              Total_DDT= sum(Result_SMH[which(CAS %in% Total_DDTs)]))%>%
    gather(PARAMETER,Result_SMH, Total_LMW_PAHs:Total_DDT) %>%
    mutate(Class="Totals",Units_SMH="µg/Kg-dry")
  
  Modified_Sed_Chem_Data = x %>%
    bind_rows(.,Totals) %>%
    left_join(ERM_PEC_PEL,by=c("PARAMETER"="ANALYTE"))%>%
    mutate(ERM_Q=Result_SMH/ERM,PEC_Q=Result_SMH/PEC,
           ERM_PEC_min_Q=Result_SMH/ERM_PEC_MIN,
           PEL_Q=Result_SMH/PEL,ERL_Q=Result_SMH/ERL)
  
  return(Modified_Sed_Chem_Data)
}          


#===========================================================


# deprecated- use summarise_at()

multiple_summarise_fun=function(x,summary_vars){
summary_vars =as_label(enquo(summary_vars))
x %>%
mutate(!!str_c(summary_vars,"_mean"):= map_dbl(data,~.x %>%
                                     pull(!!summary_vars) %>%
                                     mean(.,na.rm = T))) 

  
}



#=====================================================
#===================Quotient scores===================
  
#  Risk of benthic impact 	Mean ERM-Q	   Site Score Mean ERM Quotient
#  Low 	        ≤ 0.022	              0
#  Medium 	     > 0.022 - 0.098	    1
#  High 	      > 0.098 - 0.473	      2
#  Very High 	  > 0.473	              3
 
   
Hyland_Scores_fun=function(df,Hyland_scores){

x=ensym(Hyland_scores)


  df %>%   
mutate(!!str_c(as_label(x),"_Hyland_score"):=case_when(
  !!x<= 0.022 ~ 0,
  !!x > 0.022 & !!x <=0.098 ~ 1,
  !!x > 0.098 & !!x <=0.473 ~ 2,
  !!x > 0.473 ~ 3
),
!!str_c(as_label(x),"_Hyland_desc"):= case_when(
  !!x<= 0.022 ~"Good",
  !!x > 0.022 & !!x <=0.098 ~"Fair",
  !!x > 0.098 & !!x <=0.473 ~"High",
  !!x > 0.473 ~ "Very_High"
  ))
}
  
  
  
  
#============= Hyland and global Quotients=====================

# means >0.1 indicate 75% chance of observing toxicity           
#Hylands=c("Nickel","Dieldrin","4,4-DDD","4,4-DDT","Total_PAHs","Total_LMW_PAHs","Total_HMW_PAHs")  
Hyland_values = read_csv("App_lookup_Tables/Hyland_values.csv")
Global_values = read_csv("App_lookup_Tables/Global_values.csv")

Hyland_fun=function(x){ 
  
  x %>% 
  #filter(str_detect(ClientSampID,"PRESS"))%>%
  #filter(ClientSampID %in% "PRESS-10") %>%
  group_by(CBP_NAME,Ana_Sam_Mrs_Container_Id_Desc,YEAR) %>% 
  summarise(ERM_Q_Globalmean=mean(ERM_Q[which(PARAMETER %in% Global_values$Global)],na.rm=TRUE),
            PEC_Q_Globalmean=mean(PEC_Q[which(PARAMETER %in% Global_values$Global)],na.rm=TRUE),
            Min_Global=mean(ERM_PEC_min_Q[which(PARAMETER %in% Global_values$Global)],na.rm=TRUE),
            PEL_Q_Globalmean=mean(PEL_Q[which(PARAMETER %in% Global_values$Global)],na.rm=TRUE),
            ERM_Hyland=mean(ERM_Q[which(PARAMETER %in% Hyland_values$Hyland)],na.rm=TRUE),
            PEC_Hyland=mean(PEC_Q[which(PARAMETER %in% Hyland_values$Hyland)],na.rm=TRUE),
            ERM_PEC_min_Hyland=mean(ERM_PEC_min_Q[which(PARAMETER %in% Hyland_values$Hyland)],na.rm=TRUE),
            PEL_Hyland=mean(PEL_Q[which(PARAMETER %in% Hyland_values$Hyland)],na.rm=TRUE)                
  )}
