


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

# mg/kg or  ug/g = ppm
# ug/kg or  ng/g = ppb

#=====================load pins from rsconnect===============================
Sed_Chem=pin_get("EstProbMon_Sed_Chem_2015_2019",board="rsconnect") 

DCLS_params=pin_get("EstProbMon_DCLS_Params",board="rsconnect")

stations_summary=pin_get("EstProbMon_Stations_Summary",board="rsconnect")

#====================== Hyland Quotients======================================
#=============Evaluating quotients======================================

# Certain chemicals ar dropped during Hylands averaging of ERMs hence the lookup tables
# Dons Global averages dropped slightly less chems
# SMH Global is all chems

Global_param_drops=c("Dieldrin","Heptachlor epoxide","gamma-BHC")
Hyland_values = read_csv("App_lookup_Tables/Hyland_values.csv")
Global_values = read_csv("App_lookup_Tables/Global_values.csv")
#============================================================================

over_the_limit=function(x){
  case_when(x>=1~"Exceedance",TRUE~NA_character_)}


#====== Clean names, join chemicals to lookup table screening values =========
Sed_Chem_all=Sed_Chem %>% 
class_fun(CAS=CAS_NO,Result=RESULT,Analyte=PARAMETER,Units = UNIT)%>%
left_join(PAHs_CAS,by=c("PARAMETER"="PAH"))%>%
Full(.)%>%
mutate_at(vars(ERM_Q,ERM_PEC_min_Q,PEC_Q,PEL_Q,ERL_Q),list(Over=over_the_limit))%>%
group_by(CBP_NAME,Ana_Sam_Mrs_Container_Id_Desc,YEAR)%>%
  nest()%>%
left_join(DCLS_params %>% select(CBP_NAME,Salinity_regime,Ana_Sam_Mrs_Container_Id_Desc,TOC),by=c("CBP_NAME","Ana_Sam_Mrs_Container_Id_Desc"))%>%
mutate(Salinity_regime=
           case_when(
             str_detect(CBP_NAME,"PRESS")~"TF",
             str_detect(CBP_NAME,"VA15-0023D")~"HM",
             TRUE~Salinity_regime))%>%
mutate(
Global_SMH=map(data,~summarise_at(.x,vars(ERM_Q,ERM_PEC_min_Q,PEC_Q,PEL_Q,ERL_Q), 
                                           list(~mean(.,na.rm=T), ~max(.,na.rm=T)))%>% 
                          pivot_longer(everything(), names_to = "mean_or_max", values_to = "value")
                          ),
Global_Don=map(data,~filter(.x, !PARAMETER %in% Global_param_drops)%>%
                          summarise_at(.,vars(ERM_Q,ERM_PEC_min_Q,PEC_Q,PEL_Q,ERL_Q), 
                                           list(~mean(.,na.rm=T), ~max(.,na.rm=T)))%>%
                          pivot_longer(everything(), names_to = "mean_or_max", values_to = "value")
                          ),
Hyland=map(data,~filter(.x, PARAMETER %in% Hyland_values$Hyland)%>%
                 summarise_at(.,vars(ERM_Q,ERM_PEC_min_Q,PEC_Q,PEL_Q,ERL_Q), 
                              list(~mean(.,na.rm=T), ~max(.,na.rm=T)))%>%
                 pivot_longer(everything(), names_to = "mean_or_max", values_to = "value")%>%
                Hyland_Scores_fun2(value)))%>% 
  mutate(Quotient_score=case_when(
  Salinity_regime %in% "TF"~ map_dbl(Hyland,~pluck(.,"Hyland_score",3)),
  Salinity_regime %in% "OH"~ map_dbl(Hyland,~pluck(.,"Hyland_score",2)),
  TRUE~map_dbl(Hyland,~pluck(.,"Hyland_score",1))
))%>% 
Quotient_score_desc_fun(Quotient_score)%>% #-------start of LRM----------
mutate(LRM=map(data,~exp(.x$B0+.x$B1*log10(.x$Result_SMH))/(1+exp(.x$B0+.x$B1*log10(.x$Result_SMH))))) %>%
mutate(LRM_max=map_dbl(LRM,~max(.,na.rm=T)))%>%
mutate(LRM_Prob=map(LRM,~ifelse(.x==0,0,0.11+(0.33*.x)+(0.4*.x^2))),
       LRM_prob_max=map_dbl(LRM_Prob,~max(.,na.rm=T)),
       LRM_prob_avg=map_dbl(LRM_Prob,~mean(.,na.rm=T)))%>%
  LRM_scores_fun(LRM_prob_max,LRM_prob_avg)%>%
  mutate(ESB=map2(data,TOC, ~filter(.x, ESB_FCV > 0) %>% #-------- Start of ESB------
  mutate(ESB_ppm=Result_SMH/1000,ESB_TOC=ESB_ppm/(.y/1000),COC=ESB_TOC/ESB_FCV) %>% 
  select(PARAMETER,ESB_ppm,ESB_TOC,COC)))%>%
  mutate(ESB_sum=map_dbl(ESB,~sum(.x$COC)),
         ESB_50th=ESB_sum*1.64)%>%
  ESB_score_fun(ESB_50th)%>%
  mutate(Sed_Chem_WOE=max(Quotient_score,LRM_score,ESB_score))

pin(Sed_Chem_all,"EstProbMon_Sed_Chem_2015_2019_WOE",board="rsconnect")

Elizabeth_19= Sed_Chem_all %>% 
  filter(CBP_NAME %in% c("VA1920-0019A","VA19-0020A","VA19-0048A"))


  #VA19044B=Sed_Chem_all %>% filter(CBP_NAME %in% "VA19-0044B") %>%
  #LRM_scores_fun(LRM_prob_max,LRM_prob_avg)
  
  #VA19044B_uunnest=VA19044B %>% unnest(c(LRM_Prob,LRM,data))
  
 



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

 
x=Elizabeth_19 %>% 
  mutate(MM=map(data,~select(.x, PARAMETER,Result_SMH,ERL,ERL_Q)%>% 
                  filter(.,ERL_Q>1)))


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
class_fun=function(df,CAS,Result,Analyte,Units){
  
  CAS=enquo(CAS) 
  Result=enquo(Result)
  Analyte=enquo(Analyte)
  Units=enquo(Units)
  
  df %>%
    #mutate(!!Analyte:= gsub('[^ -~]', '', !!Analyte))%>%
    #mutate_if(!!Analyte:=gsub("<U\\+00B4>", "",!!Analyte))%>%
    mutate(Class=case_when(
      !!Analyte %in% metals_CAS$NAME~ "METAL",
      str_detect(!!Analyte,"DD")==T ~ "PESTICIDE",
      !!CAS %in% PCB_Congener_CAS$CAS  & str_detect(!!Analyte,"Surr")==F | !!Analyte %in% "2,3,4,4',5-PeCB1,2"~ "PCB",
      !!Analyte %in% Pesticides_CAS$Name | !!Analyte %in% "Hexachlorocyclohexane" ~ "PESTICIDE",
      !!Analyte %in% PAHs_CAS$PAH ~ "PAH",
      TRUE ~ NA_character_))%>%
     mutate(
       
       
       Result_SMH=case_when(Class %in% "METAL" & !!Units %in% c("ug/Kg-dry","µg/Kg-dry") ~ !!Result/1000,
                           PARAMETER=="Organic Carbon, Total" & !!Units %in% "mg/Kg-dry" ~ !!Result/1000,
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

#====================================================
#   LRM scores fun
#========================================================= 
#values taken directly from Don WOE workbooks. 
#Pavg =IF(P5=0, "-",IF(P5<0.0298,"Good",IF(P5<0.0761,"Fair",IF(P5<=0.1423,"Poor","Very Poor"))))
#Pmax =IF(N3="-","-",IF(N3<=0.5,"Good",IF(N3>=0.75,"Poor","Fair")))

LRM_scores_fun=function(x,LRM_prob_max,LRM_prob_avg){
  x %>%
    mutate(
      prob_max_desc=case_when(
        {{LRM_prob_max}}<=0.5~"Good",
        {{LRM_prob_max}}>=0.75~"Poor",
        TRUE~"Fair"),
      prob_max_score=case_when(
        {{LRM_prob_max}}<=0.5 ~ 0,
        {{LRM_prob_max}}>=0.75 ~2,
        TRUE~1),
      prob_avg_desc=case_when(
        {{LRM_prob_avg}}<=0.0298~"Good",
        {{LRM_prob_avg}}>0.0298 & LRM_prob_avg<=0.0761~"Fair",
        {{LRM_prob_avg}}>0.0761 & LRM_prob_avg<=0.14231~"Poor",
        {{LRM_prob_avg}}>0.14231 ~"Very Poor"),
      prob_avg_score=case_when(
        {{LRM_prob_avg}}<=0.0298~0,
        {{LRM_prob_avg}}>0.0298 & LRM_prob_avg<=0.0761~1,
        {{LRM_prob_avg}}>0.0761 & LRM_prob_avg<=0.14231~2,
        {{LRM_prob_avg}}>0.14231 ~3),
      LRM_score=mean(c(prob_avg_score,prob_max_score)),
      LRM_Score_desc=case_when(
        LRM_score<= 0.5 ~ "Good",
        LRM_score > 0.5 & LRM_score <=1.5 ~ "Fair",
        LRM_score > 1.5 & LRM_score <=2.5 ~ "Poor",
        LRM_score > 2.5  ~ "Very Poor"
      ))}
#=====================================================

            # ESB scores fun

# 
#=====================================================

ESB_score_fun=function(x,ESB_50th){
  
  x %>% 
    mutate(ESB_score_desc=
             case_when(
               {{ESB_50th}}<=1~"Good",
               {{ESB_50th}}>1 & {{ESB_50th}} <=3~"Fair", 
               {{ESB_50th}}>3 & {{ESB_50th}} <=5~"Poor",
               {{ESB_50th}}>5 ~"Very Poor"),
           ESB_score=case_when(
             {{ESB_50th}}<=1~0,
             {{ESB_50th}}>1 & {{ESB_50th}} <=3~1, 
             {{ESB_50th}}>3 & {{ESB_50th}} <=5~2,
             {{ESB_50th}}>5 ~3)
    )}

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
  
Hyland_Scores_fun2=function(df,Hyland_scores){
  
  x=enquo(Hyland_scores)
  
  
  df %>%   
    mutate(Hyland_score=case_when(
      !!x<= 0.022 ~ 0,
      !!x > 0.022 & !!x <=0.098 ~ 1,
      !!x > 0.098 & !!x <=0.473 ~ 2,
      !!x > 0.473 ~ 3
    ),
    Hyland_desc= case_when(
      !!x<= 0.022 ~"Good",
      !!x > 0.022 & !!x <=0.098 ~"Fair",
      !!x > 0.098 & !!x <=0.473 ~"High",
      !!x > 0.473 ~ "Very_High"
    ))
}
  

#================================

# Quotient Score desc

Quotient_score_desc_fun=function(x,Quotient_score){
   x%>% mutate(
     Quotient_risk=case_when(
       {{Quotient_score}}<=0 ~ "Low",
       {{Quotient_score}}>0 & {{Quotient_score}}<=1 ~ "Medium",
       {{Quotient_score}}>1 & {{Quotient_score}}<=1.5 ~ "Med-High",
       {{Quotient_score}}>1.5 & {{Quotient_score}}<=2 ~ "High",
       {{Quotient_score}}>2 ~ "Very High"
     )
   )}


#======================================
  
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
