

file_path="E:/DON's HD/! ! ! ! ! ! HD Reorg/1. Estuarine ProbMon/1.3. WOE/1.3.3. WOE Annual/"

All_years=data.frame(Files=list.files("E:/DON's HD/! ! ! ! ! ! HD Reorg/1. Estuarine ProbMon/1.3. WOE/1.3.3. WOE Annual")) %>% 
  slice(2:20,22)%>%
  mutate(Year=stringr::str_extract(Files, "^.{4}"))%>%
  #mutate(Number_of_Stations=sub('.*\\(','',Files)) %>%
  group_by(Year)%>%
  mutate(WOE_Files=map2(Files,file_path,~list.files(paste0(file_path,Files),full.names = T,pattern = ".xls|.xlsx")))

#=======================================================================================================================

#=======================================================================================================================

woe_only_path="E:/WOE_only/"

All_years=data.frame(Files=list.files("E:/WOE_only/")) %>% 
 # slice(2:20,23)%>%
  mutate(Year=stringr::str_extract(Files, "^.{4}"))%>%
  #mutate(Number_of_Stations=sub('.*\\(','',Files)) %>%
  group_by(Year)%>%
  mutate(WOE_Files=map2(Files,file_path,~list.files(paste0(woe_only_path,Files),full.names = T,pattern = ".xls|.xlsx")))



All_years$WOE_Files[20:21] %>% 
flatten_chr(.)%>%
map_df(.,~set_names(.x) %>% read_excel(path = .x,sheet=6,range="H1:H1"))%>%
gather("SiteName","ESB_Score")%>%
drop_na("ESB_Score")
  
map_df(.,~set_names(.x) %>% read_excel(path = .x,sheet=2,range="C5:C5"))%>%
gather("SiteName","Site_ID")%>%
drop_na("Site_ID")

map_df(.,~set_names(.x) %>% read_excel(path = .x,sheet=2,range="C5:C5"))%>%
gather("SiteName","Site_ID")%>%
drop_na("Site_ID")


map_df(.,~set_names(.x) %>% read_excel(path = .x,sheet=4,range="D2:D2")) %>%
gather("SiteName","SED_CHEM") %>%
drop_na("SED_CHEM")

 %>%
  gather("SiteName","SED_CHEM") %>%
  drop_na("SED_CHEM")
#==================================================================================================================================================

File_function= function(path,Files,Year){
  
  file_path=path
  Folder_Path=paste0(file_path,"/",Files)
  file_list=list.files(Folder_Path,pattern = ".xls|.xlsx")
  #listing_categories=c("2A","2B","3B","5A")
  #regions=c("PRO","TRO","CO","NRO")
  
  #~~~~~~~~~~~~~~Summary~~~~~~~~~~~~~~~~~~~~~~~~  
  Summary_Range=ifelse(Year<2004,"G15:G15","G19:G19")  
  
  tox_Range= function(x,y){
    
    case_when(
      y <= 2004 ~ 5,
      str_detect(x,"QA") & y > 2012 ~ 10,
      !str_detect(x,"QA") & y > 2012~7,
      str_detect(x,"QA") & y > 2004 & y <= 2012 ~ 8,
      !str_detect(x,"QA") & y > 2004 & y <= 2012 ~ 6,
      TRUE ~ 7
    )
  }
  
  
  Summ=map2(Folder_Path,file_list,~ set_names(.y) %>% read_excel(path = paste0(.x,"/",.y),sheet=3,range=Summary_Range)) %>%
    flatten() %>%
    as_tibble() %>% 
    gather("SiteName","Summary")
  
  Category=map2(Folder_Path,file_list,~ set_names(.y) %>% read_excel(path = paste0(.x,"/",.y),sheet=3,range='H19:H19')) %>%
    flatten() %>%
    as_tibble() %>% 
    gather("SiteName","Category")
  
  Cause=map2(Folder_Path,file_list,~ set_names(.y) %>% read_excel(path = paste0(.x,"/",.y),sheet=3,range='I19:I19')) %>%
    flatten() %>%
    as_tibble() %>% 
    gather("SiteName","Potential_Cause")
  
  #~~~~~~~~~~~~Coords~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Lat=map2(Folder_Path,file_list,~ set_names(.y) %>% 
             read_excel(path = paste0(.x,"/",.y),sheet=2,range="B6:B6")) %>%
    flatten() %>%
    as_tibble() %>% 
    gather("SiteName","Latitude")
  
  Lon=map2(Folder_Path,file_list,~ set_names(.y) %>% 
             read_excel(path = paste0(.x,"/",.y),sheet=2,range="C6:C6")) %>%
    flatten() %>%
    as_tibble() %>% 
    gather("SiteName","Longitude")
  
  
  Bot_Sal=map2(Folder_Path,file_list,~ set_names(.y) %>% 
                 read_excel(path = paste0(.x,"/",.y),sheet=2,range="B13:B13")) %>%
    flatten() %>%
    as_tibble() %>% 
    gather("SiteName","Bottom_Salinity") %>%
    mutate(Bottom_Salinity=as.numeric(Bottom_Salinity))
  
  
  
  remove = c("Control-corrected survivorship =", "%")
  
  tox=pmap(list(Folder_Path,file_list,Year),function(a,b,c) set_names(b) %>% 
             read_excel(path = paste0(a,"/",b),sheet=tox_Range(b,c),range="F13:F13")) %>%
    flatten() %>%
    as_tibble() %>% 
    gather("SiteName","Tox_Survival") %>%
    mutate(Tox_Survival= str_remove_all(Tox_Survival, paste(remove, collapse = "|")))
  
  
  #~~~~~~~~~~~~~~~ Region~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Region=map2(Folder_Path,file_list,~ set_names(.y) %>% 
                read_excel(path = paste0(.x,"/",.y),sheet=2,range="I7:I7")) %>%
    flatten() %>%
    as_tibble() %>% 
    gather("SiteName","Region")
  
  Sal_reg=map2(Folder_Path,file_list,~ set_names(.y) %>% 
                read_excel(path = paste0(.x,"/",.y),sheet=2,range="C13:C13")) %>%
    flatten() %>%
    as_tibble() %>% 
    gather("SiteName","Salinity")
  
  #~~~~~~~~~~~~~~~~~BIBI~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #if QA G22, sheet 11
  
  
  
  ches_bay=pmap(list(Folder_Path,file_list,Year),function(a,b,c) set_names(b) %>%  
                  read_excel(path = paste0(a,"/",b),sheet=case_when(
                    c > 2014 ~  ifelse(str_detect(b,"QA"),11,8),
                    c <= 2014 ~ ifelse(str_detect(b,"QA"),9,7)
                  ),
                  range=case_when(
                    c > 2014 ~  ifelse(str_detect(b,"QA"),"G22:G22","G20:G20"),
                    c <= 2014 ~ ifelse(str_detect(b,"QA"),"F19:F19","F19:F19")
                  )))%>%
    flatten() %>%
    as_tibble() %>% 
    gather("SiteName","ChesBay_IBI")
  
  #~~~~~~~~~~~~~~~Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  df=data.frame(plyr::join_all(list(Region,Summ,Bot_Sal,Sal_reg,tox,Category,ches_bay,Lat,Lon), by='SiteName', type='left'))
  
  return(df)
  
}


#=======================================================================================================================================================




Files_2018=list.files("T:/2017&2019_PRESS/2019_PRESS/2019 PRESS Weight-of-Evidence") %>%
#Files_2018=list.files("E:/WOE_only") %>%  
as_tibble()%>%
  slice(4:19)%>%
  mutate(Year=as.numeric(stringr::str_extract(value, "^.{4}")))%>%
  mutate(Summaries=map2(value,Year,~File_function(.x,.y))) %>%
  mutate(Final_Summaries=map2(Summaries,Year,~mutate(.x,Year=.y)))

WOE_output_2018=map_df(Files_2018$Final_Summaries,~rbind(.x)) %>%
  mutate(Site=substr(SiteName,1,11),
         DEQ_SiteName=substr(SiteName,13,23)
  )



f=Files %>% bind_rows(Summaries,.id="Year")




  
  %>%
  mutate(Bottom_Salinity=as.numeric(Bottom_Salinity))

WOE_output_2018=map_df(Files_2018$Final_Summaries,~rbind(.x)) %>%
  mutate(Site=substr(SiteName,1,11),
         DEQ_SiteName=substr(SiteName,13,23)
  )
