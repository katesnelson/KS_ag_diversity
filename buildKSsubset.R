library(dplyr)
library(sf)
library(zoo)

w<-readRDS("wwheat_panel.rds")
s<-readRDS("soy_panel.rds")
c<-readRDS("corn_panel.rds")

w_ks<-w %>% filter(.,STATE == "KS")
s_ks<-s %>% filter(.,STATE == "KS")
c_ks<-c %>% filter(.,STATE == "KS")

head(w_ks)

w_ks<-w_ks %>% select(.,-c(LSM_ENN_ALL))
s_ks<-s_ks %>% select(.,-LSM_ENN_ALL)
c_ks<-c_ks %>% select(.,-LSM_ENN_ALL)

saveRDS(w_ks, "wheat_KS_final.rds")
saveRDS(s_ks, "soy_KS_final.rds")
saveRDS(c_ks, "corn_KS_final.rds")

st_geometry(w_ks)<-NULL
st_geometry(s_ks)<-NULL
st_geometry(c_ks)<-NULL


write.csv(w_ks, "wheat_KS_final.csv")
write.csv(s_ks, "soy_KS_final.csv")
write.csv(c_ks, "corn_KS_final.csv")

########Add new data with control from Becatien 4/21/20

w2<-read.csv("Data/wheat1.csv", stringsAsFactors = F )
w2<-w2 %>% mutate(GEOID=ï..GEOID) %>% mutate(Cnty.yr=paste0(GEOID,"_",Yr)) %>% select(Cnty.yr, avg_wheat, range_wheat, volume)
w_ks<- w_ks %>% mutate(Yr=YEAR-2008) %>% mutate(Cnty.yr=paste0(GEOID,"_",Yr))
w_new<- left_join(w_ks,w2)
saveRDS(w_new, "wheat_KS_final.rds")


s2<-read.csv("Data/soy1.csv", stringsAsFactors = F)
s2<-s2 %>% mutate(GEOID=ï..GEOID) %>% mutate(Cnty.yr=paste0(GEOID,"_",Yr)) %>% select(Cnty.yr, avg_beans, range_beans, volume)
s_ks<- s_ks %>% mutate(Yr=YEAR-2008) %>% mutate(Cnty.yr=paste0(GEOID,"_",Yr))
s_new<- left_join(s_ks,s2)
saveRDS(s_new, "soy_KS_final.rds")


c2<-read.csv("Data/corn1.csv", stringsAsFactors = F)
c2<-c2 %>% mutate(GEOID=ï..GEOID) %>% mutate(Cnty.yr=paste0(GEOID,"_",Yr)) %>% select(Cnty.yr, avg_corn, range_corn, volume)
c_ks<- c_ks %>% mutate(Yr=YEAR-2008) %>% mutate(Cnty.yr=paste0(GEOID,"_",Yr))
c_new<- left_join(c_ks,c2)
saveRDS(c_new, "corn_KS_final.rds")

##########Now log transform the Yield and Water Volume variables to make them normalish 

w_new <- w_new %>% mutate(YIELD=log(YIELD), volume=log(as.numeric(volume)))
s_new <- s_new %>% mutate(YIELD=log(YIELD), volume=log(as.numeric(volume)))
c_new <- c_new %>% mutate(YIELD=log(YIELD), volume=log(as.numeric(volume)))

##########Now add the lagged, average, change, and cv variables

#wheat
w_new<- w_new %>% arrange(GEOID,Yr)

w_new<- w_new %>% group_by(GEOID) %>%  mutate(dYield = lead(YIELD, n=1) - YIELD) #change in yield dep var

w_new<- w_new %>% group_by(GEOID)  %>% #coefficient of variation for dependent (yield)
  mutate(covYield = rollapply(YIELD, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(YIELD, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))

w_new<- w_new %>% group_by(GEOID)  %>% mutate(dSIDI = lead(LSM_SIDI_ALL, n=1) - LSM_SIDI_ALL, #change in diversity ind var for change models
                                             dSHDI = lead(LSM_SHDI_ALL, n=1) - LSM_SHDI_ALL, 
                                             dRICH = lead(LSM_RICH_ALL, n=1) - LSM_RICH_ALL)

w_new<- w_new %>% group_by(GEOID)   %>% mutate(aveSIDI2yr = rollapply(LSM_SIDI_ALL, width=2, FUN = mean, align= "center", fill=NA, na.rm=T), #2 yr average diversity for change models
                                              aveSHDI2yr = rollapply(LSM_SHDI_ALL, width=2, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveRICH2yr = rollapply(LSM_RICH_ALL, width=2, FUN = mean, align= "center", fill=NA, na.rm=T))

w_new<- w_new %>% group_by(GEOID)   %>% mutate(aveSIDI3yr = rollapply(LSM_SIDI_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T), #3 yr average diversity for variability models
                                              aveSHDI3yr = rollapply(LSM_SHDI_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveRICH3yr = rollapply(LSM_RICH_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))

w_new<- w_new %>% group_by(GEOID)  %>% #coefficient of variation for diversity
  mutate(covSIDI = rollapply(LSM_SIDI_ALL, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(LSM_SIDI_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
         covSHDI = rollapply(LSM_SHDI_ALL, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(LSM_SHDI_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
         covRICH = rollapply(LSM_RICH_ALL, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(LSM_RICH_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))


w_new<- w_new %>% group_by(GEOID)  %>% mutate(lagSIDI = lag(LSM_SIDI_ALL, n=1), #lagged diversity vars
                                             lagSHDI= lag(LSM_SHDI_ALL, n=1), 
                                             lagRICH= lag(LSM_RICH_ALL, n=1)) 

w_new<- w_new %>% group_by(GEOID) %>% mutate(lagGDD = lag(GDD, n=1), #lagged climate vars
                                            lagSDD = lag(SDD, n=1), 
                                            lagTP = lag(TP, n=1), 
                                            lagvolume = lag(volume, n=1))

w_new<- w_new %>% group_by(GEOID)  %>% mutate(dGDD = lead(GDD, n=1) - GDD, #change in climate vars for change models
                                             dSDD = lead(SDD, n=1) - SDD, 
                                             dTP = lead(TP, n=1) - TP,
                                             dVol = lead(volume, n=1) - volume)

w_new<- w_new %>% group_by(GEOID)   %>% mutate(aveGDD2yr = rollapply(GDD, width=2, FUN = mean, align= "center", fill=NA, na.rm=T), #2 yr average climate vars for change models
                                              aveSDD2yr = rollapply(SDD, width=2, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveTP2yr = rollapply(TP, width=2, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveVol2yr = rollapply(volume, width=2, FUN = mean, align= "center", fill=NA, na.rm=T))

w_new<- w_new %>% group_by(GEOID)   %>% mutate(aveGDD3yr = rollapply(GDD, width=3, FUN = mean, align= "center", fill=NA, na.rm=T), #3 yr average climate vars for change models
                                              aveSDD3yr = rollapply(SDD, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveTP3yr = rollapply(TP, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveVol3yr = rollapply(volume, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))


w_new<- w_new %>% group_by(GEOID)  %>% #coefficient of variation for climate vars
  mutate(covGDD = rollapply(GDD, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(GDD, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
         covSDD = rollapply(SDD, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(SDD, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
         covTP = rollapply(TP, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(TP, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
         covVol = rollapply(volume, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(volume, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))

#soy
s_new<- s_new %>% arrange(GEOID,Yr)

s_new<- s_new %>% group_by(GEOID) %>%  mutate(dYield = lead(YIELD, n=1) - YIELD) #change in yield dep var

s_new<- s_new %>% group_by(GEOID)  %>% #coefficient of variation for dependent (yield)
  mutate(covYield = rollapply(YIELD, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(YIELD, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))

s_new<- s_new %>% group_by(GEOID)  %>% mutate(dSIDI = lead(LSM_SIDI_ALL, n=1) - LSM_SIDI_ALL, #change in diversity ind var for change models
                                             dSHDI = lead(LSM_SHDI_ALL, n=1) - LSM_SHDI_ALL, 
                                             dRICH = lead(LSM_RICH_ALL, n=1) - LSM_RICH_ALL)

s_new<- s_new %>% group_by(GEOID)   %>% mutate(aveSIDI2yr = rollapply(LSM_SIDI_ALL, width=2, FUN = mean, align= "center", fill=NA, na.rm=T), #2 yr average diversity for change models
                                              aveSHDI2yr = rollapply(LSM_SHDI_ALL, width=2, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveRICH2yr = rollapply(LSM_RICH_ALL, width=2, FUN = mean, align= "center", fill=NA, na.rm=T))

s_new<- s_new %>% group_by(GEOID)   %>% mutate(aveSIDI3yr = rollapply(LSM_SIDI_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T), #3 yr average diversity for variability models
                                              aveSHDI3yr = rollapply(LSM_SHDI_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveRICH3yr = rollapply(LSM_RICH_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))

s_new<- s_new %>% group_by(GEOID)  %>% #coefficient of variation for diversity
  mutate(covSIDI = rollapply(LSM_SIDI_ALL, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(LSM_SIDI_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
         covSHDI = rollapply(LSM_SHDI_ALL, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(LSM_SHDI_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
         covRICH = rollapply(LSM_RICH_ALL, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(LSM_RICH_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))


s_new<- s_new %>% group_by(GEOID)  %>% mutate(lagSIDI = lag(LSM_SIDI_ALL, n=1), #lagged diversity vars
                                             lagSHDI= lag(LSM_SHDI_ALL, n=1), 
                                             lagRICH= lag(LSM_RICH_ALL, n=1)) 

s_new<- s_new %>% group_by(GEOID) %>% mutate(lagGDD = lag(GDD, n=1), #lagged climate vars
                                            lagSDD = lag(SDD, n=1), 
                                            lagTP = lag(TP, n=1), 
                                            lagvolume = lag(volume, n=1))

s_new<- s_new %>% group_by(GEOID)  %>% mutate(dGDD = lead(GDD, n=1) - GDD, #change in climate vars for change models
                                             dSDD = lead(SDD, n=1) - SDD, 
                                             dTP = lead(TP, n=1) - TP,
                                             dVol = lead(volume, n=1) - volume)

s_new<- s_new %>% group_by(GEOID)   %>% mutate(aveGDD2yr = rollapply(GDD, width=2, FUN = mean, align= "center", fill=NA, na.rm=T), #2 yr average climate vars for change models
                                              aveSDD2yr = rollapply(SDD, width=2, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveTP2yr = rollapply(TP, width=2, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveVol2yr = rollapply(volume, width=2, FUN = mean, align= "center", fill=NA, na.rm=T))

s_new<- s_new %>% group_by(GEOID)   %>% mutate(aveGDD3yr = rollapply(GDD, width=3, FUN = mean, align= "center", fill=NA, na.rm=T), #3 yr average climate vars for change models
                                              aveSDD3yr = rollapply(SDD, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveTP3yr = rollapply(TP, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveVol3yr = rollapply(volume, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))


s_new<- s_new %>% group_by(GEOID)  %>% #coefficient of variation for climate vars
  mutate(covGDD = rollapply(GDD, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(GDD, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
         covSDD = rollapply(SDD, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(SDD, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
         covTP = rollapply(TP, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(TP, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
         covVol = rollapply(volume, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(volume, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))

#corn
c_new<- c_new %>% arrange(GEOID,Yr)

c_new<- c_new %>% group_by(GEOID) %>%  mutate(dYield = lead(YIELD, n=1) - YIELD) #change in yield dep var

c_new<- c_new %>% group_by(GEOID)  %>% #coefficient of variation for dependent (yield)
  mutate(covYield = rollapply(YIELD, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(YIELD, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))

c_new<- c_new %>% group_by(GEOID)  %>% mutate(dSIDI = lead(LSM_SIDI_ALL, n=1) - LSM_SIDI_ALL, #change in diversity ind var for change models
                                             dSHDI = lead(LSM_SHDI_ALL, n=1) - LSM_SHDI_ALL, 
                                             dRICH = lead(LSM_RICH_ALL, n=1) - LSM_RICH_ALL)

c_new<- c_new %>% group_by(GEOID)   %>% mutate(aveSIDI2yr = rollapply(LSM_SIDI_ALL, width=2, FUN = mean, align= "center", fill=NA, na.rm=T), #2 yr average diversity for change models
                                              aveSHDI2yr = rollapply(LSM_SHDI_ALL, width=2, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveRICH2yr = rollapply(LSM_RICH_ALL, width=2, FUN = mean, align= "center", fill=NA, na.rm=T))

c_new<- c_new %>% group_by(GEOID)   %>% mutate(aveSIDI3yr = rollapply(LSM_SIDI_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T), #3 yr average diversity for variability models
                                              aveSHDI3yr = rollapply(LSM_SHDI_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveRICH3yr = rollapply(LSM_RICH_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))

c_new<- c_new %>% group_by(GEOID)  %>% #coefficient of variation for diversity
  mutate(covSIDI = rollapply(LSM_SIDI_ALL, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(LSM_SIDI_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
         covSHDI = rollapply(LSM_SHDI_ALL, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(LSM_SHDI_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
         covRICH = rollapply(LSM_RICH_ALL, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(LSM_RICH_ALL, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))


c_new<- c_new %>% group_by(GEOID)  %>% mutate(lagSIDI = lag(LSM_SIDI_ALL, n=1), #lagged diversity vars
                                             lagSHDI= lag(LSM_SHDI_ALL, n=1), 
                                             lagRICH= lag(LSM_RICH_ALL, n=1)) 

c_new<- c_new %>% group_by(GEOID) %>% mutate(lagGDD = lag(GDD, n=1), #lagged climate vars
                                            lagSDD = lag(SDD, n=1), 
                                            lagTP = lag(TP, n=1), 
                                            lagvolume = lag(volume, n=1))

c_new<- c_new %>% group_by(GEOID)  %>% mutate(dGDD = lead(GDD, n=1) - GDD, #change in climate vars for change models
                                             dSDD = lead(SDD, n=1) - SDD, 
                                             dTP = lead(TP, n=1) - TP,
                                             dVol = lead(volume, n=1) - volume)

c_new<- c_new %>% group_by(GEOID)   %>% mutate(aveGDD2yr = rollapply(GDD, width=2, FUN = mean, align= "center", fill=NA, na.rm=T), #2 yr average climate vars for change models
                                              aveSDD2yr = rollapply(SDD, width=2, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveTP2yr = rollapply(TP, width=2, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveVol2yr = rollapply(volume, width=2, FUN = mean, align= "center", fill=NA, na.rm=T))

c_new<- c_new %>% group_by(GEOID)   %>% mutate(aveGDD3yr = rollapply(GDD, width=3, FUN = mean, align= "center", fill=NA, na.rm=T), #3 yr average climate vars for change models
                                              aveSDD3yr = rollapply(SDD, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveTP3yr = rollapply(TP, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
                                              aveVol3yr = rollapply(volume, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))


c_new<- c_new %>% group_by(GEOID)  %>% #coefficient of variation for climate vars
  mutate(covGDD = rollapply(GDD, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(GDD, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
         covSDD = rollapply(SDD, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(SDD, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
         covTP = rollapply(TP, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(TP, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
         covVol = rollapply(volume, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(volume, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))


####Save

saveRDS(w_new, "wheat_KS_final.rds")
saveRDS(s_new, "soy_KS_final.rds")
saveRDS(c_new, "corn_KS_final.rds")

write.csv(w_new, "wheat_KS_final.csv")
write.csv(s_new, "soy_KS_final.csv")
write.csv(c_new, "corn_KS_final.csv")