
############################################################################################################
###function to prep standardized data for modeling with rw fertilizer and herbicide/insecticide controls###
##########################################################################################################
prep.data.norm_nc<-function(file, savename, crop, projection){
  
  dat<-readRDS(paste0(wd,"/Data/",file,".RDS")) #read in the data
  dat$GEOID<-as.character(dat$GEOID)
  
  cnty<-readRDS(paste0(wd,"/Data/county.RDS"))
  cnty<-st_as_sf(cnty) %>% dplyr::select(.,c(5,11)) %>% st_transform(., projection)
  
  #prep for building adjacency matrix
  d_sub<-dat[!is.na(dat$YIELD),]
  names<-unique(d_sub$GEOID)
  
  #clip county file to crop producing areas only and setup indexing by CNTY
  cnty<-cnty[cnty$GEOID %in% names,] %>% arrange(., GEOID) %>% mutate (.,CNTY=seq(1,nrow(.),1)) #order dataset by county, build group index that corresponds to adjacency matrix
  
  ##Build the relational matrix for areas of interest
  neighbors<-as(cnty,"Spatial") %>% poly2nb(., queen=F)#convert county sf to spatial polygons for inla relation matrix & create neighbors list from polygon object (neighbors share one or points at boundary)
  #temp <- poly2nb(shp, queen=T)#neighbors must share more than one point at boundary
  H.adj <- nb2mat(neighbors, style ="B", zero.policy=TRUE ) #convert to a sparse matrix to reduce memory (neighbor list to binary coded neighbor weights matrix)
  H.adj <-as(H.adj, "dgTMatrix") #convert to a sparse style matrix
  saveRDS(H.adj, paste0(wd,"/Data/H.adj.",crop,".rds"))
  
  #plot spatial object
  #image(inla.graph2matrix(H.adj), xlab="", ylab="")
  
  ###add indexes, standardize and set up vars for log-linear model with rw (by binning) and quadratic (by squaring) terms
  d_sub<-left_join(d_sub,st_set_geometry(cnty[,c("GEOID","CNTY")], NULL),by="GEOID")#join new index to crop dataset
  d_sub<-d_sub %>% arrange(.,CNTY,YEAR) #Ordering dataset by county then time for modeling as is done in the Ohio example in chapter 7 of r-inla book
  
  d_sub$YEAR.id<-d_sub$Yr  #so index starts at 1
  d_sub$YEAR.id2<-(d_sub$Yr )^2
  
  d_sub<- d_sub %>% ungroup() %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))) #remove infinites from log transforms completed earlier
  
  #Add expense data
    fert<-read.csv(paste0(wd,"/fertilizer.csv"), stringsAsFactors = F)
    fert_long<-gather(fert, key=year,value=fert, X2008:X2018, factor_key=F ) %>% 
      mutate(year=substr(year,2,5)) %>% 
      mutate(Yr=as.numeric(year)-2008) %>%
      mutate(GEOID=X) %>%
      mutate(Cnty.yr=paste0(GEOID, "_",Yr))
    
    herb<-read.csv(paste0(wd,"/herbicide.csv"), stringsAsFactors = F)
    herb_long<-gather(herb, key=year,value=herb, X2008:X2018, factor_key=F ) %>% 
      mutate(year=substr(year,2,5)) %>% 
      mutate(Yr=as.numeric(year)-2008) %>%
      mutate(GEOID=X) %>%
      mutate(Cnty.yr=paste0(GEOID, "_",Yr))
    
    d_sub <-d_sub %>% left_join(., fert_long, by="Cnty.yr") %>% left_join(., herb_long, by="Cnty.yr")
  
  #Normalize water volume, fertilizer and herbicide by agricultral area
    d_sub<- d_sub %>% mutate_at(c("fert", "herb","volume"), ~./TOTAL_ACRES)
  
  #Calculate change and cv for fert and herb and recalculate using normalized water vol
      d_sub<- d_sub %>% group_by(GEOID)  %>% mutate(dVol = lead(volume, n=1) - volume, #change in climate vars for change models
                                                  dfert = lead(fert, n=1) - fert,
                                                  dherb = lead(herb, n=1) - herb)
    
    d_sub<- d_sub %>% group_by(GEOID)   %>% mutate(aveVol2yr = rollapply(volume, width=2, FUN = mean, align= "center", fill=NA, na.rm=T), #2 yr average climate vars for change models
                                                   avefert2yr = rollapply(fert, width=2, FUN = mean, align= "center", fill=NA, na.rm=T),
                                                   aveherb2yr = rollapply(herb, width=2, FUN = mean, align= "center", fill=NA, na.rm=T))
    
    d_sub<- d_sub %>% group_by(GEOID)   %>% mutate(aveVol3yr = rollapply(volume, width=3, FUN = mean, align= "center", fill=NA, na.rm=T), #3 yr average climate vars for change models
                                                   avefert3yr = rollapply(SDD, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
                                                   aveherb3yr = rollapply(TP, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))
    
    
    d_sub<- d_sub %>% group_by(GEOID)  %>% #coefficient of variation for climate vars
      mutate(covVol = rollapply(volume, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(volume, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
             covfert = rollapply(fert, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(fert, width=3, FUN = mean, align= "center", fill=NA, na.rm=T),
             covherb = rollapply(herb, width=3, FUN = sd, align= "center", fill=NA, na.rm=T)/rollapply(herb, width=3, FUN = mean, align= "center", fill=NA, na.rm=T))
    
  
  #Set up random walk (rw1) vars for model by binning
  
  d_sub<- d_sub %>% ungroup() %>% mutate_at(c("fert", "herb","volume","TP","GDD","SDD", "SOIL", "PERC_IRR","ACRES","LSM_SIDI_ALL", "LSM_SHDI_ALL", "LSM_RICH_ALL",
                                "lagvolume", "lagTP","lagGDD","lagSDD", "lagSIDI","lagSHDI","lagRICH",
                                "dSIDI", "dSHDI", "dRICH", "dGDD", "dSDD", "dTP", "dVol","dfert","dherb",
                                "aveSIDI2yr", "aveSHDI2yr", "aveRICH2yr", "aveSIDI3yr","aveSHDI3yr","aveRICH3yr",
                                "aveGDD2yr", "aveSDD2yr", "aveTP2yr", "aveVol2yr", "aveGDD3yr","aveSDD3yr","aveTP3yr","aveVol3yr","avefert2yr","aveherb2yr","avefert3yr","aveherb3yr",
                                "covSIDI", "covSHDI", "covRICH", "covGDD","covSDD","covTP","covVol","covfert","covherb"), ~inla.group(., n = 20, method = "quantile")) #0907 not sure why this is being mean right now and throuwing error  --> looks like it needs to be ungrouped
  
  
  #Scale all vars except Yield 
  
  d_sub<- d_sub %>% mutate_at(c("fert", "herb", "volume","TP","GDD","SDD", "SOIL", "PERC_IRR","ACRES","LSM_SIDI_ALL", "LSM_SHDI_ALL", "LSM_RICH_ALL",
                                "lagvolume", "lagTP","lagGDD","lagSDD", "lagSIDI","lagSHDI","lagRICH",
                                "dSIDI", "dSHDI", "dRICH", "dGDD", "dSDD", "dTP", "dVol","dfert","dherb",
                                "aveSIDI2yr", "aveSHDI2yr", "aveRICH2yr", "aveSIDI3yr","aveSHDI3yr","aveRICH3yr",
                                "aveGDD2yr", "aveSDD2yr", "aveTP2yr", "aveVol2yr", "aveGDD3yr","aveSDD3yr","aveTP3yr","aveVol3yr","avefert2yr","aveherb2yr","avefert3yr","aveherb3yr",
                                "covSIDI", "covSHDI", "covRICH", "covGDD","covSDD","covTP","covVol","covfert","covherb"), ~scale(.))
  
  saveRDS(d_sub,paste0(wd,"/Data/",savename,"norm_nc.rds"))
  
}



##########################################################
###functions to build formulas with BYM spatial effects###
##########################################################

build.formula.bym.nlnc<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" & controls == "All"){       
    formula <- paste0('YIELD ~  1 + 
    f(fert, model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
    f(herb, model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(volume, model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(TP, model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(GDD, model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(SDD, model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(CNTY, model="bym2", graph=H.adj.',crop,', scale.model=TRUE, constr = TRUE,  
           hyper=list( phi =list(prior = "pc", param = c(phi.u, phi.alpha), inital=-3), 
                       prec =list(prior = "pc.prec", param = c(u,alpha), inital = 5))) +
      f(YEAR.id)+ f(YEAR.id2) + 
      f(',independent, ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" & controls == "All"){
      formula <- paste0('YIELD ~ 1 + 
      f(fert, model="rw1", scale.model=T) +
      f(herb, model="rw1", scale.model=T) +
      f(volume, model="rw1", scale.model=T) +
      f(TP, model="rw1",scale.model=T) + 
      f(GDD,model="rw1",scale.model=T) +
      f(SDD,model="rw1",scale.model=T) +
      f(CNTY, model="bym", graph=H.adj.', crop,', scale.model=TRUE) +
      f(YEAR.id)+ f(YEAR.id2) +
      f(', independent, ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}

build.formula.bym.change1.nlnc<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc"){       
    formula <- paste0('dYield ~  1 + 
      f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',controls[3],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[4],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[5],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[6],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(CNTY, model="bym2", graph=H.adj.',crop,', scale.model=TRUE, constr = TRUE,  
           hyper=list( phi =list(prior = "pc", param = c(phi.u, phi.alpha), inital=-3), 
                       prec =list(prior = "pc.prec", param = c(u,alpha), inital = 5))) +
      f(YEAR.id)+ f(YEAR.id2)+ 
      f(',independent, ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default"){
      formula <- paste0('dYield ~ 1 + 
      f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T) + 
      f(',controls[3],',model="rw1",scale.model=T) +
      f(',controls[4],',model="rw1",scale.model=T) +
      f(',controls[5],',model="rw1",scale.model=T) +
      f(',controls[6],',model="rw1",scale.model=T) +
      f(CNTY, model="bym", graph=H.adj.', crop,', scale.model=TRUE) +
      f(YEAR.id)+ f(YEAR.id2) +
      f(', independent, ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}
build.formula.bym.change2.nlnc<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" ){       
    formula <- paste0('dYield ~  1 +
      f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',controls[3],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[4],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[5],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[6],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(CNTY, model="bym2", graph=H.adj.',crop,', scale.model=TRUE, constr = TRUE,  
           hyper=list( phi =list(prior = "pc", param = c(phi.u, phi.alpha), inital=-3), 
                       prec =list(prior = "pc.prec", param = c(u,alpha), inital = 5))) + 
      f(YEAR.id)+ f(YEAR.id2)+ 
      f(',independent[1], ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',independent[2], ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" ){
      formula <- paste0('dYield ~ 1 + 
     f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T) + 
      f(',controls[3],',model="rw1",scale.model=T) +
      f(',controls[4],',model="rw1",scale.model=T) +
      f(',controls[5],',model="rw1",scale.model=T) +
      f(',controls[6],',model="rw1",scale.model=T) +
      f(CNTY, model="bym", graph=H.adj.', crop,', scale.model=TRUE) +
      f(YEAR.id)+ f(YEAR.id2) +
      f(', independent[1], ', model="rw1", scale.model=T)+
      f(', independent[2], ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}
build.formula.bym.cov1.nlnc<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" ){       
    formula <- paste0('covYield ~  1 +
     f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',controls[3],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[4],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
       f(',controls[5],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
        f(',controls[6],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(CNTY, model="bym2", graph=H.adj.',crop,', scale.model=TRUE, constr = TRUE,  
           hyper=list( phi =list(prior = "pc", param = c(phi.u, phi.alpha), inital=-3), 
                       prec =list(prior = "pc.prec", param = c(u,alpha), inital = 5))) +
      f(YEAR.id)+ f(YEAR.id2)+ 
      f(',independent, ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" ){
      formula <- paste0('covYield ~ 1 + 
      f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T) + 
      f(',controls[3],',model="rw1",scale.model=T) +
      f(',controls[4],',model="rw1",scale.model=T) +
      f(',controls[5],',model="rw1",scale.model=T) +
      f(',controls[6],',model="rw1",scale.model=T) +
      f(CNTY, model="bym", graph=H.adj.', crop,', scale.model=TRUE) +
      f(YEAR.id)+ f(YEAR.id2) +
      f(', independent, ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}

build.formula.bym.cov2.nlnc<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" ){       
    formula <- paste0('covYield ~  1 +
     f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',controls[3],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[4],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(',controls[5],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(',controls[6],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(CNTY, model="bym2", graph=H.adj.',crop,', scale.model=TRUE, constr = TRUE,  
           hyper=list( phi =list(prior = "pc", param = c(phi.u, phi.alpha), inital=-3), 
                       prec =list(prior = "pc.prec", param = c(u,alpha), inital = 5))) +
      f(YEAR.id)+ f(YEAR.id2)+ 
      f(',independent[1], ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',independent[2], ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" ){
      formula <- paste0('covYield ~ 1 + 
      f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T) + 
      f(',controls[3],',model="rw1",scale.model=T) +
      f(',controls[4],',model="rw1",scale.model=T) +
      f(',controls[5],',model="rw1",scale.model=T) +
      f(',controls[6],',model="rw1",scale.model=T) +
      f(CNTY, model="bym", graph=H.adj.', crop,', scale.model=TRUE) +
      f(YEAR.id)+ f(YEAR.id2) +
      f(', independent[1], ', model="rw1", scale.model=T)+
      f(', independent[2], ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}

build.formula.bym.cov3.nlnc<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" ){       
    formula <- paste0('covYield ~  1 +
     f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',controls[3],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[4],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(',controls[5],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(',controls[6],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(',controls[7],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(',controls[8],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(',controls[9],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(',controls[10],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(CNTY, model="bym2", graph=H.adj.',crop,', scale.model=TRUE, constr = TRUE,  
           hyper=list( phi =list(prior = "pc", param = c(phi.u, phi.alpha), inital=-3), 
                       prec =list(prior = "pc.prec", param = c(u,alpha), inital = 5))) +
      f(YEAR.id)+ f(YEAR.id2)+ 
      f(',independent, ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" ){
      formula <- paste0('covYield ~ 1 + 
      f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T) + 
      f(',controls[3],',model="rw1",scale.model=T) +
      f(',controls[4],',model="rw1",scale.model=T) +
      f(',controls[5],',model="rw1",scale.model=T) +
      f(',controls[6],',model="rw1",scale.model=T) +
      f(',controls[7],',model="rw1",scale.model=T) +
      f(',controls[8],',model="rw1",scale.model=T) +
      f(',controls[9],',model="rw1",scale.model=T) +
      f(',controls[10],',model="rw1",scale.model=T) +
      f(CNTY, model="bym", graph=H.adj.', crop,', scale.model=TRUE) +
      f(YEAR.id)+ f(YEAR.id2) +
      f(', independent, ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}

build.formula.bym.cov4.nlnc<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" ){       
    formula <- paste0('covYield ~  1 +
     f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',controls[3],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[4],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +  
      f(',controls[5],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[6],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[7],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[8],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[9],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[10],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[11],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[12],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(CNTY, model="bym2", graph=H.adj.',crop,', scale.model=TRUE, constr = TRUE,  
           hyper=list( phi =list(prior = "pc", param = c(phi.u, phi.alpha), inital=-3), 
                       prec =list(prior = "pc.prec", param = c(u,alpha), inital = 5))) +
      f(YEAR.id)+ f(YEAR.id2)+ 
      f(',independent[1], ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',independent[2], ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" ){
      formula <- paste0('covYield ~ 1 + 
      f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T) + 
      f(',controls[3],',model="rw1",scale.model=T) +
      f(',controls[4],',model="rw1",scale.model=T) +
      f(',controls[5],',model="rw1",scale.model=T) +
      f(',controls[6],',model="rw1",scale.model=T) +
      f(',controls[7],',model="rw1",scale.model=T) +
      f(',controls[8],',model="rw1",scale.model=T) +
      f(',controls[9],',model="rw1",scale.model=T) +
      f(',controls[10],',model="rw1",scale.model=T) +
      f(',controls[11],',model="rw1",scale.model=T) +
      f(',controls[12],',model="rw1",scale.model=T) +
      f(CNTY, model="bym", graph=H.adj.', crop,', scale.model=TRUE) +
      f(YEAR.id)+ f(YEAR.id2) +
      f(', independent[1], ', model="rw1", scale.model=T)+
      f(', independent[2], ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}


############################################################
###functions to build formulas with County Fixed Effects###
###########################################################


build.formula.fe.nlnc<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" & controls == "All"){       
    formula <- paste0('YIELD ~  1 + 
    f(fert, model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
    f(herb, model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(volume, model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(TP, model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(GDD, model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(SDD, model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      as.factor(CNTY) + 
      f(YEAR.id)+ f(YEAR.id2)+ 
      f(',independent, ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" & controls == "All"){
      formula <- paste0('YIELD ~ 1 + 
      f(fert, model="rw1", scale.model=T) +
      f(herb, model="rw1", scale.model=T) +
      f(volume, model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(TP, model="rw1",scale.model=T) + 
      f(GDD,model="rw1",scale.model=T) +
      f(SDD,model="rw1",scale.model=T) +
      as.factor(CNTY) +
      f(YEAR.id)+ f(YEAR.id2) +
      f(', independent, ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}

build.formula.fe.change1.nlnc<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc"){       
    formula <- paste0('dYield ~  1 +
      f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',controls[3],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[4],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(',controls[5],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(',controls[6],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      as.factor(CNTY) + 
      f(YEAR.id)+ f(YEAR.id2)+ 
      f(',independent, ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default"){
      formula <- paste0('dYield ~ 1 + 
      f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T) + 
      f(',controls[3],',model="rw1",scale.model=T) +
      f(',controls[4],',model="rw1",scale.model=T) +
      f(',controls[5],',model="rw1",scale.model=T) +
      f(',controls[6],',model="rw1",scale.model=T) +
      as.factor(CNTY) +
      f(YEAR.id)+ f(YEAR.id2) +
      f(', independent, ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}

build.formula.fe.change2.nlnc<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" ){       
    formula <- paste0('dYield ~  1 +
      f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',controls[3],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[4],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[5],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(',controls[6],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      as.factor(CNTY) + 
      f(YEAR.id)+ f(YEAR.id2)+ 
      f(',independent[1], ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',independent[2], ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" ){
      formula <- paste0('dYield ~ 1 + 
     f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T) + 
      f(',controls[3],',model="rw1",scale.model=T) +
      f(',controls[4],',model="rw1",scale.model=T) +
      f(',controls[5],',model="rw1",scale.model=T) +
      f(',controls[6],',model="rw1",scale.model=T) +
      as.factor(CNTY) +
      f(YEAR.id)+ f(YEAR.id2) +
      f(', independent[1], ', model="rw1", scale.model=T)+
      f(', independent[2], ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}

build.formula.fe.cov1.nlnc<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" ){       
    formula <- paste0('covYield ~  1 + 
     f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',controls[3],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[4],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +  
      f(',controls[5],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +  
      f(',controls[6],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +  
      as.factor(CNTY) + 
      f(YEAR.id)+ f(YEAR.id2)+ 
      f(',independent, ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" ){
      formula <- paste0('covYield ~ 1 + 
      f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T) + 
      f(',controls[3],',model="rw1",scale.model=T) +
      f(',controls[4],',model="rw1",scale.model=T) +
      f(',controls[5],',model="rw1",scale.model=T) +
      f(',controls[6],',model="rw1",scale.model=T) +
      as.factor(CNTY) +
      f(YEAR.id)+ f(YEAR.id2) +
      f(', independent, ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}

build.formula.fe.cov2.nlnc<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" ){       
    formula <- paste0('covYield ~  1 +
     f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',controls[3],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[4],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[5],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[6],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      as.factor(CNTY) + 
      f(YEAR.id)+ f(YEAR.id2)+ 
      f(',independent[1], ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',independent[2], ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" ){
      formula <- paste0('covYield ~ 1 + 
      f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T) + 
      f(',controls[3],',model="rw1",scale.model=T) +
      f(',controls[4],',model="rw1",scale.model=T) +
      f(',controls[5],',model="rw1",scale.model=T) +
      f(',controls[6],',model="rw1",scale.model=T) +
      as.factor(CNTY) +
      f(YEAR.id)+ f(YEAR.id2) +
      f(', independent[1], ', model="rw1", scale.model=T)+
      f(', independent[2], ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}
build.formula.fe.cov3.nlnc<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" ){       
    formula <- paste0('covYield ~  1 +
     f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',controls[3],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[4],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[5],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[6],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[7],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[8],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(',controls[9],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(',controls[10],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      as.factor(CNTY) + 
      f(YEAR.id)+ f(YEAR.id2)+ 
      f(',independent, ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" ){
      formula <- paste0('covYield ~ 1 + 
      f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T) + 
      f(',controls[3],',model="rw1",scale.model=T) +
      f(',controls[4],',model="rw1",scale.model=T) +
      f(',controls[5],',model="rw1",scale.model=T) +
      f(',controls[6],',model="rw1",scale.model=T) +
      f(',controls[7],',model="rw1",scale.model=T) +
      f(',controls[8],',model="rw1",scale.model=T) +
      f(',controls[9],',model="rw1",scale.model=T) +
      f(',controls[10],',model="rw1",scale.model=T) +
      as.factor(CNTY) +
      f(YEAR.id)+ f(YEAR.id2) +
      f(', independent, ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}
build.formula.fe.cov4.nlnc<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" ){       
    formula <- paste0('covYield ~  1 +
     f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',controls[3],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[4],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[5],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[6],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[7],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[8],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[9],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(',controls[10],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(',controls[11],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      f(',controls[12],', model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  +
      as.factor(CNTY) + 
      f(YEAR.id)+ f(YEAR.id2)+ 
      f(',independent[1], ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(',independent[2], ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" ){
      formula <- paste0('covYield ~ 1 + 
      f(',controls[1],', model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(',controls[2],', model="rw1",scale.model=T) + 
      f(',controls[3],',model="rw1",scale.model=T) +
      f(',controls[4],',model="rw1",scale.model=T) +
      f(',controls[5],',model="rw1",scale.model=T) +
      f(',controls[6],',model="rw1",scale.model=T) +
      f(',controls[7],',model="rw1",scale.model=T) +
      f(',controls[8],',model="rw1",scale.model=T) +
      f(',controls[9],',model="rw1",scale.model=T) +
      f(',controls[10],',model="rw1",scale.model=T) +
      f(',controls[11],',model="rw1",scale.model=T) +
      f(',controls[12],',model="rw1",scale.model=T) +
      as.factor(CNTY) +
      f(YEAR.id)+ f(YEAR.id2) +
      f(', independent[1], ', model="rw1", scale.model=T)+
      f(', independent[2], ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}








###########################################
###function to run basic gaussian model ###
###########################################
run.model<-function(formula, data, crop, priors ="pc", name){
  
  df<-readRDS(paste0(wd,"/Data/",data,".RDS")) #read in the data #define the dataset for the crop of interest
  h<-readRDS(paste0(wd,"/Data/H.adj.", crop, ".RDS"))
  assign(paste0("H.adj.", crop), h)
   v<-NA
   n<-NA
   Q<-NA
   u<-NA
   alpha<-NA
   phi.u<-NA
   phi.alpha<-NA
   phi.prior<-NA
  if (priors == "pc"){ #setup pc priors
    v=sd(df$YIELD,na.rm=TRUE)/0.31
    n=dim(df)[1] 
    Q = INLA:::inla.pc.bym.Q(h)
    Q = INLA:::inla.scale.model(Q, constr=list(A=matrix(1, 1, n), e=0))             
    u = 0.2/0.31
    alpha = 0.01
    phi.u = 0.5
    phi.alpha = 2/3 ## prob(phi < phi.u) = phi.alpha   
    phi.prior = INLA:::inla.pc.bym.phi(Q=Q, u= phi.u, alpha = phi.alpha)
  }
  OUT<-inla(formula, data=df, family = "gaussian",     #run the inla model
            control.predictor = list( compute=TRUE), 
            control.compute=list(dic=TRUE, cpo=TRUE), 
            control.fixed = list(prec = 0.0000001)) 
  saveRDS(OUT, file=paste0(wd,"/Output/",name ,".rds")) #save the model output
 
} 



############################################################
###function to return basic diagnostics for yield models ###
############################################################
model.results<-function(model_results_name, data="d", crop="corn", diagnostics = TRUE){
  results<-readRDS(paste0(wd,"/Output/",model_results_name,".rds"))
  data<-readRDS(paste0(wd,"/Data/",data,".RDS"))
 # st_geometry(data)<-NULL
  summary<-summary(results)
  if (diagnostics == TRUE){ #if diagnostics turned on return the dic, cpo, and pit
    #R_INLA Diagnostics
    DIC<-(results$dic$dic)
    CPO<-(sum(log(results$cpo$cpo), na.rm=T)) 
    h<-ggplot(as.data.frame(results$cpo$pit), aes(x=results$cpo$pit)) + geom_histogram() + ggtitle( "Model PIT")
    #PPC Distribution Checks
    p_val<-c()
    n<-length(data$YIELD)
    for(i in (1:n)){
      p_val[i]<-inla.pmarginal(q=data$YIELD[i],
                               marginal=results$marginals.fitted.values[[i]])
    }
    p2<-ggplot(as.data.frame(p_val), aes(x=p_val)) + geom_histogram() + ggtitle("Posterior Predictive p-values")
    
    data<-cbind(data,results$summary.fitted.values$mean)
    colnames(data)[(dim(data)[2])]<-c("FittedVals")
    p1<-ggplot(data[!is.na(data$YIELD),], aes(x=YIELD,y=FittedVals)) + geom_point(size=0.2, color="blue") +
      ggtitle("Scatter Plot of Predicted and Observed Values")  + xlab("Observed Value") + ylab("Mean of Post. Pred. Distr.")
    
    
    #PPC Summary Metrics
    sq_dif<-(data$YIELD-data$FittedVals)^2
    MSE<-1/n*(sum(sq_dif,na.rm=T))
    
    pred_res2<-(data$FittedVals[!is.na(data$YIELD)] - mean(data$YIELD, na.rm=T)) ^2
    obs_res2<-(data$YIELD[!is.na(data$YIELD)] - mean(data$YIELD, na.rm=T))^2
    R2<-sum(pred_res2, na.rm=T)/sum(obs_res2, na.rm=T)
    
    #OUPUTS
    summtable<-summary$fixed[,1:5]
    summtable.pt2<-as.data.frame(summary$hyperpar[,1:5])
    summtable<-rbind(summtable,summtable.pt2)
    summtable<-kable(summtable, caption="Summary Table of Model Estimates")
    
    Diagnostics<-as.data.frame(t(c(DIC,CPO,MSE,R2)))
    colnames(Diagnostics)<-c("DIC","CPO","MSE","R2")
    Diagnostics<-kable(Diagnostics, caption="Model Diagnostic Metrics")
    
    figure<-ggarrange(h,p1,p2,ncol=3,nrow=1)
    annotate_figure(figure, top= text_grob("Model Diagnostic Plots"))
    
    my_list<-list(figure, Diagnostics)
    return(my_list)
  }  
  my_list<-list(summary)
  return(my_list)
}

#######################################################
###function to return diagnotics for change models ###
######################################################
model.results.change<-function(model_results_name, data="d", crop="corn", diagnostics = TRUE){
  results<-readRDS(paste0(wd,"/Output/",model_results_name,".rds"))
  data<-readRDS(paste0(wd,"/Data/",data,".RDS"))
  # st_geometry(data)<-NULL
  summary<-summary(results)
  if (diagnostics == TRUE){ #if diagnostics turned on return the dic, cpo, and pit
    #R_INLA Diagnostics
    DIC<-(results$dic$dic)
    CPO<-(sum(log(results$cpo$cpo), na.rm=T)) 
    h<-ggplot(as.data.frame(results$cpo$pit), aes(x=results$cpo$pit)) + geom_histogram() + ggtitle( "Model PIT")
    #PPC Distribution Checks
    p_val<-c()
    n<-length(data$dYield)
    for(i in (1:n)){
      p_val[i]<-inla.pmarginal(q=data$dYield[i],
                               marginal=results$marginals.fitted.values[[i]])
    }
    p2<-ggplot(as.data.frame(p_val), aes(x=p_val)) + geom_histogram() + ggtitle("Posterior Predictive p-values")
    
    data<-cbind(data,results$summary.fitted.values$mean)
    colnames(data)[(dim(data)[2])]<-c("FittedVals")
    p1<-ggplot(data[!is.na(data$YIELD),], aes(x=dYield,y=FittedVals)) + geom_point(size=0.2, color="blue") +
      ggtitle("Scatter Plot of Predicted and Observed Values")  + xlab("Observed Value") + ylab("Mean of Post. Pred. Distr.")
    
    
    #PPC Summary Metrics
    sq_dif<-(data$dYield-data$FittedVals)^2
    MSE<-1/n*(sum(sq_dif,na.rm=T))
    
    pred_res2<-(data$FittedVals[!is.na(data$dYield)] - mean(data$dYield, na.rm=T)) ^2
    obs_res2<-(data$dYield[!is.na(data$dYield)] - mean(data$dYield, na.rm=T))^2
    R2<-sum(pred_res2, na.rm=T)/sum(obs_res2, na.rm=T)
    
    #OUPUTS
    summtable<-summary$fixed[,1:5]
    summtable.pt2<-as.data.frame(summary$hyperpar[,1:5])
    summtable<-rbind(summtable,summtable.pt2)
    summtable<-kable(summtable, caption="Summary Table of Model Estimates")
    
    Diagnostics<-as.data.frame(t(c(DIC,CPO,MSE,R2)))
    colnames(Diagnostics)<-c("DIC","CPO","MSE","R2")
    Diagnostics<-kable(Diagnostics, caption="Model Diagnostic Metrics")
    
    figure<-ggarrange(h,p1,p2,ncol=3,nrow=1)
    annotate_figure(figure, top= text_grob("Model Diagnostic Plots"))
    
    my_list<-list(figure, Diagnostics)
    return(my_list)
  }  
  my_list<-list(summary)
  return(my_list)
}

##########################################################
###function to return diagnostics for stability models ###
##########################################################
model.results.cov<-function(model_results_name, data="d", crop="corn", diagnostics = TRUE){
  results<-readRDS(paste0(wd,"/Output/",model_results_name,".rds"))
  data<-readRDS(paste0(wd,"/Data/",data,".RDS"))
  # st_geometry(data)<-NULL
  summary<-summary(results)
  if (diagnostics == TRUE){ #if diagnostics turned on return the dic, cpo, and pit
    #R_INLA Diagnostics
    DIC<-(results$dic$dic)
    CPO<-(sum(log(results$cpo$cpo), na.rm=T)) 
    h<-ggplot(as.data.frame(results$cpo$pit), aes(x=results$cpo$pit)) + geom_histogram() + ggtitle( "Model PIT")
    #PPC Distribution Checks
    p_val<-c()
    n<-length(data$covYield)
    for(i in (1:n)){
      p_val[i]<-inla.pmarginal(q=data$covYield[i],
                               marginal=results$marginals.fitted.values[[i]])
    }
    p2<-ggplot(as.data.frame(p_val), aes(x=p_val)) + geom_histogram() + ggtitle("Posterior Predictive p-values")
    
    data<-cbind(data,results$summary.fitted.values$mean)
    colnames(data)[(dim(data)[2])]<-c("FittedVals")
    p1<-ggplot(data[!is.na(data$YIELD),], aes(x=covYield,y=FittedVals)) + geom_point(size=0.2, color="blue") +
      ggtitle("Scatter Plot of Predicted and Observed Values")  + xlab("Observed Value") + ylab("Mean of Post. Pred. Distr.")
    
    
    #PPC Summary Metrics
    sq_dif<-(data$covYield-data$FittedVals)^2
    MSE<-1/n*(sum(sq_dif,na.rm=T))
    
    pred_res2<-(data$FittedVals[!is.na(data$covYield)] - mean(data$covYield, na.rm=T)) ^2
    obs_res2<-(data$covYield[!is.na(data$covYield)] - mean(data$covYield, na.rm=T))^2
    R2<-sum(pred_res2, na.rm=T)/sum(obs_res2, na.rm=T)
    
    #OUPUTS
    summtable<-summary$fixed[,1:5]
    summtable.pt2<-as.data.frame(summary$hyperpar[,1:5])
    summtable<-rbind(summtable,summtable.pt2)
    summtable<-kable(summtable, caption="Summary Table of Model Estimates")
    
    Diagnostics<-as.data.frame(t(c(DIC,CPO,MSE,R2)))
    colnames(Diagnostics)<-c("DIC","CPO","MSE","R2")
    Diagnostics<-kable(Diagnostics, caption="Model Diagnostic Metrics")
    
    figure<-ggarrange(h,p1,p2,ncol=3,nrow=1)
    annotate_figure(figure, top= text_grob("Model Diagnostic Plots"))
    
    my_list<-list( figure, Diagnostics)
    return(my_list)
  }  
  my_list<-list(summary)
  return(my_list)
}

#################################################
###function to return maps of spatial effects ###
#################################################
plot.spatialeffect<-function(model_results_name, data, crop, type ="BYM", scale="CNTY", projection){
  data<-readRDS(paste0(wd,"/Data/",data,".RDS"))
  results<-readRDS(paste0(wd,"/Output/",model_results_name,".rds"))
  
  if (type =="BYM"){ #Spatial effect for bym models
    if (scale =="CNTY"){ #County spatial effects
      Nareas<-length(unique(data$CNTY)) #number of unique spatial locations in dataset
      asr<-results$summary.random$CNTY[1:Nareas,c(1,2,3)]#extract the area-specific residuals (this is Besag + iid)
      
      cnty<-readRDS(paste0(wd,"/Data/county.RDS"))
      cnty<-st_as_sf(cnty) %>% dplyr::filter(STATEFP == 20) %>% dplyr::select(.,c(5,11)) %>% st_transform(., projection)
      
      data2<- left_join(cnty, data, by="GEOID")
      
      #Create spatial dataframe with information for map and plot
      map.asr<-left_join(data2[,c("CNTY","GEOID","geometry")],asr,by=c("CNTY"="ID"))
      colnames(map.asr)<-c("CNTY","GEOID","Mean","Standard_Deviation","geometry")
      
      saveRDS(map.asr,paste0(wd,"/output/",model_results_name,"_asr.rds"))
      
      p1 <- ggplot(data=map.asr)+
        geom_sf(aes(fill= Mean), col="white", size=0.5, alpha=0.7) +
        labs(fill="Mean of Area Specific Residuals") +
        theme(legend.position = "bottom") + scale_fill_viridis_c(option="magma")
      
      p2 <- ggplot(data=map.asr)+
        geom_sf(aes(fill= Standard_Deviation), col="white", size=0.5, alpha=0.7) +
        labs(fill="Standard Deviation of Area Specific Residuals") +
        theme(legend.position = "bottom") + scale_fill_viridis_c(option="magma")
      
      my_list<-list(p1,p2)
      
      return (my_list)
    }
    
  } else {
    
    
    if (type=="FE"){
      Nareas<-length(unique(data$CNTY)) #number of unique spatial locations in dataset
      fe<-results$summary.fixed[2:Nareas,c(1,2)]#extract the area-specific effect
      fe$ID<-c(2:Nareas)
      
      cnty<-readRDS(paste0(wd,"/Data/county.RDS"))
      cnty<-st_as_sf(cnty) %>% dplyr::filter(STATEFP == 20) %>% dplyr::select(.,c(5,11)) %>% st_transform(., projection)
      
      data<- left_join(cnty, data, by="GEOID")
      
      #Create spatial dataframe with information for map and plot
      map.fe<-left_join(data[,c("CNTY","GEOID","geometry")],fe,by=c("CNTY"="ID"))
      colnames(map.fe)<-c("CNTY","GEOID","Mean","Standard_Deviation","geometry")
      
      saveRDS(map.fe,paste0(wd,"/output/",model_results_name,"_asr.rds"))
      
      p1 <- ggplot(data=map.fe)+
        geom_sf(aes(fill= Mean), col="white", size=0.5, alpha=0.7) +
        labs(fill="Mean of Fixed Effects") +
        theme(legend.position = "bottom") + scale_fill_viridis_c(option="magma")
      
      p2 <- ggplot(data=map.fe)+
        geom_sf(aes(fill= Standard_Deviation), col="white", size=0.5, alpha=0.7) +
        labs(fill="Standard Deviation of Fixed Effects") +
        theme(legend.position = "bottom") + scale_fill_viridis_c(option="magma")
      
      my_list<-list(p1,p2)
      
      return (my_list)
    }
  }
  
  }





##########################################################
###function to plot smooth effects of independent vars ###
##########################################################
plot.smootheffect<-function(model_results_name, var="SDI", dependent="log(Yield)"){
  
  results<-readRDS(paste0(wd,"/Output/",model_results_name,".rds"))
  
    p2<-ggplot(data=results$summary.random[[var]][ ,c(1,4,5,6)], aes (x=ID, y=`0.5quant`)) + geom_point() +
      geom_ribbon(aes(x=ID,ymin = `0.025quant`, ymax = `0.975quant`, fill="red"), alpha=0.25) + 
      theme(legend.position="none")+ xlab(var) + ylab(paste0("Estimated Effect on ",dependent))
    return(p2)
    
}



#########################################################################################################
###function to create summary plots combining all three crops effect estimates for each dependent var###
########################################################################################################
plot.lm.summ<-function(file_end, dependent, vars, crops, y.min, y.max, x.min, x.max){
  
  dependent.labels<-dependent %>% str_replace_all(c("yield"="YIELD","dyield"="dYield","covyield"="covYield"))
  vars.labels<-vars %>% str_replace_all(c("LSM_SIDI_ALL"="Diversity",
                                          "aveSIDI2yr" = "Average Diversity (2 years)",
                                          "dfert"="Change in Fertilizer",
                                          "dherb"="Change in Herbicide & Insecticide",
                                          "dVol"="Change in Irrigation Water Volume",
                                          "dTP"="Change in Total Precipitation",
                                          "dGDD"="Change in Growing Degree Days",
                                          "dSDD"="Change in Stress Degree Days",
                                          "dSIDI"="Change in Diversity",
                                          "aveSIDI3yr"="Average Diversity (3 years)",
                                          "covSIDI"="Coefficient of Variation \n for Diversity",
                                          "covfert"="Coefficient of Variation \n for Fertlizer",
                                          "covherb"="Coefficient of Variation \n for Herbicide & Insecticide",
                                          "covVol"="Coefficient of Variation \n for Irrigation Water Volume",
                                          "covTP"="Coefficient of Variation \n for Total Precipiation",
                                          "covGDD"="Coefficient of Variation \n for Growing Degree Days",
                                          "covSDD"="Coefficient of Variation \n for Stress Degree Days",
                                          "fert"="Fertilizer Expenses",
                                          "herb"="Herbicide & Insecticide Expenses",
                                          "volume"="Irrigation Water Volume",
                                          "TP"="Total Precipitation",
                                          "SDD"="Stress Degree Days",
                                          "GDD"="Growing Degree Days"))
  crops.labels<-crops
  crops.labels[crops.labels == "wheat"] <- "Winter wheat"
  crops.labels[crops.labels == "corn"] <- "Corn"
  crops.labels[crops.labels == "soy"] <- "Soy"
  
  
  my.plots<-c()
  
  foreach (i=1:length(vars))%do% {
    
    obs<-NA
    
    foreach (j=1:length(crops))%do% {
      #pull the model output by metric for ag diversity
      if (file.exists(paste0(wd,"/Output/", crops[j], "_",dependent[1],file_end,".rds"))){
        results_1<-readRDS(paste0(wd,"/Output/", crops[j], "_",dependent[1],file_end,".rds"))
        results_1<-get(vars[i],results_1$summary.random)[ ,c(1,4,5,6)] #pull the results for the effect of the metric for the crop
        results_1$Crop<-crops.labels[j]#add crop column
        obs<-c(obs, list(results_1))
      }
    }
    
    
    #combine data for a single plot
    library(data.table)
    
    results_tot<-as.data.frame(rbindlist(obs[2:4]))
    mean_lm<-results_tot %>% group_by (Crop) %>% summarise_at('ID', mean)
    sd_lm<-results_tot %>% group_by (Crop) %>% summarise_at('ID', sd)
    sd_lm$ID<-sd_lm$ID+mean_lm$ID
    sdlow_lm<-results_tot %>% group_by (Crop) %>% summarise_at('ID', sd)
    sdlow_lm$ID<-mean_lm$ID-sdlow_lm$ID
    
    #build a plot
     library("viridis")
    p1<-ggplot(data=results_tot) + geom_line(aes (x=ID, y=`0.5quant`, group = Crop, col=Crop), size=0.5) +
      xlim(x.min[1], x.max[1]) + ylim(y.min[i],y.max[i]) +
      xlab(paste0(vars.labels[i])) + ylab(paste0("Effect on ", dependent)) +    #+ ggtitle (paste0(subplot.labels[i])) 
      scale_x_continuous(breaks=c(-3,-2,-1, 0, 1,2,3), limits=c(x.min[i], x.max[i])) + #this overwrites the min-max limit which we want to be the same for all plots
      geom_ribbon(aes(x=ID, ymin = `0.025quant`, ymax = `0.975quant`, group=Crop, fill=Crop), alpha=0.25) +
      scale_color_viridis(discrete = TRUE, option = "D")+
      scale_fill_viridis(discrete = TRUE, option = "D") +
      theme_light() + theme(panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            panel.ontop = FALSE,
                            axis.text=element_text(size=10),
                            axis.title = element_text(size = 10))    
            
    
    my.plots[[i]]<-p1
    
  }
  return(my.plots)
}


###########################################################################################################################################
###function to create summary plots combining all climate response curves and input response curves for each crop-landscape metric model###
###########################################################################################################################################
plot.lm.summ.cl<-function(metrics, crops, dependent= "log (Yield)", y.min, y.max, x.min, x.max){
  
  metrics.labels<-toupper(metrics)
  crops.labels<-crops
  crops.labels[crops.labels == "wwheat"] <- "Winter wheat"
  crops.labels[crops.labels == "corn"] <- "Corn"
  crops.labels[crops.labels == "soy"] <- "Soy"
  control.labels<-c("TP","GDD","SDD","SOIL","fert","herb","vol")
  
  all.plots<-c()
  my.plots<-c()
  
  foreach (i=1:length(crops))%do% {
    
    foreach(j=1:length(control.labels))%do%{
    
        obs<-NA
    
        foreach (k=1:length(metrics))%do% {
        #pull the model output by metric for ag diversity
            if (file.exists(paste0(wd,"/Output/", crops[i], "_",metrics[k],".rds"))){
              results_1<-readRDS(paste0(wd,"/Output/", crops[i], "_",metrics[k],".rds"))
              results_1<-get(control.labels[j],results_1$summary.random)[ ,c(1,4,5,6)] #pull the results for the effect of the metric for the crop
              results_1$Metric<-metrics.labels[k]#add metric column
              results_1$control<-control.labels[j]
              results_1$Crop<-crops.labels[i]
              obs<-c(obs, list(results_1))
            }
        }
    
    
    #combine data for a single plot
    library(data.table)
    
    results_tot<-as.data.frame(rbindlist(obs[2:(length(metrics)+1)]))
    
    #build a plot
    p1<-ggplot(data=results_tot) + geom_line(aes (x=ID, y=`0.5quant`, group = Metric, col= Metric), size=0.5) +
      xlim(x.min[j,i], x.max[j,i]) + ylim(y.min[j,i],y.max[j,i]) +
      xlab(paste0(control.labels[j])) + ylab(paste0("Effect on", dependent)) +  ggtitle (paste0(crops.labels[i])) +
      geom_ribbon(aes(x=ID, ymin = `0.025quant`, ymax = `0.975quant`, group=Metric), alpha=0.025) +
      theme_classic()
    
    my.plots[[j]]<-p1
    
    }
    all.plots[[i]]<-my.plots
    }
  return(all.plots)
}



################################################################################
###function to return a summary table of model fit statistics for models
############################################################################
model.results.summ<-function(data, crop="corn", metrics=metrics){
  
  Diagnostics<-data.frame(matrix(ncol=6, nrow=length(metrics)))
  colnames(Diagnostics)<-c("Metric", "DIC","CPO","MSE","Mean P Value", "R2")
  
  foreach (i=1:length(metrics))%do% {
    
    results<-readRDS(paste0(wd,"/Output/",crop, "_", metrics[i],".rds"))
    d<-readRDS(paste0(wd,"/Data/",data,".RDS"))
    
    # return the dic, cpo
    #R_INLA Diagnostics
    Diagnostics$DIC[i]<-(results$dic$dic)
    Diagnostics$CPO[i]<-(sum(log(results$cpo$cpo), na.rm=T)) 
    
    #PPC Distribution Checks
    p_val<-c()
    n<-length(d$YIELD)
    for(j in (1:n)){
      p_val[j]<-inla.pmarginal(q=d$YIELD[j],
                               marginal=results$marginals.fitted.values[[j]])
    }
    Diagnostics$`Mean P Value`[i]<-mean(p_val, na.rm=T)
    
    d<-cbind(d,results$summary.fitted.values$mean)
    colnames(d)[(dim(d)[2]-1)]<-c("FittedVals")
    
    
    #PPC Summary Metrics
    sq_dif<-(d$YIELD-d$FittedVals)^2
    Diagnostics$MSE[i]<-1/n*(sum(sq_dif,na.rm=T))
    
    pred_res2<-(d$FittedVals[!is.na(d$YIELD)] - mean(d$YIELD, na.rm=T)) ^2
    obs_res2<-(d$YIELD[!is.na(d$YIELD)] - mean(d$YIELD, na.rm=T))^2
    Diagnostics$R2[i]<-sum(pred_res2, na.rm=T)/sum(obs_res2, na.rm=T)
    Diagnostics$Metric[i]<-metrics[i]
    #OUPUTS
    
  }
  
  
  Diag<-regulartable(Diagnostics)
  print(Diag, preview="docx")
  return(Diagnostics)
}



