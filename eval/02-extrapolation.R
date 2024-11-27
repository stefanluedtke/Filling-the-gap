setwd("~/projects/OTC_Thuenen/Filling-the-gap/")

source("code/models.R")
source("code/utils.R")

require(data.table)
require(xtable)



evalMuchMissing.all = function(){
  
  set.seed(42)
  
  fns = c("bias","bass","herring")
  res.compl = do.call("rbind",lapply(fns,function(fn){
    
    print(fn)
    data = loadData(fn)
    data = preprocData.simple(data)
    
    ###############
    # project coordinates to local grid
    sf_wgs84 <- st_as_sf(data, coords = c("WEST","SOUTH"), crs = 4326) #correct order 
    sf_etrs89_projected <- st_transform(sf_wgs84, crs = 3035)
    #sf_etrs89_projected <- st_transform(sf_wgs84, crs = 32634) #UTM Zone 34N
    coords <- st_coordinates(sf_etrs89_projected)
    orig.south=data$SOUTH
    data$WEST = coords[,1]
    data$SOUTH = coords[,2]
    ##############
    
    years = unique(data$Year)
    
  
    
    #filter some years that do not work for unknown reasons...
    if(fn=="bass"){
      years = years[years!=1999 & years != 2000]
    }
    
    xx = function(tr,te) traintestXGB(tr,te,nrounds=20,max.depth = 3)
    en = function(tr,te) doTargetEncoding.TrainTest.vy.vr(tr,te)
    res.xgb.vy.vr = do.call("rbind",lapply(years,function(y){
      print(y)
      evalMuchMissing(data,xx,testyear = y,featurefun = en,orig.south=orig.south)
    }))
    res.xgb.vy.vr$model="XGB-noSDInteraction"
    

    res.lmm= do.call("rbind",lapply(years,function(y){
      print(y)
      evalMuchMissing(data,function(tr,te) traintestLMER2(tr,te),testyear = y,featurefun = doNoFilter,orig.south=orig.south)
    }))
    res.lmm$model="LMM-noSDEffect"
    
    res.gam.m.a = do.call("rbind",lapply(years,function(y){
      print(y)
      evalMuchMissing(data,function(tr,te) trainTestGAM2(tr,te,k=15),testyear = y,featurefun = doNoFilter,orig.south=orig.south) 
    }))
    res.gam.m.a$model="GAM-M"
    
    
    
    res.compl = rbind(res.xgb.vy.vr,res.lmm,res.gam.m.a) 
    
    res.compl$type=fn
    res.compl
    
  }))

  write.table(res.compl,"results/extrapolation.csv",sep=",",row.names = F,col.names = T)
  
  #summarise results
  res.compl %>%
    group_by(model,type,variable) %>%
    summarise(rmse=sqrt(mean((pred-truth)^2)),r2=cor(pred,truth)^2) %>%
    ungroup %>%
    group_by(model,type) %>%
    summarise(rmse=mean(rmse),r2=mean(r2))-> tp
  tp
    
}


summarizeResults = function(){
  
  res.compl=fread("results/extrapolation.csv")
  
  fns = c("bass","bias","herring")
  #scaling factor: need all truth data first
  sa = do.call("rbind",lapply(fns,function(fn){
    data = loadData(fn)
    data = preprocData.simple(data)
    data$type=fn
    
    sa = data %>%
      group_by(variable,type) %>% 
      summarise(scale=sd(value))
    sa
  }))
  
  res.compl = merge(res.compl,sa,by=c("variable","type"))
  
  
  res.compl %>%
    group_by(model,type,variable) %>%
    summarise(rmse=sqrt(mean((pred-truth)^2)),r2=cor(pred,truth)^2,scale=mean(scale)) %>%
    ungroup %>%
    group_by(model,type) %>%
    summarise(rmse=mean(rmse)/mean(scale),r2=mean(r2))-> tp
  tp
  
  
  #plot them
  p = ggplot(tp)+
    geom_point(aes(x=model,y=r2))+
    facet_wrap(.~type,scales = "free_y")+
    theme_light()
  p
  ggsave(paste0("figures/extrapolation/extrapolation_r2.png"),p,width=6,heigh=3)
  
  p = ggplot(tp)+
    geom_point(aes(x=model,y=rmse))+
    facet_wrap(.~type,scales = "free_y")+
    theme_light()
  p
  ggsave(paste0("figures/extrapolation/extrapolation_rmse.png"),p,width=6,heigh=3)
  
  
  tp$ms = paste0(tp$model,"-",tp$type)
  tp = tp[order(tp$type),]
  xtable(t(tp[,-5]))
  
}



##################
#helper functions

evalMuchMissing = function(data.agewise.coord,traintestfun,testyear=2020,featurefun=doNoFilter,orig.south=NA){
  

  train = data.agewise.coord[!(data.agewise.coord$Year == testyear & orig.south < 57.5),]
  test = data.agewise.coord[(data.agewise.coord$Year == testyear & orig.south < 57.5),]
  
  tt = featurefun(train,test)
  train.enc = tt[[1]]
  test.enc = tt[[2]]

  
  res = traintestfun(train.enc,test.enc)
  
  res$Year=test$Year
  res$Sub_Div = test$Sub_Div
  res$RECT=test$RECT
  res$variable=test$variable
  
  res
  
  
}
