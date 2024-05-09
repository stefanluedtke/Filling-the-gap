setwd("~/projects/OTC_Thuenen/SpratSurvey-SL/code/toPublish/")

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
    
    years = unique(data$Year)
    
  
    
    #filter some years that do not work for unknown reasons...
    if(fn=="bass"){
      years = years[years!=1999 & years != 2000]
    }
    
    xx = function(tr,te) traintestXGB(tr,te,nrounds=20,max.depth = 3)
    en = function(tr,te) doTargetEncoding.TrainTest.vy.vr(tr,te)
    res.xgb.vy.vr = do.call("rbind",lapply(years,function(y){
      print(y)
      evalMuchMissing(data,xx,testyear = y,featurefun = en)
    }))
    res.xgb.vy.vr$model="xgb"
    
    
    # res.gam.m= do.call("rbind",lapply(years,function(y){
    #   print(y)
    #   evalMuchMissing(data,function(tr,te) trainTestGAM2(tr,te,k=15),testyear = y,featurefun = doNoFilter) 
    # }))
    # res.gam.m$model="gam2"
    
    res.lmm= do.call("rbind",lapply(years,function(y){
      print(y)
      evalMuchMissing(data,function(tr,te) traintestLMER2(tr,te),testyear = y,featurefun = doNoFilter)
    }))
    res.lmm$model="lmm2"
    
    res.gam.m.a = do.call("rbind",lapply(years,function(y){
      print(y)
      evalMuchMissing(data,function(tr,te) trainTestGAM3(tr,te,k=15),testyear = y,featurefun = doNoFilter) 
    }))
    res.gam.m.a$model="gam3"
    
    
    
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
  
  res.compl %>%
    group_by(model,type,variable) %>%
    summarise(rmse=sqrt(mean((pred-truth)^2)),r2=cor(pred,truth)^2) %>%
    ungroup %>%
    group_by(model,type) %>%
    summarise(rmse=mean(rmse),r2=mean(r2))-> tp
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

evalMuchMissing = function(data.agewise.coord,traintestfun,testyear=2020,featurefun=doNoFilter){
  

  train = data.agewise.coord[!(data.agewise.coord$Year == testyear & data.agewise.coord$SOUTH < 57.5),]
  test = data.agewise.coord[(data.agewise.coord$Year == testyear & data.agewise.coord$SOUTH < 57.5),]
  
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
