setwd("~/projects/OTC_Thuenen/Filling-the-gap/")

source("code/models.R")
source("code/utils.R")

require(data.table)
require(ggplot2)

#assign colors to the different methods
group.colors <- c(GAM1 = "#6baed6", GAM2 = "#08519c", LMM1 = "#74c476", LMM2 = "#006d2c", XGB1 = "#fd8d3c", XGB2 = "#a63603", Baseline = "#252525", GAM2.tweedie ="#a50f15")



doAllCV = function(){
  
  set.seed(42)
  
  fns = c("bass","bias","herring")
  
  for(fn in fns){
    
    data = loadData(fn)
    data = preprocData.simple(data)
    
    ###############
    # project coordinates to local grid
    sf_wgs84 <- st_as_sf(data, coords = c("SOUTH","WEST"), crs = 4326) 
    sf_etrs89_projected <- st_transform(sf_wgs84, crs = 3035)
    #sf_etrs89_projected <- st_transform(sf_wgs84, crs = 32634) #UTM Zone 34N
    coords <- st_coordinates(sf_etrs89_projected)
    data$WEST = coords[,1]
    data$SOUTH = coords[,2]
    ##############
    
    res.complete = do.call("rbind",lapply(1:10,function(it){
      
      shuffle = sample(nrow(data),nrow(data),replace = F)
      data = data[shuffle,]
      
      print(paste("iteration",it))
      
      res.complete = do.call("rbind",lapply(c(0.1,0.25,0.5,0.75,0.99),function(testfraction){ #
        
        print(paste("frac ",testfraction))
        
        res.baseline = doCV.RectPercent(data,traintestBaseline,doNoFilter,testfraction = testfraction,getPredicion = TRUE,doPrint = F)
        res.baseline$model="Baseline"
        
        res.lmer1 = doCV.RectPercent(data,traintestLMER1,doNoFilter,testfraction = testfraction,getPredicion = TRUE,doPrint = F)
        res.lmer1$model="LMM-SDEffect"

        res.lmer2 = doCV.RectPercent(data,traintestLMER2,doNoFilter,testfraction = testfraction,getPredicion = TRUE,doPrint = F)
        res.lmer2$model="LMM-noSDEffect"

        res.gam1 = doCV.RectPercent(data,trainTestGAM1,doNoFilter,testfraction = testfraction,getPredicion = TRUE,doPrint = F)
        res.gam1$model="GAM-noM"

        res.gam2 = doCV.RectPercent(data,trainTestGAM2,doNoFilter,testfraction = testfraction,getPredicion = TRUE,doPrint = F)
        res.gam2$model="GAM-M"

        res.xgb1 = doCV.RectPercent(data,traintestXGB,doTargetEncoding.TrainTest.vy.vr,testfraction = testfraction,getPredicion = TRUE,doPrint = F)
        res.xgb1$model="XGB-noSDInteraction"
        
        res.xgb2 = doCV.RectPercent(data,traintestXGB,doTargetEncoding.TrainTest.vsy.vr,testfraction = testfraction,getPredicion = TRUE,doPrint = F)
        res.xgb2$model="XGB-SDInteraction"
        
        #todo other models
        res = rbind(res.baseline,res.lmer1,res.lmer2,res.gam2,res.gam2,res.xgb1,res.xgb2)

        res$testfraction = testfraction
        res$iteration=it
        
        res
        
      }))
    }))
    
    write.table(res.complete,
                paste0("results/imputation_cv_",fn,".csv"),
                sep=",",col.names = TRUE,row.names = FALSE)
    
  }
}



plotCV.grouped = function(){
  
  for(fn in c("bass","bias","herring")){
    
    
    res.complete = fread(paste0("results/imputation_cv_",fn,".csv"))
    tp = res.complete
    
    
    #change model names for plotting
    tp$model[tp$model=="baseline"] = "Baseline"
    tp$model[tp$model=="gam1"] = "GAM1"
    tp$model[tp$model=="gam2"] = "GAM2"
    tp$model[tp$model=="lmer1"] = "LMM1"
    tp$model[tp$model=="lmer2"] = "LMM2"
    tp$model[tp$model=="xgb1"] = "XGB1"
    tp$model[tp$model=="xgb2"] = "XGB2"
    
    tp$mtype = "Baseline"
    tp$mtype[tp$model %in% c("GAM1","GAM2")] = "GAM"
    tp$mtype[tp$model %in% c("LMM1","LMM2")] = "LMM"
    tp$mtype[tp$model %in% c("XGB1","XGB2")] = "XGB"
    
    tp$ftype = "sd"
    tp$ftype[tp$model %in% c("GAM1","GAM2","LMM2","XGB1")] = "no_sd"
    
    
    #scaling factor: need all truth data first
    data = loadData(fn)
    data = preprocData.simple(data)
    #scale per age group
    sa = data %>%
      group_by(variable) %>% 
      summarise(scale=sd(value))
    tp = merge(tp,sa,by="variable")
    

    tp %>%
      dplyr::group_by(model,variable,testfraction,iteration,mtype,ftype) %>%
      dplyr::summarise(r2=cor(pred,truth)^2,rmse=sqrt(mean((pred-truth)^2))/mean(scale)) %>%
      ungroup() %>%
      group_by(model,variable,testfraction,mtype,ftype) %>%
      dplyr::summarise(r2m=mean(r2),r2l=t.test(r2)$conf.int[1],r2u=t.test(r2)$conf.int[2],
                       rmsem=mean(rmse),rmsel=t.test(rmse)$conf.int[1],rmseu=t.test(rmse)$conf.int[2])-> tp2
    
    p = ggplot(tp2)+
      geom_line(aes(x=testfraction,y=r2m,color=model,linetype=ftype,group=model))+
      facet_wrap(~variable,scales="free_y")+
      #geom_ribbon(aes(x=testfraction,ymin = r2l,ymax = r2u, fill = model), alpha = 0.1, show.legend = F)+
      theme(legend.position = "right")+
      labs(x="p",y="R2",color="Model",linetype="Features")+
      theme_light()+      
      scale_color_manual(values=group.colors)
    
    p
    
    ggsave(paste0("figures/cv/all-impute-agewise-r2-",fn,".png"),p,width=8,height=6)
    
    p = ggplot(tp2)+
      geom_line(aes(x=testfraction,y=rmsem,color=model,linetype=ftype,group=model))+
      #geom_ribbon(aes(x=testfraction,ymin = rmsel,ymax = rmseu, fill = mod), alpha = 0.1, show.legend = F)+
      facet_wrap(~variable,scales = "free_y")+ #
      theme(legend.position = "right")+
      labs(x="p",y="NRMSE",color="Model",linetype="Features")+
      theme_light()+
      scale_color_manual(values=group.colors)
    
    p
    
    ggsave(paste0("figures/cv/all-impute-agewise-rmse-",fn,".png"),p,width=8,height=6)
    
    
    tp %>%
      dplyr::group_by(model,variable,testfraction,iteration,ftype,mtype) %>%
      dplyr::summarise(r2=cor(pred,truth)^2,rmse=sqrt(mean((pred-truth)^2)),scale=mean(scale)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(model,testfraction,iteration,ftype,mtype) %>%
      dplyr::summarise(r2=mean(r2),rmse=mean(rmse)/mean(scale)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(model,testfraction,ftype,mtype) %>%
      dplyr::summarise(r2m=mean(r2),r2l=t.test(r2)$conf.int[1],r2u=t.test(r2)$conf.int[2],
                       rmsem=mean(rmse),rmsel=t.test(rmse)$conf.int[1],rmseu=t.test(rmse)$conf.int[2])-> tp2

    
    p = ggplot(tp2)+
      geom_line(aes(x=testfraction,y=r2m,color=model,linetype=ftype,group=model))+
      geom_point(aes(x=testfraction,y=r2m,color=model,group=model))+
      theme(legend.position = "right")+
      labs(x="p",y="R2",color="Model",linetype="Features")+
      theme_light()+
      scale_color_manual(values=group.colors)
    p
    
    ggsave(paste0("figures/cv/all-impute-grouped-r2-",fn,".png"),p,width=6,height=4)
    
    
    p = ggplot(tp2)+
      geom_line(aes(x=testfraction,y=rmsem,color=model,linetype=ftype,group=model))+
      geom_point(aes(x=testfraction,y=rmsem,color=model,group=model))+
      theme(legend.position = "right")+
      labs(x="p",y="NRMSE",color="Model",linetype="Features")+
      theme_light()+
      scale_color_manual(values=group.colors)
    p
    
    ggsave(paste0("figures/cv/all-impute-grouped-rmse-",fn,".png"),p,width=6,height=4)
    
  
    
  }
  
}



doCV.tweedie = function(){
  
  set.seed(42)
  
  fns = c("bias")
  
  for(fn in fns){
    
    data = loadData(fn)
    data = preprocData.simple(data)
    
    res.complete = do.call("rbind",lapply(1:5,function(it){
      
      shuffle = sample(nrow(data),nrow(data),replace = F)
      data = data[shuffle,]
      
      print(paste("iteration",it))
      
      res.complete = do.call("rbind",lapply(c(0.25,0.75,0.99),function(testfraction){ #
        
        print(paste("frac ",testfraction))
        
        res.gam2 = doCV.RectPercent(data,trainTestGAM2,doNoFilter,testfraction = testfraction,getPredicion = TRUE,doPrint = T)
        res.gam2$model="gam2"
        
        res.gam2.tweedie = doCV.RectPercent(data,trainTestGAM2.tweedie,doNoFilter,testfraction = testfraction,getPredicion = TRUE,doPrint = T)
        res.gam2.tweedie$model="gam2.tweedie"
        
        #todo other models
        res = rbind(res.gam2,res.gam2.tweedie)
        
        res$testfraction = testfraction
        res$iteration=it
        
        res
        
      }))
    }))
    
    write.table(res.complete,
                paste0("results/imputation_cv_tweedie_",fn,".csv"),
                sep=",",col.names = TRUE,row.names = FALSE)
    
  }
}


plotTweedie = function(){
  
  
  res = fread("results/imputation_cv_tweedie_bias.csv")
  
  
  #scaling factor: need all truth data first
  data = loadData("bias")
  data = preprocData.simple(data)
  #scale per age group
  sa = data %>%
    group_by(variable) %>% 
    summarise(scale=sd(value))
  res = merge(res,sa,by="variable")
  
  res %>%
    dplyr::group_by(model,variable,testfraction,iteration) %>%
    dplyr::summarise(r2=cor(pred,truth)^2,rmse=sqrt(mean((pred-truth)^2)),scale=mean(scale)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(model,testfraction,iteration) %>%
    dplyr::summarise(r2=mean(r2),rmse=mean(rmse)/mean(scale)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(model,testfraction) %>%
    dplyr::summarise(r2m=mean(r2),r2l=t.test(r2)$conf.int[1],r2u=t.test(r2)$conf.int[2],
                     rmsem=mean(rmse),rmsel=t.test(rmse)$conf.int[1],rmseu=t.test(rmse)$conf.int[2])-> tp2
  
  tp2$model[tp2$model=="gam2"] = "GAM2"
  tp2$model[tp2$model=="gam2.tweedie"] = "GAM2.tweedie"
  
  p = ggplot(tp2)+
    geom_line(aes(x=testfraction,y=r2m,color=model))+
    geom_point(aes(x=testfraction,y=r2m,color=model))+
    theme(legend.position = "right")+
    labs(x="p",y="R2",color="Model")+
    theme_light()+
    scale_color_manual(values=group.colors)
    
  p
  
  ggsave(paste0("figures/cv/tweedie-impute-r2-","bias",".png"),p,width=6,height=4)
  
  
  p = ggplot(tp2)+
    geom_line(aes(x=testfraction,y=rmsem,color=model))+
    geom_point(aes(x=testfraction,y=rmsem,color=model))+
    theme(legend.position = "right")+
    labs(x="p",y="NRMSE",color="Model")+
    theme_light()+
    scale_color_manual(values=group.colors)
  
  p
  
  ggsave(paste0("figures/cv/tweedie-impute-rmse-","bias",".png"),p,width=6,height=4)
  
  
  
}


#########
#helper functions

doCV.RectPercent = function(data.agewise.pca,traintestfun,featurefun = doSimpleFilter,testfraction=0.5,getPredicion=FALSE,getAgewiseError=FALSE,doPrint=FALSE){
  
  #to get an idea of number of rects per SD: between 5 and 12 (sometimes just 2)
  # nn = data.agewise.pca %>% group_by(Sub_Div,year) %>% summarise(n()/8) 
  foldcol =  paste0(data.agewise.pca$Year,"-",data.agewise.pca$Sub_Div)
  folds = unique(foldcol)
  
  
  cvres = do.call("rbind",lapply(1:length(folds),function(foldi){
    
    fold = folds[foldi]
    #print(fold)
    
    data.sd = data.agewise.pca[foldcol==fold,]
    #unique rects in that sd
    rects = unique(data.sd$RECT)
    #number of rects in test
    ntest = floor(length(rects)*testfraction)
    if(ntest == length(rects)){
      ntest = ntest-1
    }
    if(ntest == 0){
      ntest = 1
    }
    
    if(doPrint){
      print(paste0(foldi,"/",length(folds),", test on ",ntest,"/",length(rects)))
    }
    
    testrects = rects[1:ntest]
    trainrects = rects[(ntest+1):length(rects)]
    
    data.sd.train = data.sd[data.sd$RECT %in% trainrects,]
    data.sd.test = data.sd[data.sd$RECT %in% testrects,]
    
    #complete train data is a combination data not in this SD:year and the test data of this SD:year
    train.d = (data.agewise.pca[foldcol!=fold,])
    train.d = rbind(train.d,data.sd.train)
    
    test.d = data.sd.test
    
    trte  = featurefun(train.d,test.d)
    train = trte[[1]]
    test = trte[[2]]
    
    res = traintestfun(train,test)
    res$variable = test.d$variable # data.agewise.pca$variable[foldcol==fold]
    res$RECT = test.d$RECT # data.agewise.pca$RECT[foldcol==fold]
    res$Year = fold
    res
  }))
  
  if(getPredicion){
    return(cvres)
  }
  if(getAgewiseError){
    cvres %>%
      group_by(variable) %>%
      summarise(rmse=sqrt(mean((truth-pred)^2)),r2=cor(truth,pred)^2) -> res.agewise
    return(res.agewise)
  }
  #cvres = cvres[!is.na(cvres$pred),]
  rmse = sqrt(mean((cvres$truth-cvres$pred)^2)) 
  r2 = cor(cvres$truth,cvres$pred)^2
  
  data.frame(rmse,r2)
  
  
}
