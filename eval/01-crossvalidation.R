setwd("~/projects/OTC_Thuenen/SpratSurvey-SL/code/toPublish/")

source("code/models.R")
source("code/utils.R")

require(data.table)
require(ggplot2)

doAllCV = function(){
  
  
  fns = c("bass","bias","herring")
  
  for(fn in fns){
    
    data = loadData(fn)
    data = preprocData.simple(data)
    
    res.complete = do.call("rbind",lapply(1:10,function(it){
      
      shuffle = sample(nrow(data),nrow(data),replace = F)
      data = data[shuffle,]
      
      print(paste("iteration",it))
      
      res.complete = do.call("rbind",lapply(c(0.1,0.25,0.5,0.75,0.99),function(testfraction){ #
        
        print(paste("frac ",testfraction))
        
        res.baseline = doCV.RectPercent(data,traintestBaseline,doNoFilter,testfraction = testfraction,getPredicion = TRUE,doPrint = F)
        res.baseline$model="baseline"
        
        res.lmer1 = doCV.RectPercent(data,traintestLMER1,doNoFilter,testfraction = testfraction,getPredicion = TRUE,doPrint = F)
        res.lmer1$model="lmer1"

        res.lmer2 = doCV.RectPercent(data,traintestLMER2,doNoFilter,testfraction = testfraction,getPredicion = TRUE,doPrint = F)
        res.lmer2$model="lmer2"

        # res.gam1 = doCV.RectPercent(data,trainTestGAM1,doNoFilter,testfraction = testfraction,getPredicion = TRUE,doPrint = T)
        # res.gam1$model="gam1"

        res.gam2 = doCV.RectPercent(data,trainTestGAM2,doNoFilter,testfraction = testfraction,getPredicion = TRUE,doPrint = F)
        res.gam2$model="gam2"

        res.gam3 = doCV.RectPercent(data,trainTestGAM3,doNoFilter,testfraction = testfraction,getPredicion = TRUE,doPrint = F)
        res.gam3$model="gam3"

        res.xgb1 = doCV.RectPercent(data,traintestXGB,doTargetEncoding.TrainTest.vy.vr,testfraction = testfraction,getPredicion = TRUE,doPrint = F)
        res.xgb1$model="xgb1"
        
        res.xgb2 = doCV.RectPercent(data,traintestXGB,doTargetEncoding.TrainTest.vsy.vr,testfraction = testfraction,getPredicion = TRUE,doPrint = F)
        res.xgb2$model="xgb2"
        
        #todo other models
        res = rbind(res.baseline,res.lmer1,res.lmer2,res.gam2,res.gam3,res.xgb1,res.xgb2)

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


plotCV = function(){
  
  fn = "herring"
  
  res.complete = fread(paste0("results/imputation_cv_",fn,".csv"))
  tp = res.complete

  tp %>%
    dplyr::group_by(model,variable,testfraction,iteration) %>%
    dplyr::summarise(r2=cor(pred,truth)^2,rmse=sqrt(mean((pred-truth)^2))) %>%
    ungroup() %>%
    group_by(model,variable,testfraction) %>%
    dplyr::summarise(r2m=mean(r2),r2l=t.test(r2)$conf.int[1],r2u=t.test(r2)$conf.int[2],
                     rmsem=mean(rmse),rmsel=t.test(rmse)$conf.int[1],rmseu=t.test(rmse)$conf.int[2])-> tp2
    #dplyr::summarise(r2m=mean(r2),rmsem=mean(rmse))-> tp2
  
  p = ggplot(tp2)+
    geom_line(aes(x=testfraction,y=r2m,color=model))+
    facet_wrap(~variable,scales = "free_y")+
    #geom_ribbon(aes(x=testfraction,ymin = r2l,ymax = r2u, fill = mod), alpha = 0.1, show.legend = F)+
    theme(legend.position = "right")+
    labs(x="p",y="R2",color="Model")+
    theme_light()
  p
  
  #ggsave(paste0("figures/cv/best-impute-agewise-r2-",fn,".png"),p,width=8,height=6)
  
  p = ggplot(tp2)+
    geom_line(aes(x=testfraction,y=rmsem,color=model))+
    #geom_ribbon(aes(x=testfraction,ymin = rmsel,ymax = rmseu, fill = mod), alpha = 0.1, show.legend = F)+
    facet_wrap(~variable,scales = "free_y")+
    theme(legend.position = "right")+
    labs(x="p",y="RMSE",color="Model")+
    theme_light()
  p
  
  #ggsave(paste0("figures/cv/best-impute-agewise-rmse-",fn,".png"),p,width=8,height=6)
  
  
  
  tp %>%
    dplyr::group_by(model,variable,testfraction,iteration) %>%
    dplyr::summarise(r2=cor(pred,truth)^2,rmse=sqrt(mean((pred-truth)^2))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(model,testfraction,iteration) %>%
    dplyr::summarise(r2=mean(r2),rmse=mean(rmse)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(model,testfraction) %>%
    dplyr::summarise(r2m=mean(r2),r2l=t.test(r2)$conf.int[1],r2u=t.test(r2)$conf.int[2],
                     rmsem=mean(rmse),rmsel=t.test(rmse)$conf.int[1],rmseu=t.test(rmse)$conf.int[2])-> tp2
    #dplyr::summarise(r2m=mean(r2),rmsem=mean(rmse))-> tp2
  
  
  p = ggplot(tp2)+
    geom_line(aes(x=testfraction,y=r2m,color=model))+
    #geom_ribbon(aes(x=testfraction,ymin = r2l,ymax = r2u, fill = mod), alpha = 0.1, show.legend = F)+
    geom_point(aes(x=testfraction,y=r2m,color=model))+
    theme(legend.position = "right")+
    labs(x="p",y="R2",color="Model")+
    theme_light()
  p
  
  #ggsave(paste0("figures/cv/best-impute-r2-",fn,".png"),p,width=6,height=4)
  
  
  p = ggplot(tp2)+
    geom_line(aes(x=testfraction,y=rmsem,color=model))+
    geom_point(aes(x=testfraction,y=rmsem,color=model))+
    #geom_ribbon(aes(x=testfraction,ymin = rmsel,ymax = rmseu, fill = mod), alpha = 0.1, show.legend = F)+
    theme(legend.position = "right")+
    labs(x="p",y="RMSE",color="Model")+
    theme_light()
  p
  
  #ggsave(paste0("figures/cv/best-impute-rmse-",fn,".png"),p,width=6,height=4)
  
  
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
