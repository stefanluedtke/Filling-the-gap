setwd("~/projects/OTC_Thuenen/Filling-the-gap/")

source("code/models.R")
source("code/utils.R")

require(ggplot2)
# require(cowplot)
library(ggpubr)


#assign colors to the different methods
group.colors <- c( "GAM" = "#08519c", #"GAM-year" = "#032c57","GAM-noM" = "#6baed6",
                   "LMM-noFE" = "#fd8d3c", "LMM-noSD" = "#d1117b","LMM-full" = "#a63603",
                   "XGB-noSD" = "#a9f5ab", "XGB-SD" = "#018c39", 
                   "Baseline" = "#121212")



evalIndexStability.all = function(){
  
  set.seed(42)
  
  for(type in c("bias","bass","herring")){
    
    data = loadData(type)
    data = preprocData.simple(data)
    
    ###############
    # project coordinates to local grid
    sf_wgs84 <- st_as_sf(data, coords = c("WEST","SOUTH"), crs = 4326) #correct order 
    sf_etrs89_projected <- st_transform(sf_wgs84, crs = 3035)
    #sf_etrs89_projected <- st_transform(sf_wgs84, crs = 32634) #UTM Zone 34N
    coords <- st_coordinates(sf_etrs89_projected)
    data$WEST = coords[,1]
    data$SOUTH = coords[,2]
    ##############
    
    if(type=="bias"){
      sds = c("22","23","24","25","26","27","28_2","29","32")
    }else if(type=="bass"){
      sds = c("24","25","26","28_2")
    }else if(type=="herring"){
      sds = c("25","26","27","28_2","29")
    }
    
    stab.all = do.call("rbind",lapply(c(0.1,0.25,0.5,0.75,0.99),function(removeFrac){
      print(removeFrac)
      
      if(type=="bass"){
        skipyears = c(1991:2001,2016)
      }else if(type=="bias"){
        skipyears = c(1991:2008)
      }else if(type=="herring"){
        skipyears = c(1991:2001,2008)
      }
      
      #baseline
      stab.baseline = evalIndexStability(data,sds,doNoFilter,traintestBaseline,removeFrac,skipyears = skipyears)
      stab.baseline$removeFrac = removeFrac
      stab.baseline$type = "Baseline"
      
      #complex LMER
      stab.lmerComplex = evalIndexStability(data,sds,doNoFilter,traintestLMER3,removeFrac,skipyears = skipyears)
      stab.lmerComplex$removeFrac = removeFrac
      stab.lmerComplex$type = "LMM-SDEffect-year"
      
      #XGB (the one with the rect interaction effect)
      stab.xgb = evalIndexStability(data,sds,doTargetEncoding.TrainTest.vsy.vr,traintestXGB,removeFrac,skipyears = skipyears)
      stab.xgb$removeFrac = removeFrac
      stab.xgb$type = "XGB-SDInteraction"
      
      #GAM
      stab.gam = evalIndexStability(data,sds,doNoFilter,trainTestGAM2,removeFrac,skipyears = skipyears)
      stab.gam$removeFrac = removeFrac
      stab.gam$type = "GAM-M"
      
      stab.all = rbind(stab.baseline,stab.lmerComplex,stab.xgb,stab.gam)
    }))
    
    write.table(stab.all,paste0("results/index_leaveOut_",type,".csv"),sep=",",row.names = F,col.names = T)
    
  }
  
}



plotIndexStability = function(){
  
  for(type in c("bass","bias","herring")){
    
    stab.all = read.csv(paste0("results/index_leaveOut_",type,".csv"))
    iterations = 10
    
    stab.all$type[stab.all$type=="GAM-M"] = "GAM"
    stab.all$type[stab.all$type=="LMM-SDEffect-year"] = "LMM-full"
    stab.all$type[stab.all$type=="XGB-SDInteraction"] = "XGB-SD"
    
    
    #scale: sd of true index of that age group
    
    tpa = do.call("rbind",lapply(unique(stab.all$variable),function(v){
      do.call("rbind",lapply(unique(stab.all$type),function(t){
        do.call("rbind",lapply(unique(stab.all$removeFrac),function(r){
          index.baseline.leaveout = stab.all[stab.all$variable == v & stab.all$type==t & stab.all$removeFrac==r & stab.all$it != "truth",]
          index.baseline = stab.all[stab.all$variable == v & stab.all$type==t & stab.all$removeFrac==r & stab.all$it == "truth",]
          index.baseline = index.baseline[order(index.baseline$variable,index.baseline$Year),]
          index.baseline.leaveout = index.baseline.leaveout[order(index.baseline.leaveout$variable,index.baseline.leaveout$Year),]
          rmses = do.call("rbind",lapply(1:iterations,function(it){
            ii = index.baseline.leaveout[index.baseline.leaveout$it==it,]
            ii = ii[order(ii$variable,ii$Year),]
            sqrt(mean((ii$index-index.baseline$index)^2))
          }))
          r2s = do.call("rbind",lapply(1:iterations,function(it){
            ii = index.baseline.leaveout[index.baseline.leaveout$it==it,]
            ii = ii[order(ii$variable,ii$Year),]
            cor(ii$index,index.baseline$index)^2
          }))
          data.frame(variable=v,type=t,removeFrac=r,
                     rmsem=mean(rmses)/sd(index.baseline$index),rmsel=t.test(rmses)$conf.int[1],rmseu=t.test(rmses)$conf.int[2],
                     r2m=mean(r2s),r2l=t.test(r2s)$conf.int[1],r2u=t.test(r2s)$conf.int[2])
        }))
      }))
    }))
    
    
    p = ggplot(tpa)+
      geom_line(aes(x=removeFrac,y=r2m,color=type))+
      geom_ribbon(aes(x=removeFrac,ymin = r2l,ymax = r2u, fill = type), alpha = 0.1, show.legend = F)+
      facet_wrap(.~variable)+
      labs(x="p",y="R2",color="Model")+
      theme_light()+
      scale_color_manual(values=group.colors)
    p
    
    p = ggplot(tpa)+
      geom_line(aes(x=removeFrac,y=rmsem,color=type))+
      #geom_ribbon(aes(x=removeFrac,ymin = rmsel,ymax = rmseu, fill = type), alpha = 0.1, show.legend = F)+
      facet_wrap(.~variable)+
      labs(x="p",y="RMSE",color="Model")+
      theme_light()+
      scale_color_manual(values=group.colors)
    p
    
    
    #tp2 = tpa %>% group_by(removeFrac,type) %>% summarise(rmse=mean(rmse),r2=mean(r2))
    tpa = do.call("rbind",lapply(unique(stab.all$variable),function(v){
      do.call("rbind",lapply(unique(stab.all$type),function(t){
        do.call("rbind",lapply(unique(stab.all$removeFrac),function(r){
          index.baseline.leaveout = stab.all[stab.all$variable == v & stab.all$type==t & stab.all$removeFrac==r & stab.all$it != "truth",]
          index.baseline = stab.all[stab.all$variable == v & stab.all$type==t & stab.all$removeFrac==r & stab.all$it == "truth",]
          index.baseline = index.baseline[order(index.baseline$variable,index.baseline$Year),]
          index.baseline.leaveout = index.baseline.leaveout[order(index.baseline.leaveout$variable,index.baseline.leaveout$Year),]
          rmses = do.call("rbind",lapply(1:iterations,function(it){
            ii = index.baseline.leaveout[index.baseline.leaveout$it==it,]
            ii = ii[order(ii$variable,ii$Year),]
            sqrt(mean((ii$index-index.baseline$index)^2))
          }))
          r2s = do.call("rbind",lapply(1:iterations,function(it){
            ii = index.baseline.leaveout[index.baseline.leaveout$it==it,]
            ii = ii[order(ii$variable,ii$Year),]
            cor(ii$index,index.baseline$index)^2
          }))
          data.frame(variable=v,type=t,removeFrac=r,scale=sd(index.baseline$index),
                     rmses=rmses,
                     r2s=r2s,it=1:10)
        }))
      }))
    }))
    #now mean over all ages, then mean and conf int again
    tpa %>% group_by(removeFrac,type,it) %>% 
      summarise(rmse=mean(rmses)/mean(scale),r2=mean(r2s)) %>%
      ungroup() %>%
      group_by(removeFrac,type) %>%
      summarise(r2m=mean(r2),r2l=t.test(r2)$conf.int[1],r2u=t.test(r2)$conf.int[2],
                rmsem=mean(rmse),rmsel=t.test(rmse)$conf.int[1],rmseu=t.test(rmse)$conf.int[2])-> tp2
    
    p = ggplot(tp2)+
      geom_line(aes(x=removeFrac,y=r2m,color=type))+
      geom_point(aes(x=removeFrac,y=r2m,color=type))+
      #geom_ribbon(aes(x=removeFrac,ymin = r2l,ymax = r2u, fill = type), alpha = 0.1, show.legend = F)+
      labs(x="p",y="R2",color="Model")+
      theme_light()+
      scale_color_manual(values=group.colors)+
      theme(legend.position = "none")
    p
    
    ggsave(paste0("figures/index/index-r2-",type,".png"),p,width=6,height=4)
    
    
    p = ggplot(tp2)+
      geom_line(aes(x=removeFrac,y=rmsem,color=type))+
      geom_point(aes(x=removeFrac,y=rmsem,color=type))+
      #geom_ribbon(aes(x=removeFrac,ymin = rmsel,ymax = rmseu, fill = type), alpha = 0.1, show.legend = F)+
      labs(x="p",y="NRMSE",color="Model")+
      theme_light()+
      scale_color_manual(values=group.colors)+
      theme(legend.position = "none")
    
    p
    
    ggsave(paste0("figures/index/index-rmse-",type,".png"),p,width=6,height=4)
    
    #plot legend
    p = ggplot(tp2)+
      geom_line(aes(x=removeFrac,y=rmsem,color=type))+
      geom_point(aes(x=removeFrac,y=rmsem,color=type))+
      #geom_ribbon(aes(x=removeFrac,ymin = rmsel,ymax = rmseu, fill = type), alpha = 0.1, show.legend = F)+
      labs(x="p",y="NRMSE",color="Model")+
      theme_light()+
      scale_color_manual(values=group.colors)+
      theme(legend.position="bottom")
    p
    
    leg = get_legend(p)
    as_ggplot(leg)
    
    
    ggsave("figures/index/legend_index.png",leg,width=8,height=0.6)
    
    
  }
  
}



evalIndexStability = function(data,sds,encodeFunction,traintestFunction,removeFrac,skipyears = 2016){
  
  index.baseline = computeIndex(data,sds,encodeFunction,traintestFunction,skipyears = skipyears)
  
  iterations = 10
  
  index.baseline.leaveout = do.call("rbind",lapply(1:iterations,function(it){
    print(it)
    data.shuff = data[sample(1:nrow(data),nrow(data),replace = FALSE),]
    ind  = computeIndex.leaveOut(data.shuff,sds,encodeFunction,traintestFunction,skipyears = skipyears,removeFrac = removeFrac)
    ind$it = it
    ind
  }))

index.baseline$it  = "truth"
index.baseline.leaveout$it = as.character(index.baseline.leaveout$it)
index.baseline.all = rbind(index.baseline,index.baseline.leaveout)

}



#interpolate, then index
computeIndex = function(data,sds,encodeFunction,traintestFunction,skipyears=2016){
  
  years = unique(data$Year)
  years = years[!(years %in% skipyears)]
  
  res = doImputation(data,sds,encodeFunction,traintestFunction,years=years)
  
  res %>%
    group_by(variable,Year) %>%
    summarise(index = sum(value)) -> ret
  
}





computeIndex.leaveOut = function(data.agewise,sds,encodeFunction,traintestFunction,skipyears=2016,removeFrac = 0.5){
  
 
  ys = unique(data.agewise$Year)
  ys = ys[!(ys %in% skipyears)] 
  
  #iterate over years
  impInd = do.call("rbind",lapply(ys,function(y){
    
    otheryears = data.agewise[data.agewise$Year != y,]
    thisyear = data.agewise[data.agewise$Year == y,]
    
    # for each sd, select which rects to keep
    thisyear.train =do.call("rbind",lapply(unique(data.agewise$Sub_Div),function(sd){
      
      ty.sd = thisyear[thisyear$Sub_Div==sd,]
      rects.sd = unique(ty.sd$RECT)
      
      #number of train samples: at least one, at most n-1, round to removeFrac
      ntrain = round(length(rects.sd) * (1-removeFrac))
      if(ntrain == 0){
        ntrain = 1
      }
      if(ntrain == length(rects.sd)){
        ntrain = length(rects.sd) - 1
      }
      
      keeprects = rects.sd[1:ntrain]
      tr = ty.sd[ty.sd$RECT %in% keeprects,]
    }))
    
    train.all = rbind(thisyear.train,otheryears)
    

    sy = ys
    sy = sy[sy!=y]
    sy = c(sy,skipyears)
    impIndAll = computeIndex(train.all,sds,encodeFunction,traintestFunction,sy)
    impIndAll[impIndAll$Year == y,]
    
    
  }))
  
}
  
