setwd("~/projects/OTC_Thuenen/Filling-the-gap/")

library(grid)
library(gridExtra) 

source("code/models.R")
source("code/utils.R")
source("eval/02-extrapolation.R")

group.colors <- c( "GAM" = "#08519c", #"GAM-year" = "#032c57","GAM-noM" = "#6baed6",
                   "LMM-noFE" = "#fd8d3c", "LMM-noSD" = "#d1117b","LMM-full" = "#a63603",
                   "XGB-noSD" = "#a9f5ab", "XGB-SD" = "#018c39", 
                   "Baseline" = "#121212")


plotExtrapolationExampleMap = function(){
  
  type = "bias"
  data = loadData(type)
  
  age = "Age2"
  year= 2019
  
  ###############
  # project coordinates to local grid
  sf_wgs84 <- st_as_sf(data, coords = c("WEST","SOUTH"), crs = 4326) #correct order 
  sf_etrs89_projected <- st_transform(sf_wgs84, crs = 3035)
  #sf_etrs89_projected <- st_transform(sf_wgs84, crs = 32634) #UTM Zone 34N
  coords <- st_coordinates(sf_etrs89_projected)
  data$orig.south=data$SOUTH
  data$WEST = coords[,1]
  data$SOUTH = coords[,2]
  ##############
  
  data.imp = evalMuchMissing(data,function(tr,te) trainTestGAM2(tr,te,k=15),testyear = year,
                             featurefun = doNoFilter,orig.south = data$orig.south) 
  data2 = merge(data,data.imp,by=c("RECT","Year","Sub_Div","variable"),all.x = T)
  data2$value[!is.na(data2$pred)] = data2$pred[!is.na(data2$pred)]
  data2 = data2[,-c(10:11)]
  
  
  data.imp = evalMuchMissing(data,function(tr,te) traintestLMER2(tr,te),testyear = year,featurefun = doNoFilter,
                             orig.south = data$orig.south) 
  data3 = merge(data,data.imp,by=c("RECT","Year","Sub_Div","variable"),all.x = T)
  data3$value[!is.na(data3$pred)] = data3$pred[!is.na(data3$pred)]
  data3 = data3[,-c(10:11)]
  
  
  data.imp = evalMuchMissing(data,function(tr,te) traintestXGB(tr,te),testyear = year,
                             featurefun = doTargetEncoding.TrainTest.vy.vr,orig.south = data$orig.south) 
  data4 = merge(data,data.imp,by=c("RECT","Year","Sub_Div","variable"),all.x = T)
  data4$value[!is.na(data4$pred)] = data4$pred[!is.na(data4$pred)]
  data4 = data4[,-c(10:11)]
  
  #data$wasImp = FALSE

  data$model = "Truth"
  data2$model = "GAM"
  data3$model = "LMM-noSD"
  data4$model = "XGB-noSD"
  data.imp = rbind(data,data2,data3,data4)
  data.imp$wasImp = FALSE
  data.imp$wasImp[data.imp$orig.south<57.5] = TRUE
  data.imp$wasImp[data.imp$model == "Truth"] = FALSE
  
  data.imp$model = factor(data.imp$model,levels=c("Truth","LMM-noSD","XGB-noSD","GAM"))
  
  
  worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                           returnclass = 'sf')
  

  dd.coord = data.imp[data.imp$Year==year,]
  dd.coord = dd.coord[dd.coord$variable == age,]
  
  #merge geometry
  ices <- st_read("data/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp")
  ices.b = dplyr::filter(ices,Ecoregion=="Baltic Sea")
  dd = merge(ices.b,dd.coord,by.x="ICESNAME",by.y="RECT")
  
  dd<- dd %>%
    arrange(wasImp)
  
  maxv=1500
  dd$value[dd$value<0] = 0
  dd$value[dd$value>maxv] = maxv
  
  p = ggplot(dd)+
    geom_sf(aes(fill=value, color = wasImp), lwd = 0.75)+
    geom_sf(data = worldmap)+
    coord_sf(xlim = c(12, 24), ylim = c(54, 60))+
    theme_light()+
    scale_color_manual(values=c("transparent","red"))+
    scale_fill_viridis(limits=c(0,maxv))+ #
    ggtitle(paste0(year,", ",age))+
    theme(axis.text.x = element_text(angle = 90))+
    labs(fill="abundance")+
    facet_wrap(.~model,nrow = 1)
  p
  
  ggsave(paste0("figures/extrapolate_map_",type,"_",year,"_",age,".png"),p,width=12,height=3.5)
  
}


plotIndices = function(){
  
  
  data = loadData("bass") 
  
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
  
  
  data.imp.gam = doImputation(data,c("24","25","26","28_2"),doNoFilter,trainTestGAM2,years=2001:2021)
  data.imp.gam$model="GAM-M"
  
  data.imp.baseline = doImputation(data,c("24","25","26","28_2"),doNoFilter,traintestBaseline,years=c(2001:2015,2017:2021))
  data.imp.baseline$model="Baseline"
  
  data.imp.all = rbind(data.imp.gam,data.imp.baseline)
  
  data.imp.all %>%
    group_by(Year,model) %>%
    summarise(index = sum(value)) -> ind      
  ind = rbind(ind,data.frame(Year=2016,model="Baseline",index=NA))
  
  ind$model[ind$model=="GAM-M"] = "GAM"
  
  ind$model=factor(ind$model,levels=c("GAM","Baseline"))
  
  p= ggplot(ind)+
    geom_line(aes(x=Year,y=index,color=model,linetype=model))+
    geom_point(aes(x=Year,y=index,color=model))+
    labs(x="year",y="abundance")+
    theme_light()+
    scale_color_manual(values=group.colors)
  p
  ggsave("figures/repaired_bass_both.png",p,width=6,height=4)
  
  
  data.imp.gam %>%
    group_by(Year,variable,model) %>%
    summarise(index = sum(value)) -> ind.age   
  
  p = ggplot(ind.age)+
    geom_point(aes(x=Year,y=variable,size=index),color="#08519c")+ 
    scale_size_area("abundance",max_size=10)+
    theme_light()+
    theme(axis.text.x = element_text(angle=90),legend.position = "right")+
    labs(y="")
  p
  
  ggsave("figures/repaired_bass_bubble.png",p,width=8,height=4)
  
}




plotInterpolationExamples = function(){
  
  #plot all models for two years: 2016 (much missing) and ...
  
  type = "bass"
  data = loadData(type)
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
  
  sds = c("24","25","26","28_2") #bass
  
  data = data[data$Sub_Div %in% sds,]
  
  ####################
  #2016
  ####################
  yy = 2016
  
  #data.imp.baseline = doImputation(data,sds,doNoFilter,traintestBaseline,years=yy)
  #data.imp.baseline$model = "Baseline"
  #data.imp.lmm1 = doImputation(data,sds,doNoFilter,traintestLMER1,years=yy)
  #data.imp.lmm1$model = "LMM-noFE"
  data.imp.lmm2 = doImputation(data,sds,doNoFilter,traintestLMER2,years=yy)
  data.imp.lmm2$model = "LMM-noSD"
  #data.imp.lmm3 = doImputation(data,sds,doNoFilter,traintestLMER3,years=yy)
  #data.imp.lmm3$model = "LMM-full"
  data.imp.gam = doImputation(data,sds,doNoFilter,trainTestGAM2,years=yy)
  data.imp.gam$model = "GAM"
  data.imp.xgb1 = doImputation(data,sds,doTargetEncoding.TrainTest.vy.vr,traintestXGB,years=yy)
  data.imp.xgb1$model = "XGB-noSDInteraction"
  #data.imp.xgb2 = doImputation(data,sds,doTargetEncoding.TrainTest.vsy.vr,traintestXGB,years=yy)
  #data.imp.xgb2$model = "XGB-SDInteraction"
  
  
  data.imp = rbind(data.imp.lmm2,
                   data.imp.gam,data.imp.xgb1)
  
  
  worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                           returnclass = 'sf')
  
  
  dd.coord = data.imp[data.imp$Year==yy,]
  
  dd.coord = dd.coord %>% group_by(RECT,SOUTH,WEST,Year,Sub_Div,model) %>% 
    summarise(value=sum(value),wasImp=all(wasImp))
  
  #dd.coord = dd.coord[dd.coord$variable == age,]
  
  #merge geometry
  ices <- st_read("data/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp")
  ices.b = dplyr::filter(ices,Ecoregion=="Baltic Sea")
  dd = merge(ices.b,dd.coord,by.x="ICESNAME",by.y="RECT")
  
  dd<- dd %>%
    arrange(wasImp)
  
  #maxv=1500
  #dd$value[dd$value<0] = 0
  #dd$value[dd$value>maxv] = maxv
  
  p = ggplot(dd)+
    geom_sf(aes(fill=value, color = wasImp), lwd = 0.75)+
    geom_sf(data = worldmap)+
    coord_sf(xlim = c(12, 24), ylim = c(54, 60))+
    theme_light()+
    scale_color_manual(values=c("transparent","red"))+
    #scale_fill_viridis(limits=c(0,maxv))+ #
    scale_fill_viridis()+
    ggtitle(yy)+
    theme(axis.text.x = element_text(angle = 90))+
    labs(fill="abundance")+
    facet_wrap(.~model,nrow = 1)
  p
  
  ggsave("figures/index/maps_all_2016.png",p,width=12,height=4)
  
  
  ####################
  #2020
  ####################
  yy=2021
  
  data.imp.baseline = doImputation(data,sds,doNoFilter,traintestBaseline,years=yy)
  data.imp.baseline$model = "Baseline"
  data.imp.lmm1 = doImputation(data,sds,doNoFilter,traintestLMER1,years=yy)
  data.imp.lmm1$model = "LMM-noFE"
  data.imp.lmm2 = doImputation(data,sds,doNoFilter,traintestLMER2,years=yy)
  data.imp.lmm2$model = "LMM-noSD"
  data.imp.lmm3 = doImputation(data,sds,doNoFilter,traintestLMER3,years=yy)
  data.imp.lmm3$model = "LMM-full"
  data.imp.gam = doImputation(data,sds,doNoFilter,trainTestGAM2,years=yy)
  data.imp.gam$model = "GAM"
  data.imp.xgb1 = doImputation(data,sds,doTargetEncoding.TrainTest.vy.vr,traintestXGB,years=yy)
  data.imp.xgb1$model = "XGB-noSDInteraction"
  data.imp.xgb2 = doImputation(data,sds,doTargetEncoding.TrainTest.vsy.vr,traintestXGB,years=yy)
  data.imp.xgb2$model = "XGB-SDInteraction"
  
  
  data.imp = rbind(data.imp.baseline,data.imp.lmm1,data.imp.lmm2,data.imp.lmm3,
                   data.imp.gam,data.imp.xgb1,data.imp.xgb2)
  
  
  worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                           returnclass = 'sf')
  
  
  dd.coord = data.imp[data.imp$Year==yy,]
  
  dd.coord = dd.coord %>% group_by(RECT,SOUTH,WEST,Year,Sub_Div,model) %>% 
    summarise(value=sum(value),wasImp=all(wasImp))
  
  #dd.coord = dd.coord[dd.coord$variable == age,]
  
  #merge geometry
  ices <- st_read("data/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp")
  ices.b = dplyr::filter(ices,Ecoregion=="Baltic Sea")
  dd = merge(ices.b,dd.coord,by.x="ICESNAME",by.y="RECT")
  
  dd<- dd %>%
    arrange(wasImp)
  
  #maxv=1500
  #dd$value[dd$value<0] = 0
  #dd$value[dd$value>maxv] = maxv
  
  p = ggplot(dd)+
    geom_sf(aes(fill=value, color = wasImp), lwd = 0.75)+
    geom_sf(data = worldmap)+
    coord_sf(xlim = c(12, 24), ylim = c(54, 60))+
    theme_light()+
    scale_color_manual(values=c("transparent","red"))+
    #scale_fill_viridis(limits=c(0,maxv))+ #
    scale_fill_viridis()+
    ggtitle(yy)+
    theme(axis.text.x = element_text(angle = 90))+
    labs(fill="abundance")+
    facet_wrap(.~model,nrow = 3)
  p
  
  ggsave("figures/index/maps_all_2021.png",p,width=12,height=10)
}



plot2016 = function(){
  
  type = "bass"
  data = loadData(type)
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
  
  sds = c("24","25","26","28_2") #bass
  
  data = data[data$Sub_Div %in% sds,]
  

  data.imp = doImputation(data,sds,doNoFilter,trainTestGAM2,years=2016)
  
  p = plotMap.imputed(data.imp,2016,"all","GAM")
  p
  
  ggsave("figures/index/map2016_all.png",p,width=6,height=4)
}



plotIndexConsistencyExample = function(){
  
  type="bias"
  stab.all = read.csv(paste0("results/index_leaveOut_",type,".csv"))
  
  stab = stab.all[stab.all$removeFrac == 0.75,]
  
  stab %>%
    group_by(Year,it,removeFrac,type) %>%
    summarise(index=sum(index)) -> stab.m
  
  p = ggplot(stab) +
    geom_line(aes(x=Year,y=index,color=it))+
    #facet_grid(rows=vars(type),cols = vars(variable),scales = "free_y")
    ggh4x::facet_grid2(type~variable, scales = "free_y", independent = "y")+
    theme_light()+
    theme(axis.text = element_blank(),axis.ticks=element_blank(),legend.position = "none")
  p
  
  ggsave("figures/index/consistencyExample.png",p,width=14,height=6)
  
  
}
