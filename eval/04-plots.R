setwd("~/projects/OTC_Thuenen/SpratSurvey-SL/code/toPublish/")

source("code/models.R")
source("code/utils.R")
source("eval/02-extrapolation.R")


plotExtrapolationExampleMap = function(){
  
  type = "bias"
  data = loadData(type)
  
  age = "Age2"
  year= 2019
  
  data.imp = evalMuchMissing(data,function(tr,te) trainTestGAM3(tr,te,k=15),testyear = year,featurefun = doNoFilter) 
  data2 = merge(data,data.imp,by=c("RECT","Year","Sub_Div","variable"),all.x = T)
  data2$value[!is.na(data2$pred)] = data2$pred[!is.na(data2$pred)]
  data2 = data2[,-c(9:10)]
  
  
  data.imp = evalMuchMissing(data,function(tr,te) traintestLMER2(tr,te),testyear = year,featurefun = doNoFilter) 
  data3 = merge(data,data.imp,by=c("RECT","Year","Sub_Div","variable"),all.x = T)
  data3$value[!is.na(data3$pred)] = data3$pred[!is.na(data3$pred)]
  data3 = data3[,-c(9:10)]
  
  
  data.imp = evalMuchMissing(data,function(tr,te) traintestXGB(tr,te),testyear = year,featurefun = doTargetEncoding.TrainTest.vy.vr) 
  data4 = merge(data,data.imp,by=c("RECT","Year","Sub_Div","variable"),all.x = T)
  data4$value[!is.na(data4$pred)] = data4$pred[!is.na(data4$pred)]
  data4 = data4[,-c(9:10)]
  
  #data$wasImp = FALSE

  data$model = "Truth"
  data2$model = "GAM3"
  data3$model = "LMM2"
  data4$model = "XGB1"
  data.imp = rbind(data,data2,data3,data4)
  data.imp$wasImp = FALSE
  data.imp$wasImp[data.imp$SOUTH<57.5] = TRUE
  data.imp$wasImp[data.imp$model == "Truth"] = FALSE
  
  data.imp$model = factor(data.imp$model,levels=c("Truth","LMM2","XGB1","GAM3"))
  
  
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
  
  
  data = loadData("bass") #load again because we use all data for training
  data.imp.gam = doImputation(data,c("24","25","26","28_2"),doNoFilter,trainTestGAM3,years=2001:2021)
  data.imp.gam$model="gam"
  
  data.imp.baseline = doImputation(data,c("24","25","26","28_2"),doNoFilter,traintestBaseline,years=c(2001:2015,2017:2021))
  data.imp.baseline$model="baseline"
  
  data.imp.all = rbind(data.imp.gam,data.imp.baseline)
  data.imp.all %>%
    group_by(Year,model) %>%
    summarise(index = sum(value)) -> ind      
  ind = rbind(ind,data.frame(Year=2016,model="baseline",index=NA))
  
  
  p= ggplot(ind)+
    geom_line(aes(x=Year,y=index,color=model))+
    geom_point(aes(x=Year,y=index,color=model))+
    labs(x="year",y="abundance")+
    theme_light()
  p
  ggsave("figures/repaired_bass_both.png",p,width=6,height=4)
  
  
  data.imp.gam %>%
    group_by(Year,variable) %>%
    summarise(index = sum(value)) -> ind.age   
  
  p = ggplot(ind.age)+
    geom_point(aes(x=Year,y=variable,size=index))+ 
    scale_size_area(max_size=10)+
    theme_light()+
    theme(axis.text.x = element_text(angle=90),legend.position = "none")+
    labs(y="")
  p
  
  ggsave("figures/repaired_bass_bubble.png",p,width=6,height=4)
  
}



plot2016 = function(){
  
  type = "bass"
  data = loadData(type)
  sds = c("24","25","26","28_2") #bass
  
  data = data[data$Sub_Div %in% sds,]
  

  data.imp = doImputation(data,sds,doNoFilter,trainTestGAM3,years=2016)
  
  p = plotMap.imputed(data.imp,2016,"Age2")
  
  
  ggsave("figures/index/map2016age2.png",p,width=6,height=4)
}



plotIndexConsistencyExample = function(){
  
  type="bias"
  stab.all = read.csv(paste0("results/index_leaveOut_",type,".csv"))
  
  # stab = stab.all[stab.all$type %in% c("Baseline","GAM3") & stab.all$removeFrac == 0.75,]
  stab = stab.all[stab.all$removeFrac == 0.75,]
  
  stab %>%
    group_by(Year,it,removeFrac,type) %>%
    summarise(index=sum(index)) -> stab.m
  
  #stab.m = stab[stab$variable=="Age2",]
  
  p = ggplot(stab) +
    geom_line(aes(x=Year,y=index,color=it))+
    #facet_grid(rows=vars(type),cols = vars(variable),scales = "free_y")
    ggh4x::facet_grid2(type~variable, scales = "free_y", independent = "y")+
    theme(axis.text = element_blank(),axis.ticks=element_blank(),legend.position = "none")
  p
  
  ggsave("figures/index/consistencyExample.png",p,width=14,height=6)
  
  
}
