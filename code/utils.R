
require(sf)
require(dplyr)
require(reshape2)
require(rnaturalearth)
require(viridis)
require(ggplot2)

doNoFilter = function(train,test) list(train,test)



loadData = function(type="bias"){
  
  
  #load sprat data (to know which rectangles and years to use) and spatial data (to know rectangle def)
  ices <- st_read("data/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp")
  ices.b = dplyr::filter(ices,Ecoregion=="Baltic Sea")
  
  if(type=="bias"){
    fn0 = "data/surveys/SPR_BIAS_N_RECT.txt"
    sprat = read.table(fn0,sep=";",header = T)
    sprat = sprat[,colnames(sprat)!="SPECIES"]
  }
  if(type=="bass"){
    fn0 = "data/surveys/SPR_BASS_N_RECT.txt"
    sprat = read.table(fn0,sep=";",header = T)
  }
  if(type=="herring"){
    fn0 = "data/surveys/HER_BIAS_RECT.txt"
    sprat = read.table(fn0,sep=";",header = T,dec = ",")
    colnames(sprat)= c("Year","Sub_Div", "RECT",  "Area","Age0", "Age1", "Age2", "Age3", "Age4",
                       "Age5", "Age6", "Age7", "Age8.")
  }
  
  
  spratl = melt(sprat,id.vars=c("Year","Sub_Div","RECT","Area"))

  ices.sprat = merge(ices.b,spratl,by.x="ICESNAME",by.y="RECT")
  colnames(ices.sprat)[colnames(ices.sprat)=="ICESNAME"] = "RECT"
  ices.sprat = ices.sprat[,c("RECT","SOUTH","WEST","Year","Sub_Div","Area","variable","value")]
  ices.sprat=data.frame(ices.sprat)
  ices.sprat=ices.sprat[,colnames(ices.sprat)!="geometry"]
  ices.sprat
  
  
}



#some of the analyses and models cannot handle rectangles that appear only once in the dataset. 
#This function removes them 
preprocData.simple = function(data.agewise.pca,removeSingleRects=TRUE,removeRects = c("43G8","43G6","45G7","37G4","42G6","44G8","46G7","39G4")){
  
  if(removeSingleRects){
    data.agewise.pca %>%
      group_by(Year,Sub_Div,variable) %>%
      filter(n()>1) ->
      data.agewise.pca.correct
    nums =  table(data.agewise.pca.correct$RECT)
    rrs = names(nums)[nums<=9]
    data.agewise.pca.correct = data.agewise.pca.correct[!(data.agewise.pca.correct$RECT %in% rrs),]
    data.agewise.pca = data.frame(data.agewise.pca.correct)
  }
  
  data.agewise.pca = data.agewise.pca[!(data.agewise.pca$RECT %in% removeRects),]
  data.agewise.pca
  
}


plotMap = function(data,year,age){
  
  worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                           returnclass = 'sf')
  
  dd.coord = data[data$Year==year,]
  dd.coord = dd.coord[dd.coord$variable == age,]
  
  #merge geometry
  ices <- st_read("data/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp")
  ices.b = dplyr::filter(ices,Ecoregion=="Baltic Sea")
  dd = merge(ices.b,dd.coord,by.x="ICESNAME",by.y="RECT")

  p = ggplot(dd)+
    geom_sf(aes(fill=value),color="transparent")+
    geom_sf(data = worldmap)+
    coord_sf(xlim = c(12, 24), ylim = c(54, 60))+
    theme_light()+
    scale_fill_viridis()+ #
    ggtitle(paste0(year,", ",age))+
    theme(axis.text.x = element_text(angle = 90))+
    labs(fill="abundance")
  p
  
  
}

plotMap.imputed = function(data,year,age){
  
  worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                           returnclass = 'sf')
  
  dd.coord = data[data$Year==year,]
  dd.coord = dd.coord[dd.coord$variable == age,]
  
  #merge geometry
  ices <- st_read("data/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp")
  ices.b = dplyr::filter(ices,Ecoregion=="Baltic Sea")
  dd = merge(ices.b,dd.coord,by.x="ICESNAME",by.y="RECT")
  
  dd<- dd %>%
    arrange(wasImp)
  
  p = ggplot(dd)+
    geom_sf(aes(fill=value, color = wasImp), lwd = 0.75)+
    geom_sf(data = worldmap)+
    coord_sf(xlim = c(12, 24), ylim = c(54, 60))+
    theme_light()+
    scale_color_manual(values=c("transparent","red"))+
    scale_fill_viridis()+ #
    ggtitle(paste0(year,", ",age))+
    theme(axis.text.x = element_text(angle = 90))+
    labs(fill="abundance")
  p
  
}


doImputation = function(data,sds,encodeFunction,traintestFunction,years){
  
  data %>% distinct(Sub_Div, RECT, .keep_all = TRUE) -> su
  su = su[,c("RECT","Sub_Div","Area","SOUTH","WEST")]
  
  su = su[su$Sub_Div %in% sds,]
  su %>%
    group_by(Sub_Div) %>%
    summarise(Area = sum(Area)) -> sd.area
  
  ages=unique(data$variable)
  allrects.agewise = do.call("rbind",lapply(ages,function(a){
    
    sua = su
    sua$variable=a
    sua
    
  }))


  all = do.call("rbind",lapply(years,function(y){
    
    datay = data[data$Year == y,]
    alldaty = merge(datay,allrects.agewise,by.x=c("RECT","Sub_Div","variable","SOUTH","WEST","Area"),by.y=c("RECT","Sub_Div","variable","SOUTH","WEST","Area"),all.y=TRUE)
    
    #for each of the SDs: merge with list of all rects in that sd. impute the missing ones
    toImp = alldaty[is.na(alldaty$value),]
    notImp = alldaty[!is.na(alldaty$value),]

    if(nrow(toImp)==0){
      notImp$wasImp=FALSE
      dat.afterimp=notImp
      
    }else{
      toImp$Year = y
      # Feature encoding...
      tt = encodeFunction(data,toImp)
      tr = tt[[1]]
      te = tt[[2]]
      imp = traintestFunction(tr,te)
      toImp$value = imp$pred
      
      notImp$wasImp=FALSE
      toImp$wasImp=TRUE
      
      dat.afterimp = rbind(notImp,toImp)
    }
    
    dat.afterimp
    
  }))
}



