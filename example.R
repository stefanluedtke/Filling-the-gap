setwd("~/projects/OTC_Thuenen/SpratSurvey-SL/code/toPublish/")

source("code/models.R")
source("code/utils.R")

set.seed(42)

#load a dataset, e.g. BASS
type = "bass"
data = loadData(type)

#only keep the subdivisions that are eventually relevant for the index
sds = c("24","25","26","28_2") #bass
#sds = c("22","23","24","25","26","27","28_2","29","32") #bias
#sds = c("25","26","27","28_2","29") #herring

data = data[data$Sub_Div %in% sds,]

#plot the data
plotMap(data,2009,"Age2")

#run an imputation model
#requires selecting a feature encoding function and a model. Here, we use the XGB2 model, which is 
#encoding: doTargetEncoding.TrainTest.vy.vr 
#model: traintestXGB
data.imp = doImputation(data,sds,doTargetEncoding.TrainTest.vy.vr,traintestXGB,years=2009)

#plot again
plotMap.imputed(data.imp,2009,"Age2")

#try the baseline imputation model instead and plot again
data.imp2 = doImputation(data,sds,doNoFilter,traintestBaseline,years=2009)
plotMap.imputed(data.imp2,2009,"Age2")


#try the GAM imputation model instead and plot again
data.imp2 = doImputation(data,sds,doNoFilter,trainTestGAM2,years=2009)
plotMap.imputed(data.imp2,2009,"Age2")


#compute index: do imputation for all years, then plot
data = loadData("bass") #load again because we use all data for training
data.imp.gam = doImputation(data,sds,doNoFilter,trainTestGAM2,years=2001:2021)
data.imp.gam$model="gam"

data.imp.baseline = doImputation(data,sds,doNoFilter,traintestBaseline,years=c(2001:2015,2017:2021))
data.imp.baseline$model="baseline"

data.imp.all = rbind(data.imp.gam,data.imp.baseline)
data.imp.all %>%
  group_by(Year,model) %>%
  summarise(index = sum(value)) -> ind      
ind = rbind(ind,data.frame(Year=2016,model="baseline",index=NA))

ggplot(ind)+
  geom_line(aes(x=Year,y=index,color=model))+
  geom_point(aes(x=Year,y=index,color=model))+
  labs(x="year",y="abundance")+
  theme_light()


