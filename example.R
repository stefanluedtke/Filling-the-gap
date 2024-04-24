setwd("~/projects/OTC_Thuenen/SpratSurvey-SL/code/toPublish/")

source("code/models.R")
source("code/utils.R")

#load the BIAS dataset
data = loadData("bias")

#some of the analyses and models cannot handle rectangles that appear only once in the dataset. 
#This function removes them (they don't matter for the index anyways, as they are not in the subdivisions used for computing the index)
data = preprocData.simple(data)

#plot the data
plotMap(data,2009,"Age2")

#run an imputation model
data.imp = doImputation(data,c("24","25","26","28_2"),doTargetEncoding.TrainTest.vy.vr,traintestXGB,years=2009)

#plot again
plotMap.imputed(data.imp,2009,"Age2")

#compute index
data.imp %>%
  group_by(Year,variable) %>%
  summarise(index = sum(value)) -> ind

ggplot(ind)+geom_line(aes(x=Year,y=index))+facet_wrap(.~variable)

#try the baseline imputation model instead and plot again
data.imp2 = doImputation(data,c("24","25","26","28_2"),doNoFilter,traintestBaseline,years=2009)
plotMap.imputed(data.imp2,2009,"Age2")

