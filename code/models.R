require(xgboost)
require(lme4)
require(mgcv)

##############################################################
#Linear mixed models
##############################################################

traintestLMER = function(train,test,form = value ~ n  + (1|variable) + (1|RECT) + (1|year),allow.new.levels=FALSE){

  mod = lmer(form,data=train)
  pred = predict(mod,test,allow.new.levels=allow.new.levels)
  
  r = data.frame(pred=pred,truth=test$value)
  colnames(r) = c("pred","truth")
  r
}

traintestBaseline = function(...) traintestLMER(...,form= value ~ (Area+0|Year:Sub_Div:variable))
traintestLMER1 = function(...) traintestLMER(...,form =  value ~ (Area+0|Year:Sub_Div:variable) + (1|RECT:variable))
traintestLMER2 = function(...) traintestLMER(...,form = value ~ (Area+0|Year:variable) + (1|RECT:variable))


##############################################################
#GAMs
##############################################################

trainTestGAM1 = function(train,test,k=15){
  
  #group test samples by variable
  uyv = unique(test$variable)
  test$pred = NA
  
  
  res = do.call("rbind",lapply(uyv,function(yvi){
    #print(yvi)
    dat = train[train$variable == yvi,]
    te = test[test$variable == yvi,]
    
    gmod = gam(value ~ s(SOUTH,WEST,k=k),data=dat,method="REML")
    pred = predict(gmod,te,type="response")
    test$pred[test$variable==yvi] <<- pred
    NA
  }))
  
  res = data.frame(pred=test$pred,truth=test$value)
  
}

#https://link.springer.com/article/10.1007/s11160-023-09795-2#Tab2
trainTestGAM2 = function(train,test,k=15){
  
  #group test samples by variable
  uyv = unique(test$variable)
  test$pred = NA
  
  mobs = train %>% group_by(Year,variable) %>% summarise(mobs = mean(value),.groups = "drop")
  train = merge(train,mobs,by=c("Year","variable"))
  test = merge(test,mobs,by=c("Year","variable"))
  
  res = do.call("rbind",lapply(uyv,function(yvi){
    #print(yvi)
    dat = train[train$variable == yvi,]
    te = test[test$variable == yvi,]
    
    gmod = gam(value ~ s(SOUTH,WEST,k=k,by=mobs),data=dat,method="REML")
    pred = predict(gmod,te,type="response")
    test$pred[test$variable==yvi] <<- pred
    NA
  }))
  
  res = data.frame(pred=test$pred,truth=test$value)
  
}

trainTestGAM3 = function(train,test,k=15){
  
  #group test samples by variable
  uyv = unique(test$variable)
  test$pred = NA
  
  mobs = train %>% group_by(Year,variable) %>% summarise(mobs = mean(value),.groups = "drop")
  train = merge(train,mobs,by=c("Year","variable"))
  test = merge(test,mobs,by=c("Year","variable"))
  
  res = do.call("rbind",lapply(uyv,function(yvi){
    #print(yvi)
    dat = train[train$variable == yvi,]
    te = test[test$variable == yvi,]
    
    gmod = gam(value ~ s(SOUTH,WEST,k=k,by=mobs*Area),data=dat,method="REML")
    pred = predict(gmod,te,type="response")
    test$pred[test$variable==yvi] <<- pred
    NA
  }))
  
  res = data.frame(pred=test$pred,truth=test$value)
  
}

##############################################################
#Boosted trees
##############################################################

traintestXGB = function(train,test,max.depth=2,nrounds=30){
  
  
  train = train[ , order(colnames(train))]
  test = test[ , order(colnames(test))]
  
  train = train[,colnames(train)!="variable"]
  test = test[,colnames(test)!="variable"]
 
  
  mtrain = as.matrix(train[,colnames(train) != "value"])
  mtest = as.matrix(test[,colnames(test) != "value"])
  
  mod <- xgboost(data = mtrain, label = train$value,
                 max.depth = max.depth,  nrounds =nrounds, objective = "reg:squarederror",verbose = 0) #2/30
  
  
  pred = predict(mod,mtest)
  r = data.frame(pred=pred,truth=test$value)
  colnames(r) = c("pred","truth")
  r
}

#the different feature encodings for xgb
#############################

#encoding for XGB1
doTargetEncoding.TrainTest.vy.vr = function(data.train,data.test,removeExtra=""){
  
  data.train$var.sd.year = paste0(data.train$variable,".",data.train$Year)
  data.test$var.sd.year = paste0(data.test$variable,".",data.test$Year)
  
  dvar = encode_target_getD(data.train$var.sd.year,data.train$value)
  data.train$viy_enc = encode_target_makeTargetEnc(data.train$var.sd.year,dvar)
  data.test$viy_enc = encode_target_makeTargetEnc(data.test$var.sd.year,dvar)
  
  data.train$var.rect = paste0(data.train$variable,".",data.train$RECT)
  data.test$var.rect = paste0(data.test$variable,".",data.test$RECT)
  #print(data.test$var.rect[(!data.test$var.rect %in% data.train$var.rect )])
  
  dvar = encode_target_getD(data.train$var.rect,data.train$value)
  data.train$vr_enc = encode_target_makeTargetEnc(data.train$var.rect,dvar)
  data.test$vr_enc = encode_target_makeTargetEnc(data.test$var.rect,dvar)
  
  #feats = target.enc[,c(2:71,77:79)]
  feats.train = data.train[, (colnames(data.train) != "RECT")
                           & (colnames(data.train) != "Year") &(colnames(data.train) != "var.rect")
                           &(colnames(data.train) != "Sub_Div") &(colnames(data.train) != "var.sd.year")] 
  feats.test = data.test[, (colnames(data.test) != "RECT")
                         & (colnames(data.test) != "Year") &(colnames(data.test) != "var.rect")
                         &(colnames(data.test) != "Sub_Div") &(colnames(data.test) != "var.sd.year")]
  
  feats.train = feats.train[!(colnames(feats.train) %in% removeExtra)]
  feats.test = feats.test[!(colnames(feats.test) %in% removeExtra)]
  
  list(feats.train,feats.test)
}

#encoding for XGB1
doTargetEncoding.TrainTest.vsy.vr = function(data.train,data.test,removeExtra=""){
  
  data.train$var.sd.year = paste0(data.train$variable,".",data.train$Sub_Div,".",data.train$Year)
  data.test$var.sd.year = paste0(data.test$variable,".",data.test$Sub_Div,".",data.test$Year)
  
  dvar = encode_target_getD(data.train$var.sd.year,data.train$value)
  data.train$viy_enc = encode_target_makeTargetEnc(data.train$var.sd.year,dvar)
  data.test$viy_enc = encode_target_makeTargetEnc(data.test$var.sd.year,dvar)
  
  data.train$var.rect = paste0(data.train$variable,".",data.train$RECT)
  data.test$var.rect = paste0(data.test$variable,".",data.test$RECT)
  
  dvar = encode_target_getD(data.train$var.rect,data.train$value)
  data.train$vr_enc = encode_target_makeTargetEnc(data.train$var.rect,dvar)
  data.test$vr_enc = encode_target_makeTargetEnc(data.test$var.rect,dvar)
  
  #feats = target.enc[,c(2:71,77:79)]
  feats.train = data.train[, (colnames(data.train) != "RECT")
                           & (colnames(data.train) != "Year") &(colnames(data.train) != "var.rect")
                           &(colnames(data.train) != "Sub_Div") &(colnames(data.train) != "var.sd.year")] 
  feats.test = data.test[, (colnames(data.test) != "RECT")
                         & (colnames(data.test) != "Year") &(colnames(data.train) != "var.rect")
                         &(colnames(data.test) != "Sub_Div") &(colnames(data.test) != "var.sd.year")]
  
  feats.train = feats.train[!(colnames(feats.train) %in% removeExtra)]
  feats.test = feats.test[!(colnames(feats.test) %in% removeExtra)]
  
  list(feats.train,feats.test)
}

encode_target_getD <- function(x, y) {
  d <- aggregate(y, list(factor(x, exclude = NULL)), mean, na.rm = TRUE)
}

encode_target_makeTargetEnc <- function(x, d) {
  m <- d[is.na(as.character(d[, 1])), 2]
  l <- d[, 2]
  names(l) <- d[, 1]
  l <- l[x]
  l[is.na(l)] <- m
  l
}



traintestCatboost =  function(train,test){
  
  areatrain = train$Area
  areatest = test$Area
  train$value = train$value/train$Area
  test$value = test$value/test$Area
  train = train[,colnames(train) != "Area"]
  test = test[,colnames(test) != "Area"]
  
  X_train = train[,colnames(train) != "value"]
  X_test = test[,colnames(test) != "value"]
  y_train = train$value
  y_test = test$value
  

  
  train_pool <- catboost.load_pool(data = X_train, label = y_train)
  test_pool = catboost.load_pool(data = X_test, label = y_test) #y_test

  model <- catboost.train(train_pool,  NULL,
                          params = list(loss_function = 'RMSE', # learning_rate = 0.3,
                                        iterations = 100,depth=4, metric_period=10,logging_level="Silent",
                                        use_best_model=FALSE))
  
  prediction <- catboost.predict(model, test_pool)
  
  r = data.frame(pred=prediction*areatest,truth=test$value*areatest)
  colnames(r) = c("pred","truth")
  r
}


