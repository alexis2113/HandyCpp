library(quantmod)
library(PerformanceAnalytics)
mydata=readRDS("stockdata.rds")

#define 5 days as short period, 15 days as long period

mv5=SMA(mydata$Close,5)
mv15=SMA(mydata$Close,15)

#when short-term profit higher than long-term, short the asset

signal=ifelse(mv5<mv15,0,1)
sig=lag(signal)

#calculate overall profit

roc=ROC(type='discrete',mydata$Close)
return=roc * sig

#plot overall profit

charts.PerformanceSummary(return)

