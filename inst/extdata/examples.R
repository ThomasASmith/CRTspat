test1 = read.csv(file='data-raw/BagamoyoTrial.csv')
test2=Convert_LatLong(test1)
test3=Anonymise_TrialSite(test2)
test4=Randomize_CRT(test2)
test5=Simulate_CRT(trial=test4,efficacy=0.4,initialPrevalence=0.5,sd=0.4)
GEEresults=Analyse_CRT(trial=test5,method='GEE')
GEEresults$PointEstimates$efficacy
GEEresults$IntervalEstimates$efficacy
GEEresults$PointEstimates$ICC
MCMCresults=Analyse_CRT(TD=test5,method='MCMC03')
MCMCresults$PointEstimates$efficacy
MCMCresults$IntervalEstimates$efficacy
MCMCresults$PointEstimates$contaminationRange
MCMCresults$IntervalEstimates$contaminationRange

test6=Specify_CRTbuffer(trial = test5, bufferWidth = 0.5)
trial=test_Simulate_CRT
method='GEE'
requireBootstrap=FALSE
confidence_level=0.2
excludeBuffer=FALSE

test1 = read.csv(file='data-raw/Solarmal_baseline.csv')
test2=Convert_LatLong(df=test1)
test2$num=test2$RDT_test_result
test2=with(test2,test2[,c('x','y','num')])
test3=Aggregate_CRT(test2)
Plot_CRTmap(test3,maskbuffer=0.5)
test4 = DefineClusters(trial=test3, h=50, algo='NN')
Plot_CRTmap(test4,maskbuffer=0.5)
test5 = Randomize_CRT(test4)
Plot_CRTmap(test5,maskbuffer=0.5)

ICC=c()
clusterSize = seq(1:10)*20
for(size in clusterSize){
  test4a=DefineClusters(trial=test3, h=size, algo='NN')
  test4b= Randomize_CRT(test4a)
  baselineAnalysis = calculateICC_baseline(trial=test4b)
  ICC=c(ICC,baselineAnalysis$estimates$ICC)
}
plot(clusterSize,ICC)

estimates=methods=c()
all_methods = c('piecewise_linear','logit','sigmoid','empirical',
                'GEE','MCMC01','MCMC02','MCMC03')
for(method in all_methods){
  print(method)
  res = Analyse_CRT(trial=test_Simulate_CRT,method=method)
  estimates = c(estimates,res$PointEstimates$efficacy)
}
names(estimates)=all_methods
estimates


library(ggplot2)
ggplot2::ggplot(data=order_df(test4), aes(y=y, x=x, colour=cluster)) + geom_point()+
  scale_color_gradientn(colours = rainbow(6))

order_df= function(df){
  #calculate squared distance to corner
  #distmin = with(df,(x - min(x))^2 + (y - min(y))^2)
  df= df[order(df$x),]
  #df= df[order(distmin),]
return(df)}

junk <- test1[1:20,]

ggplot(data=junk, aes(x=long, y=lat)) +
  geom_point(size=4,color='white') +
  geom_point(size=1,alpha=1) +
theme(panel.background = element_rect(fill = "grey"))

