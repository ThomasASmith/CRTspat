test1 = read.csv(file='data-raw/BagamoyoTrial.csv')
test2=Convert_LatLong(test1)
test3=Anonymise_TrialSite(trial=test2)
test4=Randomize_CRT(test2)
test5=Simulate_CRT(trial=test4,efficacy=0.4,initialPrevalence=0.5,sd=0.4)
GEEresults=Analyse_CRT(trial=test5,method='GEE')
GEEresults$PointEstimates$efficacy
GEEresults$IntervalEstimates$efficacy
GEEresults$PointEstimates$ICC
INLAresults=Analyse_CRT(trial=test5,method='LR')
INLAresults$PointEstimates$efficacy
INLAresults$IntervalEstimates$efficacy
INLAresults$PointEstimates$contamination$contaminationRange

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


