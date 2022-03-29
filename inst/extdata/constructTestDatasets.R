# R code for creating the test datasets for package CRTspillover

# Coordinates from the Avecnet trial
load("C:/git_repos/trialdesign/RCode_DesignCRTs/ForPublication/Simulated_Trials.RData")
AvecNet_coordinates= trial[[7]][,c(1,2)]
save(AvecNet_coordinates,file='data/AvecNet_coordinates.RData')

# Test of Design_CRT
set.seed(1234)
test_AvecNet = Design_CRT()
save(test_AvecNet,file='data/test_AvecNet.RData')

# Test of Simulate_TrialSite()
set.seed(1234)
test_site = Simulate_TrialSite()
save(test_site,file='data/test_site.RData')

# Test of Simulate_CRT()
set.seed(1234)
test_Simulate_CRT =  Simulate_CRT()
save(test_Simulate_CRT,file='data/test_Simulate_CRT.RData')

#to create test results for Analyse_CRT():
set.seed(1234)
test_Analyse_CRT=Analyse_CRT(trial=test_Simulate_CRT,method='GEE',
                        requireBootstrap=FALSE,
                        alpha=0.2)
save(test_Analyse_CRT,file='data/test_Analyse_CRT.RData')

#to create test results for testing chain of construction of trials,
set.seed(1234)
Solarmal_baseline = read.csv(file='inst/extdata/Solarmal_baseline.csv')
save(Solarmal_baseline,file='data/Solarmal_baseline.RData')
testLocationsLatLong = Solarmal_baseline[, c('lat','long')]
testLocationsxy=Convert_LatLong(testLocationsLatLong) #test_site is simulated
testAnonymisedLocations=Anonymise_TrialSite(testLocationsxy)
testClusters=DefineClusters(testAnonymisedLocations)
save(testClusters,file='data/testClusters.RData')
testArms=Randomize_CRT(testClusters)
save(testArms,file='data/testArms.RData')
testBuffer=Specify_CRTbuffer(trial = testArms, bufferWidth = 0.1)
save(testBuffer,file='data/testBuffer.RData')
testOutcomes=Simulate_CRT(trial=testBuffer,efficacy=0.4,initialPrevalence=0.5,sd=0.4)
save(testOutcomes,file='data/testOutcomes.RData')

# compress files
tools::resaveRdaFiles('data')
