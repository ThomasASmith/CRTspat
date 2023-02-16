# R code for creating the test datasets for package CRTspillover

# Coordinates from the Avecnet trial
load("C:/git_repos/trialdesign/RCode_DesignCRTs/ForPublication/Simulated_Trials.RData")
AvecNet_coordinates= trial[[7]][,c(1,2)]
write.csv(AvecNet_coordinates,file='AvecNet_coordinates.csv')

library(CRTspillover)
# All the test files are in tests/testthat
setwd('c:/git_repos/CRTspillover/tests/testthat')

# TEST 1 Test of Design_CRT
set.seed(1234)
get_test1 = function(){
    Solarmal_baseline <- read.csv(file = "Solarmal_baseline.csv")
    testLocationsxy <- Convert_LatLong(Solarmal_baseline) #test_site is simulated
    coordinates <- Solarmal_baseline[]
    test_design <- Design_CRT(alpha = 0.05,
                              desiredPower = 0.8,
                              effect = 0.6,
                              ICC = 0.175,
                              pC = 0.4,
                              postulatedContaminationRange = 0.05,
                              coordinates = testLocationsxy,
                              h = 80)
    df <- test_design$descriptionCoreTrial$trial
    df$cluster <- as.numeric(df$cluster)
    df$arm <- as.character(df$arm)
    row.names(df) <- NULL
    return(df)
}
write.csv(get_test1(),file='test_design.csv',row.names=FALSE)

# TEST 2:  Test of Simulate_TrialSite()
set.seed(1234)
test_site <- Simulate_TrialSite()
write.csv(test_site,file='test_site.csv',row.names=FALSE)

# TEST 3: Test of Simulate_CRT()
set.seed(1234)
get_test3 = function(){
    Solarmal_baseline <- read.csv(file = "Solarmal_baseline.csv")
    testLocationsxy <- Convert_LatLong(Solarmal_baseline) #test_site is simulated
    testClusters1 <- DefineClusters(testLocationsxy,h = 50)
    testArms1 <- Randomize_CRT(trial = testClusters,matchedPair = TRUE)
    test_Simulate_CRT <- Simulate_CRT(trial = testArms1,
                                      theta_inp = 1.2,initialPrevalence = 0.4,
                                      ICC_inp = 0.05,efficacy = 0.2,tol = 0.05)
    test_Simulate_CRT$cluster <- as.numeric(test_Simulate_CRT$cluster)
    test_Simulate_CRT$arm <- as.character(test_Simulate_CRT$arm)
    rownames(test_Simulate_CRT) <- NULL
return(test_Simulate_CRT)}
test_Simulate_CRT <- get_test3()
write.csv(test_Simulate_CRT,file='test_Simulate_CRT.csv',row.names = FALSE)

# TEST 4: test of Analyse_CRT():
set.seed(1234)
get_test4 = function(){
  trial <- read.csv(file = "test_Simulate_CRT.csv")
  testEstimates <- Analyse_CRT(trial = trial,
                               method = 'GEE',excludeBuffer = FALSE,
                               requireBootstrap = FALSE,alpha = 0.2)
  value <- testEstimates$contamination$data$positives[4]
  return(value)
}
# return value is an integer

#TEST 5
set.seed(1234)
get_test5 = function(){
  Solarmal_baseline <- read.csv(file = "Solarmal_baseline.csv")
  testLocationsLatLong <- Solarmal_baseline[, c('lat','long')]
  testLocationsxy <- Convert_LatLong(testLocationsLatLong) #test_site is simulated
  testAnonymisedLocations <- Anonymise_TrialSite(testLocationsxy)
  testClusters <- DefineClusters(testAnonymisedLocations,h = 50)
  testArms <- Randomize_CRT(trial = testClusters,matchedPair = FALSE)
  testBuffer <- Specify_CRTbuffer(trial = testArms, bufferWidth = 0.1)
  testBuffer$cluster <- as.numeric(testBuffer$cluster)
  testBuffer$arm <- as.character(testBuffer$arm)
return(testBuffer)}

testBuffer <- get_test5()
write.csv(testBuffer,file='testBuffer.csv',row.names=FALSE)




# compress files
#tools::resaveRdaFiles('data')
