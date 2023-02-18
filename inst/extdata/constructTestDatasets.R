# R code for creating the test datasets for package CRTspillover
library(CRTspillover)

extdata <- system.file("extdata",package = 'CRTspillover')

# the test code is in tests/testthat
# the test data is in inst/extdata

# Coordinates from the Avecnet trial
#load("C:/git_repos/trialdesign/RCode_DesignCRTs/ForPublication/Simulated_Trials.RData")
#AvecNet_coordinates= trial[[7]][,c(1,2)]
#write.csv(AvecNet_coordinates,file=paste0(extdata,'AvecNet_coordinates.csv')

setwd('c:/git_repos/CRTspillover/inst/extdata')

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
    test_Clusters1 <- DefineClusters(test_Locationsxy,h = 50)
    test_Arms1 <- Randomize_CRT(trial = test_Clusters,matchedPair = TRUE)
    test_Simulate_CRT <- Simulate_CRT(trial = test_Arms1,
                                      theta_inp = 1.2,initialPrevalence = 0.4,
                                      ICC_inp = 0.05,efficacy = 0.4,tol = 0.05)
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
  test_Estimates <- Analyse_CRT(trial = trial,
                               method = 'GEE',excludeBuffer = FALSE,
                               requireBootstrap = FALSE,alpha = 0.2)
  # remove the model object before storing as text as this cannot be reloaded and is large
  test_Estimates$model.object <- NULL
  dput(test_Estimates, file = 'test_Analyse_CRT.txt')
  value <- test_Estimates$contamination$data$positives[4]
  return(value)
}
# return value of test is an integer but the output file may be required for other purposes
test_Analyse_CRT = dget(file = 'test_Analyse_CRT.txt')

# TEST 4 DOES NOT REQUIRE ANY ADDITIONAL FILES

set.seed(1234)
get_test5 = function(extdata){
  Solarmal_baseline <- read.csv(file = paste0(extdata,"/Solarmal_baseline.csv"))
  test_LocationsLatLong <- Solarmal_baseline[, c('lat','long')]
  test_Locationsxy <- Convert_LatLong(test_LocationsLatLong) #test_site is simulated
  test_AnonymisedLocations <- Anonymise_TrialSite(test_Locationsxy)
  test_Clusters <- DefineClusters(test_AnonymisedLocations,h = 50)
  test_Arms <- Randomize_CRT(trial = test_Clusters,matchedPair = FALSE)
  test_Buffer <- Specify_CRTbuffer(trial = test_Arms, bufferWidth = 0.1)
  test_Buffer$cluster <- as.numeric(test_Buffer$cluster)
  test_Buffer$arm <- as.character(test_Buffer$arm)
  return(test_Buffer)}

test_Buffer <- get_test5(extdata)
write.csv(test_Buffer,file='test_Buffer.csv',row.names=FALSE)

# TESTS 6 - 9 USE THE SAME DATASETS AS EARLIER TESTS
