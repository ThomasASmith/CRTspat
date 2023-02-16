extdata <- system.file("extdata",package = 'CRTspillover')
set.seed(1234)
test_that("Toy GEE analysis creates correct output", {
   get_test4 = function(extdata){
      trial <- read.csv(file = paste0(extdata,"/test_Simulate_CRT.csv"))
      testEstimates <- Analyse_CRT(trial = trial,
                                   method = 'GEE',excludeBuffer = FALSE,
                                   requireBootstrap = FALSE,alpha = 0.2)
      value <- testEstimates$contamination$data$positives[4]
      return(value)
   }
   expect_equal(get_test4(extdata), 151)
})

set.seed(1234)
get_test5 = function(extdata){
   Solarmal_baseline <- read.csv(file = paste0(extdata,"/Solarmal_baseline.csv"))
   testLocationsLatLong <- Solarmal_baseline[, c('lat','long')]
   testLocationsxy <- Convert_LatLong(testLocationsLatLong) #test_site is simulated
   testAnonymisedLocations <- Anonymise_TrialSite(testLocationsxy)
   testClusters <- DefineClusters(testAnonymisedLocations,h = 50)
   testArms <- Randomize_CRT(trial = testClusters,matchedPair = FALSE)
   testBuffer <- Specify_CRTbuffer(trial = testArms, bufferWidth = 0.1)
   testBuffer$cluster <- as.numeric(testBuffer$cluster)
   testBuffer$arm <- as.character(testBuffer$arm)
   return(testBuffer)}
test_that("Anonymisation, randomization, and creation of buffer produces expected trial", {
   expect_equal(get_test5(extdata), read.csv(file = paste0(extdata,"/testBuffer.csv")))
})

test_that("Analysis using T option gives expected efficacy", {
   get_test6 = function(extdata){
      trial <- read.csv(file = paste0(extdata,"/test_Simulate_CRT.csv"))
      analysis <- Analyse_CRT(trial=trial,method = 'T',link='identity')
      value <- round(as.numeric(10000 * analysis$pt.ests$efficacy))
      return(value)}
   expect_equal(get_test6(extdata), 1998)
})

test_that("Analysis using GEE option gives expected ICC", {
   get_test7 = function(extdata){
      trial <- read.csv(file = paste0(extdata,"/test_Simulate_CRT.csv"))
      analysis <- Analyse_CRT(trial=trial,method = 'GEE',link='log')
      value <- round(as.numeric(10000 * analysis$pt.ests$ICC))
      return(value)}
   expect_equal(get_test7(extdata), 464)
})

test_that("Analysis using ML option gives expected contamination range", {
   get_test8 = function(extdata){
      trial <- read.csv(file = paste0(extdata,"/test_Simulate_CRT.csv"))
      analysis <- Analyse_CRT(trial=trial,method = 'ML',link='logit',cfunc='P')
      value <- round(as.numeric(10000 * analysis$pt.ests$contaminationRange))
      return(value)}
   expect_equal(get_test8(extdata), 23563)
})

