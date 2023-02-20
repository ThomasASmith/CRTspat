extdata <- system.file("extdata",package = 'CRTspillover')
set.seed(1234)
test_that("Toy GEE analysis creates correct output", {
   get_test4 = function(extdata){
      trial <- read.csv(file = paste0(extdata,"/test_Simulate_CRT.csv"))
      test_Estimates <- Analyse_CRT(trial = trial,
                                   method = 'GEE',excludeBuffer = FALSE,
                                   alpha = 0.2)
      value <- test_Estimates$contamination$data$positives[4]
      return(value)
   }
   expect_equal(get_test4(extdata), 151)
})

set.seed(1234)
get_test5 = function(extdata){
   Solarmal_baseline <- read.csv(file = paste0(extdata,"/Solarmal_baseline.csv"))
   test_LocationsLatLong <- Solarmal_baseline[, c('lat','long')]
   test_Locationsxy <- Convert_LatLong(test_LocationsLatLong) #test_site is simulated
   test_AnonymisedLocations <- Anonymise_TrialSite(test_Locationsxy)
   test_clusters <- DefineClusters(test_AnonymisedLocations,h = 50)
   test_arms <- Randomize_CRT(trial = test_clusters,matchedPair = FALSE)
   test_buffer <- Specify_CRTbuffer(trial = test_arms, bufferWidth = 0.1)
   test_buffer$cluster <- as.numeric(test_buffer$cluster)
   test_buffer$arm <- as.character(test_buffer$arm)
   return(test_buffer)}
test_that("Anonymisation, randomization, and creation of buffer produces expected trial", {
   expect_equal(get_test5(extdata), read.csv(file = paste0(extdata,"/test_Buffer.csv")))
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

