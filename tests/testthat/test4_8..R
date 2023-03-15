set.seed(1234)
test_that("Toy GEE analysis creates correct output", {
   get_test4 = function(){
      trial <- readdata("test_Simulate_CRT.csv")
      test_Estimates <- analyseCRT(trial = trial,
                                   method = 'GEE',excludeBuffer = FALSE,
                                   alpha = 0.2)
      value <- test_Estimates$contamination$data$positives[4]
      return(value)
   }
   expect_equal(get_test4(), 151)
})

set.seed(1234)
get_test5 = function(){
   Solarmal_baseline <- readdata("Solarmal_baseline.csv")
   test.locationsLatLong <- Solarmal_baseline[, c('lat','long')]
   test.locationsxy <- latlong_as_xy(test.locationsLatLong) #test_site is simulated
   test.anonymizedlocations <- anonymize.site(test.locationsxy)
   test.clusters <- specify.clusters(test.anonymizedlocations,h = 50)
   test.arms <- randomizeCRT(trial = test.clusters,matchedPair = FALSE)
   test.buffer <- specify.buffer(trial = test.arms, buffer.width = 0.1)
   test.buffer$cluster <- as.numeric(test.buffer$cluster)
   test.buffer$arm <- as.character(test.buffer$arm)
   test.buffer$CRT.spat.full <- NULL
   test.buffer$CRT.spat.core <- NULL
   test.buffer$design <- NULL
   class(test.buffer) <- "data.frame"
   # To recreate test file
   # write.csv(test.buffer, file = "inst/extdata/test.buffer.csv", row.names = FALSE)
   return(test.buffer)}
test_that("Anonymisation, randomization, and creation of buffer produces expected trial", {
   expect_equal(get_test5(), readdata("test.buffer.csv"))
})

extdata <- system.file("extdata",package = 'CRTspillover')
test_that("Analysis using T option gives expected efficacy", {
   get_test6 = function(extdata){
      trial <- read.csv(file = paste0(extdata,"/test_Simulate_CRT.csv"))
      analysis <- analyseCRT(trial=trial,method = 'T',link='identity')
      value <- round(as.numeric(10000 * analysis$pt.ests$effect.size))
      return(value)}
   expect_equal(get_test6(extdata), 709)
})

test_that("Analysis using GEE option gives expected ICC", {
   get_test7 = function(extdata){
      trial <- read.csv(file = paste0(extdata,"/test_Simulate_CRT.csv"))
      analysis <- analyseCRT(trial=trial,method = 'GEE',link='log')
      value <- round(as.numeric(10000 * analysis$pt.ests$ICC))
      return(value)}
   expect_equal(get_test7(extdata), 464)
})

