set.seed(1234)
test_that("Toy GEE analysis creates correct output", {
   get_test4 = function(){
      # To recreate input file from vignette 2
      # write.csv(example2a$trial, file = "inst/extdata/exampleCRT.txt", row.names = FALSE)
      trial <- readdata("exampleCRT.txt")
      test_Estimates <- CRTanalysis(trial = trial,
                                   method = 'GEE',excludeBuffer = FALSE,
                                   alpha = 0.2)
      value <- round(test_Estimates$pt_ests$effect_size*1000)
      return(value)
   }
   expect_equal(get_test4(), 388)
})

set.seed(1234)
get_test5 = function(){
   test_locationsLatLong <- readdata("example_latlong.csv")
   test_locationsxy <- latlong_as_xy(test_locationsLatLong) #test_site is simulated
   test_anonymizedlocations <- anonymize_site(test_locationsxy)
   test_clusters <- specify_clusters(test_anonymizedlocations,h = 5)
   test_arms <- randomizeCRT(trial = test_clusters,matchedPair = FALSE)
   test_buffer <- specify_buffer(trial = test_arms, buffer_width = 0.1)
   trial <- test_buffer$trial
   trial$cluster <- as.numeric(trial$cluster)
   trial$arm <- as.character(trial$arm)
   # To recreate test file
   # write.csv(trial, file = "inst/extdata/example_buffer.csv", row.names = FALSE)
   return(trial)}
test_that("Anonymisation, randomization, and creation of buffer produces expected trial", {
   expect_equal(get_test5(), readdata("example_buffer.csv"))
})

test_that("Analysis using T option gives expected efficacy", {
   get_test6 = function(extdata){
      analysis <- CRTanalysis(readdata("exampleCRT.txt"),method = 'T',link='identity')
      value <- round(as.numeric(10000 * analysis$pt_ests$effect_size))
      return(value)}
   expect_equal(get_test6(extdata), 1446)
})

test_that("Analysis using GEE option gives expected coefficient of variation", {
   get_test7 = function(extdata){
      analysis <- CRTanalysis(readdata("exampleCRT.txt"),method = 'GEE',link='log')
      value <- round(as.numeric(10 * analysis$description$cv_percent))
      return(value)}
   expect_equal(get_test7(extdata), 484)
})

