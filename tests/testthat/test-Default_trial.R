set.seed(1234)
test_that("Design_CRT() creates the default trial", {
  expect_identical(Design_CRT(), test_AvecNet)
})
set.seed(1234)
test_that("Simulate_TrialSite() creates the default site", {
  expect_identical(Simulate_TrialSite(), test_site)
})
set.seed(1234)
test_that("Simulate_CRT() creates the default simulation", {
  expect_identical(Simulate_CRT(), test_Simulate_CRT)
})
set.seed(1234)
test_that("Toy GEE analysis creates correct output", {
  expect_equal(Analyse_CRT(trial=test_Simulate_CRT,method='GEE',excludeBuffer = FALSE,
                               requireBootstrap=FALSE,
                               alpha=0.2), test_Analyse_CRT)
})
set.seed(1234)
get_test = function(){
  testLocationsLatLong = Solarmal_baseline[, c('lat','long')]
  testLocationsxy=Convert_LatLong(testLocationsLatLong) #test_site is simulated
  testAnonymisedLocations=Anonymise_TrialSite(testLocationsxy)
  testClusters=DefineClusters(testAnonymisedLocations,h=50)
  testArms=Randomize_CRT(testClusters)
  testBuffer=Specify_CRTbuffer(trial = testArms, bufferWidth = 0.1)
  return(testBuffer)
}
test_that("Anonymisation, randomization, and creation of buffer produces expected trial",
          {
  expect_identical(get_test(), testBuffer)
})
