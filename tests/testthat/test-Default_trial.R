set.seed(1234)
test_that("Design_CRT() creates the default trial", {
  expect_identical(Design_CRT(alpha = 0.05,
                              desiredPower = 0.8,
                              effect = 0.6,
                              ICC = 0.175,
                              pC = 0.4,
                              postulatedContaminationRange = 0.05,
                              coordinates=CRTspillover::test_coordinates,
                              h = 80)
                            , test_design)
})
set.seed(1234)
test_that("Simulate_TrialSite() creates the default site", {
  expect_identical(Simulate_TrialSite(), test_site)
})
set.seed(1234)
test_that("Simulate_CRT() creates the default simulation", {
  expect_identical(Simulate_CRT(trial=CRTspillover::test_design$descriptionFullTrial$trial,
                  theta_inp=1.2,initialPrevalence=0.4,
                  ICC_inp=0.05,efficacy=0.2,tol=0.05),
                  test_Simulate_CRT)
})
set.seed(1234)
test_that("Toy GEE analysis creates correct output", {
  get_test4 = function(){
    testEstimates= Analyse_CRT(trial=test_Simulate_CRT,method='GEE',excludeBuffer = FALSE,
                requireBootstrap=FALSE,alpha=0.2)
    testEstimates$ModelObject=NULL
    return(testEstimates)
  }
  expect_equal(get_test4(), test_Analyse_CRT)
})
set.seed(1234)
get_test5 = function(){
  testLocationsLatLong = Solarmal_baseline[, c('lat','long')]
  testLocationsxy=Convert_LatLong(testLocationsLatLong) #test_site is simulated
  testAnonymisedLocations=Anonymise_TrialSite(testLocationsxy)
  testClusters=DefineClusters(testAnonymisedLocations,h=50)
  testArms=Randomize_CRT(trial=testClusters,matchedPair = FALSE)
  testBuffer=Specify_CRTbuffer(trial = testArms, bufferWidth = 0.1)
  return(testBuffer)
}
test_that("Anonymisation, randomization, and creation of buffer produces expected trial",      {
  expect_identical(get_test5(), testBuffer)
})
