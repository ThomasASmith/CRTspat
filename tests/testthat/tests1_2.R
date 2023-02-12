set.seed(1234)
test_that("Design_CRT() creates the default trial", {
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
  expect_equal(get_test1(),read.csv(file = "test_design.csv"))
})

set.seed(1234)
test_that("Simulate_TrialSite() creates the default site", {
  expect_equal(Simulate_TrialSite(),read.csv(file = "test_site.csv"))
})
