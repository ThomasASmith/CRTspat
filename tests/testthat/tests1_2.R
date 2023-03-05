set.seed(1234)
test_that("Design_CRT() creates the default trial", {
  get_test1 = function(){
    Solarmal_baseline <- readdata("Solarmal_baseline.csv")
    testLocationsxy <- convert.latlong.xy(Solarmal_baseline) #test_site is simulated
    test_design <- Design_CRT(alpha = 0.05,
                              desiredPower = 0.8,
                              effect = 0.6,
                              ICC = 0.175,
                              pC = 0.4,
                              buffer.width = 0.05,
                              trial = testLocationsxy,
                              h = 80,
                              outcome.type ="Dichotomous")
    test_design$CRT.design.full <- NULL
    test_design$CRT.design.core <- NULL
    test_design$input.parameters <- NULL
    class(test_design) <- "data.frame"
    test_design$cluster <- as.numeric(test_design$cluster)
    test_design$arm <- as.character(test_design$arm)
    df <- test_design[test_design$buffer == FALSE,]
    row.names(df) <- NULL
    # To recreate test file
    # write.csv(df, file = "inst/extdata/test_design.csv", row.names = FALSE)
    return(df)
  }
  expect_equal(get_test1(),readdata("test_design.csv"))
})

set.seed(1234)
test_that("Simulate_TrialSite() creates the default site", {
  get_test2 = function(){
    trial <- Simulate_TrialSite()
    class(trial) <- "data.frame"
    return(trial)
  }
  expect_equal(get_test2(),readdata("test_site.csv"))
})
