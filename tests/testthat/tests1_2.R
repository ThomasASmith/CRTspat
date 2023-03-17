set.seed(1234)
test_that("designCRT() creates the default trial", {
  get_test1 = function(){
    Solarmal_baseline <- readdata("Solarmal_baseline.csv")
    testLocationsxy <- latlong_as_xy(Solarmal_baseline) #test_site is simulated
    test_design <- designCRT(alpha = 0.05,
                              desiredPower = 0.8,
                              effect = 0.6,
                              ICC = 0.175,
                              yC = 0.4,
                              buffer.width = 0.05,
                              trial = testLocationsxy,
                              h = 80,
                              outcome.type ="d")
    test_design <- test_design$trial
    test_design$cluster <- as.numeric(test_design$cluster)
    test_design$arm <- as.character(test_design$arm)
    test_design <- test_design[test_design$buffer == FALSE,]
    row.names(test_design) <- NULL
    # To recreate test file
    # write.csv(df, file = "inst/extdata/test_design.csv", row.names = FALSE)
    return(test_design)
  }
  expect_equal(get_test1(),readdata("test_design.csv"))
})

set.seed(1234)
test_that("simulateSite() creates the default site", {
  get_test2 = function(){
    CRT <- simulateSite(geoscale = 0.5,
                        locations = 2500,
                        kappa = 4,
                        mu = 50)
    trial <- CRT$trial
    # To recreate test file
    # write.csv(trial, file = "inst/extdata/test_site.csv", row.names = FALSE)
    return(trial)
  }
  expect_equal(get_test2(),readdata("test_site.csv"))
})
