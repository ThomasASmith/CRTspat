set.seed(1234)
test_that("designCRT() creates the default trial", {
  get_test1 = function(){
    Solarmal_baseline <- readdata("Solarmal_baseline.csv")
    testLocationsxy <- latlong_as_xy(Solarmal_baseline) #test_site is simulated
    test_buffered <- as_CRTspat(testLocationsxy) %>%
      specify_clusters(h = 50, algorithm = 'NN') %>%
      randomizeCRT(matchedPair = FALSE) %>%
      specify_buffer(buffer.width = 0.05)
    test_design <- CRTpower(trial = test_buffered,
                            alpha = 0.05,
                            desiredPower = 0.8,
                            effect = 0.6,
                            ICC = 0.175,
                            yC = 0.4,
                            k = 20,
                            outcome.type ="d")

    # To recreate test file
    # design <- test_design$design
    # dump(c("design"), file = "inst/extdata/analysis_design.txt")
    return(test_design$design)
  }
  expect_equal(get_test1(),readdata("analysis_design.txt"))
})

set.seed(1234)
test_that("simulate_site() creates the default site", {
  get_test2 = function(){
    CRT <- simulate_site(geoscale = 0.5,
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
