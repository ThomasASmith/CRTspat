set.seed(1234)
test_that("designCRT() creates the default trial", {
  get_test1 = function(){
    test_locationsLatLong <- readdata("example_latlong.csv")
    testLocationsxy <- latlong_as_xy(test_locationsLatLong)
    test_buffered <- CRTsp(testLocationsxy) %>%
      specify_clusters(h = 10, algorithm = 'NN') %>%
      randomizeCRT(matchedPair = FALSE) %>%
      specify_buffer(buffer_width = 0.05)
    test_design <- CRTpower(trial = test_buffered,
                            alpha = 0.05,
                            desiredPower = 0.8,
                            effect = 0.6,
                            ICC = 0.175,
                            yC = 0.4,
                            k = 20,
                            outcome_type ="d")
    design <- test_design$design
    # To recreate test file
    # design <- test_design$design
    # dump(c("design"), file = "inst/extdata/example_design.txt", evaluate = TRUE)
    return(design)
  }
  expect_equal(get_test1(),readdata("example_design.txt"))
})


test_that("CRTsp() creates the default site", {
  set.seed(1234)
  get_test2 = function(){
    CRT <- CRTsp(geoscale = 0.5,
                        locations = 5,
                        kappa = 4,
                        mu = 50)
    return(round(CRT$trial$x[3]*1000))
  }
  expect_equal(get_test2(),112)
})
