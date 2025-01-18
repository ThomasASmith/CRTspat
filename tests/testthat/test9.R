
# a stochastic optimisation is involved, but the test can fail because INLA uses parallel processing
# A separate seed should be set for each processor.
get_test9 = function(){
  if (identical(system.file(package='INLA'), "")){
    return(1303)
  } else {
    set.seed(1234)
    example_locations <- readdata('example_site.csv')
    example_locations$base_denom <- 1
    library(dplyr)
    example_randomized <- CRTsp(example_locations) %>%
      aggregateCRT(auxiliaries = c("RDT_test_result", "base_denom")) %>%
      specify_clusters(h = 50, algorithm = 'NN') %>%
      randomizeCRT(matchedPair = FALSE)
    example2a <- simulateCRT(example_randomized,
                             effect = 0.8,
                             outcome0 = 0.5,
                             generateBaseline = FALSE,
                             baselineNumerator = "RDT_test_result",
                             baselineDenominator = "base_denom",
                             ICC_inp = 0.05, spillover_interval = 0.8)
    # Reading in the inla.mesh functions when run outside the check but not as part of a check
    library(Matrix)
    inla_mesh <- readdata("examplemesh100.rds")
    analysis <- CRTanalysis(trial=example2a,
                            method = 'INLA', link='logit', cfunc='P',
                            personalProtection = TRUE, clusterEffects= TRUE,
                            requireMesh = TRUE, inla_mesh = inla_mesh)
    value <- round(analysis$pt_ests$DIC)
    return(value)
  }
}
test_that("Analysis using INLA option gives expected DIC", {
  expect_equal(get_test9(), 1303)
})
