set.seed(1234)
# a stochastic optimisation is involved, but the test can fail because INLA uses parallel processing
# A separate seed should be set for each processor.
get_test9 = function(){
  CRT <- readdata("test.CRT.csv")
  analysis <- analyseCRT(trial=trial,
                          method = 'INLA', link='logit', cfunc='P',
                          localisedEffects = TRUE, clusterEffects= TRUE,
                          requireMesh = TRUE, inla.mesh = inlaMesh100)
  value <- round(10 * analysis$model.object$dic$dic)
  return(value)}
test_that("Analysis using INLA option gives expected DIc", {
  expect_equal(get_test9(), 50925)
})
