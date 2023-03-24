set.seed(1234)
# a stochastic optimisation is involved, but the test can fail because INLA uses parallel processing
# A separate seed should be set for each processor.
get_test9 = function(){
  CRT <- readdata("testCRT.csv")
  # Reading in the inla.mesh functions when run outside the check but not as part of a check
  library(Matrix)
  inla_mesh <- readdata("testmesh100.txt")
  analysis <- CRTanalysis(trial=CRT,
                          method = 'INLA', link='logit', cfunc='P',
                          localisedEffects = TRUE, clusterEffects= TRUE,
                          requireMesh = TRUE, inla.mesh = inla_mesh)
  value <- round(analysis$model.object$dic$dic)
  return(value)}
test_that("Analysis using INLA option gives expected DIc", {
  expect_equal(get_test9(), 3852)
})
