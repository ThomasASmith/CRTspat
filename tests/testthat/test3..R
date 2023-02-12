
set.seed(1234)
test_that("Simulate_CRT() creates the default simulation", {
 get_test3 = function(){
   Solarmal_baseline <- read.csv(file = "Solarmal_baseline.csv")
   testLocationsxy <- Convert_LatLong(Solarmal_baseline) #test_site is simulated
   testClusters1 <- DefineClusters(testLocationsxy,h = 50)
   testArms1 <- Randomize_CRT(trial = testClusters,matchedPair = TRUE)
   test_Simulate_CRT <- Simulate_CRT(trial = testArms1,
                        theta_inp = 1.2,initialPrevalence = 0.4,
                        ICC_inp = 0.05,efficacy = 0.2,tol = 0.05)
   test_Simulate_CRT$cluster <- as.numeric(test_Simulate_CRT$cluster)
   test_Simulate_CRT$arm <- as.character(test_Simulate_CRT$arm)
   rownames(test_Simulate_CRT) <- NULL
 return(test_Simulate_CRT)}
 expect_equal(get_test3(),read.csv(file = "test_Simulate_CRT.csv"))
})


