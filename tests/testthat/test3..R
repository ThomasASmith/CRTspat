
extdata <- system.file("extdata",package = 'CRTspillover')
set.seed(1234)
test_that("Simulate_CRT() creates the default simulation", {
 get_test3 = function(extdata){
   Solarmal_baseline <- read.csv(file = paste0(extdata,"/Solarmal_baseline.csv"))
   test_Locationsxy <- Convert_LatLong(Solarmal_baseline) #test_site is simulated
   test_Clusters <- DefineClusters(test_Locationsxy,h = 50)
   test_Arms <- Randomize_CRT(trial = test_Clusters,matchedPair = TRUE)
   test_Simulate_CRT <- Simulate_CRT(trial = test_Arms,
                        theta_inp = 1.2,initialPrevalence = 0.4,
                        ICC_inp = 0.05,efficacy = 0.4,tol = 0.05)
   test_Simulate_CRT$cluster <- as.numeric(test_Simulate_CRT$cluster)
   test_Simulate_CRT$arm <- as.character(test_Simulate_CRT$arm)
   test_Simulate_CRT$RDT_test_result <- NULL
   test_Simulate_CRT$X <- NULL
   rownames(test_Simulate_CRT) <- NULL
 return(test_Simulate_CRT)}
 expect_equal(get_test3(extdata),read.csv(file = paste0(extdata,"/test_Simulate_CRT1.csv")))
})

