
set.seed(1234)
test_that("Simulate_CRT() creates the default simulation", {
 get_test3 = function(extdata){
   Solarmal_baseline <- readdata("Solarmal_baseline.csv")
   test_Locationsxy <- convert.latlong.xy(Solarmal_baseline) #test_site is simulated
   test_Clusters <- specify.clusters(test_Locationsxy,h = 50)
   test_Arms <- Randomize_CRT(trial = test_Clusters,matchedPair = TRUE)
   test.simulateCRT <- Simulate_CRT(trial = test_Arms,
                        theta_inp = 1.2,initialPrevalence = 0.4,
                        ICC_inp = 0.05,efficacy = 0.4,tol = 0.05)
   test.simulateCRT$cluster <- as.numeric(test.simulateCRT$cluster)
   test.simulateCRT$arm <- as.character(test.simulateCRT$arm)
   test.simulateCRT$RDT_test_result <- NULL
   test.simulateCRT$X <- NULL
   rownames(test.simulateCRT) <- NULL
   class(test.simulateCRT) <- "data.frame"
   # write.csv(test.simulateCRT,file = "inst/extdata/test.simulateCRT1.csv", row.names = FALSE)
 return(test.simulateCRT)}
 expect_equal(get_test3(extdata),readdata("test.simulateCRT1.csv"))
})

