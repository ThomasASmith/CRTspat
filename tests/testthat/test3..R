
set.seed(1234)
test_that("Simulate_CRT() creates the default simulation", {
 get_test3 = function(){
   Solarmal_baseline <- readdata("Solarmal_baseline.csv")
   test_Locationsxy <- latlong_as_xy(Solarmal_baseline) #test_site is simulated
   test_Clusters <- specify.clusters(test_Locationsxy,h = 50)
   test_Arms <- randomizeCRT(trial = test_Clusters,matchedPair = TRUE)
   test.simulateCRT <- simulateCRT(trial = test_Arms,
                        theta_inp = 1.2,outcome0 = 0.4,
                        ICC_inp = 0.05,efficacy = 0.4,tol = 0.05)
   test.simulateCRT$cluster <- as.numeric(test.simulateCRT$cluster)
   test.simulateCRT$arm <- as.character(test.simulateCRT$arm)
   test.simulateCRT$RDT_test_result <- NULL
   test.simulateCRT$X <- NULL
   rownames(test.simulateCRT) <- NULL
   test.simulateCRT <- CRT_as_data.frame(test.simulateCRT)
   # write.csv(test.simulateCRT,file = "inst/extdata/test.simulateCRT1.csv", row.names = FALSE)
 return(test.simulateCRT)}
 expect_equal(get_test3(),readdata("test.simulateCRT1.csv"))
})

