
set.seed(1234)
test_that("Simulate_CRT() creates the default simulation", {
 get_test3 = function(){
   Solarmal_baseline <- readdata("Solarmal_baseline.csv")
   test_Locationsxy <- latlong_as_xy(Solarmal_baseline) #test_site is simulated
   test_Clusters <- specify_clusters(test_Locationsxy,h = 50)
   test_Arms <- randomizeCRT(trial = test_Clusters,matchedPair = TRUE)
   test.simulateCRT <- simulateCRT(trial = test_Arms,
                        theta_inp = 1.2,outcome0 = 0.4,
                        ICC_inp = 0.05 ,effect = 0.4, tol = 0.05)
   trial <- test.simulateCRT$trial
   trial$cluster <- as.numeric(trial$cluster)
   trial$arm <- as.character(trial$arm)
   trial$RDT_test_result <- NULL
   trial$X <- NULL
   rownames(trial) <- NULL
   # write.csv(trial,file = "inst/extdata/test.simulateCRT1.csv", row.names = FALSE)
 return(trial)}
 expect_equal(get_test3(),readdata("test_df.csv"))
})

