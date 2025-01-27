
set.seed(1234)
test_that("Simulate_CRT() creates the default simulation", {
 get_test3 = function(){
    test_Clusters <- specify_clusters(readdata("example_site.csv"),h = 50)
    test_Arms <- randomizeCRT(trial = test_Clusters,matchedPair = TRUE)
    test_simulateCRT <- simulateCRT(trial = test_Arms,
                                    spillover_interval = 1.2, outcome0 = 0.4,
                                    ICC_inp = 0.05 ,effect = 0.4, tol = 0.05)
    value <- round(test_simulateCRT$geom_full$sd_h * 100)
 return(value)}
 expect_equal(get_test3(),350)
})

