set.seed(1234)
test_that("Design_CRT() creates the default trial", {
  get_test1 = function(){
    Solarmal_baseline <- read.csv(file = "Solarmal_baseline.csv")
    testLocationsxy <- Convert_LatLong(Solarmal_baseline) #test_site is simulated
    coordinates <- Solarmal_baseline[]
    test_design <- Design_CRT(alpha = 0.05,
                    desiredPower = 0.8,
                    effect = 0.6,
                    ICC = 0.175,
                    pC = 0.4,
                    postulatedContaminationRange = 0.05,
                    coordinates = testLocationsxy,
                    h = 80)
    df <- test_design$descriptionCoreTrial$trial
    df$cluster <- as.numeric(df$cluster)
    df$arm <- as.character(df$arm)
    row.names(df) <- NULL
  return(df)
  }
expect_equal(get_test1(),read.csv(file = "test_design.csv"))
})

set.seed(1234)
test_that("Simulate_TrialSite() creates the default site", {
 expect_equal(Simulate_TrialSite(),read.csv(file = "test_site.csv"))
})

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

set.seed(1234)
test_that("Toy GEE analysis creates correct output", {
 get_test4 = function(){
             trial <- read.csv(file = "test_Simulate_CRT.csv")
             testEstimates <- Analyse_CRT(trial = trial,
             method = 'GEE',excludeBuffer = FALSE,
             requireBootstrap = FALSE,alpha = 0.2)
             value <- testEstimates$contamination$data$positives[4]
 return(value)
 }
 expect_equal(get_test4(), 116)
})

set.seed(1234)
get_test5 = function(){
   Solarmal_baseline <- read.csv(file = "Solarmal_baseline.csv")
   testLocationsLatLong <- Solarmal_baseline[, c('lat','long')]
   testLocationsxy <- Convert_LatLong(testLocationsLatLong) #test_site is simulated
   testAnonymisedLocations <- Anonymise_TrialSite(testLocationsxy)
   testClusters <- DefineClusters(testAnonymisedLocations,h = 50)
   testArms <- Randomize_CRT(trial = testClusters,matchedPair = FALSE)
   testBuffer <- Specify_CRTbuffer(trial = testArms, bufferWidth = 0.1)
   testBuffer$cluster <- as.numeric(testBuffer$cluster)
   testBuffer$arm <- as.character(testBuffer$arm)
return(testBuffer)}
test_that("Anonymisation, randomization, and creation of buffer produces expected trial", {
 expect_equal(get_test5(), read.csv(file = "testBuffer.csv"))
})

test_that("Analysis using T option gives expected efficacy", {
get_test6 = function(){
   analysis <- Analyse_CRT(trial=test_Simulate_CRT,method = 'T',link='identity')
   value <- round(as.numeric(10000 * analysis$pt.ests$efficacy))
return(value)}
   expect_equal(get_test6(), 306)
})

test_that("Analysis using GEE option gives expected ICC", {
   get_test7 = function(){
      analysis <- Analyse_CRT(trial=test_Simulate_CRT,method = 'GEE',link='log')
      value <- round(as.numeric(10000 * analysis$pt.ests$ICC))
      return(value)}
   expect_equal(get_test7(), 443)
})

test_that("Analysis using ML option gives expected contamination range", {
   get_test8 = function(){
      analysis <- Analyse_CRT(trial=test_Simulate_CRT,method = 'ML',link='logit',cfunc='P')
      value <- round(as.numeric(10000 * analysis$pt.ests$contaminationRange))
      return(value)}
   expect_equal(get_test8(), 45316)
})

set.seed(1234) # a stochastic optimisation is involved
get_test9 = function(){
   analysis <- Analyse_CRT(trial=test_Simulate_CRT,
               method = 'INLA', link='logit', cfunc='P',
               clusterEffects= TRUE, inla.mesh = inlaMesh100)
   value <- round(1000 * analysis$model.object$dic$dic)
   return(value)}
test_that("Analysis using INLA option gives expected DIc", {
   expect_equal(get_test9(), 5092287)
})
