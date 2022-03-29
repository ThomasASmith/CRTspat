#' Analysis of cluster randomized trial with contamination
#'
#' \code{Analyse_CRT} returns outputs from a statistical analysis of a cluster randomized trial (CRT).
#' @param trial trial dataframe including locations, clusters, arms, and binary outcomes
#' @param method statistical method used to analyse trial. Options are 'piecewise_linear','logit','sigmoid','empirical','GEE','MCMC01','MCMC02','MCMC03'
#' @param excludeBuffer exclude any buffer zone (records with buffer=TRUE) from the analysis
#' @param requireBootstrap logical indicator of whether bootstrap confidence intervals are required
#' @param alpha confidence level for confidence intervals and credible intervals
#' @param resamples number of bootstrap samples
#' @param burnin number of burnin interations of MCMC algorithm (for MCMC methods)
#' @return list containing the following results of the analysis
#' \itemize{
#' \item \code{description}: Description of the trial dataset
#' \item \code{method}: statistical method
#' \item \code{PointEstimates}: point estimates
#' \item \code{IntervalEstimates}: interval estimates
#' \item \code{ModelObject}: object returned by the fitting function
#' }
#' @importFrom grDevices rainbow
#' @importFrom stats binomial dist kmeans median na.omit qlogis qnorm quantile rbinom rnorm runif simulate
#' @importFrom utils head read.csv
#' @export
#'
#' @examples
#' # Standard GEE analysis of test dataset
#' exampleGEE=Analyse_CRT(trial=test_Simulate_CRT,method='GEE')

Analyse_CRT <- function(trial,
                        method='sigmoid',
                        excludeBuffer=FALSE,
                        requireBootstrap=FALSE,
                        alpha =0.05,
                        resamples=1000,
                        burnin=500){

  ##############################################################################
  # MAIN FUNCTION CODE STARTS HERE

  if("buffer" %in% colnames(trial) & excludeBuffer)
  {
    trial = trial[!trial$buffer,]
  }
  #trial needs to be ordered for MCMC analyses and for bootstrap of an08
  trial <- trial[order(trial$cluster),]

  # if nearestDiscord is not provided augment the trial data frame with distance to nearest discordant coordinate
  # (specifyBuffer assigns a buffer only if a buffer width is > 0 is input)
  if(is.null(trial$nearestDiscord)) {trial <- Specify_CRTbuffer(trial=trial,bufferWidth=0)}

  PointEstimates=IntervalEstimates=list()
  ModelObject=NA
  eta = alpha
  sd = 0.5/(qnorm(1-eta)*sqrt(2)) #initial value used in bootstrap calculations
  trial$neg=trial$denom - trial$num  #count of negatives for use in geeglm formulae

  if(method=='piecewise_linear'){
    #an01: piecewise_linear
    initialTheta <- qnorm(1-eta,sd = sqrt(2*sd^2))/(1-2*eta)*c(1/4,4)
    PointEstimates <- EstimateEffectiveness(trial, FUN=CalculateExpectP01,initialTheta)
    PointEstimates$contaminationRange <- TransformPW(PointEstimates$GAsolution3,eta)
    if(requireBootstrap){
      es = unlist(PointEstimates[1:4])
      ml01 <- c(logit(es[2]),logit(es[1]) - logit(es[2]),es[4])
      boot_an01 <- boot::boot(data=trial, statistic=BootSingleTrialAnalysis01,
                        num=resamples, sim="parametric", ran.gen=rgen_an01, mle=ml01)
      PointEstimates$bootstrapMean_efficacy = mean(boot_an01$t[,1])
      PointEstimates$bootstrapMean_contaminationRange = mean(boot_an01$t[,2])
      IntervalEstimates =
        addBootstrapIntervals(IntervalEstimates=IntervalEstimates,boot=boot_an01)
    }
  }
  if(method=='logit'){
    #an02: logit_model
    max_d <- max(abs(trial$nearestDiscord)) + 1e-5
    eta_d <- (trial$nearestDiscord + max_d)/(2 * max_d)
    min_beta_d <- min(log(eta_d/(1-eta_d)))
    eta_d_bound <- (qnorm(1-eta,sd = sqrt(2*sd^2))*c(1/4) + max_d)/(2 * max_d)
    beta_d <- log(eta_d_bound/(1-eta_d_bound))
    d <- ((min_beta_d-beta_d)/min_beta_d) #other extreme is that cont is at max_d => d = 0 or 2
    initialTheta <- c(0,log((1-eta)/eta)/log(d) - 2)
    PointEstimates <- EstimateEffectiveness(trial, FUN=CalculateExpectP02,initialTheta)
    if(requireBootstrap){
      es = unlist(PointEstimates[1:4])
      ml02 <- c(logit(es[2]),logit(es[1]) - logit(es[2]),es[4])
      boot_an02 <- boot::boot(data=trial, statistic=BootSingleTrialAnalysis02,
                        num=resamples, sim="parametric", ran.gen=rgen_an02, mle=ml02)
      PointEstimates$bootstrapMean_efficacy = mean(boot_an02$t[,1])
      PointEstimates$bootstrapMean_contaminationRange = mean(boot_an02$t[,2])
      IntervalEstimates =
        addBootstrapIntervals(IntervalEstimates=IntervalEstimates,boot=boot_an02)
    }
    PointEstimates <- PointEstimates[names(PointEstimates) != "GAsolution3"]
  }
  if(method=='sigmoid'){
    #an03: sigmoid_function
    initialTheta <- log((1-eta)/eta)/qnorm(1-eta,sd = sqrt(2*sd^2))*c(1/4,4)
    #note: this is equivalent to qlogis(0.95,scale = 1/qnorm(1-eta,sd = sqrt(2*sd^2)))
    PointEstimates <- EstimateEffectiveness(trial, FUN=CalculateExpectP03,initialTheta)
    if(requireBootstrap){
      es = unlist(PointEstimates[1:4])
      ml03 <- c(logit(es[2]),logit(es[1]) - logit(es[2]),es[4])
      boot_an03 <- boot::boot(data=trial, statistic=BootSingleTrialAnalysis03,
                        num=resamples, sim="parametric", ran.gen=rgen_an03, mle=ml03)
      PointEstimates$bootstrapMean_efficacy = mean(boot_an03$t[,1])
      PointEstimates$bootstrapMean_contaminationRange = mean(boot_an03$t[,2])
      IntervalEstimates =
        addBootstrapIntervals(IntervalEstimates=IntervalEstimates,boot=boot_an03)
    }
    PointEstimates <- PointEstimates[names(PointEstimates) != "GAsolution3"]
  }
  if(method=='empirical'){
    #an05: empirical analysis
    PointEstimates <- EmpiricalAnalysis(trial)
    if(requireBootstrap){
      es = unlist(PointEstimates)
      boot_an05 <- boot::boot(data=trial, statistic=BootEmpiricalAnalysis,
                      num=resamples, sim="parametric", ran.gen=rgen_an05, mle=es[3])
      PointEstimates$bootstrapMean_efficacy = mean(boot_an05$t)
      IntervalEstimates$efficacy <- namedCL(quantile(boot_an05$t,c(eta/2,1-eta/2)),eta=eta)
    }
  }
  if(method=='GEE'){
    #an08: GEE analysis
    GEEestimates <- GEEAnalysis(trial,eta=eta,resamples=resamples)
    PointEstimates <- GEEestimates$PointEstimates
    IntervalEstimates <- GEEestimates$IntervalEstimates
    if(requireBootstrap){
      ml <- geepack::geeglm(cbind(num,neg) ~ arm, id = cluster, corstr = "exchangeable", data=trial, family=binomial(link="logit"))

      boot_an08 <- boot::boot(data=trial, statistic=BootGEEAnalysis,
                        num=resamples, sim="parametric", ran.gen=rgen_an08, mle=ml)

      PointEstimates$bootstrapMean_efficacy = mean(boot_an08$t)
      IntervalEstimates$bootstrapEfficacy <- namedCL(quantile(boot_an08$t,c(eta/2,1-eta/2)),eta=eta)
    }
  }
  if (method %in% c('MCMC01','MCMC02','MCMC03')){
    if (method=='MCMC01'){
      # mildly informative prior
      initialTheta <- qnorm(1-eta,sd = sqrt(2*sd^2))/(1-2*eta)*c(1/4,4)
      datajags_an01<-with(trial,list(interven=as.numeric(trial$arm)-1,cluster=cluster,num=num,denom=denom,x=nearestDiscord,N=nrow(trial),ncluster=max(cluster),initialTheta=initialTheta))

      jagsout = jagsUI::autojags(data=datajags_an01, inits=NULL,
                         parameters.to.save=c("Es","beta3"), model.file=textConnection(get_model_an01()),
                         n.chains= 4, n.adapt=NULL, iter.increment=1000, n.burnin=burnin, n.thin=1)

      PointEstimates$ModelObject<-jagsout$sims.list
      es <- c(jagsout$q2.5$Es,jagsout$q50$Es,jagsout$q97.5$Es)
      cont <- TransformPW(c(jagsout$q2.5$beta3,jagsout$q50$beta3,jagsout$q97.5$beta3),eta)
    }
    if (method=='MCMC02'){
      nD <- trial$nearestDiscord/max(abs(trial$nearestDiscord))
      max_d <- max(abs(nD)) + 1e-5
      prop_d <- (nD + max_d)/(2 * max_d)
      logit_d <- log(prop_d/(1-prop_d))
      trial$d <- ((min(logit_d)- logit_d)/min(logit_d))

      eta_d_bound <- (qnorm(1-eta,sd = sqrt(2*sd^2))*c(1/4) + max_d)/(2 * max_d)
      beta_d <- log(eta_d_bound/(1-eta_d_bound))
      d <- ((min(logit_d)-beta_d)/min(logit_d)) #other extreme is that cont is at max_d => d = 0 or 2

      initialTheta <- c(0,qlogis(0.95)/log(d) - 2)

      datajags_an02<-with(trial,list(interven=as.numeric(trial$arm)-1,cluster=cluster,num=num,denom=denom,N=nrow(trial),ncluster=max(cluster),d=d,initialTheta=initialTheta))

      jagsout = jagsUI::autojags(data=datajags_an02, inits=NULL,
                                 parameters.to.save=c("Es","beta3"), model.file=textConnection(get_model_an02()),
                                 n.chains= 4, n.adapt=NULL, iter.increment=1000, n.burnin=burnin, n.thin=1)

      PointEstimates$ModelObject<-jagsout$sims.list
      es <- c(jagsout$q2.5$Es,jagsout$q50$Es,jagsout$q97.5$Es)
      cont <- TransformPW(c(jagsout$q2.5$beta3,jagsout$q50$beta3,jagsout$q97.5$beta3),eta)
    }
    if (method=='MCMC03'){
      #take a mildly informative prior (same as for models before)
      initialTheta <- qlogis(0.95)/qnorm(1-eta,sd = sqrt(2*sd^2))*c(1/4,4)
      datajags_an03<-with(trial,list(interven=as.numeric(trial$arm)-1,cluster=cluster,num=num,denom=denom,x=nearestDiscord,N=nrow(trial),ncluster=max(cluster),initialTheta=initialTheta))

      jagsout = jagsUI::autojags(data=datajags_an03, inits=NULL,
                                 parameters.to.save=c("Es","beta3"), model.file=textConnection(get_model_an03()),
                                 n.chains= 4, n.adapt=NULL, iter.increment=1000, n.burnin=burnin, n.thin=1)

      PointEstimates$ModelObject<-jagsout$sims.list
      es <- c(jagsout$q2.5$Es,jagsout$q50$Es,jagsout$q97.5$Es)
      cont <- TransformPW(c(jagsout$q2.5$beta3,jagsout$q50$beta3,jagsout$q97.5$beta3),eta)
    }
    PointEstimates$efficacy=es[2]
    IntervalEstimates$efficacy=namedCL(c(es[1],es[3]),eta=eta)
    PointEstimates$contaminationRange=cont[2]
    IntervalEstimates$contaminationRange=namedCL(c(cont[1],cont[3]),eta=eta)
  }
  description= get_description(trial)

# tidy up and consolidate the list of results
  ModelObject=PointEstimates$ModelObject
  PointEstimates <- PointEstimates[names(PointEstimates) != "GAsolution3"]
  PointEstimates <- PointEstimates[names(PointEstimates) != "ModelObject"]
  results = list(description=description,
               method=method,
               PointEstimates=PointEstimates,
               IntervalEstimates=IntervalEstimates,
               ModelObject=ModelObject)

return(results)}


#' Simple description and estimation of intra-cluster correction (ICC) from baseline data
#'
#' @param trial trial dataframe containing cluster assignments (variable cluster), numerators (num), and denominators (denom)
#' @param baselineNumerator name of numerator variable for baseline data (if present)
#' @param baselineDenominator name of denominator variable for baseline data (if present)
#' @param method method for estimating ICC (uses package 'ICCbin')
#' @param ci.type method for estimating confidence intervals for ICC (uses package 'ICCbin')
#' @return list containing calculation of average proportion and output from package 'ICCbin'
#' @export
Analyse_baseline = function(trial,
                            baselineNumerator='base_num',
                            baselineDenominator='base_denom',
                            method='aovs',
                            ci.type='aov'){
  # If numerators are not provided, assign a value of 1 to all records
  if(is.null(trial[[baselineDenominator]])) { trial[[baselineDenominator]] = 1 }
  expand= do.call(rbind.data.frame, lapply(1:nrow(trial), function(i) {
    neg=trial[[baselineDenominator]][i]-trial[[baselineNumerator]][i]
    df = rbind(data.frame(cluster = trial$cluster[i], y=1, z=0:trial[[baselineNumerator]][i]),
               data.frame(cluster = trial$cluster[i], y=0, z=0:neg))
    return(df[df$z > 0, c('cluster','y')])}))

  ests = ICCbin::iccbin(cid=as.factor(cluster), y=y, data=expand,
                                method=method, ci.type=ci.type, alpha = 0.05)
  ests$description = noquote(prettyNum(c(nrow(trial),sum(expand$y),length(expand$y),
                           round(mean(expand$y),digits=3))))
  names(ests$description) = c('locations','events','denominator','proportion')
return(ests)}


#add labels to confidence limits
namedCL=function(limits,eta=eta){
  names(limits)=c(paste0(100*eta/2,'%'),paste0(100-100*eta/2,'%'))
  return(limits)}

#inverse logit transformation
ilogit = function(logitp){return(1/(1+exp(-logitp)))}

# Minimal data description and crude efficacy estimate
get_description = function(trial){
  positives = tapply(trial$num, trial$arm, FUN=sum)
  totals = tapply(trial$denom, trial$arm, FUN=sum)
  ratios = positives/totals
  efficacy = 1 - ratios[2]/ratios[1]
  description = list(positives=positives,totals=totals,ratios=ratios,efficacy=efficacy)
  return(description)}

# Log Likelihood to be maximized
pseudoLogLikelihood <- function(par, FUN=CalculateExpectP,trial) {
  logitexpectP <- FUN(par,trial)
  transf <- 1/(1+exp(-logitexpectP)) #inverse logit transformation

  # FOR BERNOULLI CASE
  #pseudoLogLikelihood <- sum(log(ifelse(trial$num > 0 , transf,(1 - transf))))

  #alternative:
  #pseudoLogLikelihood <- sum(trial$num*logitexpectP - log(1+exp(logitexpectP)))

  # FOR BINOMIAL
  pseudoLogLikelihood <- sum(trial$num*log(transf) + trial$neg*log(1-transf))


  return(pseudoLogLikelihood)
}

##############################################################################

EstimateEffectiveness <- function(trial, FUN,initialTheta){

  #initial range for Effectiveness: in [0,1]
  GA <- GA::ga("real-valued", fitness = pseudoLogLikelihood, FUN, trial,
               lower = c(-10,-10,initialTheta[1]), upper = c(10,10,initialTheta[2]),
               maxiter = 500, run = 50, optim = TRUE,monitor = FALSE)
  #summary(GA)
  #plot(GA)

  # transform the parameters into interpretable functions
  pIhat <- 1/(1+exp(-GA@solution[1]))
  pChat <- 1/(1+exp(-(GA@solution[2] + GA@solution[1])))
  Eshat <- (pChat - pIhat)/pChat
  PointEstimates=list(controlP=pChat,
                      interventionP=pIhat,
                      efficacy=Eshat,
                      GAsolution3=GA@solution[3],
                      ModelObject=GA)
  return(PointEstimates)
}

##############################################################################
#  Different functions to estimate effectiveness
##############################################################################

# piecewise linear model
CalculateExpectP01 <- function(par,trial){

  pIhat <- par[1]
  tau <- par[2]
  theta <- par[3]

  expectP <- pIhat + (tau/(2*theta))*(theta - trial$nearestDiscord)
  expectP <- ifelse((trial$nearestDiscord < -theta),tau+pIhat,expectP)
  expectP <- ifelse((trial$nearestDiscord > theta),pIhat,expectP)
  return(expectP)
}

##############################################################################

# sigmoidal function based on logit distance
CalculateExpectP02 <- function(par,trial){

  max_d <- max(abs(trial$nearestDiscord)) + 1e-5
  prop_d <- (trial$nearestDiscord + max_d)/(2 * max_d)
  logit_d <- log(prop_d/(1-prop_d))
  d <- ((min(logit_d)- logit_d)/min(logit_d))

  pIhat <- par[1]
  tau <- par[2]
  theta <- par[3]

  expectP <- pIhat + tau/(1 + d^(theta+2))
  return(expectP)
}

##############################################################################

# sigmoid function model
CalculateExpectP03 <- function(par,trial){

  pIhat <- par[1]
  tau <- par[2]
  theta <- par[3]

  expectP <- pIhat + (tau)/(1 + exp(theta*(trial$nearestDiscord)))

  return(expectP)
}

##############################################################################

# standard non-model based analyses
EmpiricalAnalysis <- function(trial){
  description <- get_description(trial)
  PointEstimates=list(controlP=description$ratios[1],
                      interventionP=description$ratios[2],
                      efficacy=description$efficacy)

  return(PointEstimates)
}

##############################################################################

GEEAnalysis <- function(trial,eta=eta, resamples=resamples){

  fit <- geepack::geeglm(cbind(num,neg) ~ arm, id = cluster, corstr = "exchangeable", data=trial, family=binomial(link="logit"))
  summary_fit = summary(fit)
  logitpC = summary_fit$coefficients[1,1]
  logitpI = summary_fit$coefficients[1,1] + summary_fit$coefficients[2,1]
  se_logitpC = summary_fit$coefficients[1,2]
  se_logitpI = sqrt(fit$geese$vbeta[1,1] + fit$geese$vbeta[2,2] + 2*fit$geese$vbeta[1,2])
  z=-qnorm(eta/2) #standard deviation score for calculating confidence intervals
  CL_pC = namedCL(ilogit(c(logitpC-z*se_logitpC,logitpC+z*se_logitpC)),eta=eta)
  CL_pI = namedCL(ilogit(c(logitpI-z*se_logitpI,logitpI+z*se_logitpI)),eta=eta)
  CL_eff = estimateCLEfficacy(mu=summary_fit$coefficients[,1], Sigma=fit$geese$vbeta ,eta=eta, resamples=resamples)

  # Intracluster correlation
  ICC = noLabels(summary_fit$corr[1]) #with corstr = "exchangeable", alpha is the ICC
  se_ICC = noLabels(summary_fit$corr[2])
  CL_ICC = namedCL(noLabels(c(ICC-z*se_ICC,ICC+z*se_ICC)),eta=eta)

  clusterSize = nrow(trial)/nlevels(as.factor(trial$cluster))
  DesignEffect = 1 + (clusterSize - 1)*ICC #Design Effect
  CL_DesignEffect = 1 + (clusterSize - 1)*CL_ICC

  PointEstimates=list(controlP=ilogit(logitpC),
                      interventionP=ilogit(logitpI),
                      efficacy=(1 - ilogit(logitpI)/ilogit(logitpC)),
                      ICC=ICC,
                      DesignEffect=DesignEffect,
                      ModelObject=fit)
  IntervalEstimates=list(controlP=CL_pC,
                         interventionP=CL_pI,
                         efficacy=CL_eff,
                         ICC=CL_ICC,
                         DesignEffect=CL_DesignEffect)
  estimates=list(PointEstimates=PointEstimates,IntervalEstimates=IntervalEstimates)
  return(estimates)}

noLabels=function(x){
  xclean = as.matrix(x)
  dimnames(xclean) = NULL
  xclean = as.vector(xclean)
  return(xclean)}

estimateCLEfficacy = function(mu, Sigma,eta=eta, resamples=resamples){
  # Use resampling approach to avoid need for complicated Taylor approximation
  # use at least 10000 samples (this is very cheap)
  resamples1 = max(resamples,10000,na.rm = TRUE)
  samples = MASS::mvrnorm(n = resamples1,mu=mu,Sigma=Sigma)
  pC = ilogit(samples[,1])
  pI = ilogit(samples[,1] + samples[,2])
  eff = 1 - pI/pC
  CL = quantile(eff, probs = c(eta/2, 1-eta/2))
return(CL)}

#functions for an01, an02 and an03

rgen_an01<-function(data,mle){
  out<-data

  #simulate data for numerator num
  modelp <- CalculateExpectP01(mle,out)
  transf <- 1/(1+exp(-modelp))
  out$num <- rbinom(length(transf),out$denom,transf) #simulate from bernoulli distribution

  return(out)
}

rgen_an02<-function(data,mle){
  out<-data

  #simulate data for numerator num
  modelp <- CalculateExpectP02(mle,out)
  transf <- 1/(1+exp(-modelp))
  out$num <- rbinom(length(transf),out$denom,transf) #simulate from bernoulli distribution

  return(out)
}

rgen_an03<-function(data,mle){
  out<-data

  #simulate data for numerator num
  modelp <- CalculateExpectP03(mle,out)
  transf <- 1/(1+exp(-modelp))
  out$num <- rbinom(length(transf),out$denom,transf) #simulate from bernoulli distribution

  return(out)
}

BootSingleTrialAnalysis01 <- function(resampledData,eta=eta,sd=sd) {

  initialTheta <- qnorm(1-eta,sd = sqrt(2*sd^2))/(1-2*eta)*c(1/4,4)

  an01 <- unlist(EstimateEffectiveness(resampledData, FUN=CalculateExpectP01,initialTheta)[1:4])

  an01[4] <- TransformPW(an01[4],eta)

  return(c(an01[c(3,4)]))
}

BootSingleTrialAnalysis02 <- function(resampledData,eta=eta,sd=sd) {

  max_d <- max(abs(resampledData$nearestDiscord)) + 1e-5
  eta_d <- (resampledData$nearestDiscord + max_d)/(2 * max_d)
  min_beta_d <- min(log(eta_d/(1-eta_d)))

  eta_d_bound <- (qnorm(1-eta,sd = sqrt(2*sd^2))*c(1/4) + max_d)/(2 * max_d)
  beta_d <- log(eta_d_bound/(1-eta_d_bound))
  d <- ((min_beta_d-beta_d)/min_beta_d) #other extreme is that cont is at max_d => d = 0 or 2

  initialTheta <- c(0,log((1-eta)/eta)/log(d) - 2)

  an02 <- unlist(EstimateEffectiveness(resampledData, FUN=CalculateExpectP02,initialTheta)[1:4])

  an02[4] <- TransformLogitDist(resampledData$nearestDiscord,an02[4],eta)

  return(c(an02[c(3,4)]))
}

BootSingleTrialAnalysis03 <- function(resampledData,eta=eta,sd=sd) {

  initialTheta <- log((1-eta)/eta)/qnorm(1-eta,sd = sqrt(2*sd^2))*c(1/4,4)

  an03 <- unlist(EstimateEffectiveness(resampledData, FUN=CalculateExpectP03,initialTheta)[1:4])

  an03[4] <- TransformSigm(an03[4],eta)

  return(c(an03[c(3,4)]))
}



addBootstrapIntervals = function(IntervalEstimates,boot){
  IntervalEstimates$efficacy <- namedCL(quantile(boot$t[,1],c(eta/2,1-eta/2)),eta=eta)
  IntervalEstimates$contaminationRange <- namedCL(quantile(boot$t[,2],c(eta/2,1-eta/2)),eta=eta)
  return(IntervalEstimates)}

##############################################################################
#functions for an05

rgen_an05 <- function(data,mle){

  description <- describeBy(data$num/data$denom, group=data$arm)
  pChat <- description$control$mean
  pIhat <- description$intervention$mean

  #simulate data for numerator num
  modelp <- ifelse(as.numeric(data$arm) - 1 > 0, pIhat,pChat)
  data$num <- rbinom(length(modelp),data$denom,modelp) #simulate from Binomial distribution

  return(data)
}

# standard non-model based analyses
BootEmpiricalAnalysis <- function(resampledData){
  description <- describeBy(resampledData$num/resampledData$denom, group=resampledData$arm)
  #reports summary statistic by a grouping variable
  pChat <- description$control$mean
  pIhat <- description$intervention$mean
  Eshat <- 1 -  pIhat/pChat

  return(Eshat)
}

##############################################################################
#functions for an08

rgen_an08 <- function(data,mle){
  out<-data
  out$num <- unlist(simulate(mle))
  return(out)
}

BootGEEAnalysis <- function(resampledData){

  resampledData <- resampledData[order(resampledData$cluster),]
  fit <- geepack::geeglm(cbind(num,neg) ~ arm, id = cluster, corstr = "exchangeable", data=resampledData, family=binomial(link="logit"))

  pChat <- exp(fit$coefficients[1])/(1+exp(fit$coefficients[1]))
  pIhat <- exp(sum(fit$coefficients))/(1+exp(sum(fit$coefficients)))
  Eshat <- 1 - pIhat/pChat

  return(Eshat)
}

##############################################################################
#Backtransform for contamination range

TransformPW <- function(cont,eta){
  ContTrans <- cont*(1 - 2*eta)
  return(ContTrans)
}

TransformLogitDist <- function(nearestDiscord,cont,eta){

  #logit model
  max_d <- max(abs(nearestDiscord)) + 1e-5
  min_beta_d <- min(log((nearestDiscord + max_d)/(2 * max_d)/(1-(nearestDiscord + max_d)/(2 * max_d))))
  if (sign(nearestDiscord[which.max(abs(nearestDiscord))]) < 0){
    d <- exp(log((1-eta)/eta)/(cont+2)) #the sign changes
  } else{
    d <- exp(log(eta/(1-eta))/(cont+2))
  }
  beta_d <- (1-d)*min_beta_d
  eta_d <- exp(beta_d)/(1+exp(beta_d))

  ContTrans <- abs(2*max_d*eta_d - max_d)

  return(ContTrans)
}

TransformSigm <- function(cont,eta){

  ContTrans <- (log((1-eta)/eta))/cont
  return(ContTrans)
}


##############################################################################
#model formulations for MCMCan01, MCMCan02 and MCMCan03,
#interven makes it much more stable

# some reformulation st random effects can be assigned to part wout intervention
get_model_an01 = function(){ textstr=
"model{
  for(i in 1:N){
    num[i] ~ dbin(p[i],denom[i]) #every datapoint is drawn from a Binomial distribution with probability p
    p[i] <- 1/(1 + exp(-logitp[i])) #transformation of p
    logitp[i] <- ifelse(x[i] < -beta3,beta1[cluster[i]],ifelse(x[i] > beta3,beta1[cluster[i]] + beta2,beta1[cluster[i]] + beta2*interven[i]*(beta3 + x[i])/(2*beta3)))
  }

  # intercepts and random effects
  for(ic in 1:ncluster) {
    beta1[ic] ~ dnorm(intercept, tau)
  }

  # priors
  beta2 ~ dnorm(0, 0.001)
  intercept ~ dnorm(0, 0.001)
  tau <- 1/(sigma*sigma) #precision
  sigma ~ dunif(0,3) #uniform prior for the standard deviation of the random effects
  beta3 ~ dunif(initialTheta[1],initialTheta[2])
  pC <- 1/(1+exp(-intercept))
  pI <- 1/(1+exp(-intercept-beta2))
  Es <- 1 - pI/pC
  } "
return(textstr)}

##############################################################################

# some reformulation st random effects can be assigned to part without intervention
get_model_an02 = function(){ textstr=
"model{
  for(i in 1:N){
    num[i] ~ dbin(p[i],denom[i]) #every datapoint is drawn from a Binomial distribution with probability p
    p[i] <- 1/(1 + exp(-logitp[i])) #transformation of p
    logitp[i] <- beta1[cluster[i]] + beta2*interven[i]*(d[i]^(beta3 + 2)/(1 + d[i]^(beta3+2)))
  }

  # intercepts and random effects
  for(ic in 1:ncluster) {
    beta1[ic] ~ dnorm(intercept, tau)
  }

  # priors
  beta2 ~ dnorm(0, 0.001)
  intercept ~ dnorm(0, 0.001)
  tau <- 1/(sigma*sigma) #precision
  sigma ~ dunif(0,3) #from Tom
  beta3 ~ dunif(initialTheta[1], initialTheta[2])
  pC <- 1/(1+exp(-intercept))
  pI <- 1/(1+exp(-intercept-beta2))
  Es <- 1 - pI/pC
  } "
return(textstr)}

##############################################################################

# some reformulation st random effects can be assigned to part wout intervention
get_model_an03 = function(){ textstr=
"model{
  for(i in 1:N){
    num[i] ~ dbin(p[i],denom[i]) #every datapoint is drawn from a Binomial distribution with probability p
    p[i] <- 1/(1 + exp(-logitp[i])) #transformation of p
    logitp[i] <- beta1[cluster[i]] + beta2*interven[i]*(1/(1 + exp(-beta3*x[i])))
  }

  # intercepts and random effects
  for(ic in 1:ncluster) {
    beta1[ic] ~ dnorm(intercept, tau)
  }

  # priors
  beta2 ~ dnorm(0, 0.001)
  intercept ~ dnorm(0, 0.001)
  tau <- 1/(sigma*sigma) #precision
  sigma ~ dunif(0,3) #from Tom
  beta3 ~ dunif(initialTheta[1],initialTheta[2])
  pC <- 1/(1+exp(-intercept))
  pI <- 1/(1+exp(-intercept-beta2))
  Es <- 1 - pI/pC} "
return(textstr)}

