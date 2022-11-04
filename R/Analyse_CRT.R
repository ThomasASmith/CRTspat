#' Analysis of cluster randomized trial with contamination
#'
#' \code{Analyse_CRT} returns outputs from a statistical analysis of a cluster randomized trial (CRT).
#' @param trial trial dataframe including locations, clusters, arms, and binary outcomes
#' @param method statistical method used to analyse trial. Options are 'L0','L1','L2','L3','GEE','M1','M2','M3'
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
                        method='L3',
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
  #trial needs to be ordered for MCMC analyses and for bootstrap of GEE
  trial <- trial[order(trial$cluster),]

  # if nearestDiscord is not provided augment the trial data frame with distance to nearest discordant coordinate
  # (specifyBuffer assigns a buffer only if a buffer width is > 0 is input)
  if(is.null(trial$nearestDiscord)) {trial <- Specify_CRTbuffer(trial=trial,bufferWidth=0)}

  PointEstimates=IntervalEstimates=list()
  ModelObject=NA
  eta = alpha
  sd = 0.5/(qnorm(1-eta)*sqrt(2)) #initial value used in bootstrap calculations
  trial$neg=trial$denom - trial$num  #count of negatives for use in geeglm formulae

############### Analyses that ignore contamination #################

  if(method=='L0'){
    # empirical analysis
    PointEstimates <- EmpiricalAnalysis(trial)
    if(requireBootstrap){
      boot_an05 <- boot::boot(data=trial, statistic=BootEmpiricalAnalysis,
                              R=resamples, sim="parametric", ran.gen=rgen_an05, mle=PointEstimates)
      PointEstimates$bootstrapMean_efficacy = mean(boot_an05$t)
      IntervalEstimates$efficacy <- namedCL(quantile(boot_an05$t,c(eta/2,1-eta/2)),eta=eta)
    }
  }

  if(method=='GEE'){
    #GEE analysis with cluster effects
    GEEestimates <- GEEAnalysis(trial,eta=eta,resamples=resamples)
    PointEstimates <- GEEestimates$PointEstimates
    IntervalEstimates <- GEEestimates$IntervalEstimates
    if(requireBootstrap){
      ml <- geepack::geeglm(cbind(num,neg) ~ arm, id = cluster, corstr = "exchangeable", data=trial, family=binomial(link="logit"))

      boot_output <- boot::boot(data=trial, statistic=BootGEEAnalysis,
                                R=resamples, sim="parametric", ran.gen=rgen_GEE, mle=ml)

      PointEstimates$bootstrapMean_efficacy = mean(boot_output$t)
      IntervalEstimates$bootstrapEfficacy <- namedCL(quantile(boot_output$t,c(eta/2,1-eta/2)),eta=eta)
    }
  }

############### ML Methods with contamination functions and logistic link #################
  if(method %in% c('L1','L2','L3')){
    i = which (method == c('L1','L2','L3'))
    LPMethods = c('CalculateLinearPredictor01',
                  'CalculateLinearPredictor02',
                  'CalculateLinearPredictor03')
    FUN2 <- FUN1 <- eval(parse(text=LPMethods[i]))
    par = SingleTrialAnalysis(trial=trial,FUN2=FUN2)
    PointEstimates <- FittingResults(trial, par=par,FUN1=FUN1)
    if(requireBootstrap){
      mle=list(par=par,FUN1=FUN1,link='logit')
      boot_estimates = data.frame(V1=c(),V2=c(),V3=c(),V4=c())
      # resampling can crash because of resamples containing data from only one arm
      # to prevent this crashing the whole program, bootstrapping is done in batches
      # of 5 resamples using 'try' to avoid crashing out
      resamples1=5
      while(nrow(boot_estimates) < resamples){
        # sample a random value each time through the loop so the seed is change
        boot_output <- try(expr= boot::boot(data=trial, statistic=SingleTrialAnalysis,
                                            R=resamples1, sim="parametric", ran.gen=rgen,
                                            mle=mle, FUN2=FUN2),silent=TRUE)
        if(!substr(boot_output[1],1,5)=="Error"){
          new_estimates=as.data.frame(t(matrix(data= unlist(apply(as.data.frame(boot_output$t), 1,
                                                                  FittingResults,trial=trial,
                                                                  FUN1=FUN1)),ncol=resamples1,nrow=length(PointEstimates))))
          boot_estimates = rbind(boot_estimates,new_estimates)
        }
      }
      colnames(boot_estimates) = names(PointEstimates)
      IntervalEstimates =  computeIntervals(df=boot_estimates,eta=eta)
    }
  }

  if (method %in% c('M1','M2','M3')){
    # mildly informative prior
    initialTheta <- qnorm(1-eta,sd = sqrt(2*sd^2))/(1-2*eta)*c(1/4,4)
    datajags<-with(trial,list(eta=eta,cluster=cluster,num=num,denom=denom,d=nearestDiscord,N=nrow(trial),ncluster=max(cluster),initialTheta=initialTheta))

    jagsout = jagsUI::autojags(data=datajags, inits=NULL,
                       parameters.to.save=c("Es","cont"), model.file=textConnection(get_model_an02()),
                       n.chains= n.chains, n.adapt=NULL, iter.increment=200, n.burnin=burnin, n.thin=1)

    PointEstimates$ModelObject<-jagsout$sims.list
    es <- c(jagsout$q2.5$Es,jagsout$q50$Es,jagsout$q97.5$Es)
    cont <- c(jagsout$q2.5$cont,jagsout$q50$cont,jagsout$q97.5$cont)
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

#logit transformation
logit = function(p){return(log(p/(1-p)))}

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
pseudoLogLikelihood <- function(par, FUN=FUN ,trial) {
  logitexpectP <- FUN(par,trial)
  transf <- 1/(1+exp(-logitexpectP)) #inverse logit transformation

  # FOR BINOMIAL
  pseudoLogLikelihood <- sum(trial$num*log(transf) + trial$neg*log(1-transf))


  return(pseudoLogLikelihood)
}

##############################################################################

FittingResults <- function(trial, FUN1, par){


  # transform the parameters into interpretable functions
  pIhat <- 1/(1+exp(-par[1]))
  pChat <- 1/(1+exp(-(par[2] + par[1])))
  Eshat <- (pChat - pIhat)/pChat

  #estimate contamination range
  #The absolute value of deltaP is used so that a positive range is obtained even with negative efficacy
  trial$deltaP <- abs(unlist(ifelse(trial$arm=='control', pChat - 1/(1+exp(-FUN1(trial=trial,par=par))),
                   1/(1+exp(-FUN1(trial=trial,par=par)))- pIhat)))
  #the maximum contamination should be where the distance is zero.
  max_deltaP <- max(trial$deltaP)
  trial <- trial[order(trial$nearestDiscord),]
  thetaL <- trial$nearestDiscord[which(trial$deltaP > 0.05*max_deltaP)[1]]
  trial <- trial[order(-trial$nearestDiscord),]
  thetaU <- trial$nearestDiscord[which(trial$deltaP > 0.05*max_deltaP)[1]]
  #contamination range
  theta <- thetaU - thetaL

  PointEstimates=list(controlP=pChat,
                      interventionP=pIhat,
                      efficacy=Eshat,
                      contaminationRange = theta)
  return(PointEstimates)
}

##############################################################################
#  Different functions to estimate effectiveness
##############################################################################

# piecewise linear model (on the logit scale) for contamination function

CalculateLinearPredictor01 <- function(par,trial){

  theta <- par[3]

  lp <- par[1] + (par[2]/(2*theta))*(theta - trial$nearestDiscord)
  lp <- ifelse((trial$nearestDiscord < -theta),par[2]+par[1],lp)
  lp <- ifelse((trial$nearestDiscord > theta),par[1],lp)
  return(lp)
}

##############################################################################

# sigmoidal function based on logit distance (on the logit scale) for contamination function
CalculateLinearPredictor02 <- function(par,trial){
  scaled_logit = get_scaled_logit(trial)
  lp <- par[1] + par[2]/(1 +  exp(-scaled_logit*par[3]))
  return(lp)
}

get_scaled_logit = function(trial){
  # add a small constant to the calculated maximum of d to avoid division by 0 at the limits
  max_d <- max(abs(trial$nearestDiscord)) + 1e-5
  # prop_d is d as a proporton of the maximum in the dataset
  prop_d <- (trial$nearestDiscord + max_d)/(2 * max_d)
  logit_d <- -log(prop_d/(1-prop_d))
  # scale the logit by its maximum value to avoid taking powers of large numbers
  scale_factor <- max(abs(logit_d))
  scaled_logit <- logit_d/scale_factor
return(logit_d)}

##############################################################################

# sigmoid function (on the logit scale) for contamination function
CalculateLinearPredictor03 <- function(par,trial){

  lp <- par[1] + par[2]/(1 + exp(par[3]*(trial$nearestDiscord)))
  return(lp)
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

#functions for analysis of L1, L2, L3

rgen<-function(data,mle){
  par=mle$par
  FUN1=mle$FUN1
  #simulate data for numerator num
  modelp <- FUN1(par=par,trial=data)
  if(mle$link == 'logit') {
    transf = ilogit(modelp)
  } else {
    transf = modelp
  }
  data$num <- rbinom(length(transf),data$denom,transf) #simulate from binomial distribution
  return(data)
}

SingleTrialAnalysis <- function(trial, FUN2=FUN2) {

  GA <- GA::ga("real-valued", fitness = pseudoLogLikelihood, FUN=FUN2,
               trial=trial,
               lower = c(-10,-10,-100), upper = c(10,10,100),
               maxiter = 500, run = 50, optim = TRUE,monitor = FALSE)
  result <- GA@solution

  return(result)
}

computeIntervals = function(df,eta){
  varnames= colnames(df)
  IntervalEstimates = as.list(varnames)
  for(i in 1:length(varnames)){
    IntervalEstimates[[i]]=namedCL(quantile(df[,varnames[i]],c(eta/2,1-eta/2)),eta=eta)
  }
  names(IntervalEstimates)=varnames
return(IntervalEstimates)}

##############################################################################
#functions for Empirical Analysis

# standard descriptive analyses
EmpiricalAnalysis <- function(trial){
  description <- get_description(trial)
  PointEstimates=list(controlP=unname(description$ratios[1]),
                      interventionP=unname(description$ratios[2]),
                      efficacy=unname(description$efficacy),
                      contaminationRange = NA)
  return(PointEstimates)
}


rgen_an05 <- function(data,mle){

  description <- psych::describeBy(data$num/data$denom, group=data$arm)
  pChat <- description$control$mean
  pIhat <- description$intervention$mean

  #simulate data for numerator num
  modelp <- ifelse(as.numeric(data$arm) - 1 > 0, pIhat,pChat)
  data$num <- rbinom(length(modelp),data$denom,modelp) #simulate from Binomial distribution

  return(data)
}

# standard non-model based analyses
BootEmpiricalAnalysis <- function(resampledData){
  description <- psych::describeBy(resampledData$num/resampledData$denom, group=resampledData$arm)
  #reports summary statistic by a grouping variable
  pChat <- description$control$mean
  pIhat <- description$intervention$mean
  Eshat <- 1 -  pIhat/pChat

  return(Eshat)
}

##############################################################################
#functions for bootstrap resampling of GEE

rgen_GEE <- function(data,mle){
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

# MCMC model with cluster random effect and piecewise linear model (on the logit scale) for contamination function
get_model_an01 = function(){ textstr=
"model{
  for(i in 1:N){
    num[i] ~ dbin(p[i],denom[i]) #every datapoint is drawn from a Binomial distribution with probability p
    fixedp[i] <- ifelse(d[i] < -beta,pC,ifelse(d[i] > beta,pI,pI + (pC - pI) * (beta - d[i])/(2*beta)))
    logitp[i] <- logit(fixedp[i]) + alpha[cluster[i]] #cluster random effect on logit scale
    p[i] <- 1/(1 + exp(-logitp[i])) #back transformation to linear scale
  }

# random effects
for(ic in 1:ncluster) {
  alpha[ic] ~ dnorm(0, tau)
}

# priors
pC ~ dunif(0,1)
pI ~ dunif(0,1)
tau <- 1/(sigma*sigma) #precision
sigma ~ dunif(0,3)
beta ~ dunif(initialTheta[1],initialTheta[2])

# derived quantities
Es <- 1 - pI/pC
cont <- (1 - eta)*2*beta # beta is half-width of the contaminated zone
}
"
return(textstr)}

##############################################################################

# MCMC model with cluster random effect and cumulative normal contamination function
get_model_an02 = function(){ textstr=
  "model{
    for(i in 1:N){
      num[i] ~ dbin(p[i],denom[i]) #every datapoint is drawn from a Binomial distribution with probability p
      fixedp[i] <- pI + (pC-pI)*pnorm(d[i],0,beta)  #model excluding random effect
      logitp[i] <- logit(fixedp[i]) + alpha[cluster[i]] #cluster random effect on logit scale
      p[i] <- 1/(1 + exp(-logitp[i])) #back transformation to linear scale
    }

    # random effects
    for(ic in 1:ncluster) {
      alpha[ic] ~ dnorm(0, tau)
    }

    # priors
    pC ~ dunif(0,1)
    pI ~ dunif(0,1)
    tau <- 1/(sigma*sigma) #precision
    sigma ~ dunif(0,3)
    beta ~ dunif(initialTheta[1],initialTheta[2])

    # derived quantities
    Es <- 1 - pI/pC
    cont <- 2 * qnorm(1 - 0.5*eta,0,beta) # beta is the precision of the contamination function
  }
"
return(textstr)}

##############################################################################
# MCMC model with cluster random effect and sigmoid tablemodel for contamination function
get_model_an03 = function(){ textstr=
  "model{
    for(i in 1:N){
      num[i] ~ dbin(p[i],denom[i]) #every datapoint is drawn from a Binomial distribution with probability p
      fixedp[i] <- pI + (pC-pI)/(1 + exp(-beta*d[i]))  #model excluding random effect
      logitp[i] <- logit(fixedp[i]) + alpha[cluster[i]] #cluster random effect on logit scale
      p[i] <- 1/(1 + exp(-logitp[i])) #back transformation to linear scale
    }

    # random effects
    for(ic in 1:ncluster) {
      alpha[ic] ~ dnorm(0, tau)
    }

    # priors
    pC ~ dunif(0,1)
    pI ~ dunif(0,1)
    tau <- 1/(sigma*sigma) #precision
    sigma ~ dunif(0,3)
    beta ~ dunif(initialTheta[1],initialTheta[2])

    # derived quantities
    Es <- 1 - pI/pC
    cont <- 2* log((1-0.5*eta)/(0.5*eta))/beta #contamination diameter
  }
"
return(textstr)}
