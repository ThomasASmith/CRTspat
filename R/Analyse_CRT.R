#' Analysis of cluster randomized trial with contamination
#'
#' \code{Analyse_CRT} returns outputs from a statistical analysis of a cluster randomized trial (CRT).
#' @param trial trial dataframe including locations, clusters, arms, and outcomes
#' @param method statistical method used to analyse trial.
#' options are:
#' 'EMP'  : empirical,
#' 'ML'   : maximum likelihood,
#' 'GEE'  : generalised estimating equations
#' 'INLA' : INLA
#' @param cfunc transformation defining the contamination function
#' options are:
#' 'S': piecewise linear (slope),
#' 'L': inverse logistic (sigmoid),
#' 'P': inverse probit,
#' 'X': contamination not modelled'
#' @param numerator name of numerator variable for efficacy data (if present)
#' @param denominator name of denominator variable for efficacy data (if present)
#' @param excludeBuffer exclude any buffer zone (records with buffer=TRUE) from the analysis
#' @param alpha confidence level for confidence intervals and credible intervals
#' @param requireBootstrap logical indicator of whether bootstrap confidence intervals are required
#' @param baselineOnly logical: indicator of whether required analysis is of efficacy or of baseline only
#' @param baselineNumerator name of numerator variable for baseline data (if present)
#' @param baselineDenominator name of denominator variable for baseline data (if present)
#' @param localisedEffects logical: indicator of whether the model includes local effects with no contamination
#' @param clusterEffects logical: indicator of whether the model includes cluster random effects
#' @param spatialEffects logical: indicator of whether the model includes spatial random effects
#' @param resamples number of bootstrap samples
#' @param inlaMesh name of pre-existing INLA input object created by CreateMesh()
#' @return list containing the following results of the analysis
#' \itemize{
#' \item \code{description}: Description of the trial dataset
#' \item \code{method}: statistical method
#' \item \code{PointEstimates}: point estimates
#' \item \code{IntervalEstimates}: interval estimates
#' \item \code{ModelObject}: object returned by the fitting function
#' \item \code{contamination}: function values and statistics describing the estimated contamination
#' }
#' @importFrom grDevices rainbow
#' @importFrom stats binomial dist kmeans median na.omit qlogis qnorm quantile rbinom rnorm runif simulate
#' @importFrom utils head read.csv
#' @export
#'
#' @examples
#' # Standard GEE analysis of test dataset ignoring contamination
#' exampleGEE=Analyse_CRT(trial=test_Simulate_CRT,method='GEE')

Analyse_CRT <- function(trial,
                        method='ML',
                        cfunc='L',
                        numerator='num',
                        denominator='denom',
                        excludeBuffer=FALSE,
                        alpha = 0.05,
                        requireBootstrap=FALSE,
                        baselineOnly=FALSE,
                        baselineNumerator='base_num',
                        baselineDenominator='base_denom',
                        localisedEffects=FALSE,
                        clusterEffects=TRUE,
                        spatialEffects=FALSE,
                        resamples=1000,
                        inlaMesh=NULL){

  ##############################################################################
  # MAIN FUNCTION CODE STARTS HERE
  cat('\n=====================    ANALYSIS OF CLUSTER RANDOMISED TRIAL    =================\n')
  #store the input trial so that any changes to this can be removed before output
  input_trial = trial
  cluster=NULL

  if("buffer" %in% colnames(trial) & excludeBuffer)
  {
    trial = trial[!trial$buffer,]
  }

  #trial needs to be ordered for some analyses
  trial <- trial[order(trial$cluster),]

  if(method=='EMP' | method=='GEE') cfunc = 'X'
  if(baselineOnly){
    if(method %in% c('EMP','ML','GEE')){
      method='GEE'
      cat('Analysis of baseline only, using method GEE\n')
    } else if(method %in% c('LR','CRE','SPDE','SPCRE')) {
      method='CRE'
      cat('Analysis of baseline only, using method CRE\n')
    }
    # cfunc='Z' is used to remove the estimation of effect size from the model
    cfunc= 'Z'
    trial$y1=trial[[baselineNumerator]]
    trial$y0=trial[[baselineDenominator]]-trial[[baselineNumerator]]
    trial$y_off=trial[[baselineDenominator]]

  } else {
    trial$y1=trial[[numerator]]
    trial$y0=trial[[denominator]]-trial[[numerator]]
    trial$y_off=trial[[denominator]]

    # if nearestDiscord is not provided augment the trial data frame with distance to nearest discordant coordinate
    # (specifyBuffer assigns a buffer only if a buffer width is > 0 is input)
    if(is.null(trial$nearestDiscord)) {trial <- Specify_CRTbuffer(trial=trial,bufferWidth=0)}
  }
  PointEstimates=ModelObject=list()
  IntervalEstimates=list(controlP=NA,interventionP=NA,efficacy=NA)
  sd = 0.5/(qnorm(1-alpha)*sqrt(2)) #initial value used in bootstrap calculations
  description= get_description(trial)
  # Specify the function used for calculating the linear predictor
  LPfunction = c('CalculateNoEffect',
                 'CalculateNoContaminationFunction',
                 'CalculatePiecewiseLinearFunction',
                 'CalculateLogisticFunction',
                 'CalculateProbitFunction')[which (cfunc == c('Z','X','S','L','P'))]
  FUN2 <- FUN1 <- eval(parse(text=LPfunction))

  if(method=='EMP'){
    # empirical analysis that ignores contamination
    PointEstimates <- EmpiricalAnalysis(trial)
    PointEstimates$contaminationParameter = NA #contamination is not estimated
    if(requireBootstrap){
      boot_emp <- boot::boot(data=trial, statistic=BootEmpiricalAnalysis,
                             R=resamples, sim="parametric", ran.gen=rgen_emp, mle=PointEstimates)
      PointEstimates$bootstrapMean_efficacy = mean(boot_emp$t)
      IntervalEstimates$efficacy <- namedCL(quantile(boot_emp$t,c(alpha/2,1-alpha/2)),alpha=alpha)
    }
  }
  if(method=='GEE'){
    #GEE analysis of cluster effects
    # create model formula
    fterms = ifelse(cfunc == 'Z', 'cbind(y1,y0) ~ 1', 'cbind(y1,y0) ~ arm')
    formula <- stats::as.formula(paste(fterms, collapse = " + "))

    fit <- geepack::geeglm(formula=formula, id = cluster, corstr = "exchangeable", data=trial, family=binomial(link="logit"))
    summary_fit = summary(fit)

    z=-qnorm(alpha/2) #standard deviation score for calculating confidence intervals
    logitpC = summary_fit$coefficients[1,1]
    se_logitpC = summary_fit$coefficients[1,2]
    CL_pC = namedCL(invlogit(c(logitpC-z*se_logitpC,logitpC+z*se_logitpC)),alpha=alpha)

    # Intracluster correlation
    ICC = noLabels(summary_fit$corr[1]) #with corstr = "exchangeable", alpha is the ICC
    se_ICC = noLabels(summary_fit$corr[2])
    CL_ICC = namedCL(noLabels(c(ICC-z*se_ICC,ICC+z*se_ICC)),alpha=alpha)

    clusterSize = nrow(trial)/nlevels(as.factor(trial$cluster))
    DesignEffect = 1 + (clusterSize - 1)*ICC #Design Effect
    CL_DesignEffect = 1 + (clusterSize - 1)*CL_ICC

    # remove the temporary objects from the dataframe
    fit$ModelObject$data$y1=fit$ModelObject$data$y0=fit$ModelObject$data$y_off=NULL
    PointEstimates=list(controlP=invlogit(logitpC),
                        ICC=ICC,
                        DesignEffect=DesignEffect,
                        ModelObject=fit)
    IntervalEstimates=list(controlP=CL_pC,
                           ICC=CL_ICC,
                           DesignEffect=CL_DesignEffect)

    # Estimation of efficacy does not apply if analysis is of baseline only (cfunc='Z')
    if(cfunc=='X'){
      logitpI = summary_fit$coefficients[1,1] + summary_fit$coefficients[2,1]
      se_logitpI = sqrt(fit$geese$vbeta[1,1] + fit$geese$vbeta[2,2] + 2*fit$geese$vbeta[1,2])

      CL_pI = namedCL(invlogit(c(logitpI-z*se_logitpI,logitpI+z*se_logitpI)),alpha=alpha)
      CL_eff = estimateCLEfficacy(mu=summary_fit$coefficients[,1], Sigma=fit$geese$vbeta ,alpha=alpha, resamples=resamples)

      PointEstimates$interventionP=invlogit(logitpI)
      PointEstimates$efficacy=(1 - invlogit(logitpI)/invlogit(logitpC))
      IntervalEstimates$interventionP=CL_pI
      IntervalEstimates$efficacy=CL_eff
    }
    PointEstimates$contaminationParameter = NA #contamination is not estimated

  } else if(method == 'ML'){
    ############### ML Methods with contamination functions and logistic link #################

    par = SingleTrialAnalysis(trial=trial,FUN2=FUN2)
    PointEstimates <- FittingResults(trial, par=par,FUN1=FUN1)
    if(requireBootstrap){
      mle=list(par=par,FUN1=FUN1,link='logit')
      boot_estimates = data.frame(V1=c(),V2=c(),V3=c(),V4=c())
      # resampling can crash because of resamples containing data from only one arm
      # to prevent this crashing the whole program, bootstrapping is done in batches
      # of 5 resamples using 'try' to avoid crashing out
      resamples1=5
      tries=0
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
        tries=tries+5
        cat('\r',nrow(boot_estimates),' bootstrap samples analysed, out of',tries,' tries    \r' )
      }
      colnames(boot_estimates) = names(PointEstimates)
      IntervalEstimates = as.data.frame(t(apply(boot_estimates,2,function(x){
        namedCL(quantile(x,c(alpha/2,1-alpha/2)),alpha=alpha)})))
      IntervalEstimates <- setNames(split(IntervalEstimates,
                                          seq(nrow(IntervalEstimates))),
                                    rownames(IntervalEstimates))

    }
  } else if(method == 'INLA'){
    if(is.null(inlaMesh)){
      inlaMesh = createMesh(trial=trial,
                            offset = -0.1,
                            max.edge = 0.25,
                            inla.alpha = 2,
                            maskbuffer = 0.5,
                            ncells= 50)
    }

    y_off=NULL
    # specify functional form of sigmoid in distance from boundary
    # 'L' inverse logit; 'P' inverse probit; 'X' or 'Z' do not model contamination
    FUN = switch(cfunc, 'L' = "invlogit(x)", 'P' = "stats::pnorm(x)", 'X' = NULL, 'Z' = NULL)

    # create model formula
    fterms = switch(cfunc,
                    'Z' = 'y ~ 1',
                    'X' = 'y ~ 0 + b0',
                    'L' ='y ~ 0 + b0 + pvar',
                    'P' ='y ~ 0 + b0 + pvar')
    if(localisedEffects) fterms = c(fterms, 'b1')
    if(clusterEffects) fterms = c(fterms, 'f(cluster, model = "iid")')
    if(spatialEffects) fterms = c(fterms, 'f(s, model = spde)')
    formula_as_text = paste(fterms, collapse = " + ")
    formula <- stats::as.formula(formula_as_text)

    spde = inlaMesh$spde
    cat('Estimating scale parameter for contamination range','\n')

    effectse = list(df=data.frame(b0=rep(1, nrow(trial)),
                                  b1=ifelse(trial$arm=='intervention',-1,0),
                                  cluster = trial$cluster),
                    s = inlaMesh$indexs)
    effectsp = list(df=data.frame(b0=rep(1, nrow(inlaMesh$prediction)),
                                  b1=ifelse(inlaMesh$prediction$arm=='intervention',-1,0),
                                  cluster = inlaMesh$prediction$cluster),
                    s = inlaMesh$indexs)
    lc=NULL
    beta2=NA
    if(cfunc %in% c('L','P')){
      ### FOR DEBUGGING SKIP THE ESTIMATION OF beta2
      beta2 = 1

      #      beta2 =  stats::optimize(f=estimateContamination,
      #                       interval=c(0.01,10),
      #                       trial=trial,
      #                       FUN=FUN,
      #                       inlaMesh=inlaMesh,
      #                       formula=formula,
      #                       tol = 0.1)$minimum
      x = trial$nearestDiscord*exp(beta2)
      trial$pvar = eval(parse(text = FUN))
      effectse$df$pvar = trial$pvar
      x = inlaMesh$prediction$nearestDiscord*exp(beta2)
      inlaMesh$prediction$pvar = ifelse(cfunc == 'X', rep(NA,nrow(inlaMesh$prediction)), eval(parse(text = FUN)))
      effectsp$df$pvar = inlaMesh$prediction$pvar
      # set up linear contrasts (not required for cfunc='X' or 'Z')
      lc=NULL
      if(grepl('pvar', formula_as_text, fixed = TRUE)){
        lc <- INLA::inla.make.lincomb(b0 = 1, pvar = 1)
        if(grepl('b1', formula_as_text, fixed = TRUE)){
          lc <- INLA::inla.make.lincomb(b0 = 1, pvar = 1, b1 = 1)
        }
      } else if(grepl('b1', formula_as_text, fixed = TRUE)){
          lc <- INLA::inla.make.lincomb(b0 = 1, b1 = 1)
      }
    }
    # stack for estimation stk.e
    stk.e <- INLA::inla.stack(
      tag = "est",
      data = list(y = trial$y1, y_off=trial$y_off),
      A = list(1, A=inlaMesh$A),
      effects = effectse
    )

    # stack for prediction stk.p
    stk.p <- INLA::inla.stack(
      tag = "pred",
      data = list(y = NA, y_off = NA),
      A = list(1, inlaMesh$Ap),
      effects = effectsp
    )

    # stk.full comprises both stk.e and stk.p
    stk.full <- INLA::inla.stack(stk.e, stk.p)
    cat('Starting full INLA analysis                                      \n')

    inlaResult <- INLA::inla(formula,
                             family = "binomial",
                             Ntrials = y_off,
                             lincomb = lc,
                             control.family = list(link = "logit"),
                             data = INLA::inla.stack.data(stk.full),
                             control.fixed = list(correlation.matrix=TRUE),
                             control.predictor = list(
                               compute = TRUE,link = 1,
                               A = INLA::inla.stack.A(stk.full)),
                             control.compute = list(dic = TRUE),
    )
    # Augment the inla results list with application specific quantities

    index <- INLA::inla.stack.index(stack = stk.full, tag = "pred")$data
    inlaMesh$prediction$proportion <- invlogit(inlaResult$summary.linear.predictor[index,'0.5quant'])
    results = list(modelObject=inlaResult,inlaMesh = inlaMesh,PointEstimates = list(),IntervalEstimates = list(),method=method)
    cat(formula_as_text,'\n')

    #Compute sample-based confidence limits for intervened outcome and efficacy if intervention effects are estimated
    if(grepl('pvar', formula_as_text, fixed = TRUE) | grepl('b1', formula_as_text, fixed = TRUE)){
      # Specify the means of the variables
      mu <- c(b0= inlaResult$summary.fixed['b0','mean'],lc= inlaResult$summary.lincomb.derived['lc','mean'])
      # Specify the covariance matrix of the variables
      cov = inlaResult$misc$lincomb.derived.covariance.matrix[c('lc','b0'),c('lc','b0')]
      sample = as.data.frame(MASS::mvrnorm(n = 10000, mu = mu, Sigma = cov))
      sample$pC = invlogit(sample$b0)
      sample$pI = invlogit(sample$lc)
      sample$Es = 1 - sample$pI/sample$pC
      bounds = apply(sample,2,function(x){
        quantile(x,c(alpha/2,0.5, 1-alpha/2),alpha=alpha)})
      bounds <- data.frame(setNames(split(bounds,seq(ncol(bounds))),colnames(bounds)))

    } else {
      pC = unlist(invlogit(inlaResult$summary.fixed['b0',c("0.025quant","0.5quant","0.975quant")]))
      bounds=data.frame(pC=pC,pI=pC,Es=rep(0,3))
    }
    results$PointEstimates$controlP = bounds[2,'pC']
    results$PointEstimates$interventionP = bounds[2,'pI']
    results$PointEstimates$efficacy = bounds[2,'Es']
    results$PointEstimates$contaminationParameter=beta2

    # TODO: correct in case alpha <>0.05
    # Extract interval estimates
    results$IntervalEstimates$controlP = stats::setNames(bounds[c(1,3),'pC'],c('2.5%','97.5%'))
    results$IntervalEstimates$interventionP = stats::setNames(bounds[c(1,3),'pI'],c('2.5%','97.5%'))
    results$IntervalEstimates$efficacy = stats::setNames(bounds[c(1,3),'Es'],c('2.5%','97.5%'))

    results$description=description
  }
  if(method %in% c('EMP','ML','GEE')){
    # tidy up and consolidate the list of results
    ModelObject=PointEstimates$ModelObject
    PointEstimates <- PointEstimates[names(PointEstimates) != "GAsolution3"]
    PointEstimates <- PointEstimates[names(PointEstimates) != "ModelObject"]
    results = list(description=description,
                   method=method,
                   PointEstimates=PointEstimates,
                   IntervalEstimates=IntervalEstimates,
                   ModelObject=ModelObject)
  }
  results$contamination = getContaminationCurve(trial=trial,
                                                PointEstimates=results$PointEstimates,
                                                FUN1=FUN1)
  results$PointEstimates$contaminationRange =
    ifelse(cfunc=='X',NA,results$contamination$contaminationRange)
  results$contamination$contaminationRange = NULL
  cat('Analysis model: ',method,' Contamination option: ',cfunc,'\n')
  cat('Estimated Proportions-      Control: ',results$PointEstimates$controlP,' (95% CL: ',
      unlist(results$IntervalEstimates$controlP),')\n')
  cat('                       Intervention: ',results$PointEstimates$interventionP,' (95% CL: ',
      unlist(results$IntervalEstimates$interventionP),')\n')
  cat('Efficacy: ',results$PointEstimates$efficacy,' (95% CL: ',
      unlist(results$IntervalEstimates$efficacy),')\n')
  cat('Contamination Range: ',results$PointEstimates$contaminationRange,'\n')
  return(results)}

getContaminationCurve = function(trial, PointEstimates, FUN1){

  y_off=y1=cats=nearestDiscord=NULL

  range_d = max(trial$nearestDiscord)-min(trial$nearestDiscord)
  d = min(trial$nearestDiscord) + range_d*(seq(1:1001)-1)/1000
  par = with(PointEstimates,c(logit(controlP),
                              logit(interventionP)-logit(controlP),
                              contaminationParameter))
  curve = invlogit(FUN1(trial=data.frame(nearestDiscord=d),par=par))

  #estimate contamination range
  #The absolute value of deltaP is used so that a positive range is obtained even with negative efficacy
  deltaP <- abs(PointEstimates$controlP - PointEstimates$interventionP)
  thetaL = thetaU = NA
  if (abs(PointEstimates$controlP - curve[1000]) > (0.025*deltaP)) {
    thetaL <- d[min(which(abs(PointEstimates$controlP - curve) > (0.025*deltaP)))]
  }
  if (abs(PointEstimates$interventionP - curve[1000]) < (0.025*deltaP)) {
    thetaU <- d[min(which(abs(PointEstimates$interventionP - curve) < (0.025*deltaP)))]
  }
  if(is.na(thetaU)) thetaU=max(trial$nearestDiscord)
  if(is.na(thetaL)) thetaL=min(trial$nearestDiscord)
  #contamination range
  contaminatedInterval = c(thetaL,thetaU)
  if(thetaL > thetaU) contaminatedInterval= c(thetaU,thetaL)
  contaminationRange = thetaU - thetaL
  if(deltaP == 0) contaminationRange = NA
  # To remove warnings from plotting ensure that contamination interval is non-zero
  if(is.na(contaminationRange) || contaminationRange == 0){contaminatedInterval <- c(-0.0001,0.0001)}
  #categorisation of trial data for plotting
  trial$cats<- cut(trial$nearestDiscord, breaks = c(-Inf,min(d)+seq(1:9)*range_d/10,Inf), labels = FALSE)
  data = data.frame(trial %>%
                      group_by(cats) %>%
                      dplyr::summarize(positives = sum(y1),
                                       total = sum(y_off),
                                       d = median(nearestDiscord)))
  # proportions and binomial confidence intervals by category
  data$p = data$positives/data$total
  data$upper = with(data, p + 1.96*(sqrt(p*(1-p)/total)))
  data$lower = with(data, p - 1.96*(sqrt(p*(1-p)/total)))

  returnList = list(FittedCurve=data.frame(d=d,contaminationFunction=curve),
                    contaminationRange=contaminationRange,
                    contaminatedInterval=contaminatedInterval,
                    data=data)
  return(returnList)}

#add labels to confidence limits
namedCL=function(limits,alpha=alpha){
  names(limits)=c(paste0(100*alpha/2,'%'),paste0(100-100*alpha/2,'%'))
  return(limits)}

#logit transformation
logit = function(p=p){return(log(p/(1-p)))}

#inverse logit transformation
invlogit = function(x=x){return(1/(1+exp(-x)))}

# Minimal data description and crude efficacy estimate
get_description = function(trial){
  positives = tapply(trial$y1, trial$arm, FUN=sum)
  totals = tapply(trial$y_off, trial$arm, FUN=sum)
  ratios = positives/totals
  efficacy = 1 - ratios[2]/ratios[1]
  description = list(positives=positives,totals=totals,ratios=ratios,efficacy=efficacy,method='DE')
  return(description)}

# Log Likelihood to be maximized
LogLikelihood <- function(par, FUN=FUN ,trial) {
  logitexpectP <- FUN(par,trial)
  transf <- 1/(1+exp(-logitexpectP)) #inverse logit transformation

  # FOR BINOMIAL
  LogLikelihood <- sum(trial$y1*log(transf) + trial$y0*log(1-transf))


  return(LogLikelihood)
}

##############################################################################

FittingResults <- function(trial, FUN1, par){

  # transform the parameters into interpretable functions
  controlP <- invlogit(par[1])
  interventionP <- invlogit(par[2] + par[1])
  efficacy <- (controlP - interventionP)/controlP

  PointEstimates=list(controlP=controlP,
                      interventionP=interventionP,
                      efficacy=efficacy,
                      contaminationParameter=par[3])
  return(PointEstimates)
}

##############################################################################
#  Different functions for the linear predictor
##############################################################################

CalculateNoEffect <- function(par,trial){
  lp <- par[1]
  return(lp)
}

# step function for the case with no contamination

CalculateNoContaminationFunction <- function(par,trial){
  lp <- ifelse(trial$nearestDiscord < 0, par[1], par[1] + par[2])
  return(lp)
}


# piecewise linear model (on the logit scale) for contamination function

CalculatePiecewiseLinearFunction <- function(par,trial){
  # constrain the slope parameter to be positive (par[2] is positive if efficacy is negative)
  theta <- exp(par[3])
  lp <-ifelse(trial$nearestDiscord > -theta, par[1] +
                par[2]*(theta + trial$nearestDiscord)/(2*theta), par[1])
  lp <-ifelse(trial$nearestDiscord > theta, par[1] + par[2], lp)
  return(lp)
}


# sigmoid (logit) function (on the logit scale) for contamination function
CalculateLogisticFunction <- function(par,trial){
  theta <- exp(par[3])
  lp <- par[1] + par[2]*invlogit(theta*trial$nearestDiscord)
  return(lp)
}

# inverse probit function (on the logit scale) for contamination function
CalculateProbitFunction <- function(par,trial){
  theta <- exp(par[3])
  lp <- par[1] + par[2]*stats::pnorm(theta*trial$nearestDiscord)
  return(lp)
}

##############################################################################
# Functions for GEE analysis

noLabels=function(x){
  xclean = as.matrix(x)
  dimnames(xclean) = NULL
  xclean = as.vector(xclean)
  return(xclean)}

estimateCLEfficacy = function(mu, Sigma,alpha=alpha, resamples=resamples){
  # Use resampling approach to avoid need for complicated Taylor approximation
  # use at least 10000 samples (this is very cheap)
  resamples1 = max(resamples,10000,na.rm = TRUE)
  samples = MASS::mvrnorm(n = resamples1,mu=mu,Sigma=Sigma)
  pC = invlogit(samples[,1])
  pI = invlogit(samples[,1] + samples[,2])
  eff = 1 - pI/pC
  CL = quantile(eff, probs = c(alpha/2, 1-alpha/2))
  return(CL)}

##############################################################################
#functions for analysis of Maximum Likelihood models

rgen<-function(data,mle){
  par=mle$par
  FUN1=mle$FUN1
  #simulate data for numerator y1
  modelp <- FUN1(par=par,trial=data)
  if(mle$link == 'logit') {
    transf = invlogit(modelp)
  } else {
    transf = modelp
  }
  data$y1 <- rbinom(length(transf),data$y_off,transf) #simulate from binomial distribution
  return(data)
}

SingleTrialAnalysis <- function(trial, FUN2=FUN2) {

  GA <- GA::ga("real-valued", fitness = LogLikelihood, FUN=FUN2,
               trial=trial,
               lower = c(-10,-10,-100), upper = c(10,10,100),
               maxiter = 500, run = 50, optim = TRUE,monitor = FALSE)
  result <- GA@solution

  return(result)
}



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


rgen_emp <- function(data,mle){

  description <- psych::describeBy(data$y1/data$y_off, group=data$arm)
  pChat <- description$control$mean
  pIhat <- description$intervention$mean

  #simulate data for numerator num
  modelp <- ifelse(as.numeric(data$arm) - 1 > 0, pIhat,pChat)
  data$y1 <- rbinom(length(modelp),data$y_off,modelp) #simulate from Binomial distribution

  return(data)
}

# standard non-model based analyses
BootEmpiricalAnalysis <- function(resampledData){
  description <- psych::describeBy(resampledData$y1/resampledData$y_off, group=resampledData$arm)
  #reports summary statistic by a grouping variable
  pChat <- description$control$mean
  pIhat <- description$intervention$mean
  Eshat <- 1 -  pIhat/pChat

  return(Eshat)
}

########## FUNCTIONS FOR SPATIAL PARTIAL DIFFERENTIAL EQUATION MODEL IN INLA

#' \code{createMesh} Create prediction mesh and other inputs required for INLA analyis of a CRT.
#' @param trial trial dataframe including locations, clusters, arms, and binary outcomes
#' @param offset (see inla.mesh.2d documentation)
#' @param max.edge (see inla.mesh.2d documentation)
#' @param inla.alpha parameter related to the smoothness
#' @param maskbuffer (see inla.mesh.2d documentation)
#' @param ncells resolution of mesh in terms of maximum of linear dimension in pixels
#' @return list containing the mesh
#' \itemize{
#' \item \code{prediction}: Data.table containing the prediction points and covariate values
#' \item \code{A}: projection matrix from the observations to the mesh nodes.
#' \item \code{Ap}: projection matrix from the prediction points to the mesh nodes.
#' \item \code{indexs}:  index set for the SPDE model
#' \item \code{spde}: SPDE model
#' }
#' @export
#'
#' @examples
#' # low resolution mesh for test dataset
#' exampleMesh=createMesh(trial=test_Simulate_CRT,ncells=10)
createMesh = function(trial = trial,
                      offset = -0.1,
                      max.edge = 0.25,
                      inla.alpha = 2,
                      maskbuffer = 0.5,
                      ncells= 50){
  cat('Creating mesh for INLA analysis: resolution parameter= ',ncells)
  # create buffer around area of points
  trial.coords = base::matrix(c(trial$x,trial$y),ncol=2)
  sptrial = sp::SpatialPoints(trial.coords)
  buf1 <- rgeos::gBuffer(sptrial, width=maskbuffer, byid=TRUE)
  buffer <- rgeos::gUnaryUnion(buf1)

  # estimation mesh construction

  mesh <- INLA::inla.mesh.2d(
    boundary = buffer, offset = offset,
    cutoff = 0.05, max.edge = max.edge
  )

  # set up SPDE (Stochastic Partial Differential Equation) model
  spde <- INLA::inla.spde2.matern(mesh = mesh, alpha = inla.alpha, constr = TRUE)
  indexs <- INLA::inla.spde.make.index("s", spde$n.spde)
  A <- INLA::inla.spde.make.A(mesh = mesh, loc = trial.coords)

  # 8.3.6 Prediction data from https://www.paulamoraga.com/book-geospatial/sec-geostatisticaldatatheory.html
  bb <- sp::bbox(buffer)
  x <- seq(bb[1, "min"] - 1, bb[1, "max"] + 1, length.out = ncells)
  y <- seq(bb[2, "min"] - 1, bb[2, "max"] + 1, length.out = ncells)
  pred.coords <- as.matrix(expand.grid(x, y))
  buf.coords = buffer@polygons[[1]]@Polygons[[1]]@coords
  ind <- sp::point.in.polygon(
    pred.coords[, 1], pred.coords[, 2],
    buf.coords[, 1], buf.coords[, 2]
  )
  #prediction locations
  pred.coords <- pred.coords[which(ind == 1), ]

  #projection matrix for the prediction locations
  Ap <- INLA::inla.spde.make.A(mesh = mesh, loc = pred.coords)

  # Distance matrix calculations for the prediction stack
  # Create all pairwise comparisons
  pairs = tidyr::crossing(row=seq(1:nrow(pred.coords)),col=seq(1:nrow(trial)))
  # Calculate the distances
  calcdistP = function(row, col) sqrt((trial$x[col] - pred.coords[row,1])^2 + (trial$y[col] - pred.coords[row,2])^2)
  distP = apply(pairs, 1, function(y) calcdistP(y['row'],y['col']))
  distM = matrix(distP,nrow=nrow(pred.coords),ncol=nrow(trial),byrow=TRUE)
  nearestNeighbour = apply(distM,1,function(x)return(array(which.min(x))))
  armP = trial$arm[nearestNeighbour]
  clusterP = trial$cluster[nearestNeighbour]
  prediction = data.frame(x=pred.coords[,1],
                          y=pred.coords[,2],
                          nearestNeighbour=nearestNeighbour,
                          arm=armP,
                          cluster=clusterP)
  prediction <- with(prediction, prediction[order(y,x),])
  prediction$shortestDistance= apply(distM,1,min)
  calcNearestDiscord = function(x){
    discords = (trial$arm != prediction$arm[x])
    nearestDiscord = min(distM[x,discords])
    return(nearestDiscord)
  }
  rows=seq(1:nrow(prediction))
  prediction$nearestDiscord = sapply(rows,FUN=calcNearestDiscord)
  prediction$nearestDiscord = with(prediction,ifelse(arm=='control',-nearestDiscord,nearestDiscord))
  inlaMesh  = list(prediction=prediction,
                   A=A,
                   Ap=Ap,
                   indexs=indexs,
                   spde=spde)
  return(inlaMesh)}

# Use profiling to estimate beta2
estimateContamination = function(beta2 = beta2,
                                 trial=trial,
                                 FUN=FUN,
                                 inlaMesh=inlaMesh,
                                 formula=formula){
  y_off=NULL
  x = trial$nearestDiscord*exp(beta2)
  trial$pvar = -eval(parse(text = FUN))
  stk.e <- INLA::inla.stack(
    tag = "est",
    data = list(y = trial$y1, y_off=trial$y_off),
    A = list(1, A=inlaMesh$A),
    effects = list(data.frame(b0=rep(1, nrow(trial)),
                              b1=ifelse(trial$arm=='intervention',-1,0),
                              pvar = trial$pvar,
                              cluster = trial$cluster),
                   s = inlaMesh$indexs)
  )
  # run the model with just the estimation stack (no predictions needed at this stage)
  result.e <- INLA::inla(formula,
                         family = "binomial", Ntrials = y_off,
                         control.family = list(link = "logit"),
                         data = INLA::inla.stack.data(stk.e),
                         control.predictor = list(
                           compute = TRUE,link = 1,
                           A = INLA::inla.stack.A(stk.e)),
                         control.compute = list(dic = TRUE))
  loss = result.e$dic$family.dic
  cat("\rDIC: ",loss," Contamination parameter: ",beta2,"  \r")
  return(loss)}




