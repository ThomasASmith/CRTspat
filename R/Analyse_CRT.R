#' Analysis of cluster randomized trial with contamination
#'
#' \code{Analyse_CRT} returns outputs from a statistical analysis of a cluster randomized trial (CRT).
#' @param trial trial dataframe including locations, clusters, arms, and binary outcomes
#' @param method statistical method used to analyse trial.
#' options are:
#' 'EMP'  : empirical,
#' 'ML'   : maximum likelihood,
#' 'GEE'  : generalised estimating equations
#' 'LR'   : INLA model linear in contamination function.
#' 'CRE'  : INLA model linear in contamination function with cluster random effects
#' 'SPDE' : INLA model linear in contamination function with Matern spatial process
#' 'SPCRE': INLA model linear in contamination function with cluster random effects and Matern spatial process
#' @param cont transformation defining the contamination function
#' options are:
#' 'S': piecewise linear (slope),
#' 'L': inverse logistic (sigmoid),
#' 'P': inverse probit,
#' 'X': contamination not modelled'
#' @param excludeBuffer exclude any buffer zone (records with buffer=TRUE) from the analysis
#' @param requireBootstrap logical indicator of whether bootstrap confidence intervals are required
#' @param alpha confidence level for confidence intervals and credible intervals
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
#' # Standard GEE analysis of test dataset
#' exampleGEE=Analyse_CRT(trial=test_Simulate_CRT,method='GEE')

Analyse_CRT <- function(trial,
                        method='ML',
                        cont='L',
                        excludeBuffer=FALSE,
                        requireBootstrap=FALSE,
                        alpha = 0.05,
                        resamples=1000,
                        inlaMesh=NULL){

  ##############################################################################
  # MAIN FUNCTION CODE STARTS HERE
  cluster=NULL

  if("buffer" %in% colnames(trial) & excludeBuffer)
  {
    trial = trial[!trial$buffer,]
  }
  #trial needs to be ordered for some analyses
  trial <- trial[order(trial$cluster),]

  # if nearestDiscord is not provided augment the trial data frame with distance to nearest discordant coordinate
  # (specifyBuffer assigns a buffer only if a buffer width is > 0 is input)
  if(is.null(trial$nearestDiscord)) {trial <- Specify_CRTbuffer(trial=trial,bufferWidth=0)}

  PointEstimates=IntervalEstimates=ModelObject=list()
  sd = 0.5/(qnorm(1-alpha)*sqrt(2)) #initial value used in bootstrap calculations
  trial$neg=trial$denom - trial$num  #count of negatives for use in geeglm formulae

  description= get_description(trial)

  # Specify the function used for calculating the linear predictor
  if(method=='EMP' | method=='GEE') cont = 'X'
  LPfunction = c('CalculateNoContaminationFunction',
                 'CalculatePiecewiseLinearFunction',
                 'CalculateLogisticFunction',
                 'CalculateProbitFunction')[which (cont == c('X','S','L','P'))]
  FUN2 <- FUN1 <- eval(parse(text=LPfunction))

  if(method=='EMP'){
    # empirical analysis that ignores contamination
    PointEstimates <- EmpiricalAnalysis(trial)
    PointEstimates$contaminationParameter = NA #contamination is not estimated
    if(requireBootstrap){
      boot_an05 <- boot::boot(data=trial, statistic=BootEmpiricalAnalysis,
                              R=resamples, sim="parametric", ran.gen=rgen_an05, mle=PointEstimates)
      PointEstimates$bootstrapMean_efficacy = mean(boot_an05$t)
      IntervalEstimates$efficacy <- namedCL(quantile(boot_an05$t,c(alpha/2,1-alpha/2)),alpha=alpha)
    }
  } else if(method=='GEE'){
    #GEE analysis with cluster effects
    ModelObject <- GEEAnalysis(trial,alpha=alpha,resamples=resamples)
    PointEstimates <- ModelObject$PointEstimates
    PointEstimates$contaminationParameter = NA #contamination is not estimated
    IntervalEstimates <- ModelObject$IntervalEstimates
    if(requireBootstrap){
      ml <- geepack::geeglm(cbind(num,neg) ~ arm, id = cluster, corstr = "exchangeable", data=trial, family=binomial(link="logit"))

      boot_output <- boot::boot(data=trial, statistic=BootGEEAnalysis,
                                R=resamples, sim="parametric", ran.gen=rgen_GEE, mle=ml)

      PointEstimates$bootstrapMean_efficacy = mean(boot_output$t)
      IntervalEstimates$bootstrapEfficacy <- namedCL(quantile(boot_output$t,c(alpha/2,1-alpha/2)),alpha=alpha)
    }
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
      IntervalEstimates =  computeIntervals(df=boot_estimates,alpha=alpha)
    }
  } else if(method %in% c('LR','CRE','SPDE','SPCRE')){
      results = inlaModel(trial=trial,cont=cont,method=method,alpha=alpha,inlaMesh=inlaMesh)
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
  cat('\n=====================    ANALYSIS OF CLUSTER RANDOMISED TRIAL    =================\n')
  cat('Analysis model: ',method,' Contamination option: ',cont,'\n')
  cat('Estimated Proportions- Control: ',round(results$PointEstimates$controlP,digits=3),' (95% CL: ',
      round(results$IntervalEstimates$controlP,digits=3),
      ') - Intervention: ',round(results$PointEstimates$interventionP,digits=3),' (95% CL: ',
      round(results$IntervalEstimates$interventionP,digits=3),')\n')
  cat('Efficacy: ',round(results$PointEstimates$efficacy,digits=3),' (95% CL: ',
      round(results$IntervalEstimates$efficacy,digits=3),')\n')
  results$contamination = getContaminationCurve(trial=trial,
                          PointEstimates=results$PointEstimates,
                          FUN1=FUN1)
return(results)}

getContaminationCurve = function(trial, PointEstimates, FUN1){

  denom=num=cats=nearestDiscord=NULL

  range_d = max(trial$nearestDiscord)-min(trial$nearestDiscord)
  d = min(trial$nearestDiscord) + range_d*(seq(1:1001)-1)/1000
  par = with(PointEstimates,c(logit(controlP),
                               logit(interventionP)-logit(controlP),
                               contaminationParameter))
  curve = invlogit(FUN1(trial=data.frame(nearestDiscord=d),par=par))

  #estimate contamination range
  #The absolute value of deltaP is used so that a positive range is obtained even with negative efficacy
  deltaP <- abs(PointEstimates$controlP - PointEstimates$interventionP)
  thetaL <- d[which(abs(curve - PointEstimates$controlP) > 0.025*deltaP)][1]
  thetaU <- d[which(abs(curve - PointEstimates$interventionP) < 0.025*deltaP)][1]
  #contamination range
  contaminatedInterval <- c(thetaL,thetaU)
  contaminationRange = thetaU - thetaL
  # To remove warnings from plotting ensure that contamination range is non-zero
  if(!is.na(contaminationRange)){contaminatedInterval <- c(-0.0001,0.0001)}
  #categorisation of trial data for plotting
  trial$cats<- cut(trial$nearestDiscord, breaks = c(-Inf,min(d)+seq(1:9)*range_d/10,Inf), labels = FALSE)
  data = data.frame(trial %>%
    group_by(cats) %>%
    dplyr::summarize(positives = sum(num),
                     total = sum(denom),
                     d = median(nearestDiscord)))
    # proportions and binomial confidence intervals by category
    data$p = data$positives/data$total
    data$upper = with(data, p + 1.96*(sqrt(p*(1-p)/total)))
    data$lower = with(data, p - 1.96*(sqrt(p*(1-p)/total)))
  cat('Contamination Range',contaminationRange,'\n')

  returnList = list(FittedCurve=data.frame(d=d,contaminationFunction=curve),
                  contaminationRange=contaminationRange,
                  contaminatedInterval=contaminatedInterval,
                  data=data)
return(returnList)}

inlaModel = function(trial,cont,method,alpha=0.05,inlaMesh=NULL){
    if(is.null(inlaMesh)){
      inlaMesh = createMesh(trial=trial,
                            offset = -0.1,
                            max.edge = 0.25,
                            inla.alpha = 2,
                            maskbuffer = 0.5,
                            ncells= 50)
    }

    denom=NULL
    # specify functional form of sigmoid in distance from boundary
    # 'L' inverse logit; 'P' inverse probit; 'X' do not model contamination
    FUN = switch(cont, 'L' = "invlogit(-x)", 'P' = "stats::pnorm(-x)", 'X' = NULL)

    # create model formula
    fterms = ifelse(cont == 'X', 'y ~ 0 + b0', 'y ~ 0 + b0 + pvar')

    fterms = c(fterms, switch(method,
                              'CRE' = 'f(cluster, model = "iid")',
                              'SPDE' = 'f(s, model = spde)',
                              'SPCRE' = 'f(s, model = spde) + f(cluster, model = "iid")'))

    # if (method == 'LR') the formula is unchanged
    formula <- as.formula(paste(fterms, collapse = " + "))

    spde = inlaMesh$spde
    cat('Estimating scale parameter for contamination range','\n')
    effectse = list(df=data.frame(b0=rep(1, nrow(trial)),
                    cluster = trial$cluster),
                    s = inlaMesh$indexs)
    effectsp = list(df=data.frame(b0=rep(1, nrow(inlaMesh$prediction)),
                    cluster = inlaMesh$prediction$cluster),
                   s = inlaMesh$indexs)
    lc=NULL
    beta2=NA
    if(cont != 'X'){
      beta2 =  stats::optimize(f=estimateContamination,
                       interval=c(0.1,10),
                       trial=trial,
                       FUN=FUN,
                       inlaMesh=inlaMesh,
                       formula=formula,
                       tol = 0.01)$minimum
      x = trial$nearestDiscord*beta2
      trial$pvar = eval(parse(text = FUN))
      effectse$df$pvar = trial$pvar
      x = inlaMesh$prediction$nearestDiscord*beta2
      inlaMesh$prediction$pvar = ifelse(cont == 'X', rep(NA,nrow(inlaMesh$prediction)), eval(parse(text = FUN)))
      effectsp$df$pvar = inlaMesh$prediction$pvar
      # set up linear contrasts (not required for cont='X')
      lc <- INLA::inla.make.lincomb(b0 = 1, pvar = 1)
    }

    # stack for estimation stk.e
    stk.e <- INLA::inla.stack(
      tag = "est",
      data = list(y = trial$num, denom=trial$denom),
      A = list(1, A=inlaMesh$A),
      effects = effectse
    )

    # stack for prediction stk.p
    stk.p <- INLA::inla.stack(
      tag = "pred",
      data = list(y = NA, denom = NA),
      A = list(1, inlaMesh$Ap),
      effects = effectsp
    )

    # stk.full comprises both stk.e and stk.p
    stk.full <- INLA::inla.stack(stk.e, stk.p)
    cat('Starting full INLA analysis','\n')

    inlaResult <- INLA::inla(formula,
                   family = "binomial",
                   Ntrials = denom,
                   lincomb = lc,
                   control.family = list(link = "logit"),
                   data = INLA::inla.stack.data(stk.full),
                   control.predictor = list(
                     compute = TRUE,link = 1,
                     A = INLA::inla.stack.A(stk.full)),
                   control.compute = list(dic = TRUE),
                   INLA::control.inla(strategy = 'simplified.laplace', huge = TRUE) #this is to make it run faster
    )
    # Augment the inla results list with application specific quantities

    index <- INLA::inla.stack.index(stack = stk.full, tag = "pred")$data
    inlaMesh$prediction$proportion <- invlogit(inlaResult$summary.linear.predictor[index,'0.5quant'])
    results = list(modelObject=inlaResult,inlaMesh = inlaMesh,PointEstimates = list(),IntervalEstimates = list(),method=method)

    beta0= inlaResult$summary.fixed['b0',c("0.025quant","0.5quant","0.975quant")]
    beta1= ifelse (cont != 'X', inlaResult$summary.fixed['pvar',c("0.025quant","0.5quant","0.975quant")],rep(0,3))

    # This is a different parameterisation from that used for the ML models
    pC = unlist(invlogit(beta0 + beta1))
    pI = unlist(invlogit(beta0))
    Es = 1 - pI/pC

    results$PointEstimates$controlP = unname(pC[2])
    results$PointEstimates$interventionP = unname(pI[2])
    results$PointEstimates$efficacy = unname(Es[2])
    results$PointEstimates$contaminationParameter=-beta2

    # Extract interval estimates 2.5%
    # these are 95% intervals irrespective of alpha

    results$IntervalEstimates$controlP = stats::setNames(c(pC[1],pC[3]),c('2.5%','97.5%'))
    results$IntervalEstimates$interventionP = stats::setNames(c(pI[1],pI[3]),c('2.5%','97.5%'))
    results$IntervalEstimates$efficacy = stats::setNames(c(Es[1],Es[3]),c('2.5%','97.5%'))

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
  cluster=y=NULL
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
namedCL=function(limits,alpha=alpha){
  names(limits)=c(paste0(100*alpha/2,'%'),paste0(100-100*alpha/2,'%'))
  return(limits)}

#logit transformation
logit = function(p=p){return(log(p/(1-p)))}

#inverse logit transformation
invlogit = function(x=x){return(1/(1+exp(-x)))}

# Minimal data description and crude efficacy estimate
get_description = function(trial){
  positives = tapply(trial$num, trial$arm, FUN=sum)
  totals = tapply(trial$denom, trial$arm, FUN=sum)
  ratios = positives/totals
  efficacy = 1 - ratios[2]/ratios[1]
  description = list(positives=positives,totals=totals,ratios=ratios,efficacy=efficacy,method='DE')
  return(description)}

# Log Likelihood to be maximized
LogLikelihood <- function(par, FUN=FUN ,trial) {
  logitexpectP <- FUN(par,trial)
  transf <- 1/(1+exp(-logitexpectP)) #inverse logit transformation

  # FOR BINOMIAL
  LogLikelihood <- sum(trial$num*log(transf) + trial$neg*log(1-transf))


  return(LogLikelihood)
}

##############################################################################

FittingResults <- function(trial, FUN1, par){

  # transform the parameters into interpretable functions
  controlP <- 1/(1+exp(-par[1]))
  interventionP <- 1/(1+exp(-(par[2] + par[1])))
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

# step function for the case with no contamination

CalculateNoContaminationFunction <- function(par,trial){
  lp <- ifelse(trial$nearestDiscord < 0, par[1], par[1] + par[2])
  return(lp)
}


# piecewise linear model (on the logit scale) for contamination function

CalculatePiecewiseLinearFunction <- function(par,trial){

  theta <- par[3]
  lp <- par[1] + (par[2]/(2*theta))*(theta - trial$nearestDiscord)
  lp <- ifelse((trial$nearestDiscord < -theta),par[2]+par[1],lp)
  lp <- ifelse((trial$nearestDiscord > theta),par[1],lp)
  return(lp)
}


# sigmoid (logit) function (on the logit scale) for contamination function
CalculateLogisticFunction <- function(par,trial){

  lp <- par[1] + par[2]/(1 + exp(par[3]*(trial$nearestDiscord)))
  return(lp)
}

# inverse probit function (on the logit scale) for contamination function
CalculateProbitFunction <- function(par,trial){

  lp <- par[1] + par[2]*stats::pnorm(-par[3]*(trial$nearestDiscord))
  return(lp)
}

##############################################################################

GEEAnalysis <- function(trial,alpha=alpha, resamples=resamples){

  cluster=NULL
  fit <- geepack::geeglm(cbind(num,neg) ~ arm, id = cluster, corstr = "exchangeable", data=trial, family=binomial(link="logit"))
  summary_fit = summary(fit)
  logitpC = summary_fit$coefficients[1,1]
  logitpI = summary_fit$coefficients[1,1] + summary_fit$coefficients[2,1]
  se_logitpC = summary_fit$coefficients[1,2]
  se_logitpI = sqrt(fit$geese$vbeta[1,1] + fit$geese$vbeta[2,2] + 2*fit$geese$vbeta[1,2])
  z=-qnorm(alpha/2) #standard deviation score for calculating confidence intervals
  CL_pC = namedCL(invlogit(c(logitpC-z*se_logitpC,logitpC+z*se_logitpC)),alpha=alpha)
  CL_pI = namedCL(invlogit(c(logitpI-z*se_logitpI,logitpI+z*se_logitpI)),alpha=alpha)
  CL_eff = estimateCLEfficacy(mu=summary_fit$coefficients[,1], Sigma=fit$geese$vbeta ,alpha=alpha, resamples=resamples)

  # Intracluster correlation
  ICC = noLabels(summary_fit$corr[1]) #with corstr = "exchangeable", alpha is the ICC
  se_ICC = noLabels(summary_fit$corr[2])
  CL_ICC = namedCL(noLabels(c(ICC-z*se_ICC,ICC+z*se_ICC)),alpha=alpha)

  clusterSize = nrow(trial)/nlevels(as.factor(trial$cluster))
  DesignEffect = 1 + (clusterSize - 1)*ICC #Design Effect
  CL_DesignEffect = 1 + (clusterSize - 1)*CL_ICC

  PointEstimates=list(controlP=invlogit(logitpC),
                      interventionP=invlogit(logitpI),
                      efficacy=(1 - invlogit(logitpI)/invlogit(logitpC)),
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

#functions for analysis of Maximum Likelihood models

rgen<-function(data,mle){
  par=mle$par
  FUN1=mle$FUN1
  #simulate data for numerator num
  modelp <- FUN1(par=par,trial=data)
  if(mle$link == 'logit') {
    transf = invlogit(modelp)
  } else {
    transf = modelp
  }
  data$num <- rbinom(length(transf),data$denom,transf) #simulate from binomial distribution
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

computeIntervals = function(df,alpha){
  varnames= colnames(df)
  IntervalEstimates = as.list(varnames)
  for(i in 1:length(varnames)){
    IntervalEstimates[[i]]=namedCL(quantile(df[,varnames[i]],c(alpha/2,1-alpha/2)),alpha=alpha)
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

  cluster=NULL
  resampledData <- resampledData[order(resampledData$cluster),]
  fit <- geepack::geeglm(cbind(num,neg) ~ arm, id = cluster, corstr = "exchangeable", data=resampledData, family=binomial(link="logit"))

  pChat <- exp(fit$coefficients[1])/(1+exp(fit$coefficients[1]))
  pIhat <- exp(sum(fit$coefficients))/(1+exp(sum(fit$coefficients)))
  Eshat <- 1 - pIhat/pChat

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
#' exampleMesh=createMesh(trial=test_Simulate_CRT,ncells=20)
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
  denom=NULL
  x = trial$nearestDiscord*beta2
  trial$pvar = eval(parse(text = FUN))
  stk.e <- INLA::inla.stack(
    tag = "est",
    data = list(y = trial$num, denom=trial$denom),
    A = list(1, A=inlaMesh$A),
    effects = list(data.frame(b0=rep(1, nrow(trial)),
                              pvar = trial$pvar,
                              cluster = trial$cluster),
                   s = inlaMesh$indexs)
  )
  # run the model with just the estimation stack (no predictions needed at this stage)
  result.e <- INLA::inla(formula,
                   family = "binomial", Ntrials = denom,
                   control.family = list(link = "logit"),
                   data = INLA::inla.stack.data(stk.e),
                   control.predictor = list(
                     compute = TRUE,link = 1,
                     A = INLA::inla.stack.A(stk.e)),
                   control.compute = list(dic = TRUE),
                   INLA::control.inla(strategy = 'simplified.laplace', huge = TRUE) #this is to make it run faster
  )
  loss = result.e$dic$family.dic
  cat("DIC: ",loss," Contamination parameter: ",beta2,"\n")
  return(loss)}




