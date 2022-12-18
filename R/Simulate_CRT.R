#' Simulation of cluster randomized trial with contamination
#'
#' \code{Simulate_CRT} generates simulated data for a cluster randomized trial (CRT) with geographic contamination between arms and a binary outcome. Contamination is simulated as arising from a diffusion-like process.
#' For details see Multerer (PhD thesis).
#'
#' @param trial a dataframe containing locations (x,y), cluster assignments, and arm assignments
#' @param efficacy simulated efficacy (defaults to 0)
#' @param initialPrevalence prevalence in control arm (assumed equal to initial proportion)
#' @param generateBaseline logical indicator of whether baseline data should be simulated
#' @param matchedPair logical: indicator of whether pair-matching on the baseline data should be used in randomization
#' @param baselineNumerator name of numerator variable for baseline data (if present)
#' @param baselineDenominator name of denominator variable for baseline data (if present)
#' @param ICC_inp Intra Cluster Correlation, provided as input when baseline data are to be simulated
#' @param sd standard deviation of the normal kernel measuring spatial smoothing leading to contamination
#' @param theta_inp input contamination range
#' @param tol tolerance of output ICC
#' @return A trial simulation object comprising a data frame containing the following numerical quantities:
#' \itemize{
#' \item \code{x}: x-coordinates of location
#' \item \code{y}: y-coordinates of location
#' \item \code{arm}: assignment to trial arm
#' \item \code{infectiousness_proxy}: infectiousness proxy
#' \item \code{nearestDiscord}: distance to nearest discordant location
#' \item \code{base_denom}: denominator for baseline
#' \item \code{base_num}: numerator for baseline
#' \item \code{denom}: number of samples evaluated at the location
#' \item \code{num}: number of positive samples (0 or 1)
#' \item \code{index}: row number
#' }
#' @export
#'
#' @examples
#' example_simulated_CRT =  Simulate_CRT(trial=CRTspillover::testArms,,
#'                                      efficacy=0.25,
#'                                      ICC_inp=0.05,
#'                                      initialPrevalence=0.5,
#'                                      matchedPair = TRUE,
#'                                      sd=0.6,
#'                                      tol=0.05)
Simulate_CRT = function(trial=NULL,
                        efficacy=0,
                        initialPrevalence=NULL,
                        generateBaseline=TRUE,
                        matchedPair = TRUE,
                        baselineNumerator='base_num',
                        baselineDenominator='base_denom',
                        ICC_inp=NULL,
                        sd=NULL,
                        theta_inp=NULL,
                        tol=0.005){

  ##############################################################################
  #  Simulation of cluster randomized trial with contamination
  #  Written by Tom Smith, July 2017
  #  Adapted by Lea Multerer, September 2017
  ##############################################################################
  cat('\n=====================    SIMULATION OF CLUSTER RANDOMISED TRIAL    =================\n')
    bw=NULL
    # use contamination range if this is available
    if (!is.null(theta_inp)){
      sd = theta_inp/(sqrt(2)*qnorm(0.95))
    }
    if (is.null(sd)){
      print('Error: contamination range or s.d. of spatial kernel must be provided')
      trial = NULL
      return(trial)
    }

    # generate baseline data if required and exposure proxy if this is not provided
    if(!"infectiousness_proxy" %in% colnames(trial) &
       baselineNumerator %in% colnames(trial) &
       baselineDenominator %in% colnames(trial)) {
       trial$infectiousness_proxy = trial[[baselineNumerator]]/trial[[baselineDenominator]]
    } else {
        if(generateBaseline) {
          # simulate baseline data with a specified ICC

          # compute approximate diagonal of clusters
          approx_diag = sqrt((max(trial$x)-min(trial$x))^2 + (max(trial$y)-min(trial$y))^2)/sqrt(length(unique(trial$cluster)))

          # number of positives required to match the specified prevalence
          npositives = round(initialPrevalence*nrow(trial))

          # For the smoothing step compute contributions to the relative effect size
          # from other households as a function of distance to the other households

          euclid <- distance_matrix(trial$x, trial$y)
          cat("Estimating the smoothing required to achieve the target ICC of",ICC_inp,"\n")
          bw =  stats::optimize(f=ICCdeviation,interval=c(0.1,10),
                                trial=trial,ICC_inp=ICC_inp,approx_diag=approx_diag,sd=sd,
                                euclid=euclid,npositives=npositives,tol = tol)$minimum
          # overprint the output that was recording progress
          cat("\r                                                         \n")
          # create a baseline dataset using the optimized bandwidth
          trial = syntheticBaseline(bw=bw,trial=trial,sd=sd,euclid=euclid,npositives=npositives)
        }
    }
    # Assign expected proportions to each location assuming a fixed efficacy.

    # Indicator of whether the source is intervened is (as.numeric(trial$arm[i]) - 1
    # smoothedIntervened is the value of infectiousness_proxy decremented
    # by the effect of intervention and smoothed to allow for mosquito movement

    euclid <- distance_matrix(trial$x, trial$y)
    smoothedIntervened <- gauss(sd, euclid) %*% (trial$infectiousness_proxy * (1 -efficacy*(as.numeric(trial$arm) - 1)))

    # distances to nearest discordant households
    discord <- outer(trial$arm, trial$arm, "!=") #returns true & false.
    euclidd <- ifelse(discord,euclid,99999.9)

    #for the control arm return the minimal distance with a minus sign
    trial$nearestDiscord <- ifelse(trial$arm == 'control', - apply(euclidd, MARGIN=2, min), apply(euclidd, MARGIN=2, min))

    # npositives is the total number of positives to be distributed among the households
    if(!("denom" %in% colnames(trial))) trial$denom=1
    npositives = round(initialPrevalence*sum(trial$denom) * (1 - 0.5 * efficacy))

    # scale to input value of initial prevalence by assigning required number of infections with probabilities proportionate
    # to smoothedIntervened multiplied by the denominator
    expected_allocation = smoothedIntervened*trial$denom/sum(smoothedIntervened*trial$denom)

    positives = sample(x=nrow(trial),size=npositives,replace=FALSE,prob=expected_allocation)
    trial$num=0
    trial$num[positives]=1
    return(trial)
}

# Function required for optimising bandwidth
ICCdeviation = function(bw=bw,
                        trial=trial,
                        ICC_inp=ICC_inp,
                        approx_diag=approx_diag,
                        sd=sd,
                        euclid=euclid,
                        npositives=npositives){
  cluster=NULL
  trial = syntheticBaseline(bw=bw,trial=trial,sd=sd,euclid=euclid,npositives=npositives)
  trial$base_neg = trial$base_denom - trial$base_num
  fit <- geepack::geeglm(cbind(base_num,base_neg) ~ 1, id = cluster, corstr = "exchangeable", data=trial, family=binomial(link="logit"))
  summary_fit = summary(fit)
  # Intracluster correlation
  ICC = noLabels(summary_fit$corr[1]) #with corstr = "exchangeable", alpha is the ICC
  cat("\rbandwidth: ",bw*approx_diag,"  ICC=",ICC,"          ")
  loss = abs(ICC - ICC_inp)
return(loss)}

syntheticBaseline = function(bw,trial,sd,euclid,npositives){
  #assign initial pattern
  trial$infectiousness_proxy <- KDESmoother(trial$x,trial$y,kernnumber=200,bandwidth=bw,low=0.0,high=1.0)

  # Smooth the exposure proxy to allow for mosquito movement
  # Note that the s.d. in each dimension of the 2 d gaussian is sd/sqrt(2)
  # smoothedBaseline is the amount received by the each cluster from the contributions (infectiousness_proxy) of each source
  smoothedBaseline <-  gauss(sd, euclid) %*% trial$infectiousness_proxy

  # scale to input value of initial prevalence by assigning required number of infections with probabilities proportionate
  # to infectiousness_proxy

  positives = sample(x=nrow(trial),size=npositives,replace=FALSE,prob=smoothedBaseline)
  trial$base_denom=1
  trial$base_num=0
  trial$base_num[positives]=1
  return(trial)
}






##############################################################################
# compute a euclidian distance matrix
distance_matrix <- function(x,y){
  xdist<- outer(x, x, '-') #replicates x to a matrix of size length(x)xlength(x) and then applys the minus
  ydist<- outer(y, y, '-')
  euclid <- sqrt(xdist*xdist+ydist*ydist) #* is elementwise multiplication
}


##############################################################################
# add lognormal noise: not sure this function is needed
# X is the input vector comprising a sample from a smoothed distribution
# varXY is the required variance

add_noise = function(X, varXY){
  muY = 1
  varY = (varXY - var(X))/(1 + var(X))
  mu = log(muY/sqrt(1+varY/muY))
  var = log(1 + varY/muY)
  Y = stats::rlnorm(length(XY), meanlog = mu, sdlog = sqrt(var))
  XY = X*Y
return(XY)}


##############################################################################
# contribution of i to j as a function of the Gaussian process
# used in simulating contamination

gauss <- function(sd,euclid){ #definition of a gauss function
  f <- (1/(2*pi*sd^2))*exp(-(euclid^2)/(2*(sd^2)))
  totalf <- rowSums(f) #sums over rows of matrix f
  # Careful here, totalf sums over very small numbers, consider replacing
  return(f/totalf)
}


##############################################################################
# generate a random pattern of vectorial capacity with smoothing


KDESmoother <- function(x,y,kernnumber,bandwidth,low,high){

  sam <- sample(1:length(x),kernnumber,replace=F)

  xdist <- outer(x, x[sam], '-')
  ydist <- outer(y, y[sam], '-')
  euclid <- sqrt(xdist*xdist+ydist*ydist) #is not a square matrix

  f <- (1/(2*pi*bandwidth^2))*exp(-(euclid^2)/(2*(bandwidth^2)))
  totalf <- (1/ncol(euclid))*rowSums(f) #sums over rows of matrix f

  smoother <- low + totalf*((high-low)/max(totalf))

  return(smoother)
}


