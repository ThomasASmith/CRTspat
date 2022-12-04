#' Simulation of cluster randomized trial with contamination
#'
#' \code{Simulate_CRT} generates simulated data for a cluster randomized trial (CRT) with geographic contamination between arms and a binary outcome. Contamination is simulated as arising from a diffusion-like process.
#' For details see Multerer (PhD thesis).
#'
#' @param trial a dataframe containing locations (x,y), cluster assignments, and arm assignments
#' @param efficacy simulated efficacy (defaults to 0)
#' @param initialPrevalence prevalence in control arm (assumed equal to initial proportion)
#' @param generateBaseline logical indicator of whether baseline data should be simulated
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
#' example_simulated_CRT =  Simulate_CRT(trial=CRTspillover::testArms,
#'                                      efficacy=0.25,
#'                                      ICC_inp=0.05,
#'                                      initialPrevalence=0.5,
#'                                      sd=0.6,
#'                                      tol=0.05)
Simulate_CRT = function(trial=NULL,
                        efficacy=NULL,
                        initialPrevalence=NULL,
                        generateBaseline=TRUE,
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
          ##############################################################################
          # simulate baseline data with a specified ICC

            # parameters of the spatial pattern of prevalence

            initial_guess_bandwidth <- 0.0375/ICC_inp #bandwidth for KDE. The value of 0.0375 is from analysis of
            # a sample of simulations based on Rusinga geography

            # compute approximate diagonal of clusters
            approx_diag = sqrt((max(trial$x)-min(trial$x))^2 + (max(trial$y)-min(trial$y))^2)/sqrt(length(unique(trial$cluster)))

            # scale bandwidth by cluster diagonal
            scaled_initial_guess = initial_guess_bandwidth/approx_diag

            # number of positives required to match the specified prevalence
            npositives = round(initialPrevalence*nrow(trial))

            # For the smoothing step compute contributions to the relative effect size
            # from other households as a function of distance to the other households

            euclid <- distance_matrix(trial$x, trial$y)

            bw = scaled_initial_guess*c(0.25,0.5,1,2,4)
            iter = 0
            deviation = 999
            ICCs = c()
            cat("Estimating the smoothing required to achieve the target ICC","\n")
            while(abs(deviation) > tol){
              iter=iter+1
              #use a linear model fitted to log(bw) to estimate the required bandwidth
              if(iter > 5){
                lp = ICCs - ICC_inp
                logbw_i = min(2,lm(formula = log(bw) ~ lp)$coefficients[1]) #impose a maximum on bw
                bw = c(bw, exp(logbw_i))
              }
              #assign initial pattern
              trial$infectiousness_proxy <- KDESmoother(trial$x,trial$y,kernnumber=200,bandwidth=bw[iter],low=0.0,high=1.0)

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
              baseline_analysis = Analyse_baseline(trial=trial,
                                                   baselineNumerator='base_num',
                                                   baselineDenominator='base_denom',
                                                   method='aovs',
                                                   ci.type=NULL)
              ICC = baseline_analysis$estimates$ICC
              ICCs = c(ICCs,ICC)
              cat("bandwidth: ",bw[iter]*approx_diag,"  ICC=",ICC,"\n")
              deviation = ICC - ICC_inp
            }
        }
    }

    # Assign expected proportions to each location assuming a fixed efficacy.

    # Indicator of whether the source is intervened is (as.numeric(trial$arm[i]) - 1
    # smoothedIntervened is the value of infectiousness_proxy smoothed to allow for mosquito movement and decremented
    # by the effect of intervention and smoothed to allow

    euclid <- distance_matrix(trial$x, trial$y)
    smoothedIntervened <- gauss(sd, euclid) %*% (trial$infectiousness_proxy * (1 -efficacy*(as.numeric(trial$arm) - 1)))

    # distances to nearest discordant households
    discord <- outer(trial$arm, trial$arm, "!=") #returns true & false.
    euclidd <- ifelse(discord,euclid,99999.9) #ifelse(test_expression,yes,no), yes returns the elements for where test_expression is true.
    trial$nearestDiscord <- ifelse(trial$arm == 'control', - apply(euclidd, MARGIN=2, min), apply(euclidd, MARGIN=2, min))
    #for the households with intervention return the minimal distance with a minus sign

    # npositives is the total number of positives to be distributed among the households
    npositives = round(initialPrevalence*nrow(trial) * (1 - 0.5 * efficacy))

    # scale to input value of initial prevalence by assigning required number of infections with probabilities proportionate
    # to smoothedIntervened
    positives = sample(x=nrow(trial),size=npositives,replace=FALSE,prob=smoothedIntervened)
    trial$denom=1
    trial$num=0
    trial$num[positives]=1
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
  Y = rlnorm(length(XY), meanlog = mu, sdlog = sqrt(var))
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


