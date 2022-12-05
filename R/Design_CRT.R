#' Design a CRT of a malaria intervention with contamination
#'
#' \code{Design_CRT} estimates the required number of clusters and the extent of contamination between arms for a CRT with the input set of locations.
#' Outputs are:
#' (i) Estimates of the required numbers of clusters.
#' (ii) A proposal for the cluster and arm assignments to the input coordinates. (A warning is output if the number of locations is too small to allow randomisation of sufficient clusters).
#' (iii) the proportion of households in the input geography falling within the core of the clusters (i.e. outside the contamination range of locations in the opposite arm)
#'
#' @param alpha confidence level
#' @param desiredPower desired power
#' @param effect required effect size
#' @param ICC Intra-Cluster Correlation obtained from other studies
#' @param pC baseline prevalence
#' @param postulatedContamination contamination range in km, obtained from other studies
#' @param coordinates dataframe containing coordinates of households. Columns 'x' and 'y' should contain Cartesian (x,y) coordinates. Units are expected to be km.
#' @param h  proposal for the number of coordinates in each cluster
#' @param algo algorithm for cluster boundaries, choose between
#' "TSP": travelling salesman problem heuristic;
#' "NN": nearest neighbor;
#' "kmeans": kmeans algorithm
#' @param reuseTSP indicator of whether a pre-existing path should be used by the TSP algorithm
#' @return A trial design object comprising a list with the following attributes:
#' \itemize{
#' \item \code{arm}: vector of assignments to trial arms
#' \item \code{alpha}: confidence level
#' \item \code{desiredPower}: desired power
#' \item \code{effect}: Required effect size
#' \item \code{ICC}: Intra-Cluster Correlation obtained from other studies
#' \item \code{nominalDE}: calculated Design Effect
#' \item \code{pC}: baseline prevalence
#' \item \code{n_ind}: required individuals per arm in an individually randomized trial
#' \item \code{postulatedContamination}: contamination range in km, obtained from other studies
#' \item \code{h}: proposal for the number of households in each cluster
#' \item \code{algo}: algorithm used for cluster boundaries
#' \item \code{assignments}: data frame containing locations, clusters and arm assignments
#' \item \code{min_c}: minimum number of clusters required
#' }
#' @export
#'
#' @examples
#'
#' exampleDesign = Design_CRT(coordinates=CRTspillover::test_site,
#'                 ICC=0.10, effect=0.4, pC=0.35, postulatedContamination=0.5, h=100)
Design_CRT = function(  alpha = 0.05,  #Step A: confidence level
                        desiredPower = 0.8,  #Step B: power
                        effect = 0.4, #Step C: Required effect size
                        ICC = 0.175,  #Step D: ICC, obtained from other studies
                        pC = 0.4,     #Step E: baseline prevalence
                        postulatedContamination = 0.25,  #Step F: postulated contamination range in km, obtained from other studies
                        coordinates=CRTspillover::AvecNet_coordinates, # Step G		coordinates of households in study area
                        h = 80,       #Step H: proposal for the number of households in each cluster
                        #algorithm for cluster boundaries, choose between
                        #"TSP": travelling salesman problem heuristic; "NN": nearest neighbor; "kmeans": kmeans algorithm
                        algo = "kmeans",
                        reuseTSP=FALSE){

# convert power and significance level to normal deviates
Zsig = -qnorm(alpha/2)
Zpow = qnorm(desiredPower)

##############################################################################
# Calculations (see e.g. Hemming et al, 2011
# https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-11-102
# )
##############################################################################
# Step I: Calculations for the required minimum number of clusters for both arms

# required number of individuals in an individual RCT
pI = pC * (1- effect) # probability in intervened group
d = pC -pI # difference between groups
sigma2 = 1/2*(pI*(1-pI) + pC*(1-pC))
n_ind <- 2*sigma2*((Zsig + Zpow)/d)^2

# see below for calculations of design effect and minimum numbers of clusters required


##############################################################################
# Step J: specify or compute cluster boundaries
########################################################################################
trial = DefineClusters(trial=coordinates, h=h, algo=algo, reuseTSP=reuseTSP)

##############################################################################
#Step K: Random assignment of clusters to arms
trial = Randomize_CRT(trial)

# augment the trial data frame with distance to nearest discordant coordinate
# (specifyBuffer assigns a buffer only if a buffer width is > 0 is input)
trial <- Specify_CRTbuffer(trial=trial,bufferWidth=0)

output = list(pC = pC,
              alpha = alpha,
              n_ind=n_ind,
              desiredPower = desiredPower,
              inputClusterSize = h,
              algo = algo,
              postulatedContaminationRange = postulatedContamination,
              effect = effect,
              ICC = ICC,
              h = h,
              assignments = trial)

output$describeFullTrial = describeTrial(trial=trial,pC = pC, d = d, desiredPower = desiredPower,
                                         n_ind =n_ind, sigma2 = sigma2, Zsig = Zsig, ICC = ICC)
#proportion of coordinates in core
output$proportionInCore <- sum(abs(trial$nearestDiscord) >= postulatedContamination)/dim(trial)[1]

output$core_trial <- Specify_CRTbuffer(trial=trial,bufferWidth=postulatedContamination)

output$descriptionCoreTrial =  describeTrial(trial=output$core_trial, pC = pC, d = d, desiredPower = desiredPower,
                                             n_ind =n_ind, sigma2 = sigma2, Zsig = Zsig, ICC = ICC)
return(output)}

# Characteristics of a trial design
describeTrial = function(trial,pC, d, desiredPower, n_ind, sigma2, Zsig, ICC){

  arm=unique(trial[c("cluster", "arm")])[,2] #assignments

  k = length(arm) # number of clusters assigned

  ##############################################################################
  # Step L computation of characteristics of the randomization

  sd_distance = stats::sd(trial$nearestDiscord)

  #mean number of locations randomised to each arm
  mean_h = mean(table(trial$cluster))
  #standard deviation of locations randomised to each arm
  sd_h = stats::sd(table(trial$cluster))

  # coefficient of variation of the cluster sizes
  cv_h = sd_h/mean_h

  # design effect (Variance Inflation Factor) allowing for varying cluster sizes
  DE <- 1 + (cv_h^2 +1)*(mean_h-1)*ICC

  # number of individuals required per arm in CRT
  n_crt = n_ind * DE

  # minimum numbers of clusters required allowing for varying cluster sizes
  min_c <- ceiling(n_ind*DE/mean_h)

  print(paste0(k,' clusters assigned, ',min_c,' clusters required to achieve desired power of ', 100*desiredPower, '%'))

  power = 1 - stats::pnorm(sqrt(k/2*ICC)* d/sqrt(sigma2) - Zsig)   #Hemming eqn 28

  CRT_description = list(assignments=trial,
                         nominalDE = DE,
                         sd_distance = sd_distance,
                         mean_h = mean_h,
                         sd_h = sd_h,
                         requiredClusters = min_c,
                         availableClusters = length(arm),
                         power=power)

return(CRT_description)}
