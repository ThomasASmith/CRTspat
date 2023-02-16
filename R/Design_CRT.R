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
#' @param postulatedContaminationRange contamination range in km, obtained from other studies
#' @param coordinates dataframe containing coordinates of households. Columns 'x' and 'y' should contain Cartesian (x,y) coordinates. Units are expected to be km.
#' @param h  proposal for the number of coordinates in each cluster
#' @param algo algorithm for cluster boundaries, choose between
#' "TSP": travelling salesman problem heuristic;
#' "NN": nearest neighbor;
#' "kmeans": kmeans algorithm
#' @param reuseTSP indicator of whether a pre-existing path should be used by the TSP algorithm
#' @return A list comprising a list the following elements:
#' \itemize{
#' \item \code{arm}: vector of cluster assignments by numerical order of cluster to trial arms
#' \item \code{alpha}: confidence level
#' \item \code{desiredPower}: desired power
#' \item \code{effect}: Required effect size
#' \item \code{ICC}: Intra-Cluster Correlation obtained from other studies
#' \item \code{nominalDE}: calculated Design Effect
#' \item \code{pC}: baseline prevalence
#' \item \code{n_ind}: required individuals per arm in an individually randomized trial
#' \item \code{postulatedContaminationRange}: contamination range in km, obtained from other studies
#' \item \code{h}: proposed number of households in each cluster
#' \item \code{algo}: algorithm used for cluster boundaries
#' \item \code{min_c}: minimum number of clusters required
#' }
#' @export
#'
#' @examples
#'
#' exampleDesign = Design_CRT(coordinates=readdata('test_site.csv'),
#'                 ICC=0.10, effect=0.4, pC=0.35, postulatedContaminationRange=0.25, h=100)
Design_CRT = function(  alpha = 0.05,  #Step A: confidence level
                        desiredPower = 0.8,  #Step B: power
                        effect, #Step C: Required effect size
                        ICC,#Step D: ICC, obtained from other studies
                        pC, #Step E: baseline prevalence
                        postulatedContaminationRange = 0,  #Step F: postulated contamination range in km, obtained from other studies
                        coordinates, # Step G		coordinates of households in study area
                        h,  #Step H: proposal for the number of households in each cluster
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
  n_ind <- 2*sigma2*((Zsig + Zpow)/d)^2 #required individuals per arm in individually randomized trial

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
  trial <- Specify_CRTbuffer(trial=trial,bufferWidth=postulatedContaminationRange)

  output = list(pC = pC,
                alpha = alpha,
                n_ind=n_ind,
                desiredPower = desiredPower,
                inputClusterSize = h,
                algo = algo,
                postulatedContaminationRange = postulatedContaminationRange,
                effect = effect,
                ICC = ICC,
                h = h)

  cat('=====================CLUSTER RANDOMISED TRIAL DESIGN =================\n')
  cat('Significance level: ',alpha,'\n')
  cat('required effect size: ',effect,'\n')
  cat('assumed prevalence in absence of intervention ',pC,'\n')
  cat('pre-specified intra-cluster correlation: ',ICC,'\n\n')
  cat('=====================         FULL TRIAL AREA        =================\n')
  output$descriptionFullTrial = describeTrial(trial=trial,pC = pC, d = d, desiredPower = desiredPower,
                                              n_ind =n_ind, sigma2 = sigma2, Zsig = Zsig, ICC = ICC)
  if(postulatedContaminationRange > 0){

    cat('\n=====================    EXCLUDING BUFFER ZONES        =================\n')
    cat('buffer of width ',postulatedContaminationRange,' km.\n')
    output$descriptionCoreTrial =  describeTrial(trial=trial[trial$buffer==FALSE,], pC = pC, d = d, desiredPower = desiredPower,
                                                 n_ind =n_ind, sigma2 = sigma2, Zsig = Zsig, ICC = ICC)
    #proportion of coordinates in core
    output$descriptionCoreTrial$proportionInCore <- sum(abs(trial$nearestDiscord) >= postulatedContaminationRange)/dim(trial)[1]
    pbuff = sum(ifelse(trial$buffer,1,0))/nrow(trial)*100
    cat(pbuff,'% of locations in buffer zone.\n')
  }
  return(output)}

# Characteristics of a trial design
describeTrial = function(trial,pC, d, desiredPower, n_ind, sigma2, Zsig, ICC){

  arm=unique(trial[c("cluster", "arm")])[,2] #assignments

  k = length(arm)/2 # number of clusters assigned to each arm

  ##############################################################################
  # Step L computation of characteristics of the randomization

  sd_distance = stats::sd(trial$nearestDiscord)

  #mean number of locations randomized to each arm
  mean_h = mean(table(trial$cluster))

  #standard deviation of locations randomized to each arm
  sd_h = stats::sd(table(trial$cluster))

  # coefficient of variation of the cluster sizes
  cv_h = sd_h/mean_h

  # design effect (Variance Inflation Factor) allowing for varying cluster sizes(Hemming eqn 6)
  DE <- 1 + (cv_h^2 +1)*(mean_h-1)*ICC

  # number of individuals required per arm in CRT with equal cluster sizes
  n_crt = n_ind * DE

  # minimum numbers of clusters required assuming varying cluster sizes per arm (Hemming eqn 8)
  min_k = ceiling(n_ind*(1+((cv_h + 1)*mean_h - 1)*ICC)/mean_h)

  # power with k clusters per arm
  # power = stats::pnorm(sqrt(k*mean_h/(2*(1+(mean_h-1)*ICC)))* d/sqrt(sigma2) - Zsig)   #Hemming eqn 27

  power = stats::pnorm(sqrt(k*mean_h/(2*DE))* d/sqrt(sigma2) - Zsig)   #unequal cluster sizes

  cat('calculated design effect: ',DE,'\n')
  cat('Locations- total: ',nrow(trial),'. Per cluster mean: ',mean_h,'S.D.: ', sd_h,'\n')
  cat('S.D. of distance to nearest discordant location (km): ',sd_distance,'\n')
  cat('Minimum clusters (total over both arms) for power of ',desiredPower*100,'%: ',min_k,'. Available clusters: ',length(arm),'\n')
  if(min_k < length(arm)){
    cat("** Warning: insufficient clusters available to achieve required power **\n")
  }
  CRT_description = list(trial=trial,
                         nominalDE = DE,
                         sd_distance = sd_distance,
                         mean_h = mean_h,
                         sd_h = sd_h,
                         clustersRequired = 2*min_k,
                         clustersAssigned = length(arm),
                         power=power)

  return(CRT_description)}
