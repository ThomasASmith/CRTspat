#' Design a CRT of a malaria intervention with contamination
#'
#' \code{Design_CRT} estimates the required number of clusters and the extent of contamination between arms for a CRT with the input set of locations.
#' Outputs are:
#' (i) Estimates of the required numbers of clusters.
#' (ii) A proposal for the cluster and arm assignments to the input coordinates. (A warning is output if the number of locations is too small to allow randomisation of sufficient clusters).
#' (iii) the proportion of households in the input geography falling within the core of the clusters (i.e. outside the contamination range of locations in the opposite arm)
#'
#' @param alpha confidence level
#' @param power power
#' @param effect required effect size
#' @param ICC Intra-Cluster Correlation obtained from other studies
#' @param pC baseline prevalence
#' @param cont contamination range in km, obtained from other studies
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
#' \item \code{power}: power
#' \item \code{effect}: Required effect size
#' \item \code{ICC}: Intra-Cluster Correlation obtained from other studies
#' \item \code{DE}: calculated Design Effect
#' \item \code{pC}: baseline prevalence
#' \item \code{cont}: contamination range in km, obtained from other studies
#' \item \code{coordinates}: dataframe containing (x,y) coordinates of locations
#' \item \code{h}: proposal for the number of households in each cluster
#' \item \code{algo}: algorithm used for cluster boundaries
#' \item \code{assignments}: data frame containing locations, clusters and arm assignments
#' \item \code{min_c}: minimum number of clusters required
#' }
#' @export
#'
#' @examples
#'
#' exampleDesign = Design_CRT(coordinates=test_site, ICC=0.10, effect=0.4, pC=0.35, cont=0.5, h=100)
Design_CRT = function(  alpha = 0.05,  #Step A: confidence level
                        power = 0.8,  #Step B: power
                        effect = 0.4, #Step C: Required effect size
                        ICC = 0.175,  #Step D: ICC, obtained from other studies
                        pC = 0.4,     #Step E: baseline prevalence
                        cont = 0.25,  #Step F: contamination range in km, obtained from other studies
                        coordinates=AvecNet_coordinates, # Step G		coordinates of households in study area
                        h = 80,       #Step H: proposal for the number of households in each cluster
                        #algorithm for cluster boundaries, choose between
                        #"TSP": travelling salesman problem heuristic; "NN": nearest neighbor; "kmeans": kmeans algorithm
                        algo = "kmeans",
                        reuseTSP=FALSE){

#  (Code originally tested on R version 4.0)
#
#  copyright Lea Multerer and Thomas Smith, 2020
#  (lea.multerer@swisstph.ch, thomas-a.smith@swisstph.ch)

# convert power and significance level to normal deviates
Zsig = -qnorm(alpha/2)
Zpow = qnorm(power)

##############################################################################
# Calculations
##############################################################################
# Step I: Calculations for the required minimum number of clusters for both arms


n_ind <- (Zsig + Zpow)^2*(pC*(1-pC) + pC*(1-effect)*(1-pC*(1-effect)))/(pC*effect)^2
DE <- 1 + (h-1)*ICC #design effect
min_c <- ceiling((1 + n_ind*DE/h)*2)

##############################################################################
# Step J: specify or compute cluster boundaries
########################################################################################
trial = DefineClusters(trial=coordinates, h=h, algo=algo, reuseTSP=reuseTSP)

##############################################################################
#Step K: Random assignment of clusters to arms
trial = Randomize_CRT(trial)

arm=unique(trial[c("cluster", "arm")])[,2]

##############################################################################
# Step L computation of characteristics of the randomization

# augment the trial data frame with distance to nearest discordant coordinate
# (specifyBuffer assigns a buffer only if a buffer width is > 0 is input)
trial <- Specify_CRTbuffer(trial=trial,bufferWidth=0)

#proportion of coordinates in core
core <- sum(abs(trial$nearestDiscord) >= cont)/dim(trial)[1]
if (min_c > length(arm)){
  print(paste0('*** WARNING: ',length(arm),' clusters assigned, ',min_c,' clusters required to achieve desired power of ', 100*power, '%.***'))
}

output_trial = list(arm=arm, alpha = alpha, power = power,
             effect = effect, ICC = ICC, DE = DE, pC = pC, cont = cont,
             coordinate_source=coordinates, h = h, algo = algo,
             core = core, assignments=trial, min_c=min_c)
return(output_trial)}

