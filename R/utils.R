#' Aggregate multiple cluster randomized trial (CRT) records with identical location
#'
#' \code{Aggregate_CRT} aggregates data for binomial outcomes a trial dataframe contains multiple records with the same location, and outputs a dataframe with unique locations.
#' @param trial a dataframe containing locations (x,y) and optionally numerators, and corresponding denominators
#' @param auxiliaries a vector of the names of auxiliary variables to be summed across locations
#' @return A dataframe with unique locations in which the numerators and denominators are summed over the locations
#' @export
Aggregate_CRT = function(trial,auxiliaries=NULL){
  x=y=NULL
  df= trial %>% distinct(x, y, .keep_all = TRUE)
  df=df[order(df$x,df$y),]
  if(!is.null(auxiliaries)){
    for(i in 1:length(auxiliaries)){
      varName=auxiliaries[i]
      df1 = trial %>%
        group_by(x,y) %>%
        summarize(sumVar = sum(get(varName),na.rm = TRUE), .groups = 'drop')
      df1=df1[order(df1$x,df1$y),]
      df[[varName]] = df1$sumVar
    }
  }
  return(df)}


#' Randomize a two-armed cluster randomized trial
#'
#' \code{Randomize_CRT} carries out randomization of clusters for a CRT and augments the trial dataframe with assignments to arms
#'
#' @param trial dataframe containing locations and cluster assignments
#' @param matchedPair logical: indicator of whether pair-matching on the baseline data should be used in randomization
#' @param baselineNumerator name of numerator variable for baseline data (if present)
#' @param baselineDenominator name of denominator variable for baseline data (if present)
#' @return The input dataframe augmented with variable arm coded 'control' or 'intervention'
#' @export
#'
#' @examples
#' #Randomize the clusters in the example trial
#' set.seed(352)
#' exampletrial=Randomize_CRT(testClusters)
Randomize_CRT = function(trial,
    matchedPair = TRUE,
    baselineNumerator='base_num',
    baselineDenominator='base_denom'){

  cluster= base_num= base_denom= NULL
  trial$cluster = as.factor(trial$cluster)
  # Randomization, assignment to arms
  nclusters=length(unique(trial$cluster))
  #uniformly distributed numbers, take mean and boolean of that
  rand_numbers <- runif(nclusters,0,1)
  if(!is.null(trial[[baselineNumerator]]) & matchedPair){
    trial$base_num = trial[[baselineNumerator]]
    trial$base_denom = trial[[baselineDenominator]]
    cdf = data.frame(trial %>%
              group_by(cluster) %>%
              dplyr::summarize(positives = sum(base_num),
                               total = sum(base_denom)))
    cdf$p = cdf$positives/cdf$total
    cdf = cdf[order(cdf$p),]
    cdf$pair = rep(seq(1,nclusters/2),2)
    cdf$rand_numbers = rand_numbers
    cdf = cdf[with(cdf,order(pair,rand_numbers)),]
    cdf$arm = rep(c(1,0),nclusters/2)
    arm = cdf$arm[order(cdf$cluster)]
  } else {
    arm <- ifelse(rand_numbers > median(rand_numbers),1,0)
  }
  trial$arm <- factor(arm[trial$cluster[]], levels= c(0,1),labels = c("control", "intervention"))
return(trial)}


#' Specification of buffer zone in a cluster randomized trial
#'
#' \code{Specify_CRTbuffer} specifies a buffer zone in a cluster randomized trial (CRT) with geographical clustering by adding a flag to mark locations that are within a specified distance of locations in the opposite arm.
#' Locations are added incrementally to the buffer until a predetermined buffering distance is achieved.
#'
#' @param trial dataframe containing location coordinates and arm assignments
#' @param bufferWidth minimum distance between locations in the core areas of opposing arms
#' @return The input dataframe augmented with distances to the nearest discordant location and a logical indicator of whether the location is in the buffer zone
#' @export
#'
#' @examples
#' #Specify a buffer of 200m
#' exampletrial=Specify_CRTbuffer(testArms, bufferWidth=0.2)
#'
Specify_CRTbuffer= function(trial=trial,bufferWidth=0){
  # nearestDiscord: nearest coordinate in the discordant arm, for the control
  # coordinates return the minimal distance with a minus sign
  dist_trial <- as.matrix(dist(trial[,c('x','y')], method = "euclidean"))
  discord <- outer(trial$arm, trial$arm, "!=") #true & false.
  discord_dist_trial <- ifelse(discord,dist_trial,Inf)
  trial$nearestDiscord <- ifelse(trial$arm == 'control', -apply(discord_dist_trial, MARGIN=2, min), apply(discord_dist_trial, MARGIN=2, min))
  # Find the row and column numbers of the nearest discordant pair of neighbours
  discord_distance=ifelse(discord,0.2*round(abs(5*dist_trial),digits=1),Inf)
  # round off distances to nearest 20m to speed up calculation
  discord_distance =
    if(min(abs(discord_distance)) < bufferWidth){
      trial$buffer = FALSE
      pbuff=nbuff=width = 0
      while (width < bufferWidth){
        # round off distances to nearest 50m to speed up calculation
        width = min(abs(discord_distance))
        neighbours = which(abs(discord_distance) == min(abs(discord_distance)))
        trial$buffer[neighbours %% nrow(trial)] = TRUE

        # Set distances for buffered locations to inf to exclude them from further iterations
        discord_distance =
          ifelse(abs(discord_distance) == width,Inf,discord_distance)
        nbuff=sum(ifelse(trial$buffer,1,0))
        cat('\rAdded',nbuff-pbuff,'locations to buffer (',width,'km)  ')
        pbuff=nbuff
      }
      cat('\rLocations in buffer:',nbuff,'of',nrow(trial),'(',round(nbuff*100/nrow(trial),digits=1),'%)\n\n')
    }
  return(trial)}

#' Assign locations to clusters in a cluster randomized trial
#'
#' \code{DefineClusters} algorithmically assigns locations to clusters in a CRT with geographical clustering, by one of three algorithms:
#' (i)    Nearest neighbour (NN): (assigns equal numbers of locations to each cluster)
#' (ii)   kmeans clustering (kmeans) : (aims to partition locations so that each belongs to the cluster with the nearest centroid)
#' (iii)  travelling salesman problem heuristic : assigns locations sequentially along a travelling salesman path
#'
#' @param trial dataframe containing (x,y) coordinates of households
#' @param h  proposal for the number of locations in each cluster
#' @param nclusters  number of clusters
#' @param algo algorithm for cluster boundaries, options are: "TSP": travelling salesman problem heuristic; "NN": nearest neighbor; "kmeans": kmeans algorithm
#' @param reuseTSP indicator of whether a pre-existing path should be used by the TSP algorithm
#' @return The input coordinates as a dataframe augmented with the following attributes
#' \itemize{
#' \item \code{x}: x coordinates
#' \item \code{y}: y coordinates
#' \item \code{cluster}: cluster assignments
#' \item \code{h}: number of households in each cluster
#' \item \code{algo}: algorithm used for cluster boundaries
#' \item \code{path}: TSP path (if algo= 'TSP')
#' }
#' @export
#'
#' @examples
#' #Assign clusters to the test trial dataset averaging h=40 using the kmeans algorithm
#' exampletrial = DefineClusters(trial=test_site, h=40, algo='kmeans')

DefineClusters = function(trial=trial, h=NULL, nclusters=NULL, algo='NN', reuseTSP=FALSE){

  TSP_ClusterDefinition <- function(coordinates,h,nclusters,reuseTSP){

    if(!"path" %in% colnames(coordinates) | !reuseTSP){
      # Code originally from Silkey and Smith, SolarMal

      # Order the coordinates along an optimised travelling salesman path
      dist_coordinates <- dist(coordinates, method = "euclidean")
      tsp_coordinates <- TSP::TSP(dist_coordinates) # object of class TSP
      tsp_coordinates <- TSP::insert_dummy(tsp_coordinates)
      tour <- TSP::solve_TSP(tsp_coordinates,"repetitive_nn") #solves TSP, expensive
      path <- TSP::cut_tour(x=tour,cut='dummy')
      coordinates$path <- path

    }
    #order coordinates
    coordinates$order = seq(1:nrow(coordinates))
    coordinates <- coordinates[order(coordinates$path),]

    n1 = (nclusters-1)*h
    nclusters_1 = nclusters - 1
    # The last cluster may be a different size (if h is not a factor of the population size) )
    coordinates$cluster = NA
    coordinates$cluster[1:n1] <- c(rep(1:nclusters_1, each=h)) #add cluster assignment
    coordinates$cluster[which(is.na(coordinates$cluster))] = nclusters
    coordinates <- coordinates[order(coordinates$order),]
    return(coordinates)
  }

  NN_ClusterDefinition <- function(coordinates,h,nclusters){

    # algorithm is inspired by this website: ??? (comment from Lea)

    # initialize cluster, calculate euclidean distance
    dist_coordinates <- as.matrix(dist(coordinates, method = "euclidean"))
    coordinates$cluster <- NA

    nclusters_1 = nclusters -1
    for(i in 1:nclusters_1){

      #find unassigned coordinates
      cluster_unassigned <- which(is.na(coordinates$cluster))
      dist_coordinates_unassigned <- dist_coordinates[cluster_unassigned,cluster_unassigned]
      cluster_na <- rep(NA, length(cluster_unassigned))

      #find the coordinate furthest away from all the others
      index <- which.max(rowSums(dist_coordinates_unassigned))

      #find the n nearest neighbors of index
      cluster_na[head(order(dist_coordinates_unassigned[index,]), h)] <- i
      coordinates$cluster[cluster_unassigned] = cluster_na
    }
    # The last cluster may be a different size (if h is not a factor of the population size) )
    coordinates$cluster[which(is.na(coordinates$cluster))] = nclusters

    return(coordinates)
  }

  kmeans_ClusterDefinition <- function(coordinates,nclusters){

    #kmeans as implemented in R base
    km <- kmeans(x=coordinates,centers=nclusters)
    coordinates$cluster <- km$cluster

    return(coordinates)
  }

  # Local data from study area (ground survey and/or satellite images)
  coordinates = data.frame(x=as.numeric(as.character(trial$x)),
                           y=as.numeric(as.character(trial$y)))

  # the number of clusters and the target cluster size must be integers.
  # cluster size may vary slightly
  if(is.null(nclusters)){ nclusters = ceiling(nrow(coordinates)/h)}
  if(is.null(h)){ h = ceiling(nrow(coordinates)/nclusters)}
  ##############################################################################
  # derive cluster boundaries
  ##############################################################################

  if(algo=='TSP'){
    TSPoutput <- TSP_ClusterDefinition(coordinates,h,nclusters,reuseTSP)
    trial$path = TSPoutput$path
    trial$cluster = TSPoutput$cluster
  } else if(algo == 'NN'){
    trial$cluster <- NN_ClusterDefinition(coordinates,h,nclusters)$cluster
  } else if(algo == 'kmeans'){
    trial$cluster <- kmeans_ClusterDefinition(coordinates,nclusters)$cluster
  } else {
    stop('unknown method')
  }

  return(trial)}


#' Convert lat long co-ordinates to x,y
#'
#' \code{Convert_LatLong} converts co-ordinates expressed as decimal degrees into x,y using the equirectangular projection (valid for small areas).
#' Coordinates centred on the origin are returned.
#'
#' @param df data frame containing latitudes and longitudes in decimal degrees
#' @param latvar name of column containing latitudes in decimal degrees
#' @param longvar name of column containing longitudes in decimal degrees
#' @return The input dataframe with the lat-long coordinates replaced with Cartesian coordinates in units of km, centred on (0,0)
#' \itemize{
#' \item \code{x}: x co-ordinates
#' \item \code{y}: y co-ordinates
#' }
#' @export
Convert_LatLong = function(df,latvar='lat',longvar='long'){
  colnames(df)[colnames(df) == latvar] <- "lat"
  colnames(df)[colnames(df) == longvar] <- "long"
  R = 6371 # radius of the earth
  latradians = with(df,pi/180*lat)
  longradians = with(df,pi/180*long)
  meanlat = mean(latradians)
  meanlong = mean(longradians)
  drops <- c("lat","long")
  df=df[ , !(names(df) %in% drops)]
  df$y = R * (latradians - meanlat)* cos(longradians)
  df$x = R * (longradians - meanlong)
  return(df)
}


#' Anonymise locations in a trial site
#'
#' \code{Anonymise_TrialSite} carries out rotation of x,y coordinates a random angle about a random origin. Coordinates centred on the origin are returned.
#'
#' @param trial dataframe with Cartesian co-ordinates of households (columns x and y)
#' @return data.frame with modified co-ordinates of households (other columns remain unchanged)
#' @export
#'
#' @examples
#' #Rotate and reflect test site locations
#' transformedTestlocations=Anonymise_TrialSite()

Anonymise_TrialSite = function(trial=CRTspillover::test_site){

  # random rotation angle
  theta= 2 * pi * runif(n=1)
  x= trial$x
  y= trial$y
  rangex =max(x) - min(x)
  rangey =max(y) - min(y)
  translation = c(rangex * rnorm(n=1), rangey*rnorm(n=1))

  xy=t(matrix(c(x,y),ncol=2,nrow=length(x)))
  xytranslated= xy+translation

  rotation=matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),nrow=2,ncol=2)

  #Rotate
  xytrans = rotation %*% xytranslated

  # Recentre on origin
  recentred = xytrans - c(mean(xytrans[1,]),mean(xytrans[2,]))
  trial$x = recentred[1,]
  trial$y = recentred[2,]

  return(trial)}

# Local data from study area (ground survey and/or satellite images)


