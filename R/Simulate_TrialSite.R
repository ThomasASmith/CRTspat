#' Create a set of Cartesian co-ordinates for a simulated CRT
#'
#' \code{Simulate_TrialSite} Create a set of Cartesian co-ordinates using the Thomas algorithm to simulate a human settlement pattern
#'
#' @param scale standard deviation of random displacement from each settlement cluster center
#' @param households number of households in population
#' @param kappa intensity of Poisson process of settlement cluster centers
#' @param mu mean  number of points per settlement cluster
#' @return data frame with co-ordinates of households
#' @export
#'
#' @examples
#' #Generate a simulated area with 10,000 households
#' example_area = Simulate_TrialSite(scale = 2, households=10000, kappa=3, mu=40)
Simulate_TrialSite <- function(
  scale = 0.25,
  households = 2500,
  kappa = 4,
  mu = 50
){
  scaling = scale*20
  #Poisson point pattern with Thomas algorithm
  p <- spatstat.core::rThomas(kappa,scale,mu,
        win = spatstat.geom::owin(c(0,scaling),c(0,scaling)))
  #expected number of points: kappa*mu*scaling^2

  # create households and specify co-ordinates
  hhID <- c(1:households)
  x <- p$x[seq(1:households)]
  y <- p$y[seq(1:households)]
  coordinates <- data.frame(x=x,y=y)
return(coordinates)}

