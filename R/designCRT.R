#' Power and sample size calculations for a CRT
#'
#' \code{CRTpower} carries out power and sample size calculations for CRTs.
#'
#' @param trial dataframe or \code{'CRTsp'} object: optional list of locations
#' @param locations numeric: total number of units available for randomization (required if \code{trial} is not specified)
#' @param alpha numeric: confidence level
#' @param desiredPower numeric: desired power
#' @param effect numeric: required effect size
#' @param yC numeric: baseline (control) value of outcome
#' @param outcome_type character: with options -
#' \code{'y'}: continuous;
#' \code{'n'}: count;
#' \code{'e'}: event rate;
#' \code{'p'}: proportion;
#' \code{'d'}: dichotomous.
#' @param sigma2 numeric: variance of the outcome (required for \code{outcome_type = 'y'})
#' @param denominator numeric: rate multiplier (for \code{outcome_type = 'n'} or \code{outcome_type = 'e'})
#' @param N numeric: mean of the denominator for proportions (for \code{outcome_type = 'p'})
#' @param ICC numeric: Intra-cluster correlation
#' @param cv_percent numeric: Coefficient of variation of the outcome (expressed as a percentage)
#' @param c integer: number of clusters in each arm (required if \code{trial} is not specified)
#' @param sd_h standard deviation of number of units per cluster (required if \code{trial} is not specified)
#' @returns A list of class \code{'CRTsp'} object comprising the input data, cluster and arm assignments,
#' trial description and results of power calculations
#' @export
#' @details
#' Power and sample size calculations are for an unmatched two-arm trial. For counts
#' or event rate data the formula of [Hayes & Bennett, 1999](https://academic.oup.com/ije/article/28/2/319/655247) is used. This requires as an input the
#' between cluster coefficient of variation (\code{cv_percent}). For continuous outcomes and proportions the formulae of
#' [Hemming et al, 2011](https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-11-102) are used. These make use of
#' the intra-cluster correlation in the outcome (\code{ICC}) as an input. If the coefficient of variation and not the ICC is supplied then
#' the intra-cluster correlation is computed from the coefficient of variation using the formulae
#' from [Hayes & Moulton](https://www.taylorfrancis.com/books/mono/10.1201/9781584888178/cluster-randomised-trials-richard-hayes-lawrence-moulton). If incompatible values for \code{ICC} and \code{cv_percent} are supplied
#' then the value of the \code{ICC} is used.\cr\cr
#' The calculations do not consider any loss in power due to contamination, loss to follow-up etc..\cr\cr
#' If geolocations are not input then power and sample size calculations are based on the scalar input parameters.\cr\cr
#' If a trial dataframe or \code{'CRTsp'} object is input then this is used to determine the number of locations. If this input object
#' contains cluster assignments then the numbers and sizes of clusters in the input data are used to estimate the power. If buffer zones have been specified
#' then separate calculations are made for the core area and for the full site.\cr\cr
#' The output is an object of class \code{'CRTsp'} containing any input trial dataframe and values for:
#' - The required numbers of clusters to achieve the specified power.
#' - The design effect based on the input ICC.
#' - Calculations of the nominal power (ignoring any bias caused by contamination, loss to follow-up etc.)\cr
#' @examples
#' {# Power calculations for a binary outcome without input geolocations
#' examplePower1 = CRTpower(locations = 3000, ICC = 0.10, effect = 0.4, alpha = 0.05,
#'     outcome_type = 'd', desiredPower = 0.8, yC=0.35, c = 20, sd_h = 5)
#' summary(examplePower1)
#' # Power calculations for a rate outcome without input geolocations
#' examplePower2 = CRTpower(locations = 2000, cv_percent = 40, effect = 0.4, denominator = 2.5,
#'     alpha = 0.05, outcome_type = 'e', desiredPower = 0.8, yC = 0.35, c = 20, sd_h=5)
#' summary(examplePower2)
#' # Example with input geolocations and randomisation
#' examplePower3 = CRTpower(trial = readdata('example_site.csv'), desiredPower = 0.8,
#'     effect=0.4, yC=0.35, outcome_type = 'd', ICC = 0.05, c = 20)
#' summary(examplePower3)
#' }
CRTpower <- function(trial = NULL, locations = NULL, alpha = 0.05, desiredPower = 0.8,
    effect = NULL, yC = NULL, outcome_type = "d", sigma2 = NULL, denominator = 1,
    N = 1, ICC = NULL, cv_percent = NULL, c = NULL, sd_h = 0) {

    if(is.null(trial)) trial <- data.frame(x = c(), y = c())
    CRT <- CRTsp(trial)

    # populate a design list with a data about the input trial (if
    # available) and the input parameters

    # sigma2 is only required for continuous data, so its absence
    # should otherwise not crash the program
    if (!identical(outcome_type, "y")) sigma2 <- NA

    if (is.null(ICC)) ICC <- NA
    if (is.null(cv_percent)) cv_percent <- NA
    if (is.na(ICC) & is.na(cv_percent)) {
        stop("*** Value must be supplied for either ICC or cv_percent ***")
    }

    design <- ifelse(is.null(CRT$design$locations), list(), CRT$design)
    design$locations <- ifelse((nrow(CRT$trial) == 0), locations, nrow(CRT$trial))
    parnames <- c("alpha", "desiredPower", "effect", "yC", "outcome_type",
        "sigma2", "denominator", "N", "ICC", "cv_percent", "c", "sd_h")

    # Identify which variables to retrieve from the pre-existing design
    from_old <- lapply(mget(parnames), FUN = is.null)
    design[parnames] <- ifelse(from_old, design[parnames], mget(parnames))
    missing <- lapply(design[parnames], FUN = is.null)
    if (TRUE %in% missing) {
        stop("*** Value(s) must be supplied for: ", parnames[missing == TRUE],
            " ***")
    }
    CRT <- CRTsp(CRT, design = design)
}


# Characteristics of a trial design. The input is a data frame or
# CRTsp object. The output list conforms to the requirements for a
# CRTsp object
get_geom <- function(trial = NULL, design = NULL) {

    sd_distance <- clustersRequired <- DE <- power <- NULL
#   check if the power calculations need to be reconstructed from scratch
    if(!is.null(design$locations)) {
        locations <- design$locations
        c <- design$c
        mean_h <- locations/(2 * c)
        sd_h <- design$sd_h
    } else {
        outcome_type <- mean_h <- c <- sd_h <- locations <- NULL
    }
    geom <- list(locations = locations,
                 sd_h = sd_h,
                 c = c,
                 records = 0,
                 mean_h = mean_h,
                 DE = NULL,
                 power = NULL,
                 clustersRequired = NULL)

    # cluster size
    geom$c <- ifelse(is.null(geom$c), NA, round(geom$c))
    geom$mean_h <- geom$locations/(2 * geom$c)

    # overwrite values from the design with those from the data frame if these are present
    if (!is.null(trial) & nrow(trial) > 0) {
        coordinates <- data.frame(cbind(x = trial$x, y = trial$y))
        geom$records <- nrow(trial)
        geom$locations <- nrow(dplyr::distinct(coordinates))

        if (!is.null(trial$cluster)) {
            # reassign the cluster levels in case some are not
            # represented in this geom (otherwise nlevels() counts
            # clusters that are not present)
            trial$cluster <- as.factor(as.character(trial$cluster))

            geom$c <- floor(nlevels(trial$cluster)/2)
            # mean number of locations randomized in each cluster
            geom$mean_h <- mean(table(trial$cluster))

            # standard deviation of locations assigned to each cluster
            geom$sd_h <- stats::sd(table(trial$cluster))

        }
        if (!is.null(trial$arm)) {
            geom$sd_distance <- stats::sd(trial$nearestDiscord)
            arms <- unique(cbind(trial$cluster, trial$arm))[, 2]  #assignments
        }
    }

    if (!is.null(design$effect)) {
        if (is.null(geom$locations)) {
            stop("*** Number of locations is a required input ***")
        }
        if (is.null(geom$c)) {
            stop("*** Number of clusters is a required input ***")
        }
        if (identical(geom$sd_h, 0)) {
            message("*** Assuming all clusters are the same size ***")
        }
        effect <- design$effect
        yC <- design$yC

        # convert power and significance level to Zvalues
        Zsig <- -qnorm(design$alpha/2)
        Zpow <- qnorm(design$desiredPower)


        # effective cluster sizes (inflating for multiple observations
        # at the same location)
        mean_eff <- geom$mean_h * design$N
        sd_eff <- geom$sd_h * design$N

        # coefficient of variation of the cluster sizes
        cv_eff <- sd_eff/mean_eff

        # outcome in intervened group
        link <- switch(design$outcome_type, y = "identity", n = "log", e = "log",
                       p = "logit", d = "logit")
        yI <- ifelse(link == "identity", yC - design$effect, yC * (1 - design$effect))

        # difference between groups
        # d <- yC - yI

        # input value of the coefficient of variation of between cluster variation in outcome
        k <- ifelse(is.null(design$cv_percent), NA, design$cv_percent/100)

        if(is.null(design$ICC)) {
            design$ICC <- switch(link, "identity" = (k * yC)^2/sigma2,
                                "log" = NA,
                                "logit" = k^2 * yC/(1 - yC)
                          )
        }

        if(is.null(design$cv_percent)) {
            k <- switch(link, "identity" = sqrt(design$ICC * sigma2)/yC,
                                 "log" = NA,
                               "logit" = sqrt(design$ICC * (1 - yC)/yC)
            )
            design$cv_percent <- 100 * k
        }

        if(identical(link,"log")) {
            if(is.null(k)){
                stop("*** Between cluster coefficient of variation is a required input ***")
            }
            # use the formulae from Hayes & Bennett (1999) https://doi.org/10.1093/ije/28.2.319

            denom_per_cluster <- design$denominator * mean_eff
            # clusters required (both arms)
            geom$clustersRequired <- 2 * ceiling(1 + (Zsig + Zpow)^2 * ((yC + yI)/denom_per_cluster + (yC^2 + yI^2) * k^2)/((yC - yI)^2))

            # power with c clusters per arm and unequal cluster sizes
            geom$power <- stats::pnorm(sqrt((c - 1) * ((yC - yI)^2)/((yC + yI)/denom_per_cluster + (yC^2 + yI^2) * k^2)) - Zsig)

            # the design effect is the ratio of the required denominator to that required for an individually randomised trial
            required_denom_RCT <- 2 * ((Zsig + Zpow)^2 * (yC + yI)/((yC - yI)^2))
            geom$DE <- denom_per_cluster * geom$clustersRequired/required_denom_RCT

        } else {
            # design effect (Variance Inflation Factor) as a function of ICC allowing for
            # varying cluster sizes(Hemming eqn 6)
            geom$DE <- 1 + (cv_eff^2 + 1) * (mean_eff - 1) * design$ICC
            if (identical(link, "identity")) {
                # with normal models, sigma2 is an input variable
                sigma2 <- design$sigma2
            } else if (identical(link, "logit")) {
                # This is the variance for a Bernoulli. The cluster sizes are
                # inflated for the binomial case (below)
                sigma2 <- 1/2 * (yI * (1 - yI) + yC * (1 - yC))
            }
            # required individuals per arm in individually randomized trial
            n_ind <- 2 * sigma2 * ((Zsig + Zpow)/(yC - yI))^2

            # number of individuals required per arm in CRT with equal
            # cluster sizes
            n_crt <- n_ind * geom$DE

            # minimum total numbers of clusters required assuming varying
            # cluster sizes per arm (Hemming eqn 8)
            geom$clustersRequired <- 2 * ceiling(n_ind * (1 + ((cv_eff + 1) *
                                                                   mean_eff - 1) * design$ICC)/mean_eff)

            # power with c clusters per arm and unequal cluster sizes
            geom$power <- stats::pnorm(sqrt(c * mean_eff/(2 * geom$DE)) * (yC - yI)/sqrt(sigma2) - Zsig)
        }
    }
    return(geom)
}

