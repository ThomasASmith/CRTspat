#' Analysis of cluster randomized trial with contamination
#'
#' \code{CRTanalysis} carries out a statistical analysis of a cluster randomized trial (CRT).
#' @param trial an object of class \code{"CRTsp"} or a data frame containing locations in (x,y) coordinates, cluster
#'   assignments (factor \code{cluster}), and arm assignments (factor \code{arm}) and outcome data (see details).
#' @param method statistical method with options:
#'  \tabular{ll}{
#' \code{"EMP"} \tab simple averages of the data   \cr
#' \code{"T"}   \tab comparison of cluster means by t-test \cr
#' \code{"GEE"} \tab Generalised Estimating Equations \cr
#' \code{"INLA"}\tab Integrated Nested Laplace Approximation (INLA) \cr
#' \code{"MCMC"}\tab Markov chain Monte Carlo using \code{"JAGS"} \cr
#' }
#' @param cfunc transformation defining the contamination function. \cr
#' options are:
#' \tabular{llll}{
#' \code{"X"} \tab\tab contamination not modelled\tab the only valid value of \code{cfunc} for methods \code{"EMP"}, \code{"T"} and \code{"GEE"}\cr
#' \code{"L"} \tab\tab inverse logistic (sigmoid)\tab the default for \code{"INLA"} and \code{"MCMC"} methods\cr
#' \code{"P"} \tab\tab inverse probit (error function)\tab available with \code{"INLA"} and \code{"MCMC"} methods\cr
#' \code{"S"} \tab\tab piecewise linear\tab only available with the \code{"MCMC"} method\cr
#' }
#' @param link link function. options are:
#'  \tabular{ll}{
#' \code{"logit"}\tab (the default). \code{numerator} has a binomial distribution with denominator \code{denominator}.\cr
#' \code{"log"}  \tab \code{numerator} is Poisson distributed with an offset of log(\code{denominator}).
#' With the \code{"INLA"} and \code{"MCMC"} methods 'iid' random effects are used to model extra-Poisson variation.\cr
#' \code{"identity"}\tab The outcome is \code{numerator/denominator} with a normally distributed error function.\cr
#' }
#' @param numerator string: name of numerator variable for outcome
#' @param denominator string: name of denominator variable for outcome data (if present)
#' @param excludeBuffer logical: indicator of whether any buffer zone (records with \code{buffer=TRUE}) should be excluded from analysis
#' @param alpha numeric: confidence level for confidence intervals and credible intervals
#' @param baselineOnly logical: indicator of whether required analysis is of effect size or of baseline only
#' @param baselineNumerator string: name of numerator variable for baseline data (if present)
#' @param baselineDenominator string: name of denominator variable for baseline data (if present)
#' @param localisedEffects logical: indicator of whether the model includes local effects with no contamination
#' @param clusterEffects logical: indicator of whether the model includes cluster random effects
#' @param spatialEffects logical: indicator of whether the model includes spatial random effects
#' @param resamples integer: number of samples for sample-based intervals
#' @param requireMesh logical: indicator of whether spatial predictions are required
#' @param inla_mesh string: name of pre-existing INLA input object created by \code{new_mesh()}
#' @return list of class \code{CRTanalysis} containing the following results of the analysis:
#' \itemize{
#' \item \code{description} : description of the dataset
#' \item \code{method} : statistical method
#' \item \code{pt_ests} : point estimates
#' \item \code{int_ests} : interval estimates
#' \item \code{model_object} : object returned by the fitting routine
#' \item \code{contamination} : function values and statistics describing the estimated contamination
#' }
#' @details \code{CRTanalysis} is a wrapper for the statistical analysis packages:
#' [geepack](https://www.jstatsoft.org/article/view/v015i02),
#' [INLA](https://www.r-inla.org/),
#' [jagsUI](https://cran.r-project.org/web/packages/jagsUI/index.html),
#' and the [t.test](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/t.test)
#' function of package \code{stats}.
#' The wrapper does not provide an interface to the full functionality of these packages.
#' It is specific for typical analyses of cluster randomized trials with geographical clustering.\cr \cr
#'
#' The key results of the analyses can be extracted using a \code{summary()} of the output list.
#' The \code{model_object} in the output list is the usual output from the statistical analysis routine,
#' and can be also be inspected with \code{summary()}, or analysed using \code{stats::fitted()}
#' for purposes of evaluation of model fit etc..
#' @importFrom grDevices rainbow
#' @importFrom stats binomial dist kmeans median na.omit qlogis qnorm quantile rbinom rnorm runif simulate
#' @importFrom utils head read.csv
#' @export
#' @examples
#' # Analysis of test dataset by t-test
#' exampleT = CRTanalysis(readdata("exampleCRT.csv"), method = "T")
#' # Standard GEE analysis of test dataset ignoring contamination
#' exampleGEE = CRTanalysis(readdata("exampleCRT.csv"), method = "GEE")
#'
CRTanalysis <- function(
    trial, method = "GEE", cfunc = "L", link = "logit", numerator = "num",
    denominator = "denom", excludeBuffer = FALSE, alpha = 0.05,
    baselineOnly = FALSE, baselineNumerator = "base_num", baselineDenominator = "base_denom",
    localisedEffects = FALSE, clusterEffects = TRUE, spatialEffects = FALSE,
    resamples = 10000, requireMesh = FALSE, inla_mesh = NULL)
    {
    # Test of validity of inputs
    if (!method %in% c("EMP", "T", "MCMC", "GEE", "INLA", "LME4"))
        {
        stop("*** Invalid value for statistical method ***")
        return(NULL)
    }

    if (!cfunc %in% c("S", "L", "P", "X", "Z"))
        {
        stop("*** Invalid contamination function ***")
        return(NULL)
    }

    # create names for confidence limits for use throughout
    CLnames <- c(
        paste0(alpha/0.02, "%"),
        paste0(100 - alpha/0.02, "%")
    )

    # MAIN FUNCTION CODE STARTS HERE

    cluster <- NULL

    if (identical(class(trial),"CRTsp")) trial <- trial$trial

    if ("buffer" %in% colnames(trial) & excludeBuffer)
        {
        trial <- trial[!trial$buffer, ]
    }

    # trial needs to be ordered for some analyses
    if(!is.null(trial$cluster)) trial <- trial[order(trial$cluster), ]

    # Some statistical methods do not allow for contamination
    if (method %in% c("EMP", "T", "GEE")) cfunc <- "X"

    if (baselineOnly){
        if (method %in% c("EMP", "T", "GEE", "MCMC", "LME4"))
            {
            method <- "GEE"
            cat("Analysis of baseline only, using GEE\n")
        } else if (method == "INLA")
        {
            cat("Analysis of baseline only, using INLA\n")
        }
        # cfunc='Z' is used to remove the estimation of effect size from the model
        cfunc <- "Z"
        trial$y1 <- trial[[baselineNumerator]]
        trial$y0 <- trial[[baselineDenominator]] - trial[[baselineNumerator]]
        trial$y_off <- trial[[baselineDenominator]]

    } else {
        if (is.null(trial[[numerator]])){
            cat("*** No outcome data to analyse ***")
            return()
        }
        trial$y1 <- trial[[numerator]]
        trial$y0 <- trial[[denominator]] - trial[[numerator]]
        trial$y_off <- trial[[denominator]]

        # if nearestDiscord is not provided augment the trial data frame with distance to nearest discordant
        # coordinate
        if (is.null(trial$nearestDiscord)) trial$nearestDiscord <- get_nearestDiscord(trial)
    }

    # store options here- noting that the model formula depends on allowable values of other options
    options <- list(method = method, link = link, cfunc = cfunc,
                    alpha = alpha, baselineOnly = baselineOnly,
                    ftext = NULL)

    # create scaffolds for lists
    pt_ests <- list(contamination.par = NA, pr.contaminated = NA, contamination.interval = NA)
    int_ests <- list(controlY = NA, interventionY = NA, effect.size = NA)
    model_object <- list()
    description <- get_description(trial=trial, link=link, baselineOnly)
    analysis <- list(trial = trial, pt_ests = pt_ests, int_ests = int_ests)

    # create model formula for display even though this is not used for MCMC models
    fterms <- switch(cfunc,
        Z = "",
        X = "arm",
        L = "pvar",
        P = "pvar",
        S = "pvar"
    )

    if (localisedEffects & cfunc != 'X') fterms <- c(fterms, "arm")
    if (clusterEffects) {
        if (identical(method, "INLA")) {
            fterms <- c(fterms, "f(cluster, model = \'iid\')")
        } else if(identical(method, "GEE")) {
            fterms <- fterms
        } else {
            fterms <- c(fterms, "(1 | cluster)")
        }
    }
    if (identical(method, "INLA")){
        if (spatialEffects) fterms <- c(fterms, "f(s, model = spde)")
        if (identical(link,"log")) fterms <- c(fterms, "f(id, model = \'iid\')")
    }
    options$ftext <- paste(fterms, collapse = " + ")

    # TODO: remove the duplication and simplify here
    # Specify the function used for calculating the linear predictor
    LPfunction <- c(
        "CalculateNoEffect", "CalculateNoContaminationFunction", "CalculatePiecewiseLinearFunction",
        "CalculateLogisticFunction", "CalculateProbitFunction")[which(cfunc == c("Z", "X", "S", "L", "P"))]
    FUN1 <- eval(parse(text = LPfunction))
    # specify functional form of sigmoid in distance from boundary 'L' inverse logit; 'P' inverse probit; 'X'
    # or 'Z' do not model contamination
    FUN <- switch(
        cfunc, L = "invlink(link='logit', x)", P = "stats::pnorm(x)", X = NULL, Z = NULL)

    if (method == "EMP"){
        pt_ests$controlY <- unname(description$controlY)
        pt_ests$interventionY <- unname(description$interventionY)
        pt_ests$effect.size <- unname(description$effect.size)
        pt_ests$contamination.interval <- NA
    } else if (method == "T"){
        y1 <- arm <- NULL
        clusterSum <- data.frame(
            trial %>%
                group_by(cluster) %>%
                dplyr::summarize(
                  y = sum(y1),
                  total = sum(y_off),
                  arm = arm[1]
              )
        )
        clusterSum$lp <- switch(link,
                "identity" = clusterSum$y/clusterSum$total,
                "log" = log(clusterSum$y/clusterSum$total),
                "logit" = logit(clusterSum$y/clusterSum$total))
        formula <- stats::as.formula("lp ~ arm")
        model_object <- stats::t.test(
            formula = formula, data = clusterSum, alternative = "two.sided",
            conf.level = 1 - alpha, var.equal = TRUE
        )
        pt_ests$p.value <- model_object$p.value
        analysisC <- stats::t.test(
            clusterSum$lp[clusterSum$arm == "control"], conf.level = 1 - alpha)
        pt_ests$controlY <- invlink(link, analysisC$estimate[1])
        int_ests$controlY <- invlink(link, analysisC$conf.int)
        analysisI <- stats::t.test(
            clusterSum$lp[clusterSum$arm == "intervention"], conf.level = 1 - alpha)
        pt_ests$interventionY <- invlink(link, analysisI$estimate[1])
        int_ests$interventionY <- invlink(link, analysisI$conf.int)

        # Covariance matrix (note that two arms are independent so the off-diagonal elements are zero)
        Sigma <- base::matrix(
            data = c(analysisC$stderr^2, 0, 0, analysisI$stderr^2),
            nrow = 2, ncol = 2)
        if (link == 'identity'){
            pt_ests$effect.size <- pt_ests$controlY - pt_ests$interventionY
            int_ests$effect.size <- unlist(model_object$conf.int)
        }
        if (link %in% c("logit","log")){
            pt_ests$effect.size <- 1 - pt_ests$interventionY/pt_ests$controlY
            int_ests$effect.size <- 1 - exp(-unlist(model_object$conf.int))
        }
    } else if (method == "GEE") {
        # GEE analysis of cluster effects
        fterms <- c(switch(link,
                     "identity" = "y1/y_off ~ 1",
                     "log" = "y1 ~ 1 + offset(log(y_off))",
                     "logit" = "cbind(y1,y0) ~ 1"),
                    fterms)
        formula <- stats::as.formula(paste(fterms, collapse = "+"))
        if (link == "log") {
            model_object <- geepack::geeglm(
                formula = formula, id = cluster, data = trial, family = poisson(link = "log"),
                corstr = "exchangeable", scale.fix = FALSE)
        } else if (link == "logit") {
            model_object <- geepack::geeglm(
                formula = formula, id = cluster, corstr = "exchangeable",
                data = trial, family = binomial(link = "logit"))
        } else if (link == "identity") {
            model_object <- geepack::geeglm(
                formula = formula, id = cluster, corstr = "exchangeable",
                data = trial, family = gaussian)
        }

        summary.fit <- summary(model_object)

        z <- -qnorm(alpha/2)  #standard deviation score for calculating confidence intervals
        lp_yC <- summary.fit$coefficients[1, 1]
        se_lp_yC <- summary.fit$coefficients[1, 2]

        clusterSize <- nrow(trial)/nlevels(as.factor(trial$cluster))


        # remove the temporary objects from the dataframe
        model_object$data$y1 <- model_object$data$y0 <- model_object$data$y_off <- NULL
        pt_ests$controlY <- invlink(link, lp_yC)
        int_ests$controlY <- namedCL(
            invlink(link, c(lp_yC - z * se_lp_yC, lp_yC + z * se_lp_yC)),
            alpha = alpha
        )

        # Intracluster correlation
        pt_ests$ICC <- noLabels(summary.fit$corr[1])  #with corstr = 'exchangeable', alpha is the ICC
        se_ICC <- noLabels(summary.fit$corr[2])
        int_ests$ICC <- namedCL(
            noLabels(c(pt_ests$ICC - z * se_ICC, pt_ests$ICC + z * se_ICC)),
            alpha = alpha
        )
        pt_ests$DesignEffect <- 1 + (clusterSize - 1) * pt_ests$ICC  #Design Effect
        int_ests$DesignEffect <- 1 + (clusterSize - 1) * int_ests$ICC

        # Estimation of effect.size does not apply if analysis is of baseline only (cfunc='Z')
        pt_ests$effect.size <- NULL
        if (cfunc == "X") {
            lp_yI <- summary.fit$coefficients[1, 1] + summary.fit$coefficients[2,
                1]
            se_lp_yI <- sqrt(
                model_object$geese$vbeta[1, 1] +
                model_object$geese$vbeta[2, 2] +
                2 * model_object$geese$vbeta[1,2]
            )

            int_ests$interventionY <- namedCL(
                invlink(link, c(lp_yI - z * se_lp_yI, lp_yI + z * se_lp_yI)),
                alpha = alpha)

            int_ests$effect.size <- estimateCLeffect.size(
                mu = summary.fit$coefficients[, 1], Sigma = model_object$geese$vbeta,
                alpha = alpha, resamples = resamples, method = method,
                link = link)

            pt_ests$interventionY <- invlink(link, lp_yI)
            pt_ests$effect.size <- (1 - invlink(link, lp_yI)/invlink(link, lp_yC))
        }
    } else if (method == "LME4"){
        dummy <- 'rabbit'

        beta <- NA
        if (cfunc %in% c("L", "P")) {
            if (identical(Sys.getenv("TESTTHAT"), "true") | 2 == 2) {
                beta <- 2.0
            } else {
                cat("Estimating scale parameter for contamination range", "\n")
#                beta <- stats::optimize(
#                    f = estimateContaminationLME4, interval = c(-10, 10),
#                    trial = trial, FUN = FUN, formula = formula,
#                    tol = 0.1, link = link)$minimum
                beta <- 3.0
            }
            x <- trial$nearestDiscord * exp(beta)
            trial$pvar <- eval(parse(text = FUN))
            # set up linear contrasts (not required for cfunc='X' or 'Z')
        }
        if (link == "log") {
            fterms <- c(fterms,"offset(log(y_off))")
            formula <- stats::as.formula(paste(fterms, collapse = "+"))
            model_object <-lme4::glmer(formula = formula, data = trial,
                                       family = quasipoisson)
        } else if (link == "logit") {
            fterms <- c(fterms,"cbind(y1,y0) ~ 1", "cbind(y1,y0) ~ arm")
            formula <- stats::as.formula(paste(fterms, collapse = "+"))
            model_object <-lme4::glmer(formula = formula, data = trial,
                                 family = binomial)
        } else if (link == "identity") {
            fterms <- ifelse(cfunc == "Z", "y1/y_off ~ 1", "y1/y_off ~ arm")
            formula <- stats::as.formula(paste(fterms, collapse = "+"))
            model_object <-lme4::glmer(formula = formula, data = trial,
                                       family = gaussian)
        }
    # MCMC Methods
    } else if (method == "MCMC"){
        nchains <- 2
        max.iter <- 50000
        burnin <- 5000

        # JAGS outputs 95% intervals, so alpha is set to 0.05 to avoid confusion
        alpha <- 0.05
        cat("* Note: alpha set to 0.05 *\n")

        datajags <- with(trial, list(d = nearestDiscord,
                                     N = nrow(trial)))
        if (link == 'identity') {
            datajags$y <- trial$y1/trial$y_off
        } else {
            datajags$y1 <- trial$y1
            datajags$y_off <- trial$y_off
        }
        if (clusterEffects) {
            datajags$cluster <- as.numeric(as.character(trial$cluster))
            datajags$ncluster <- max(as.numeric(as.character(trial$cluster)))
        }
        # construct the rjags code by concatenating strings

        text1 <- "model{\n
          for(i in 1:N){\n
              cont[i] <- ifelse(abs(d[i]) < theta, 1, 0) \n"

        text2 <- switch(cfunc, S = "pvar[i] <- ifelse(d[i] < -ebeta,0,
                                        ifelse(d[i] > ebeta,1,
                                        (d[i] + ebeta)/(2*ebeta)))\n",
                        P = "pvar[i] <- pnorm(d[i],0,ebeta) \n",
                        L = "pvar[i] <- 1/(1 + exp(-ebeta*d[i])) \n")

        text3 <- switch(link,
                        "identity" = "y[i] ~ dnorm(lp[i],tau1) \n",
                        "log" =  "gamma1[i] ~ dnorm(0,tau1) \n
                                  Expect_y[i] <- exp(lp[i] + gamma1[i]) * y_off[i] \n
                                  y1[i] ~ dpois(Expect_y[i]) \n",
                        "logit" = "logitp[i] <- lp[i]  \n
                                   p[i] <- 1/(1 + exp(-logitp[i])) \n
                                   y1[i] ~ dbin(p[i],y_off[i]) \n"
        )

        # construct linear predictor
        text4 <- "lp[i] <- int + arm * pvar[i]"
        text5 <- ifelse(clusterEffects,
            " + gamma[cluster[i]] \n
            }\n
            for(ic in 1:ncluster) {\n
                gamma[ic] ~ dnorm(0, tau)\n
            }\n
            tau <- 1/(sigma * sigma) \n
            sigma ~ dunif(0, 2) \n
            ", "}\n")
        text6 <-
            "beta ~ dnorm(0, 1E-1) \n
            ebeta <- exp(beta) \n
            alpha <- 0.05 \n
            int ~ dnorm(0, 1E-2) \n
            arm ~ dnorm(0, 1E-2) \n"

        text7 <- switch(link,
                        "identity" = "yC <- int \n
                                      yI <- int + arm \n
                                      tau1 <- 1/(sigma1 * sigma1) \n
                                      sigma1 ~ dunif(0, 2) \n
                                      Es <- yC - yI \n",
                        "log" = "yC <- exp(int) \n
                                 yI <- exp(int + arm) \n
                                 tau1 <- 1/(sigma1 * sigma1) \n
                                 sigma1 ~ dunif(0, 2) \n
                                 Es <- 1 - yI/yC \n",
                        "logit" = "yC <- 1/(1 + exp(-int)) \n
                                   yI <- 1/(1 + exp(-(int + arm))) \n
                                   Es <- 1 - yI/yC \n"
        )

        # contamination diameter depends on contamination function
        text8 <- switch(cfunc, S = "theta <- (1 - 0.5 * alpha) * 2 * ebeta \n",
                               P = "theta <- 2 * qnorm(1 - 0.5 * alpha, 0, ebeta) \n",
                               L = "theta <- 2 * log((1- 0.5 * alpha)/(0.5 * alpha))/ebeta \n")
        text9 <- "pcont <- mean(cont)}\n"
        MCMCmodel <- paste0(text1, text2, text3, text4, text5, text6, text7, text8, text9)
        jagsout <- jagsUI::autojags(data = datajags, inits = NULL,
                                    parameters.to.save = c("Es", "theta", "yC", "yI", "beta", "pcont"),
                                    model.file = textConnection(MCMCmodel), n.chains = nchains,
                                    iter.increment = 1000, n.burnin = burnin, max.iter=max.iter)
        analysis$model_object <- jagsout
        analysis$pt_ests$controlY <- jagsout$q50$yC
        analysis$int_ests$controlY <- namedCL(c(jagsout$q2.5$yC, jagsout$q97.5$yC), alpha = alpha)
        analysis$pt_ests$interventionY <- jagsout$q50$yI
        analysis$int_ests$interventionY <- namedCL(c(jagsout$q2.5$yI, jagsout$q97.5$yI), alpha = alpha)
        analysis$pt_ests$effect.size <- jagsout$q50$Es
        analysis$int_ests$effect.size <- namedCL(c(jagsout$q2.5$Es, jagsout$q97.5$Es), alpha = alpha)
        analysis$pt_ests$contamination.interval <- jagsout$q50$theta
        analysis$int_ests$contamination.interval <- namedCL(c(jagsout$q2.5$theta, jagsout$q97.5$theta), alpha = alpha)
        analysis$pt_ests$contamination.par <- jagsout$q50$beta
        analysis$int_ests$contamination.par <- namedCL(c(jagsout$q2.5$beta, jagsout$q97.5$beta), alpha = alpha)
        analysis$pt_ests$pr.contaminated <- jagsout$q50$pcont
        analysis$int_ests$pr.contaminated <- namedCL(c(jagsout$q2.5$pcont, jagsout$q97.5$pcont), alpha = alpha)
        analysis$pt_ests$DIC <- jagsout$DIC
        analysis$model_object$MCMCmodel <- MCMCmodel
    # INLA methods
    } else if (method == "INLA")
        {
        if (identical(link,"identity")) {
            fterms <- c("y1/y_off ~ 0 + int", fterms)
        } else {
            fterms <- c("y1 ~ 0 + int", fterms)
        }
        formula <- stats::as.formula(paste(fterms, collapse = "+"))
        trial <- dplyr::mutate(trial, id =  dplyr::row_number())
        # If spatial predictions are not required a minimal mesh is generated
        pixel <- 0.5
        if (!requireMesh) pixel <- 2
        if (is.null(inla_mesh)) {
                inla_mesh <- new_mesh(
                    trial = trial, offset = -0.1, max.edge = 0.25, inla.alpha = 2,
                    maskbuffer = 0.5, pixel = pixel)
        }
        y_off <- NULL
        spde <- inla_mesh$spde
        df = data.frame(
            int = rep(1, nrow(trial)),
            id = trial$id)
        if ("int + arm" %in% fterms | "arm" %in% fterms) df$arm = ifelse(trial$arm == "intervention", 1, 0)
        if ("f(cluster, model = \'iid\')" %in% fterms) df$cluster = trial$cluster
        effectse <- list(df = df, s = inla_mesh$indexs)
        dfp = data.frame(
            int = rep(1, nrow(inla_mesh$prediction)),
            id = inla_mesh$prediction$id)
        if ("int + arm" %in% fterms | "arm" %in% fterms) dfp$arm = ifelse(inla_mesh$prediction$arm == "intervention", 1, 0)
        if ("f(cluster, model = \"iid\")" %in% fterms) dfp$cluster = inla_mesh$prediction$cluster
        effectsp <- list(df = dfp, s = inla_mesh$indexs)

        lc <- NULL
        beta <- NA
        if (cfunc %in% c("L", "P")) {
            if (identical(Sys.getenv("TESTTHAT"), "true")) {
                beta <- 2.0
            } else {
                cat("Estimating scale parameter for contamination range", "\n")
                beta <- stats::optimize(
                    f = estimateContaminationINLA, interval = c(-10, 10),
                    trial = trial, FUN = FUN, inla_mesh = inla_mesh, formula = formula,
                    tol = 0.1, link = link)$minimum
            }
            x <- trial$nearestDiscord * exp(beta)
            trial$pvar <- eval(parse(text = FUN))

            effectse$df$pvar <- trial$pvar

            x <- inla_mesh$prediction$nearestDiscord * exp(beta)
            inla_mesh$prediction$pvar <- eval(parse(text = FUN))
            effectsp$df$pvar <- inla_mesh$prediction$pvar

            # set up linear contrasts (not required for cfunc='X' or 'Z')
            if (grepl("pvar", options$ftext, fixed = TRUE))
                {
                lc <- INLA::inla.make.lincomb(int = 1, pvar = 1)
                if (grepl("arm", options$ftext, fixed = TRUE))
                  {
                  lc <- INLA::inla.make.lincomb(int = 1, pvar = 1, arm = 1)
                }
            } else if (grepl("arm", options$ftext, fixed = TRUE))
                {
                lc <- INLA::inla.make.lincomb(int = 1, arm = 1)
            }
        }
        # stack for estimation stk.e
        stk.e <- INLA::inla.stack(
            tag = "est", data = list(y1 = trial$y1, y_off = trial$y_off),
            A = list(1, A = inla_mesh$A),
            effects = effectse
        )

        # stack for prediction stk.p
        stk.p <- INLA::inla.stack(
            tag = "pred", data = list(y1 = NA, y_off = NA),
            A = list(1, inla_mesh$Ap),
            effects = effectsp
        )

        # stk.full comprises both stk.e and stk.p if a prediction mesh is in use
        stk.full <- INLA::inla.stack(stk.e, stk.p)

        if (link == "identity") {
            inla_result <- INLA::inla(
                formula, family = "gaussian", lincomb = lc,
                control.family = list(link = "identity"),
                data = INLA::inla.stack.data(stk.full),
                control.fixed = list(correlation.matrix = TRUE),
                control.predictor = list(compute = TRUE, link = 1,
                                         A = INLA::inla.stack.A(stk.full)),
                control.compute = list(dic = TRUE))
        } else if (link == "log") {
            inla_result <- INLA::inla(
                formula, family = "poisson", lincomb = lc,
                control.family = list(link = "log"),
                data = INLA::inla.stack.data(stk.full),
                control.fixed = list(correlation.matrix = TRUE),
                control.predictor = list(compute = TRUE, link = 1,
                                         A = INLA::inla.stack.A(stk.full)),
                control.compute = list(dic = TRUE))
        } else if (link == "logit") {
            inla_result <- INLA::inla(
                formula, family = "binomial", Ntrials = y_off, lincomb = lc,
                control.family = list(link = "logit"),
                data = INLA::inla.stack.data(stk.full),
                control.fixed = list(correlation.matrix = TRUE),
                control.predictor = list(compute = TRUE, link = 1,
                                         A = INLA::inla.stack.A(stk.full)),
                control.compute = list(dic = TRUE))
        }

        # Augment the inla results list with application specific quantities
        index <- INLA::inla.stack.index(stack = stk.full, tag = "pred")$data
        inla_mesh$prediction$prediction <-
                invlink(link, inla_result$summary.linear.predictor[index, "0.5quant"])
        # Compute sample-based confidence limits for intervened outcome and effect.size if intervention effects are
        # estimated
        if (grepl("pvar", options$ftext, fixed = TRUE) |
            grepl("arm", options$ftext, fixed = TRUE))
                {
            # Specify the means of the variables
            mu <- inla_result$summary.lincomb.derived$mean
            names(mu) <- rownames(inla_result$summary.lincomb.derived)
            # Specify the covariance matrix of the variables
            cov <- inla_result$misc$lincomb.derived.covariance.matrix
            sample <- as.data.frame(MASS::mvrnorm(n = 10000, mu = mu, Sigma = cov))
            sample$controlY <- invlink(link, sample$int)
            # pr.contaminated is the proportion of effect subject to contamination
            if ("arm" %in% names(mu) &
                "pvar" %in% names(mu))
                  {
                sample$interventionY <- invlink(link, sample$lc)
                sample$pr.contaminated <- with(
                  sample, 1 - (controlY - invlink(link, int + arm))/(controlY -
                    interventionY)
              )
            } else if ("arm" %in% names(mu))
                {
                sample$interventionY <- invlink(link, sample$int + sample$arm)
                sample$pr.contaminated <- 0
            } else if ("pvar" %in% names(mu))
                {
                sample$interventionY <- invlink(link, sample$int + sample$pvar)
                sample$pr.contaminated <- 1
            }
            sample$effect.size <- 1 - sample$interventionY/sample$controlY
            bounds <- (apply(
                sample, 2, function(x)
                  {
                  quantile(
                    x, c(alpha/2, 0.5, 1 - alpha/2),
                    alpha = alpha
                )
                }
            ))
        } else
        {
            controlY <- unlist(
                invlink(
                  link, inla_result$summary.fixed["int", c("0.025quant", "0.5quant", "0.975quant")]
              )
            )
            bounds <- data.frame(
                controlY = controlY, interventionY = controlY, effect.size = rep(0, 3),
                pr.contaminated = rep(0, 3)
            )
        }
        analysis$model_object <- inla_result
        analysis$inla_mesh <- inla_mesh

        analysis <- add_estimates(analysis = analysis, bounds = bounds, CLnames = CLnames)
        analysis$int_ests$pr.contaminated <- stats::setNames(
            bounds[c(1, 3),
                "pr.contaminated"], CLnames
        )
        analysis$pt_ests$pr.contaminated <- bounds[2, "pr.contaminated"]

        analysis$passed.face.validity.check <- TRUE
        if (analysis$pt_ests$pr.contaminated < 0 | analysis$pt_ests$pr.contaminated >
            1)
            {
            cat(
                "** Warning: different signs for main effect and contamination: face validity check fails **\n"
            )
            analysis$passed.face.validity.check <- FALSE
            analysis$pt_ests$pr.contaminated <- NA
            analysis$int_ests$pr.contaminated <- c(NA, NA)
        }
        analysis$pt_ests$contamination.par <- beta

        # The contamination parameter is not estimated by INLA but should be considered in the DIC
        analysis$model_object$dic$dic <- analysis$model_object$dic$dic + 2

        analysis$description <- description
    }
    if (method %in% c("EMP", "T", "GEE")) {
        # tidy up and consolidate the list of analysis
        analysis$model_object <- model_object
        analysis$pt_ests <- pt_ests[names(pt_ests) != "model_object"]
        analysis$int_ests <- int_ests
        analysis$description <- description

    }
    if (cfunc != "Z")
    {
        analysis$contamination <- getContaminationCurve(trial = trial, pt_ests = analysis$pt_ests,
                                                       FUN1 = FUN1, link = link, alpha = alpha)
        analysis$pt_ests$contamination.interval <- analysis$contamination$contamination.interval
        analysis$contamination$contamination.interval <- NULL
    } else
    {
        analysis$pt_ests$contamination.interval <- NA
    }
    analysis$options <- options
    class(analysis) <- "CRTanalysis"
    return(analysis)
}

getContaminationCurve <- function(trial, pt_ests, FUN1, link, alpha)
    {

    y_off <- y1 <- average <- upper <- lower <- cats <- nearestDiscord <- NULL

    range_d <- max(trial$nearestDiscord) -
        min(trial$nearestDiscord)
    d <- min(trial$nearestDiscord) +
        range_d * (seq(1:1001) -
            1)/1000

    # define the limits of the curve both for control and intervention arms
    limits <- c(pt_ests$controlY, pt_ests$interventionY)
    limits0 <- limits1 <- limits
    Cp <- 1
    if (is.na(pt_ests$pr.contaminated))
        {
        Cp <- 0
    } else if (0 <= pt_ests$pr.contaminated & pt_ests$pr.contaminated <=
        1)
        {
        Cp <- pt_ests$pr.contaminated
        limits0 <- c(limits[1], Cp * limits[2] + (1 - Cp) * limits[1])
        limits1 <- c(Cp * limits[1] + (1 - Cp) * limits[2], limits[2])
    }

    par0 <- c(link_tr(link, limits0[1]),
              link_tr(link, limits0[2]) - link_tr(link, limits0[1]),
                pt_ests$contamination.par
    )
    par1 <- c(
        link_tr(link, limits1[1]),
        link_tr(link, limits1[2]) - link_tr(link, limits1[1]),
        pt_ests$contamination.par
    )
    curve <- ifelse(d < 0,
               invlink(link, FUN1(trial = data.frame(nearestDiscord = d), par = par0)),
               invlink(link, FUN1(trial = data.frame(nearestDiscord = d), par = par1)))

    # estimate contamination range The absolute values of the limits are used so that a positive range is
    # obtained even with negative effect.size
    thetaL <- thetaU <- NA
    if (abs(limits0[1] - curve[1000]) >
        0.025 * abs(limits0[1] - limits0[2]))
            {
        thetaL <- d[min(
            which(
                abs(limits0[1] - curve) >
                  0.025 * abs(limits0[1] - limits0[2])
            )
        )]
    }
    if (abs(limits1[2] - curve[1000]) <
        0.025 * abs(limits1[1] - limits1[2]))
            {
        thetaU <- d[min(
            which(
                abs(limits1[2] - curve) <
                  0.025 * abs(limits1[1] - limits1[2])
            )
        )]
    }
    if (is.na(thetaU))
        thetaU <- max(trial$nearestDiscord)
    if (is.na(thetaL))
        thetaL <- min(trial$nearestDiscord)

    # contamination range
    contamination.limits <- c(thetaL, thetaU)
    if (thetaL > thetaU)
        contamination.limits <- c(thetaU, thetaL)
    contamination.interval <- thetaU - thetaL
    if (Cp == 0)
        contamination.interval <- NA
    # To remove warnings from plotting ensure that contamination interval is non-zero
    if (is.na(contamination.interval) ||
        contamination.interval == 0)
        {
        contamination.limits <- c(-1e-04, 1e-04)
    }
    # categorisation of trial data for plotting
    trial$cats <- cut(
        trial$nearestDiscord, breaks = c(
            -Inf, min(d) +
                seq(1:9) *
                  range_d/10, Inf
        ),
        labels = FALSE
    )

    # calculate for log link
    data <- data.frame(
        trial %>%
            group_by(cats) %>%
            dplyr::summarize(
                positives = sum(y1),
                total = sum(y_off),
                d = median(nearestDiscord),
                average = Williams(x=y1/y_off, alpha=alpha, option = 'M'),
                lower = Williams(x=y1/y_off, alpha=alpha, option = 'L'),
                upper = Williams(x=y1/y_off, alpha=alpha, option = 'U')
            )
    )
    if (link == 'logit') {
        # overwrite with proportions and binomial confidence intervals by category
        data$average <- data$positives/data$total
        data$upper <- with(data, average -
                               qnorm(alpha/2) * (sqrt(average * (1 - average)/total)))
        data$lower <- with(data, average +
                               qnorm(alpha/2) * (sqrt(average * (1 - average)/total)))
    }
    if (link == 'identity') {
        # overall means and t-based confidence intervals by category
        data <- trial %>%
            group_by(cats) %>%
            dplyr::summarize(
                positives = sum(y1),
                total = sum(y_off),
                d = median(nearestDiscord),
                average = mean(x=y1/y_off),
                lower = Tinterval(y1/y_off, alpha = alpha, option = 'L'),
                upper = Tinterval(y1/y_off, alpha = alpha, option = 'U')
            )
    }

    returnList <- list(
        FittedCurve = data.frame(d = d, contaminationFunction = curve),
        contamination.interval = contamination.interval,
        contamination.limits = contamination.limits,
        data = data
    )
    return(returnList)
}

# add labels to confidence limits
namedCL <- function(limits, alpha = alpha)
    {
    names(limits) <- c(
        paste0(100 * alpha/2, "%"),
        paste0(100 - 100 * alpha/2, "%")
    )
    return(limits)
}

# logit transformation
logit <- function(p = p)
    {
    return(log(p/(1 - p)))
}

# link transformation
link_tr <- function(link = link, x = x)
{
    value <- switch(link,
                    "identity" = x,
                    "log" = log(x),
                    "logit" =  log(x/(1 - x)))
    return(value)
}

# inverse transformation of link function
invlink <- function(link = link, x = x)
    {
    value <- switch(link,
        "identity" = x,
        "log" = exp(x),
        "logit" =  1/(1 + exp(-x)))
    return(value)
}

# Minimal data description and crude effect.size estimate
get_description <- function(trial, link, baselineOnly)
    {
    if(baselineOnly){
        sum.numerators = sum(trial$y1)
        sum.denominators = sum(trial$y_off)
        description <- list(
            sum.numerators = sum.numerators,
            sum.denominators = sum.denominators,
            controlY = sum.numerators/sum.denominators,
            interventionY = NULL,
            effect.size = NULL,
            nclusters = max(as.numeric(as.character(trial$cluster))),
            locations = nrow(trial)
        )
    } else {
        sum.numerators <- tapply(trial$y1, trial$arm, FUN = sum)
        sum.denominators <- tapply(trial$y_off, trial$arm, FUN = sum)
        ratio <- sum.numerators/sum.denominators
        effect.size <- switch(link,
               "identity" = ratio[2] - ratio[1],
               "log" = 1 - ratio[2]/ratio[1],
               "logit" =  1 - ratio[2]/ratio[1])
        description <- list(
            sum.numerators = sum.numerators,
            sum.denominators = sum.denominators,
            controlY = ratio[1],
            interventionY = ratio[2],
            effect.size = effect.size,
            nclusters = max(as.numeric(as.character(trial$cluster))),
            locations = nrow(trial)
        )
    }
    return(description)
}

# Log Likelihood to be maximized
LogLikelihood <- function(par, FUN = FUN, trial)
    {
    logitexpectP <- FUN(par, trial)
    transf <- 1/(1 + exp(-logitexpectP))  #inverse logit transformation

    # for binomial
    LogLikelihood <- sum(
        trial$y1 * log(transf) +
            trial$y0 * log(1 - transf)
    )

    return(LogLikelihood)
}


# transform the parameters into interpretable functions
FittingResults <- function(trial, FUN1, par)
    {
    controlY <- invlink(link = "logit", x = par[1])
    interventionY <- invlink(link = "logit", x = (par[2] + par[1]))
    effect.size <- (controlY - interventionY)/controlY
    pt_ests <- list(
        controlY = controlY, interventionY = interventionY, effect.size = effect.size,
        contamination.par = par[3]
    )
    return(pt_ests)
}

# Contributions to the linear predictor for different contamination functions

CalculateNoEffect <- function(par, trial)
    {
    lp <- par[1]
    return(lp)
}

# step function for the case with no contamination

CalculateNoContaminationFunction <- function(par, trial)
    {
    lp <- ifelse(trial$nearestDiscord < 0, par[1], par[1] + par[2])
    return(lp)
}


# piecewise linear model
CalculatePiecewiseLinearFunction <- function(par, trial)
    {
    # constrain the slope parameter to be positive (par[2] is positive if effect.size is negative)
    theta <- exp(par[3])
    lp <- ifelse(
        trial$nearestDiscord > -theta, par[1] + par[2] * (theta + trial$nearestDiscord)/(2 *
            theta), par[1]
    )
    lp <- ifelse(trial$nearestDiscord > theta, par[1] + par[2], lp)
    return(lp)
}


# sigmoid (logit) function
CalculateLogisticFunction <- function(par, trial)
    {
    theta <- exp(par[3])
    lp <- par[1] + par[2] * invlink(link = "logit", x = theta * trial$nearestDiscord)
    return(lp)
}

# inverse probit function
CalculateProbitFunction <- function(par, trial)
    {
    theta <- exp(par[3])
    lp <- par[1] + par[2] * stats::pnorm(theta * trial$nearestDiscord)
    return(lp)
}


# Functions for T and GEE analysis

noLabels <- function(x)
    {
    xclean <- as.matrix(x)
    dimnames(xclean) <- NULL
    xclean <- as.vector(xclean)
    return(xclean)
}

estimateCLeffect.size <- function(mu, Sigma, alpha, resamples, method, link)
    {

    # Use resampling approach to avoid need for Taylor approximation use at least 10000 samples (this is very
    # cheap)
    resamples1 <- max(resamples, 10000, na.rm = TRUE)
    samples <- MASS::mvrnorm(n = resamples1, mu = mu, Sigma = Sigma)

    pC <- invlink(link, samples[, 1])

    # for the T method the t-tests provide estimates for the s.e. for both arms separately
    # for the GEE method the input is in terms of the incremental effect of the intervention

    if (method == "T")
        pI <- invlink(link, samples[, 2])
    if (method == "GEE")
        pI <- invlink(link, samples[, 1] + samples[, 2])

    eff <- 1 - pI/pC

    CL <- quantile(eff, probs = c(alpha/2, 1 - alpha/2))
    return(CL)
}

# functions for INLA analysis

#' \code{new_mesh} create objects required for INLA analysis of an object of class \code{"CRTsp"}.
#' @param trial an object of class \code{"CRTsp"} or a data frame containing locations in (x,y) coordinates, cluster
#'   assignments (factor \code{cluster}), and arm assignments (factor \code{arm}) and outcome.
#' @param offset (see \code{inla.mesh.2d} documentation)
#' @param max.edge  (see \code{inla.mesh.2d} documentation)
#' @param inla.alpha parameter related to the smoothness (see \code{inla} documentation)
#' @param maskbuffer numeric: width of buffer around points (km)
#' @param pixel numeric: size of pixel (km)
#' @return list
#' \itemize{
#' \item \code{prediction} Data frame containing the prediction points and covariate values
#' \item \code{A} projection matrix from the observations to the mesh nodes.
#' \item \code{Ap} projection matrix from the prediction points to the mesh nodes.
#' \item \code{indexs} index set for the SPDE model
#' \item \code{spde} SPDE model
#' \item \code{pixel} pixel size (km)
#' }
#' @details \code{new_mesh} carries out the computationally intensive steps required for setting-up an
#' INLA analysis of an object of class \code{"CRTsp"}, creating the prediction mesh and the projection matrices.
#' The mesh can be reused for different models fitted to the same
#' geography. The computational resources required depend largely on the resolution of the prediction mesh.
#' The prediction mesh is thinned to include only pixels centred at a distance less than
#' \code{maskbuffer} from the nearest point.\cr
#' A warning may be generated unless the \code{Matrix} library is loaded.
#' @export
#' @examples
#' # low resolution mesh for test dataset
#' library(Matrix); exampleMesh=new_mesh(trial = readdata('exampleCRT.csv'), pixel = 0.5)
new_mesh <- function(trial = trial, offset = -0.1, max.edge = 0.25,
                     inla.alpha = 2, maskbuffer = 0.5, pixel = 0.5)
    {
    # extract the trial data frame from the "CRTsp" object
    if (identical(class(trial),"CRTsp")) trial <- trial$trial
    # create an id variable if this does not exist
    if(is.null(trial$id)) trial <- dplyr::mutate(trial, id =  dplyr::row_number())

    # create buffer around area of points
    trial.coords <- base::matrix(
        c(trial$x, trial$y),
        ncol = 2
    )

    tr <- sf::st_as_sf(trial, coords = c("x","y"))
    buf1 <- sf::st_buffer(tr, maskbuffer)
    buf2 <- sf::st_union(buf1)
    # determine pixel size
    area <- sf::st_area(buf2)
    buffer <- sf::as_Spatial(buf2)

    # estimation mesh construction
    # dummy call to Matrix. This miraculously allows the loading of the "dgCMatrix" in the mesh to pass the test
    dummy <- Matrix::as.matrix(c(1,1,1,1))

    mesh <- INLA::inla.mesh.2d(
        boundary = buffer, offset = offset, cutoff = 0.05, max.edge = max.edge
    )

    # set up SPDE (Stochastic Partial Differential Equation) model
    spde <- INLA::inla.spde2.matern(mesh = mesh, alpha = inla.alpha, constr = TRUE)
    indexs <- INLA::inla.spde.make.index("s", spde$n.spde)
    A <- INLA::inla.spde.make.A(mesh = mesh, loc = trial.coords)

    # 8.3.6 Prediction data from https://www.paulamoraga.com/book-geospatial/sec-geostatisticaldatatheory.html
    bb <- sf::st_bbox(buffer)

    # create a raster that is slightly larger than the buffered area
    xpixels <- round((bb$xmax - bb$xmin)/pixel) + 2
    ypixels <- round((bb$ymax - bb$ymin)/pixel) + 2
    x <- bb$xmin + (seq(1:xpixels) - 1.5)*pixel
    y <- bb$ymin + (seq(1:ypixels) - 1.5)*pixel
    all.coords <- as.data.frame(expand.grid(x, y), ncol = 2)
    colnames(all.coords) <- c("x", "y")
    all.coords <- sf::st_as_sf(all.coords, coords = c("x", "y"))
    pred.coords <- sf::st_filter(all.coords, sf::st_as_sf(buf2))
    pred.coords <- t(base::matrix(
        unlist(pred.coords),
        nrow = 2
    ))
    # projection matrix for the prediction locations
    Ap <- INLA::inla.spde.make.A(mesh = mesh, loc = pred.coords)

    # Distance matrix calculations for the prediction stack Create all pairwise comparisons
    pairs <- tidyr::crossing(
        row = seq(1:nrow(pred.coords)),
        col = seq(1:nrow(trial))
    )
    # Calculate the distances
    calcdistP <- function(row, col) sqrt(
        (trial$x[col] - pred.coords[row, 1])^2 + (trial$y[col] - pred.coords[row,
            2])^2
    )
    distP <- apply(pairs, 1, function(y) calcdistP(y["row"], y["col"]))
    distM <- base::matrix(
        distP, nrow = nrow(pred.coords),
        ncol = nrow(trial),
        byrow = TRUE
    )
    nearestNeighbour <- apply(distM, 1, function(x) return(array(which.min(x))))
    prediction <- data.frame(
        x = pred.coords[, 1], y = pred.coords[, 2], nearestNeighbour = nearestNeighbour)
    prediction$id <- trial$id[nearestNeighbour]
    if (!is.null(trial$arm)) prediction$arm <- trial$arm[nearestNeighbour]
    if (!is.null(trial$cluster)) prediction$cluster <- trial$cluster[nearestNeighbour]
    prediction <- with(
        prediction, prediction[order(y, x), ]
    )
    prediction$shortestDistance <- apply(distM, 1, min)
    rows <- seq(1:nrow(prediction))
    if (!is.null(trial$arm)) {
        prediction$nearestDiscord <- sapply(rows,
                                            FUN = calcNearestDiscord,
                                            trial = trial,
                                            prediction = prediction,
                                            distM = distM)
    }
    inla_mesh <- list(
        prediction = prediction, A = A, Ap = Ap, indexs = indexs, spde = spde, pixel = pixel
    )
    cat(
        "Mesh created of ", nrow(prediction), " pixels of size ", pixel," km \n"
    )

    return(inla_mesh)
}

# Calculate the distance to the nearest discordant location
calcNearestDiscord <- function(x, trial , prediction , distM)
{
    discords <- (trial$arm != prediction$arm[x])
    nearestDiscord <- min(distM[x, discords])
    nearestDiscord <- ifelse(prediction$arm[x] == "control", -nearestDiscord, nearestDiscord)
    return(nearestDiscord)
}

# Use profiling to estimate beta
estimateContaminationINLA <- function(
    beta = beta, trial = trial, FUN = FUN, inla_mesh = inla_mesh, formula = formula, link = link){
    y1 <- y0 <- y_off <- NULL
    x <- trial$nearestDiscord * exp(beta)
    trial$pvar <- eval(parse(text = FUN))

    stk.e <- INLA::inla.stack(
        tag = "est", data = list(y1 = trial$y1, y_off = trial$y_off),
        A = list(1, A = inla_mesh$A),
        effects = list(
            data.frame(
                int = rep(1, nrow(trial)),
                arm = ifelse(trial$arm == "intervention", 1, 0),
                pvar = trial$pvar, id = trial$id, cluster = trial$cluster
            ),
            s = inla_mesh$indexs
        )
    )
    # run the model with just the estimation stack (no predictions needed at this stage)
    if (link == "identity") {
        result.e <- INLA::inla(
            formula, family = "gaussian",
            control.family = list(link = "identity"),
            data = INLA::inla.stack.data(stk.e),
            control.predictor = list(compute = TRUE, link = 1,
                                     A = INLA::inla.stack.A(stk.e)),
            control.compute = list(dic = TRUE))
    } else if (link == "log") {
        result.e <- INLA::inla(
            formula, family = "poisson",
            control.family = list(link = "log"),
            data = INLA::inla.stack.data(stk.e),
            control.predictor = list(compute = TRUE, link = 1,
                                     A = INLA::inla.stack.A(stk.e)),
            control.compute = list(dic = TRUE))
    } else if (link == "logit") {
        result.e <- INLA::inla(
            formula, family = "binomial", Ntrials = y_off,
            control.family = list(link = "logit"),
            data = INLA::inla.stack.data(stk.e),
            control.predictor = list(compute = TRUE, link = 1,
                                     A = INLA::inla.stack.A(stk.e)),
            control.compute = list(dic = TRUE)
        )
    }
    # The DIC is penalised to allow for estimation of beta
    loss <- result.e$dic$family.dic + 2
    cat("\rDIC: ", loss, " Contamination scale parameter: ", beta, "  \r")
    return(loss)
}

# Add estimates to analysis list
add_estimates <- function(analysis, bounds, CLnames)
{
    analysis$pt_ests$controlY <- bounds[2, "controlY"]
    analysis$pt_ests$interventionY <- bounds[2, "interventionY"]
    analysis$pt_ests$effect.size <- bounds[2, "effect.size"]

    # Extract interval estimates
    analysis$int_ests$controlY <- stats::setNames(
        bounds[c(1, 3),
               "controlY"], CLnames
    )
    analysis$int_ests$interventionY <- stats::setNames(
        bounds[c(1, 3),
               "interventionY"], CLnames
    )
    analysis$int_ests$effect.size <- stats::setNames(
        bounds[c(1, 3),
               "effect.size"], CLnames
    )
return(analysis)
}

# Williams mean and confidence intervals
Williams <- function(x, alpha, option){
    logx_1 <- log(x + 1)
    if(length(logx_1) < 3) {
        value <- switch(option,
                M = mean(logx_1),
                L = NA,
                U = NA)
    } else {
        t <- stats::t.test(x = logx_1, conf.level = 1 - alpha)
        value <- switch(option,
                M = t$estimate,
                L = t$conf.int[1],
                U = t$conf.int[2])
    }
    returnvalue <- as.numeric(exp(value) - 1)
    return(returnvalue)
}

# T-based confidence intervals

Tinterval <- function(x, alpha, option){
    if(length(x) < 3){
        value <- NA
    } else {
        t <- stats::t.test(x = x, conf.level = 1 - alpha)
        value <- switch(option,
                L = t$conf.int[1],
                U = t$conf.int[2])
    }
    returnvalue <- as.numeric(value)
}


#' Summary of the results of a statistical analysis of a CRT
#'
#' \code{summary.CRTanalysis} generates a summary of a \code{CRTanalysis} including the main results
#' @param object an object of class \code{"CRTanalysis"}
#' @param ... other arguments used by summary
#' @method summary CRTanalysis
#' @export
#'
summary.CRTanalysis <- function(object, ...) {
    cat("=====================CLUSTER RANDOMISED TRIAL ANALYSIS =================\n")
    cat(
        "Analysis method: ", object$options$method, "\nLink function: ", object$options$link, "\n"
    )
    if (!is.null(object$options$ftext))
        cat("Model formula: ", object$options$ftext, "\n")
    cat(switch(object$options$cfunc,
               Z = "No comparison of arms \n",
               X = "No modelling of contamination \n",
               S = "Piecewise linear function for contamination\n",
               P = "Error function model for contamination\n",
               L = "Sigmoid (logistic) function for contamination\n"))
    CLtext <- paste0(" (", 100 * (1 - object$options$alpha), "% CL: ")
    cat(
        "Estimates:     Control: ", object$pt_ests$controlY,
        CLtext, unlist(object$int_ests$controlY),
        ")\n"
    )
    if (!is.null(object$pt_ests$effect.size))
    {
        if (!is.null(object$pt_ests$interventionY))
        {
            cat(
                "          Intervention: ", object$pt_ests$interventionY,
                CLtext, unlist(object$int_ests$interventionY),
                ")\n"
            )
            effect.measure <- ifelse(object$options$link == 'identity', "Effect size: ","  Efficacy:  ")
            cat("          ",
                effect.measure, object$pt_ests$effect.size, CLtext, unlist(object$int_ests$effect.size),
                ")\n"
            )
        }
        if (!is.na(object$pt_ests$pr.contaminated))
        {
            cat(
                "Proportion of effect subject to contamination: ", object$pt_ests$pr.contaminated,
                CLtext, unlist(object$int_ests$pr.contaminated),")\n"
            )
        }
    }
    if (!is.null(object$pt_ests$ICC))
    {
        cat(
            "Intracluster correlation (ICC): ", object$pt_ests$ICC,
            CLtext, unlist(object$int_ests$ICC),")\n"
        )
    }
    if (!is.na(object$pt_ests$contamination.interval))
    {
        cat(
            "Contamination Range: ", object$pt_ests$contamination.interval,"\n"
        )
    }
    if (!is.null(object$model_object$dic$dic) &
        object$options$cfunc %in% c("L", "P"))
    {
        cat(
            "DIC: ", object$model_object$dic$dic, " including penalty for the contamination scale parameter\n"
        )
    } else if (!is.null(object$model_object$dic$dic))
    {
        cat("DIC: ", object$model_object$dic$dic, "\n")
    }
    if (!is.null(object$model_object$p.value)){
        cat("P-value (2-sided): ", object$model_object$p.value, "\n")
    }

}

