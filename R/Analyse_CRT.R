#' Analysis of cluster randomized trial with contamination
#'
#' \\code{Analyse_CRT} carries out a statistical analysis of a cluster randomized trial (CRT).
#' @param trial dataframe: including locations, clusters, arms, and outcomes
#' @param method statistical method with options:
#' 'EMP' - empirical;
#' 'T' - comparison of cluster means by t-test;
#' 'GEE' - Generalised Estimating Equations;
#' 'INLA' - Integrated Nested Laplace Approximation (INLA);
#' 'MCMC' - Markov chain Monte Carlo using JAGS.
#' @param cfunc transformation defining the contamination function
#' options are:
#' 'X': contamination not modelled- the only valid value of 'cfunc' for methods 'EMP,'T' and 'GEE';
#' 'L': inverse logistic (sigmoid)- the default for 'INLA' and 'MCMC' methods;
#' 'P': inverse probit (error function)- available with 'INLA' and 'MCMC' methods;
#' 'S': piecewise linear- only available with the 'MCMC' method
#' @param link link function- options are:
#' 'logit': (the default). The 'numerator', has a binomial distribution with denominator 'denominator'.
#' 'log': The 'numerator' is Poisson distributed with an offset of log('denominator').
#' With the 'INLA' and 'MCMC' methods 'iid' random effects are used to model extra-Poisson variation.
#' 'identity': The outcome is 'numerator'/'denominator'. Normally distributed error function.
#' @param numerator string: name of numerator variable for outcome (if present)
#' @param denominator string: name of denominator variable for outcome data (if present)
#' @param excludeBuffer logical: indicator of whether any buffer zone (records with buffer=TRUE) should be excluded from analysis
#' @param alpha numeric: confidence level for confidence intervals and credible intervals
#' @param baselineOnly logical: indicator of whether required analysis is of effect size or of baseline only
#' @param baselineNumerator string: name of numerator variable for baseline data (if present)
#' @param baselineDenominator string: name of denominator variable for baseline data (if present)
#' @param localisedEffects logical: indicator of whether the model includes local effects with no contamination
#' @param clusterEffects logical: indicator of whether the model includes cluster random effects
#' @param spatialEffects logical: indicator of whether the model includes spatial random effects
#' @param resamples integer: number of samples for sample-based intervals
#' @param requireMesh logical: indicator of whether spatial predictions are required
#' @param inla.mesh name of pre-existing INLA input object created by createMesh()
#' @return list containing the following results of the analysis
#' \\itemize{
#' \\item \\code{description}: Description of the trial dataset
#' \\item \\code{method}: statistical method
#' \\item \\code{pt.ests}: point estimates
#' \\item \\code{int.ests}: interval estimates
#' \\item \\code{model.object}: object returned by the fitting function
#' \\item \\code{contamination}: function values and statistics describing the estimated contamination
#' }
#' @importFrom grDevices rainbow
#' @importFrom stats binomial dist kmeans median na.omit qlogis qnorm quantile rbinom rnorm runif simulate
#' @importFrom utils head read.csv
#' @export
#'
#' @examples
#' # Standard GEE analysis of test dataset ignoring contamination
#' exampleGEE=Analyse_CRT(trial=readdata('test_Simulate_CRT.csv'),method='GEE')

Analyse_CRT <- function(
    trial, method = "GEE", cfunc = "L", link = "logit", numerator = "num",
    denominator = "denom", excludeBuffer = FALSE, alpha = 0.05,
    baselineOnly = FALSE, baselineNumerator = "base_num", baselineDenominator = "base_denom",
    localisedEffects = FALSE, clusterEffects = TRUE, spatialEffects = FALSE,
    resamples = 10000, requireMesh = FALSE, inla.mesh = NULL)
    {
    # Test of validity of inputs
    if (!method %in% c("EMP", "T", "MCMC", "GEE", "INLA"))
        {
        cat("Error: Invalid value for statistical method\n")
        return(NULL)
    }
    if (!cfunc %in% c("S", "L", "P", "X"))
        {
        cat("Error: Invalid contamination function\n")
        return(NULL)
    }

    # MAIN FUNCTION CODE STARTS HERE

    cluster <- NULL

    trial <- convertCRTtodataframe(CRT = trial)

    if ("buffer" %in% colnames(trial) & excludeBuffer)
        {
        trial <- trial[!trial$buffer, ]
    }

    # create names for confidence limits for use throughout
    CLnames <- c(
        paste0(alpha/0.02, "%"),
        paste0(100 - alpha/0.02, "%")
    )

    # trial needs to be ordered for some analyses
    trial <- trial[order(trial$cluster), ]

    # Some statistical methods do not allow for contamination
    if (method %in% c("EMP", "T", "GEE")) cfunc <- "X"

    if (baselineOnly){
        if (method %in% c("EMP", "T", "GEE", "MCMC"))
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

    if (method %in% c("MCMC", "INLA")) {
        # create model formula for display even though this is only used for INLA models

        fterms <- switch(link,
                         "identity" = "y1/y_off ~ 0",
                         "log" = "y1 ~ 0",
                         "logit" = "y1 ~ 0"
        )
        fterms <- c(fterms, switch(
            cfunc, Z = "b0",
            X = "b0 + b1",
            L = "b0 + pvar",
            P = "b0 + pvar"
        ))

        if (localisedEffects & cfunc != 'X')
            fterms <- c(fterms, "b1")
        if (clusterEffects)
            fterms <- c(fterms, "f(cluster, model = \"iid\")")
        if (spatialEffects)
            fterms <- c(fterms, "f(s, model = spde)")
        if (link == 'log')
            fterms <- c(fterms, "f(id, model = \"iid\")")
        options$ftext <- paste(fterms, collapse = " + ")
    }
    model.object <- list()
    pt.ests <- list(contamination.par = NA, pr.contaminated = NA, contamination.interval = NA)
    int.ests <- list(controlY = NA, interventionY = NA, effect.size = NA)
    description <- get_description(trial=trial, link=link, baselineOnly)
    # Specify the function used for calculating the linear predictor
    LPfunction <- c(
        "CalculateNoEffect", "CalculateNoContaminationFunction", "CalculatePiecewiseLinearFunction",
        "CalculateLogisticFunction", "CalculateProbitFunction")[which(cfunc == c("Z", "X", "S", "L", "P"))]
    FUN2 <- FUN1 <- eval(parse(text = LPfunction))


    if (method == "EMP"){
        pt.ests$controlY <- unname(description$controlY)
        pt.ests$interventionY <- unname(description$interventionY)
        pt.ests$effect.size <- unname(description$effect.size)
        pt.ests$contamination.interval <- NA
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
        model.object <- stats::t.test(
            formula = formula, data = clusterSum, alternative = "two.sided",
            conf.level = 1 - alpha, var.equal = TRUE
        )
        pt.ests$p.value <- model.object$p.value
        analysisC <- stats::t.test(
            clusterSum$lp[clusterSum$arm == "control"], conf.level = 1 - alpha)
        pt.ests$controlY <- invlink(link, analysisC$estimate[1])
        int.ests$controlY <- invlink(link, analysisC$conf.int)
        analysisI <- stats::t.test(
            clusterSum$lp[clusterSum$arm == "intervention"], conf.level = 1 - alpha)
        pt.ests$interventionY <- invlink(link, analysisI$estimate[1])
        int.ests$interventionY <- invlink(link, analysisI$conf.int)

        # Covariance matrix (note that two arms are independent so the off-diagonal elements are zero)
        Sigma <- matrix(
            data = c(analysisC$stderr^2, 0, 0, analysisI$stderr^2),
            nrow = 2, ncol = 2)
        if (link == 'identity'){
            pt.ests$effect.size <- pt.ests$controlY - pt.ests$interventionY
            int.ests$effect.size <- unlist(model.object$conf.int)
        }
        if (link %in% c("logit","log")){
            pt.ests$effect.size <- 1 - pt.ests$interventionY/pt.ests$controlY
            int.ests$effect.size <- 1 - exp(-unlist(model.object$conf.int))
        }
    } else if (method == "GEE") {
        # GEE analysis of cluster effects

        if (link == "log")
        {
            fterms <- ifelse(
                cfunc == "Z", "y1 ~ 1 + offset(log(y_off))", "y1 ~ arm + offset(log(y_off))")
            formula <- stats::as.formula(fterms)
            fit <- geepack::geeglm(
                formula = formula, id = cluster, data = trial, family = poisson(link = "log"),
                corstr = "exchangeable", scale.fix = FALSE)
        } else if (link == "logit") {
            fterms <- ifelse(cfunc == "Z", "cbind(y1,y0) ~ 1", "cbind(y1,y0) ~ arm")
            formula <- stats::as.formula(fterms)
            fit <- geepack::geeglm(
                formula = formula, id = cluster, corstr = "exchangeable",
                data = trial, family = binomial(link = "logit"))
        } else if (link == "identity") {
            fterms <- ifelse(cfunc == "Z", "y1/y_off ~ 1", "y1/y_off ~ arm")
            formula <- stats::as.formula(fterms)
            fit <- geepack::geeglm(
                formula = formula, id = cluster, corstr = "exchangeable",
                data = trial, family = gaussian
            )
        }
        options$ftext <- paste(fterms, collapse = " + ")

        summary_fit <- summary(fit)

        z <- -qnorm(alpha/2)  #standard deviation score for calculating confidence intervals
        lp_yC <- summary_fit$coefficients[1, 1]
        se_lp_yC <- summary_fit$coefficients[1, 2]

        clusterSize <- nrow(trial)/nlevels(as.factor(trial$cluster))


        # remove the temporary objects from the dataframe
        fit$model.object$data$y1 <- fit$model.object$data$y0 <- fit$model.object$data$y_off <- NULL
        pt.ests$controlY <- invlink(link, lp_yC)
        int.ests$controlY <- namedCL(
            invlink(link, c(lp_yC - z * se_lp_yC, lp_yC + z * se_lp_yC)),
            alpha = alpha
        )

        pt.ests$model.object <- fit

        # Intracluster correlation
        pt.ests$ICC <- noLabels(summary_fit$corr[1])  #with corstr = 'exchangeable', alpha is the ICC
        se_ICC <- noLabels(summary_fit$corr[2])
        int.ests$ICC <- namedCL(
            noLabels(c(pt.ests$ICC - z * se_ICC, pt.ests$ICC + z * se_ICC)),
            alpha = alpha
        )
        pt.ests$DesignEffect <- 1 + (clusterSize - 1) * pt.ests$ICC  #Design Effect
        int.ests$DesignEffect <- 1 + (clusterSize - 1) * int.ests$ICC

        # Estimation of effect.size does not apply if analysis is of baseline only (cfunc='Z')
        pt.ests$effect.size <- NA
        if (cfunc == "X")
        {
            lp_yI <- summary_fit$coefficients[1, 1] + summary_fit$coefficients[2,
                1]
            se_lp_yI <- sqrt(
                fit$geese$vbeta[1, 1] + fit$geese$vbeta[2, 2] + 2 * fit$geese$vbeta[1,
                  2]
            )

            int.ests$interventionY <- namedCL(
                invlink(link, c(lp_yI - z * se_lp_yI, lp_yI + z * se_lp_yI)),
                alpha = alpha)

            int.ests$effect.size <- estimateCLeffect.size(
                mu = summary_fit$coefficients[, 1], Sigma = fit$geese$vbeta,
                alpha = alpha, resamples = resamples, method = method,
                link = link)

            pt.ests$interventionY <- invlink(link, lp_yI)
            pt.ests$effect.size <- (1 - invlink(link, lp_yI)/invlink(link, lp_yC))
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
        text4 <- "lp[i] <- b0 + b1 * pvar[i]"
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
            b0 ~ dnorm(0, 1E-2) \n
            b1 ~ dnorm(0, 1E-2) \n"

        text7 <- switch(link,
                        "identity" = "yC <- b0 \n
                                      yI <- b0 + b1 \n
                                      tau1 <- 1/(sigma1 * sigma1) \n
                                      sigma1 ~ dunif(0, 2) \n
                                      Es <- yC - yI \n",
                        "log" = "yC <- exp(b0) \n
                                 yI <- exp(b0 + b1) \n
                                 tau1 <- 1/(sigma1 * sigma1) \n
                                 sigma1 ~ dunif(0, 2) \n
                                 Es <- 1 - yI/yC \n",
                        "logit" = "yC <- 1/(1 + exp(-b0)) \n
                                   yI <- 1/(1 + exp(-(b0 + b1))) \n
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
        analysis <- list(model.object = jagsout, pt.ests = list(), int.ests = list())
        analysis$pt.ests$controlY <- jagsout$q50$yC
        analysis$int.ests$controlY <- namedCL(c(jagsout$q2.5$yC, jagsout$q97.5$yC), alpha = alpha)
        analysis$pt.ests$interventionY <- jagsout$q50$yI
        analysis$int.ests$interventionY <- namedCL(c(jagsout$q2.5$yI, jagsout$q97.5$yI), alpha = alpha)
        analysis$pt.ests$effect.size <- jagsout$q50$Es
        analysis$int.ests$effect.size <- namedCL(c(jagsout$q2.5$Es, jagsout$q97.5$Es), alpha = alpha)
        analysis$pt.ests$contamination.interval <- jagsout$q50$theta
        analysis$int.ests$contamination.interval <- namedCL(c(jagsout$q2.5$theta, jagsout$q97.5$theta), alpha = alpha)
        analysis$pt.ests$contamination.par <- jagsout$q50$beta
        analysis$int.ests$contamination.par <- namedCL(c(jagsout$q2.5$beta, jagsout$q97.5$beta), alpha = alpha)
        analysis$pt.ests$pr.contaminated <- jagsout$q50$pcont
        analysis$int.ests$pr.contaminated <- namedCL(c(jagsout$q2.5$pcont, jagsout$q97.5$pcont), alpha = alpha)
        analysis$pt.ests$DIC <- jagsout$DIC
        analysis$model.object$MCMCmodel <- MCMCmodel
    # INLA methods
    } else if (method == "INLA")
        {

        trial <- dplyr::mutate(trial, id =  dplyr::row_number())
        # If spatial predictions are not required a minimal mesh is generated
        ncells <- 50
        if (!requireMesh) ncells <- 5
        if (is.null(inla.mesh)) {
                inla.mesh <- createMesh(
                    trial = trial, offset = -0.1, max.edge = 0.25, inla.alpha = 2,
                    maskbuffer = 0.5, ncells = ncells
                )
        }
        y_off <- NULL
        # specify functional form of sigmoid in distance from boundary 'L' inverse logit; 'P' inverse probit; 'X'
        # or 'Z' do not model contamination
        FUN <- switch(
            cfunc, L = "invlink(link='logit', x)", P = "stats::pnorm(x)", X = NULL, Z = NULL)

        formula <- stats::as.formula(options$ftext)

        spde <- inla.mesh$spde

        effectse <- list(
            df = data.frame(
                b0 = rep(1, nrow(trial)),
                b1 = ifelse(trial$arm == "intervention", 1, 0),
                id = trial$id,
                cluster = trial$cluster
            ),
            s = inla.mesh$indexs
        )
        effectsp <- list(
            df = data.frame(
                b0 = rep(1, nrow(inla.mesh$prediction)),
                b1 = ifelse(inla.mesh$prediction$arm == "intervention", 1, 0),
                id = inla.mesh$prediction$id,
                cluster = inla.mesh$prediction$cluster
            ),
            s = inla.mesh$indexs
        )

        lc <- NULL
        beta <- NA
        if (cfunc %in% c("L", "P"))
            {
            cat("Estimating scale parameter for contamination range", "\n")
            beta <- stats::optimize(
                f = estimateContamination, interval = c(-10, 10),
                trial = trial, FUN = FUN, inla.mesh = inla.mesh, formula = formula,
                tol = 0.1, link = link)$minimum
            x <- trial$nearestDiscord * exp(beta)
            trial$pvar <- eval(parse(text = FUN))

            effectse$df$pvar <- trial$pvar

            x <- inla.mesh$prediction$nearestDiscord * exp(beta)
            inla.mesh$prediction$pvar <- eval(parse(text = FUN))
            effectsp$df$pvar <- inla.mesh$prediction$pvar

            # set up linear contrasts (not required for cfunc='X' or 'Z')
            if (grepl("pvar", options$ftext, fixed = TRUE))
                {
                lc <- INLA::inla.make.lincomb(b0 = 1, pvar = 1)
                if (grepl("b1", options$ftext, fixed = TRUE))
                  {
                  lc <- INLA::inla.make.lincomb(b0 = 1, pvar = 1, b1 = 1)
                }
            } else if (grepl("b1", options$ftext, fixed = TRUE))
                {
                lc <- INLA::inla.make.lincomb(b0 = 1, b1 = 1)
            }
        }
        # stack for estimation stk.e
        stk.e <- INLA::inla.stack(
            tag = "est", data = list(y1 = trial$y1, y_off = trial$y_off),
            A = list(1, A = inla.mesh$A),
            effects = effectse
        )

        # stack for prediction stk.p
        stk.p <- INLA::inla.stack(
            tag = "pred", data = list(y1 = NA, y_off = NA),
            A = list(1, inla.mesh$Ap),
            effects = effectsp
        )

        # stk.full comprises both stk.e and stk.p if a prediction mesh is in use
        stk.full <- INLA::inla.stack(stk.e, stk.p)

        if (link == "identity") {
            inla.result <- INLA::inla(
                formula, family = "gaussian", lincomb = lc,
                control.family = list(link = "identity"),
                data = INLA::inla.stack.data(stk.full),
                control.fixed = list(correlation.matrix = TRUE),
                control.predictor = list(compute = TRUE, link = 1,
                                         A = INLA::inla.stack.A(stk.full)),
                control.compute = list(dic = TRUE))
        } else if (link == "log") {
            inla.result <- INLA::inla(
                formula, family = "poisson", lincomb = lc,
                control.family = list(link = "log"),
                data = INLA::inla.stack.data(stk.full),
                control.fixed = list(correlation.matrix = TRUE),
                control.predictor = list(compute = TRUE, link = 1,
                                         A = INLA::inla.stack.A(stk.full)),
                control.compute = list(dic = TRUE))
        } else if (link == "logit") {
            inla.result <- INLA::inla(
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
        inla.mesh$prediction$prediction <-
                invlink(link, inla.result$summary.linear.predictor[index, "0.5quant"])
        # Compute sample-based confidence limits for intervened outcome and effect.size if intervention effects are
        # estimated
        if (grepl("pvar", options$ftext, fixed = TRUE) |
            grepl("b1", options$ftext, fixed = TRUE))
                {
            # Specify the means of the variables
            mu <- inla.result$summary.lincomb.derived$mean
            names(mu) <- rownames(inla.result$summary.lincomb.derived)
            # Specify the covariance matrix of the variables
            cov <- inla.result$misc$lincomb.derived.covariance.matrix
            sample <- as.data.frame(MASS::mvrnorm(n = 10000, mu = mu, Sigma = cov))
            sample$controlY <- invlink(link, sample$b0)
            # pr.contaminated is the proportion of effect subject to contamination
            if ("b1" %in% names(mu) &
                "pvar" %in% names(mu))
                  {
                sample$interventionY <- invlink(link, sample$lc)
                sample$pr.contaminated <- with(
                  sample, 1 - (controlY - invlink(link, b0 + b1))/(controlY -
                    interventionY)
              )
            } else if ("b1" %in% names(mu))
                {
                sample$interventionY <- invlink(link, sample$b0 + sample$b1)
                sample$pr.contaminated <- 0
            } else if ("pvar" %in% names(mu))
                {
                sample$interventionY <- invlink(link, sample$b0 + sample$pvar)
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
                  link, inla.result$summary.fixed["b0", c("0.025quant", "0.5quant", "0.975quant")]
              )
            )
            bounds <- data.frame(
                controlY = controlY, interventionY = controlY, effect.size = rep(0, 3),
                pr.contaminated = rep(0, 3)
            )
        }
        analysis <- list(
            model.object = inla.result, inla.mesh = inla.mesh, pt.ests = list(),
            int.ests = list())

        analysis <- add_estimates(analysis = analysis, bounds = bounds, CLnames = CLnames)
        analysis$int.ests$pr.contaminated <- stats::setNames(
            bounds[c(1, 3),
                "pr.contaminated"], CLnames
        )
        analysis$pt.ests$pr.contaminated <- bounds[2, "pr.contaminated"]

        analysis$passed.face.validity.check <- TRUE
        if (analysis$pt.ests$pr.contaminated < 0 | analysis$pt.ests$pr.contaminated >
            1)
            {
            cat(
                "** Warning: different signs for main effect and contamination: face validity check fails **\n"
            )
            analysis$passed.face.validity.check <- FALSE
            analysis$pt.ests$pr.contaminated <- NA
            analysis$int.ests$pr.contaminated <- c(NA, NA)
        }
        analysis$pt.ests$contamination.par <- beta

        # The contamination parameter is not estimated by INLA but should be considered in the DIC
        analysis$model.object$dic$dic <- analysis$model.object$dic$dic + 2

        analysis$description <- description
    }
    if (method %in% c("EMP", "T", "GEE"))
        {
        # tidy up and consolidate the list of analysis
        model.object <- pt.ests$model.object
        pt.ests <- pt.ests[names(pt.ests) !=
            "model.object"]
        analysis <- list(
            description = description, pt.ests = pt.ests,
            int.ests = int.ests, model.object = model.object
        )
    }
    if (cfunc != "Z")
    {
        analysis$contamination <- getContaminationCurve(trial = trial, pt.ests = analysis$pt.ests,
                                                       FUN1 = FUN1, link = link, alpha = alpha)
        analysis$pt.ests$contamination.interval <- analysis$contamination$contamination.interval
        analysis$contamination$contamination.interval <- NULL
    } else
    {
        analysis$pt.ests$contamination.interval <- NA
    }
    analysis$options <- options
    class(analysis) <- "CRTanalysis"
    return(analysis)
}

getContaminationCurve <- function(trial, pt.ests, FUN1, link, alpha)
    {

    y_off <- y1 <- average <- upper <- lower <- cats <- nearestDiscord <- NULL

    range_d <- max(trial$nearestDiscord) -
        min(trial$nearestDiscord)
    d <- min(trial$nearestDiscord) +
        range_d * (seq(1:1001) -
            1)/1000

    # define the limits of the curve both for control and intervention arms
    limits <- c(pt.ests$controlY, pt.ests$interventionY)
    limits0 <- limits1 <- limits
    Cp <- 1
    if (is.na(pt.ests$pr.contaminated))
        {
        Cp <- 0
    } else if (0 <= pt.ests$pr.contaminated & pt.ests$pr.contaminated <=
        1)
        {
        Cp <- pt.ests$pr.contaminated
        limits0 <- c(limits[1], Cp * limits[2] + (1 - Cp) * limits[1])
        limits1 <- c(Cp * limits[1] + (1 - Cp) * limits[2], limits[2])
    }

    par0 <- c(link_tr(link, limits0[1]),
              link_tr(link, limits0[2]) - link_tr(link, limits0[1]),
                pt.ests$contamination.par
    )
    par1 <- c(
        link_tr(link, limits1[1]),
        link_tr(link, limits1[2]) - link_tr(link, limits1[1]),
        pt.ests$contamination.par
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
        description <- list(
            sum.numerators = sum(trial$y1),
            sum.denominators = sum(trial$y_off),
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
    pt.ests <- list(
        controlY = controlY, interventionY = interventionY, effect.size = effect.size,
        contamination.par = par[3]
    )
    return(pt.ests)
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

#' \\code{createMesh} Create prediction mesh and other inputs required for INLA analyis of a CRT.
#' @param trial trial dataframe including locations, clusters, arms, and binary outcomes
#' @param offset (see inla.mesh.2d documentation)
#' @param max.edge (see inla.mesh.2d documentation)
#' @param inla.alpha parameter related to the smoothness
#' @param maskbuffer (see inla.mesh.2d documentation)
#' @param ncells resolution of mesh in terms of maximum of linear dimension in pixels
#' @return list containing the mesh
#' \\itemize{
#' \\item \\code{prediction}: Data.table containing the prediction points and covariate values
#' \\item \\code{A}: projection matrix from the observations to the mesh nodes.
#' \\item \\code{Ap}: projection matrix from the prediction points to the mesh nodes.
#' \\item \\code{indexs}:  index set for the SPDE model
#' \\item \\code{spde}: SPDE model
#' }
#' @export
#'
#' @examples
#' # low resolution mesh for test dataset
#' exampleMesh=createMesh(trial = readdata('test_Simulate_CRT.csv'), ncells = 7)
createMesh <- function(
    trial = trial, offset = -0.1, max.edge = 0.25, inla.alpha = 2, maskbuffer = 0.5,
    ncells = 50)
    {
    cat(
        "Creating mesh for INLA analysis: resolution parameter= ", ncells,
        "\n"
    )
    # create an id variable if this does not exist
    if(is.null(trial$id)) trial <- dplyr::mutate(trial, id =  dplyr::row_number())

    # create buffer around area of points
    trial.coords <- base::matrix(
        c(trial$x, trial$y),
        ncol = 2
    )
    sptrial <- sp::SpatialPoints(trial.coords)
    buf1 <- rgeos::gBuffer(sptrial, width = maskbuffer, byid = TRUE)
    buffer <- rgeos::gUnaryUnion(buf1)

    # estimation mesh construction

    mesh <- INLA::inla.mesh.2d(
        boundary = buffer, offset = offset, cutoff = 0.05, max.edge = max.edge
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
    buf.coords <- buffer@polygons[[1]]@Polygons[[1]]@coords
    ind <- sp::point.in.polygon(
        pred.coords[, 1], pred.coords[, 2], buf.coords[, 1], buf.coords[,
            2]
    )
    # prediction locations
    pred.coords <- pred.coords[which(ind == 1),
        ]

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
    distM <- matrix(
        distP, nrow = nrow(pred.coords),
        ncol = nrow(trial),
        byrow = TRUE
    )
    nearestNeighbour <- apply(distM, 1, function(x) return(array(which.min(x))))
    armP <- trial$arm[nearestNeighbour]
    clusterP <- trial$cluster[nearestNeighbour]
    idP <- trial$id[nearestNeighbour]
    prediction <- data.frame(
        x = pred.coords[, 1], y = pred.coords[, 2], nearestNeighbour = nearestNeighbour,
        arm = armP, id = idP, cluster = clusterP
    )
    prediction <- with(
        prediction, prediction[order(y, x),
            ]
    )
    prediction$shortestDistance <- apply(distM, 1, min)
    rows <- seq(1:nrow(prediction))
    prediction$nearestDiscord <- sapply(rows,
                                        FUN = calcNearestDiscord,
                                        trial = trial,
                                        prediction = prediction,
                                        distM = distM)
    inla.mesh <- list(
        prediction = prediction, A = A, Ap = Ap, indexs = indexs, spde = spde
    )
    return(inla.mesh)
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
estimateContamination <- function(
    beta = beta, trial = trial, FUN = FUN, inla.mesh = inla.mesh, formula = formula, link = link){
    y1 <- y0 <- y_off <- NULL
    x <- trial$nearestDiscord * exp(beta)
    trial$pvar <- eval(parse(text = FUN))

    stk.e <- INLA::inla.stack(
        tag = "est", data = list(y1 = trial$y1, y_off = trial$y_off),
        A = list(1, A = inla.mesh$A),
        effects = list(
            data.frame(
                b0 = rep(1, nrow(trial)),
                b1 = ifelse(trial$arm == "intervention", 1, 0),
                pvar = trial$pvar, id = trial$id, cluster = trial$cluster
            ),
            s = inla.mesh$indexs
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
    analysis$pt.ests$controlY <- bounds[2, "controlY"]
    analysis$pt.ests$interventionY <- bounds[2, "interventionY"]
    analysis$pt.ests$effect.size <- bounds[2, "effect.size"]

    # Extract interval estimates
    analysis$int.ests$controlY <- stats::setNames(
        bounds[c(1, 3),
               "controlY"], CLnames
    )
    analysis$int.ests$interventionY <- stats::setNames(
        bounds[c(1, 3),
               "interventionY"], CLnames
    )
    analysis$int.ests$effect.size <- stats::setNames(
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


#' Summarise analysis
#'
#' \code{summary.CRTanalysis} generates a summary description of an analysis of a CRT
#' @param ... other arguments
#' @param object name of analysis
#' @export
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
        "Estimates:     Control: ", object$pt.ests$controlY,
        CLtext, unlist(object$int.ests$controlY),
        ")\n"
    )
    if (object$pt.ests$effect.size != 0)
    {
        if (!is.null(object$pt.ests$interventionY))
        {
            cat(
                "          Intervention: ", object$pt.ests$interventionY,
                CLtext, unlist(object$int.ests$interventionY),
                ")\n"
            )
            effect.measure <- ifelse(object$options$link == 'identity', "Effect size: ","  Efficacy:  ")
            cat("          ",
                effect.measure, object$pt.ests$effect.size, CLtext, unlist(object$int.ests$effect.size),
                ")\n"
            )
        }
        if (!is.na(object$pt.ests$pr.contaminated))
        {
            cat(
                "Proportion of effect subject to contamination: ", object$pt.ests$pr.contaminated,
                CLtext, unlist(object$int.ests$pr.contaminated),")\n"
            )
        }
    }
    if (!is.null(object$pt.ests$ICC))
    {
        cat(
            "Intracluster correlation (ICC): ", object$pt.ests$ICC,
            CLtext, unlist(object$int.ests$ICC),")\n"
        )
    }
    if (!is.na(object$pt.ests$contamination.interval))
    {
        cat(
            "Contamination Range: ", object$pt.ests$contamination.interval,"\n"
        )
    }
    if (!is.null(object$model.object$dic$dic) &
        object$options$cfunc %in% c("L", "P"))
    {
        cat(
            "DIC: ", object$model.object$dic$dic, " including penalty for the contamination scale parameter\n"
        )
    } else if (!is.null(object$model.object$dic$dic))
    {
        cat("DIC: ", object$model.object$dic$dic, "\n")
    }
}

