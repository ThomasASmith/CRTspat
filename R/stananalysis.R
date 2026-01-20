stananalysis <- function(analysis){

  trial <- analysis$trial
  link <- analysis$options$link
  cfunc <- analysis$options$cfunc
  alpha <- analysis$options$alpha
  fterms <- analysis$options$fterms
  linearity <- analysis$options$linearity
  personalProtection <- analysis$options$personalProtection
  distance <- analysis$options$distance
  scale_par <- analysis$options$scale_par
  log_sp_prior <- analysis$options$log_sp_prior
  clusterEffects <- analysis$options$clusterEffects
  spatialEffects <- analysis$options$spatialEffects
  pixel <- analysis$options$pixel
  control <- analysis$options$control
  iter <- ifelse(is.null(control$iter), 2000, control$iter)
  # thinning is needed if the number of iterations is large to avoid memory problems
  thin <- ceiling(iter/2000)
  control$iter <- NULL
  FUN <- get_FUN(cfunc)

  # fix the maximum of the intercept parameter to correspond to the maximum of the data
  max_intercept  <- ifelse(link %in% c("logit", "cloglog"),
                           link_tr(link = link, x = 0.9999),
                           link_tr(link = link, x = with(analysis$trial, max(y1/y_off))))

  datastan <- list(N = nrow(trial), max_intercept = max_intercept)

  # construct the stan code by concatenating strings
  cb <- "
    }
    " # new line and close brace are needed repeatedly
  functionblock <- ""
  datablock <- "data{
      int<lower=0> N;
      real max_intercept;"
  parameterblock <- "parameters{
      real<upper = max_intercept> intercept;"
  transformedparameterblock <-"transformed parameters {
      vector[N] lp;"
  transformedparameterblock1 <-""
  transformedparameterblock2 <-"
      for(i in 1:N){
         lp[i] = intercept"
  transformedparameterblock3 <- ""
  modelblock <- "model {
      for(i in 1:N){"
  generatedquantitiesblock <-"generated quantities {
      vector[N] log_lik;
      for(i in 1:N){"

  if(identical(link,"identity")){
    datastan$y <- trial$y1/trial$y_off
    datablock <- paste0(datablock,"
        vector[N] y;")
    transformedparameterblock3 <- paste0(transformedparameterblock3,"
            lp[i] = lp[i] * (1 + pr[i] * effect);")
    modelblock <- paste0(modelblock,"
         y[i] ~ normal(lp[i], sigma1);
      }")
    generatedquantitiesblock <- paste0(generatedquantitiesblock,"
         log_lik[i] = normal_lpdf(y1[i] | lp[i], sigma1);")
  } else if(identical(link, 'log')){
    datastan$y1 <- trial$y1
    datastan$y_off <- trial$y_off
    datablock <- paste0(datablock,"
      array[N] int y1;
      vector[N] y_off;")
    modelblock <- paste0(modelblock,"
         y1[i] ~ poisson(Expect_y[i]);
      }")
    transformedparameterblock <- paste0(transformedparameterblock,"
      vector[N] Expect_y;")
    if (!identical(cfunc, "D")) {
      transformedparameterblock3 <- paste0(transformedparameterblock3,"
            Expect_y[i] = exp(lp[i]) * y_off[i];")
    } else {
      transformedparameterblock1 <- paste0(transformedparameterblock1,"
      real<lower=0, upper=1> efficacy;
      efficacy = 1 - exp(effect);")
      transformedparameterblock3 <- paste0(transformedparameterblock3,"
         Expect_y[i] = exp(lp[i]) * y_off[i];
         Expect_y[i] = Expect_y[i] * (1 - pr[i] * efficacy);")
    }
    generatedquantitiesblock <- paste0(generatedquantitiesblock,"
         log_lik[i] = poisson_lpmf(y1[i] | Expect_y[i]);")
  } else if(identical(link, 'logit')){
    datastan$y1 <- trial$y1
    datastan$y_off <- trial$y_off
    datablock <- paste0(datablock,"
      array[N] int y1;
      array[N] int y_off;")
    modelblock <- paste0(modelblock,"
         y1[i] ~ binomial(y_off[i], p[i]);
      }")
    transformedparameterblock <- paste0(transformedparameterblock,"
      vector[N] p;")
    if (!identical(cfunc, "D")) {
      transformedparameterblock3 <- paste0(transformedparameterblock3,"
         p[i] = 1/(1 + exp(-lp[i]));")
    } else {
      transformedparameterblock1 <- paste0(transformedparameterblock1,"
      real<lower=0, upper=1> efficacy;
      efficacy = (exp(-effect)-1)/(exp(intercept) + exp(-effect));")
      transformedparameterblock3 <- paste0(transformedparameterblock3,"
         p[i] = 1/(1 + exp(-lp[i]));
         p[i] = p[i] * (1 - pr[i] * efficacy);")
    }
    generatedquantitiesblock <- paste0(generatedquantitiesblock,"
         log_lik[i] = log((p[i] * y1[i]) + (1 - p[i])*(y_off[i] - y1[i]));")
  } else if(identical(link, 'cloglog')){
    datastan$y1 <- trial$y1
    datastan$y_off <- trial$y_off
    datablock <- paste0(datablock,"
      array[N] int y1;
      array[N] int y_off;")
    parameterblock <- paste0(parameterblock,"
      vector[N] gamma1;
      real<lower=0, upper=2> sigma1;")
    transformedparameterblock1 <- paste0(transformedparameterblock1,"
      vector[N] p;")
    transformedparameterblock3 <- paste0(transformedparameterblock3,"
         p[i] = 1 - exp(-exp(lp[i] + gamma1[i]) * y_off[i]);")
    modelblock <- paste0(modelblock,"
         gamma1[i] ~ normal(0, sigma1);
         y1[i] ~ bernoulli(p[i]);
      }")
    if (identical(cfunc, "D")) {
      transformedparameterblock1 <- paste0(transformedparameterblock1,"
      real<lower=0, upper=1> efficacy;
      efficacy = (exp(-effect)-1)/(exp(intercept) + exp(-effect));")
      transformedparameterblock3 <- paste0(transformedparameterblock3,"
         p[i] = p[i] * (1 - pr[i]*efficacy);")
    }
    generatedquantitiesblock <- paste0(generatedquantitiesblock,"
         log_lik[i] = log1m_exp(-exp(p[i])) - exp(p[i]);")
  }
  if ("arm" %in% fterms) {
    datastan$intervened <- ifelse(trial$arm == "intervention", 1, 0)
    datablock <- paste0(datablock,"
      vector[N] intervened;")
    parameterblock <- paste0(parameterblock,"
      real arm;")
    transformedparameterblock2 <- paste0(transformedparameterblock2,
                                         " + arm * intervened[i]")
  }

  if (clusterEffects) {
    datastan$cluster <- as.numeric(as.character(trial$cluster))
    datastan$ncluster <- max(as.numeric(as.character(trial$cluster)))
    datablock <- paste0(datablock,"
      array[N] int cluster;
      int ncluster;")
    parameterblock <- paste0(parameterblock,"
      vector[ncluster] gamma;
      real<lower=0, upper=2> sigma;")
    transformedparameterblock2 <- paste0(transformedparameterblock2,
                                         " + gamma[cluster[i]]")
    modelblock <- paste0(modelblock, "
      for(ic in 1:ncluster) {
        gamma[ic] ~ normal(0, sigma);
      }")
  }

  if (spatialEffects) {
    # Distance matrix calculations for the ICAR analysis

    # Calculate euclidean distance to discordant arm if this does not yet exist
    trial <- compute_distance(trial, distance = "nearestDiscord")$trial

    # Create all pairwise distances
    pred_coords <- get_pred_coords(trial = trial, maskbuffer = pixel/2, pixel = pixel)
    geodata <- assemble_geodata(trial = trial, pred_coords = pred_coords)
    adjacency_matrix <- adjacency_matrix(dist = geodata$dist)
    message(paste0('CAR term based on grid of ', geodata$N1, ' pixels, each ',pixel*1000, 'm square'))
    datastan$N3 <- geodata$N1
    datastan$pixel <- geodata$pixel
    datastan$N_edges <- nrow(adjacency_matrix)
    datastan$node1 <- adjacency_matrix$node1
    datastan$node2 <- adjacency_matrix$node2
    # This is based on Poisson regression with an ICAR component from:
    # https://mc-stan.org/learn-stan/case-studies/icar_stan.html
    functionblock <- "
    functions {
      real icar_normal_lpdf(vector phi, int N3, array[] int node1,
                            array[] int node2) {
        return -0.5 * dot_self(phi[node1] - phi[node2]);
      }
    }
  "
    datablock <- paste0(datablock, "
      int<lower=1> N3; // total number of pixels
      int pixel[N]; // mapping from samples to pixels
      int<lower=0> N_edges;
      array[N_edges] int<lower=1, upper=N3> node1; // node1[i] adjacent to node2[i]
      array[N_edges] int<lower=1, upper=N3> node2; // and node1[i] < node2[i] ")

    transformedparameterblock2 <- paste0(transformedparameterblock2,
                                         " + phi[pixel[i]] * sigma")
    parameterblock <- paste0(parameterblock,"
      vector[N3] phi; // spatial effects ")
    modelblock <- paste0(modelblock,"
      phi ~ icar_normal(N3, node1, node2);
      // soft sum-to-zero constraint on phi
      // more efficient than mean(phi) ~ normal(0, 0.001)
      sum(phi) ~ normal(0, 0.001 * N3);")
  }

  if (identical(linearity," with estimation of scale parameter. ")) {
    # Create vector of candidate values of scale_par
    # by dividing the prior (on log scale) into equal bins
    # log_sp is the central value of each bin
    nbins <- 20
    binsize <- (log_sp_prior[2] - log_sp_prior[1])/(nbins - 1)
    log_sp <- log_sp_prior[1] + c(0, seq(1:(nbins - 1))) * binsize
    # calculate effect corresponding to first value of sp
    Pr <- compute_effect(trial = trial, distance = distance,
                         scale_par = exp(log_sp[1]), FUN = FUN)
    for(i in 1:(nbins - 1)){
      Pri <- compute_effect(trial = trial,
                            distance = distance, scale_par = exp(log_sp[1 + i]), FUN = FUN)
      Pr <- data.frame(cbind(Pr, Pri))
    }
    log_sp1 <- c(log_sp + binsize/2, log_sp[nbins - 1] + binsize/2)
    datastan$Pr <- as.matrix(Pr)
    datastan$nbins <- nbins
    datastan$log_sp <- log_sp
    datablock <- paste0(datablock,"
      int nbins;
      matrix[N, nbins] Pr;
      vector[nbins] log_sp;")
    parameterblock <- paste0(parameterblock,"
      real log_scale_par;")
    modelblock <- paste0(modelblock,"
      log_scale_par ~ normal(0, 2);")
    transformedparameterblock <- paste0(transformedparameterblock,"
      vector[N] pr;
      real scale_par;
      real wt;
      real increment;")
    transformedparameterblock1 <- paste0(transformedparameterblock1,"
      scale_par = exp(log_scale_par);
      wt = -9.0;
      increment = (log_sp[nbins] - log_sp[1])/(nbins - 1);
      if (log_scale_par < log_sp[1]){
         wt = 9;
         pr = Pr[, 1];
      } else {
         for (j in 1:(nbins - 1)) {
             if(log_scale_par > log_sp[j] && log_scale_par < log_sp[j + 1]){
                 wt = (log_scale_par - log_sp[j])/increment;
                 pr = (1 - wt) * Pr[,j] + wt * Pr[,j + 1];
             }
         }
         if (wt < 0){
             pr = Pr[, nbins];
         }
      }")
    if (!identical(cfunc, "D")) {
      transformedparameterblock2 <- paste0(transformedparameterblock2,
                                           " + pr[i] * effect")
    }
    cfunc <- "O"
  }
  if (identical(cfunc, "R")) {
    trial <- compute_distance(trial, distance = distance,
                              scale_par = scale_par)$trial
    datastan$d <- trial[[distance]]
    datastan$mind <- min(trial[[distance]])
    datastan$maxd <- max(trial[[distance]])
    datablock <- paste0(datablock,"
      vector[N] d;
      real mind;
      real maxd;")
    transformedparameterblock2 <- paste0(transformedparameterblock2,
                                         " + effect*(d[i]-mind)/(maxd-mind)")
  }

  if ("effect" %in% fterms) {
    parameterblock <- paste0(parameterblock,"
      real<upper=2> log_effect;
     ")
    transformedparameterblock <- paste0(transformedparameterblock,"
      real effect;
      effect = -exp(log_effect);")
    modelblock <- paste0(modelblock,"
      log_effect ~ normal(0, 2);")
  }
  generatedquantitiesblock <- paste0(generatedquantitiesblock,"
      }
    }")

  stancode <- paste0(
    functionblock,
    datablock, cb,
    parameterblock, cb,
    transformedparameterblock,
    transformedparameterblock1,
    transformedparameterblock2, ";",
    transformedparameterblock3, cb, cb,
    modelblock, cb,
    generatedquantitiesblock)

  if(analysis$options$verbose) message(cat(stancode))
  message(paste0("\n", "*** Fitting stan model***\n"))
  options(mc.cores = parallel::detectCores())
  fit <- rstan::stan(model_code = stancode,
                     model_name = 'test4',
                     data = datastan,
                     iter = iter,
                     thin = thin,
                     control = control)

  if (identical(cfunc, "E")) cfunc = "ES"
  parameters_to_save <- switch(cfunc,
                               O = c("intercept", "effect", "scale_par"),
                               X = c("intercept"),
                               Z = c("intercept"),
                               R = c("intercept", "effect"))
  if ("arm" %in% fterms) parameters_to_save <- c(parameters_to_save, "arm")
  message(paste0("\n", "*** Calculating goodness-of-fit of stan model***\n"))
  sample <- data.frame(rstan::extract(fit, pars = parameters_to_save, permuted = TRUE))
  # int is a reserved word in stan, so intercept was spelled out
  sample$int <- sample$intercept
  # Use of loo package to compare fit of stan models
  log_lik_1 <- loo::extract_log_lik(fit, merge_chains = FALSE)
  r_eff <- loo::relative_eff(exp(log_lik_1), cores = getOption("mc.cores", 1))
  analysis$model_object <- fit
  analysis$loo <- loo::loo(log_lik_1, r_eff = r_eff, cores = getOption("mc.cores", 1))
  analysis$trial <- trial
  analysis <- extractEstimates(analysis = analysis, sample = sample)
  # distance must be re-computed in the case of surrounds with estimated scale parameter
  analysis$trial <- compute_distance(trial, distance = distance,
                                     scale_par = analysis$options$scale_par)$trial
  analysis$pt_ests$elpd <- analysis$loo$estimates['elpd_loo', 'Estimate']
  return(analysis)
}

get_pred_coords <- function(trial, maskbuffer, pixel) {
  # create buffer around area of points
  tr <- sf::st_as_sf(trial, coords = c("x","y"))
  buf1 <- sf::st_buffer(tr, maskbuffer)
  buf2 <- sf::st_union(buf1)
  # determine pixel size
  area <- sf::st_area(buf2)
  buffer <- sf::as_Spatial(buf2)
  bb <- sf::st_bbox(buffer)

  # create a raster that is slightly larger than the buffered area
  xpixels <- round((bb$xmax - bb$xmin)/pixel) + 2
  ypixels <- round((bb$ymax - bb$ymin)/pixel) + 2
  x <- bb$xmin + (seq(1:xpixels) - 1.5)*pixel
  y <- bb$ymin + (seq(1:ypixels) - 1.5)*pixel
  all_coords <- as.data.frame(expand.grid(x, y), ncol = 2)
  colnames(all_coords) <- c("x", "y")
  all_coords <- sf::st_as_sf(all_coords, coords = c("x", "y"))
  pred_coords <- sf::st_filter(all_coords, sf::st_as_sf(buf2))
  pred_coords <- t(base::matrix(unlist(pred_coords), nrow = 2))
  pred_coords <- data.frame(pred_coords)
  names(pred_coords) <- c("x", "y")
  return(pred_coords)}

assemble_geodata <- function(trial, pred_coords, sampled_only = TRUE){
  all_coords <- rbind(trial[, c("x", "y")], pred_coords[, c("x", "y")])
  all_dist <- as.matrix(dist(all_coords,method = 'euclidean'))
  N <- nrow(trial)
  pred_coords$sampled <- FALSE
  for(i in seq(1:N)){
    nearest_pixel <- which(all_dist[i, N + seq(1:nrow(pred_coords))]
                           == min(all_dist[i, N + seq(1:nrow(pred_coords))]))
    pred_coords$sampled[nearest_pixel] <- TRUE
  }
  ordered_coords <- pred_coords[order(pred_coords$sampled, decreasing = TRUE), ]
  if(sampled_only) ordered_coords <- ordered_coords[ordered_coords$sampled == TRUE, ]
  N1 <- sum(ordered_coords$sampled)
  all_coords_ordered <- rbind(trial[, c("x", "y")], ordered_coords[, c("x", "y")])
  all_dist_ordered <- as.matrix(dist(all_coords_ordered,method = 'euclidean'))
  trial$pixel <- NA
  for(i in seq(1:N)){
    nearest_pixels <- which(all_dist_ordered[i, N + seq(1:nrow(ordered_coords))]
                            == min(all_dist_ordered[i, N + seq(1:nrow(ordered_coords))]))
    trial$pixel[i] <- nearest_pixels[1]
  }
  dist <- as.matrix(dist(ordered_coords, method = 'euclidean'))
  geodata <- list(N = N, N1 = N1, N2 = nrow(ordered_coords) - N1,
                  pixel = trial$pixel, dist = dist)
  return(geodata)
}

adjacency_matrix <- function(dist){
  # The distances are between points in a regular grid, so the
  # minimum is the horizontal or vertical distance and the diagonal is sqrt(2)*min
  # Adjacent points are those with distance <= sqrt(2) * min
  mindist <- min(dist[dist > 0])
  adjacency <- dist < 1.9 * mindist # i.e. < 2
  ivector <- jvector <- c()
  for(i in seq(1:nrow(dist))) {
    jvector <- c(jvector, which(adjacency[i, ]))
    ivector <- c(ivector, rep(i, times = length(which(adjacency[i, ]))))
  }
  adjacency_df <- data.frame(node1 = ivector, node2= jvector)
  adjacency_df <- adjacency_df[adjacency_df$node1 < adjacency_df$node2, ]
  return(adjacency_df)
}
