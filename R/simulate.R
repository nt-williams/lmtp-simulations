#' Custom SuperLearner Function for the Highly Adaptive Lasso
#' 
#' Consult the SuperLearner package for more information.
SL.hal90011 <- function(Y, X, newX = NULL, max_degree = NULL,
                        fit_type = "glmnet", n_folds = 3,
                        use_min = TRUE, family = stats::gaussian(),
                        obsWeights = rep(1, length(Y)), id = id, ...) {
  n <- length(Y)
  if (!is.matrix(X)) {
    X_in <- as.matrix(X)
  } else {
    X_in <- X
  }
  if (!is.null(newX) & !is.matrix(newX)) {
    newX_in <- as.matrix(newX)
  } else {
    newX_in <- newX
  }
  hal_out <- hal9001::fit_hal(Y = Y, X = X_in, max_degree = max_degree, id = id,
                              fit_type = fit_type, n_folds = n_folds, use_min = use_min,
                              lambda = seq(1 / n^2, 1 / sqrt(n), length.out = 50),
                              family = "gaussian", weights = obsWeights, yolo = FALSE)
  if (!is.null(newX)) {
    pred <- stats::predict(hal_out, new_data = newX_in)
  } else {
    pred <- stats::predict(hal_out, new_data = X_in)
  }
  fit <- list(object = hal_out)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- "SL.hal9001"
  return(out)
}

#' Simulation shift function
#'
#' @param data A dataframe.
#' @param trt Name of treatment variable to be shifted.
mtp <- function(data, trt) {         
  a <- data[[trt]]
  (a - 1) * (a - 1 >= 1) + a * (a - 1 < 1)
}

#' Run a single LMTP Simulation
#'
#' @param seed A seed to set for reproducibility.
#' @param n The number of observations to be used for estimation.
#' @param type An integer between 1 and 4 indicating the type of model misspecification.
#' @param estimator The estimator to simulate ("tml", "sdr", "ipw", "sub").
lmtp_simulate <- function(seed, n, type, estimator) {
  set.seed(seed)
  sim_df <- suppressMessages(datagen(dag, n))
  
  a <- c("A_1", "A_2", "A_3", "A_4")                # names of treatment variables
  w <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4")) # names of covariates
  y <- "Y"                                          # final outcome
  
  wrong <- "SL.mean"
  right <- c("SL.glm", "SL.mean", "SL.hal90011")
  
  if (n == 200) folds <- 10
  if (n > 200 && n < 1000) folds <- 5
  if (n > 1000) folds <- 2
  
  switch(
    type,
    `1` = {
      lrnrs_trt <- lrnrs_out <- list(right, right, right, right)
    },
    `2` = {
      lrnrs_trt <- list(right, right, wrong, wrong)
      lrnrs_out <- list(wrong, wrong, right, right)
    },
    `3` = {
      lrnrs_trt <- list(wrong, wrong, wrong, right)
      lrnrs_out <- list(right, right, right, wrong)
    },
    `4` = {
      lrnrs_trt <- lrnrs_out <- list(wrong, wrong, wrong, wrong)
    }
  )
  
  if (estimator == "sub") {
    out <- lmtp::lmtp_sub(data = sim_df, trt = a, outcome = y,
                          time_vary = w, k = Inf, shift = mtp, 
                          learners = lrnrs_out, folds = folds)
    return(out)
  }
  
  if (estimator == "ipw") {
    out <- lmtp::lmtp_ipw(data = sim_df, trt = a, outcome = y,
                          time_vary = w, k = Inf, shift = mtp,
                          learners = lrnrs_trt, folds = folds)
    return(out)
  }
  
  if (estimator == "tml") {
    out <- lmtp::lmtp_tmle(data = sim_df, trt = a, outcome = y,
                           time_vary = w, k = Inf, shift = mtp,
                           learners_outcome = lrnrs_out,
                           learners_trt = lrnrs_trt,
                           folds = folds)
    return(out)
  }
  
  if (estimator == "sdr") {
    out <- lmtp::lmtp_sdr(data = sim_df, trt = a, outcome = y,
                          time_vary = w, k = Inf, shift = mtp,
                          learners_outcome = lrnrs_out,
                          learners_trt = lrnrs_trt,
                          folds = folds)
    return(out)
  }
}

globals <- ls()

#' Compute an instance of the LMTP simulation
#'
#' @param tasks A data frame containing task definitions.
#' @param instance An integer representing a Slurm instance 
#'   to run a partition of the lmtp simulation task list.
#' @param machines The number of machines to be used for the simulation.
#' @param save Path to a folder where results will be saved.
#' @param estimator The estimator to simulate ("tml", "sdr", "ipw", "sub").
partition <- function(tasks, instance, machines, save, estimator) {
  index <-
    (1:nrow(tasks))[(1:nrow(tasks) - 1) %% machines + 1 == instance]
  
  out <- list()
  for (i in 1:length(index)) {
    row <- index[i]
    seed <- tasks[row, "seed"]
    type <- tasks[row, "type"]
    n <- tasks[row, "n"]
    out[[i]] <- try(lmtp_simulate(seed, n, type, estimator))
    
    saveRDS(out[[i]], file.path(save, paste0(estimator, "_", index[i], ".rds")))
  }
}
