SL.caretRF <- function(Y, X, newX, family, obsWeights, id, ...) {
  index <- origami::make_folds(n = length(Y), cluster_ids = id, V = 5)
  index <- lapply(index, function(x) x$training_set)
  control <- caret::trainControl(method = "cv", search = 'random', index = index,
                                 verboseIter = TRUE, classProbs = TRUE)
  SL.caret(Y, X, newX, family, obsWeights, method = 'ranger', tuneLength = 100,
           trControl = control, ...)
}

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
#' @param data 
#' @param trt 
#'
#' @return
#' @export
#'
#' @examples
mtp <- function(data, trt) {         
  a <- data[[trt]]
  (a - 1) * (a - 1 >= 1) + a * (a - 1 < 1)
}

#' Simulate data and estimate MTP effect
#'
#' @param tasks
#' @param row the row of data frame containing the task definitions to use
#' @param estimator the estimator to use for estimation ("tml", "sdr", "sub", "ipw")
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

#' Compute an instance of the lmtp simulation
#'
#' @param tasks
#' @param instance an integer between 1-1000 representing an instance 
#'   to run a partition of the lmtp simulation task list
#' @param machines
#' @param save
#' @param estimator
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
