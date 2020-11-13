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
  
  wrong <- "SL.mean"                        # a insufficient learner stack
  right <- c("SL.glm", "SL.mean", "SL.hal") # a sufficient learner stack
  
  if (n == 200) folds <- 10
  if (n > 200 && n < 1000) folds <- 5
  if (n > 1000) n <- 2
  
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
