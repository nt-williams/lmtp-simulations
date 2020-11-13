SL.hal <- function(Y, X, newX = NULL, id = NULL,  max_degree = NULL,
                   fit_type = c("glmnet", "lassi"), n_folds = ifelse(length(Y) < 500, 10, 5),
                   lambda = seq(1 / length(Y)^2, 1 / sqrt(length(Y)), length.out = 50), 
                   use_min = TRUE, family = stats::gaussian(),
                   obsWeights = rep(1, length(Y)), ...) {
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
  
  if (family$family == "gaussian") {
    hal_out <- hal9001::fit_hal(
      Y = Y, X = X_in, max_degree = max_degree, fit_type = fit_type, id = id, 
      n_folds = n_folds, use_min = use_min, family = "gaussian",
      weights = obsWeights, yolo = FALSE
    )
  }
  
  if (family$family == "binomial") {
    hal_out <- hal9001::fit_hal(
      Y = Y, X = X_in, max_degree = max_degree, fit_type = fit_type, lambda = lambda,
      n_folds = n_folds, use_min = use_min, family = "binomial", id = id,
      weights = obsWeights, yolo = FALSE
    )
  }

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

predict.SL.hal9001 <- function(object, newdata, ...) {
  if (!is.matrix(newdata)) {
    newdata_in <- as.matrix(newdata)
  } else {
    newdata_in <- newdata
  }
  
  pred <- stats::predict(object$object, new_data = newdata_in)
  return(pred)
}
