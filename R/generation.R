#' Calculate L probabilities
#'
#' @param A 
#' @param L 
#'
#' @return
#' @export
#'
#' @examples
prob_L  <- function(A, L) {
  plogis(-0.3 * L + 0.5 * A) 
}

prob_A <- function(t, prev_A, L) {
  (t < 4) * plogis(-2 + 1 / (1 + 2 * L + prev_A)) + 
    (t == 4) * plogis(1 - 3 * prev_A + L)
}

#' Calculate Y probabilities
#'
#' @param A 
#' @param L 
#'
#' @return
#' @export
#'
#' @examples
prob_Y <- function(A, L) {
  plogis(-2 + 1 / (1 - 1.2 * A - 0.3 * L)) 
}

tau <- 4

#' The DGM for the simulation
dag <-
  simcausal::DAG.empty() +
  simcausal::node("L", t = 1, distr = "rcat.b1",
                  probs = c(0.5, 0.25, 0.25)) +
  simcausal::node("A", t = 1, distr = "rbinom", size = 5,
                  prob = (L[1] > 1) * 0.5 + (L[1] > 2) * 0.1) +
  simcausal::node("L", t = 2:tau, distr = "rbern",
                  prob = prob_L(A[t - 1], L[t - 1])) +
  simcausal::node("A", t = 2:tau, distr = "rbinom", size = 5,
                  prob = prob_A(t, A[t - 1], L[t])) +
  simcausal::node("Y", t = (tau + 1), distr = "rbern",
                  prob = prob_Y(A[tau], L[tau]), EFU = TRUE)

#' Generate data from the DGM
#'
#' @param DAG
#' @param n number of observations to draw
#' @param vecfun
datagen <- function(DAG = dag, n, vecfun = c("prob_A", "prob_L", "prob_Y")) {
  dag <- simcausal::set.DAG(DAG, vecfun = vecfun)
  data <- suppressWarnings(simcausal::sim(dag, n = as.integer(n)))
  names(data)[substr(names(data), 1, 1) == 'Y'] <- 'Y'
  return(data)
}
