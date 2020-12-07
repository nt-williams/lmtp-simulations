#' Truth calculation MTP.
#'
#' @param trt A vector of observed treatments.
tmtp <- function(trt) {                      
  (trt - 1) * (trt - 1 >= 1) + trt * (trt - 1 < 1)
}

#' Compute Density, g
#'
#' @param t The current time point.
#' @param A The observed treatment at time t.
#' @param L The observed covariate at time t.
#' @param prev_A The observed treatment at time t-1.
density_natural <- function(t, A, L, prev_A) {
  (t == 1) * dbinom(A, 5, prob = (L > 1) * 0.5 + (L > 2) * 0.1) + 
    (t > 1) * dbinom(A, 5, prob = prob_A(t, prev_A, L))
}

#' Compute Density, g_d
#'
#' @param t The current time point.
#' @param A The observed treatment at time t.
#' @param L The observed covariate at time t.
#' @param prev_A The observed treatment at time t-1.
density_mtp <- function(t, A, L, prev_A) {
  density_natural(t, A, L, prev_A) * (A - 1 < 1) + 
    density_natural(t, A + 1, L, prev_A) * (A >= 1)
}

#' Compute density ratio
#'
#' @param t The current time point.
#' @param A The observed treatment at time t.
#' @param prev_A The observed treatment at time t-1.
#' @param L The observed covariate at time t.
density_ratio <- function(t, A, L, prev_A) {
  density_mtp(t, A, L, prev_A) / density_natural(t, A, L, prev_A)
}

#' Compute joint density of A and L
#'
#' @param t The current time point
#' @param A The observed treatment at time t.
#' @param L The observed covariate at time t.
#' @param prev_A The observed treatment at time t-1.
#' @param prev_L The observed covariate at time t-1.
prob_AL <- function(t, A, L, prev_A, prev_L) {
  density_natural(t, A, L, prev_A) * (L * prob_L(prev_A, prev_L) + 
    (1 - L) * (1 - prob_L(prev_A, prev_L)))
}

#' Recursive expectation of Y3
#'
#' @param A3 Vector of treatment observations at time 3.
#' @param L3 Vector of covariate observations at time 3.
prob_Y3 <- function(A3, L3) {
  prob_Y(tmtp(0), 1) * prob_AL(4, 0, 1, A3, L3) +
    prob_Y(tmtp(1), 1) * prob_AL(4, 1, 1, A3, L3) +
    prob_Y(tmtp(2), 1) * prob_AL(4, 2, 1, A3, L3) +
    prob_Y(tmtp(3), 1) * prob_AL(4, 3, 1, A3, L3) +
    prob_Y(tmtp(4), 1) * prob_AL(4, 4, 1, A3, L3) +
    prob_Y(tmtp(5), 1) * prob_AL(4, 5, 1, A3, L3) +
    prob_Y(tmtp(0), 0) * prob_AL(4, 0, 0, A3, L3) +
    prob_Y(tmtp(1), 0) * prob_AL(4, 1, 0, A3, L3) +
    prob_Y(tmtp(2), 0) * prob_AL(4, 2, 0, A3, L3) +
    prob_Y(tmtp(3), 0) * prob_AL(4, 3, 0, A3, L3) +
    prob_Y(tmtp(4), 0) * prob_AL(4, 4, 0, A3, L3) +
    prob_Y(tmtp(5), 0) * prob_AL(4, 5, 0, A3, L3)
}

#' Recursive expectation of Y2
#'
#' @param A3 Vector of treatment observations at time 2.
#' @param L3 Vector of covariate observations at time 2.
prob_Y2 <- function(A2, L2) {
  prob_Y3(tmtp(0), 1) * prob_AL(3, 0, 1, A2, L2) +
    prob_Y3(tmtp(1), 1) * prob_AL(3, 1, 1, A2, L2) +
    prob_Y3(tmtp(2), 1) * prob_AL(3, 2, 1, A2, L2) +
    prob_Y3(tmtp(3), 1) * prob_AL(3, 3, 1, A2, L2) +
    prob_Y3(tmtp(4), 1) * prob_AL(3, 4, 1, A2, L2) +
    prob_Y3(tmtp(5), 1) * prob_AL(3, 5, 1, A2, L2) +
    prob_Y3(tmtp(0), 0) * prob_AL(3, 0, 0, A2, L2) +
    prob_Y3(tmtp(1), 0) * prob_AL(3, 1, 0, A2, L2) +
    prob_Y3(tmtp(2), 0) * prob_AL(3, 2, 0, A2, L2) +
    prob_Y3(tmtp(3), 0) * prob_AL(3, 3, 0, A2, L2) +
    prob_Y3(tmtp(4), 0) * prob_AL(3, 4, 0, A2, L2) +
    prob_Y3(tmtp(5), 0) * prob_AL(3, 5, 0, A2, L2)
}

#' Recursive expectation of Y1
#'
#' @param A3 Vector of treatment observations at time 1.
#' @param L3 Vector of covariate observations at time 1.
prob_Y1 <- function(A1, L1) {
  prob_Y2(tmtp(0), 1) * prob_AL(2, 0, 1, A1, L1) +
    prob_Y2(tmtp(1), 1) * prob_AL(2, 1, 1, A1, L1) +
    prob_Y2(tmtp(2), 1) * prob_AL(2, 2, 1, A1, L1) +
    prob_Y2(tmtp(3), 1) * prob_AL(2, 3, 1, A1, L1) +
    prob_Y2(tmtp(4), 1) * prob_AL(2, 4, 1, A1, L1) +
    prob_Y2(tmtp(5), 1) * prob_AL(2, 5, 1, A1, L1) +
    prob_Y2(tmtp(0), 0) * prob_AL(2, 0, 0, A1, L1) +
    prob_Y2(tmtp(1), 0) * prob_AL(2, 1, 0, A1, L1) +
    prob_Y2(tmtp(2), 0) * prob_AL(2, 2, 0, A1, L1) +
    prob_Y2(tmtp(3), 0) * prob_AL(2, 3, 0, A1, L1) +
    prob_Y2(tmtp(4), 0) * prob_AL(2, 4, 0, A1, L1) +
    prob_Y2(tmtp(5), 0) * prob_AL(2, 5, 0, A1, L1)
}

#' Calculate the true value of Y under the MTP and the efficiency bound
#'
#' @param n Sample size to use for truth calculation.
#' @param seed Seed for data generation.
#' @param DAG The DGM definition.
true <- function(n, seed, DAG) {
  set.seed(seed)
  sim_data <- datagen(DAG, n)
  
  rat_4 <- density_ratio(4, sim_data$A_4, sim_data$L_4, sim_data$A_3)
  rat_3 <- density_ratio(3, sim_data$A_3, sim_data$L_3, sim_data$A_2)
  rat_2 <- density_ratio(2, sim_data$A_2, sim_data$L_2, sim_data$A_1)
  rat_1 <- density_ratio(1, sim_data$A_1, sim_data$L_1, 0)
  
  Zn <- t(apply(data.frame(rat_1, rat_2, rat_3, rat_4), 1, cumprod))
  
  Y_d <- data.frame(prob_Y1(tmtp(sim_data$A_1), sim_data$L_1),
                    prob_Y2(tmtp(sim_data$A_2), sim_data$L_2),
                    prob_Y3(tmtp(sim_data$A_3), sim_data$L_3),
                    prob_Y(tmtp(sim_data$A_4), sim_data$L_4), 
                    sim_data$Y)
  
  Y_n <- data.frame(prob_Y1(sim_data$A_1, sim_data$L_1),
                    prob_Y2(sim_data$A_2, sim_data$L_2),
                    prob_Y3(sim_data$A_3, sim_data$L_3),
                    prob_Y(sim_data$A_4, sim_data$L_4))
  
  list(true = mean(Y_d[, 1]), 
       bound = var(rowSums(Zn * (Y_d[, -1] - Y_n)) + Y_d[, 1]))
}
