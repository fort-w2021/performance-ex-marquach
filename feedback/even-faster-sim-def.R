# simulate and estimate linear model with t-distributed errors for a fixed design
# inputs: reps: how many replications?
#         seed: RNG seed
#         data: data.frame containing all and only numeric covariates
#         true_coefs: coefficient vector to use for simulating new responses
#         df: degrees of freedom of the residual error t-distribution
# output: a matrix of coefficient vectors (each column is one replicate), with
#         attribute "seed" = RNG seed of the generating call for reproducibility.
simulate_even_faster <- function(reps, seed, data, true_coef = 0:ncol(data), df = 4L) {
  check_simulate_inputs(reps, seed, data, true_coef, df)

  set.seed(seed)
  coefs <- matrix(0, nrow = length(true_coef), ncol = reps)

  design <- model.matrix(~., data = data)
  expected <- design %*% true_coef
  # change: pre-compute factor (X'X)^-1 X' for normal equation
  #           beta_hat = (X'X)^-1 X'y
  # (also: never use solve(matrix) * vector, always use solve(matrix, vector):
  # numerically much more stable & slightly faster)
  coefficient_matrix <- solve(crossprod(design), t(design))

  for (rep in seq_len(reps)) {
    coefs[, rep] <- simulate_once_even_faster(expected, coefficient_matrix, df)
  }
  return(structure(coefs, seed = seed))
}

simulate_once_even_faster <- function(expected, coefficient_matrix, df) {
  response <- simulate_response_faster(expected, df)
  # change: get coefficients with simple matrix-vector-multiplication
  coefficient_matrix %*% response
}
