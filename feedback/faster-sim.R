# simulate and estimate linear model with t-distributed errors for a fixed design
# inputs: reps: how many replications?
#         seed: RNG seed
#         data: data.frame containing all and only numeric covariates
#         true_coefs: coefficient vector to use for simulating new responses
#         df: degrees of freedom of the residual error t-distribution
# output: a matrix of coefficient vectors (each column is one replicate), with
#         attribute "seed" = RNG seed of the generating call for reproducibility.
simulate_faster <- function(reps, seed, data, true_coef = 0:ncol(data), df = 4L) {
  # change 0: input checks
  check_simulate_inputs(reps, seed, data, true_coef, df)

  set.seed(seed)
  # change 1: pre-allocate container for results
  coefs <- matrix(0, nrow = length(true_coef), ncol = reps)

  # change 2a: pre-compute design X, expected (X*beta):
  design <- model.matrix(~., data = data)
  expected <- design %*% true_coef

  for (rep in seq_len(reps)) {
    coefs[, rep] <- simulate_once_faster(expected, design, df)
  }
  return(structure(coefs, seed = seed))
}

simulate_once_faster <- function(expected, design, df) {
  response <- simulate_response_faster(expected, df)
  estimate_coef_faster(response, design)
}

simulate_response_faster <- function(expected, df) {
  # change 2b: reuse <expected>, don't re-compute
  expected + rt(length(expected), df = df)
}

estimate_coef_faster <- function(response, design) {
  # change 3: .lm.fit instead of lm to avoid unnecessary input checks
  model <- .lm.fit(y = response, x = design)
  coef(model)
}

# change 0: input checks
check_simulate_inputs <- function(reps, seed, data, true_coef, df) {
  checkmate::assert_count(reps)
  checkmate::assert_integer(seed, lower = 0)
  checkmate::assert_data_frame(data, types = "numeric")
  checkmate::assert_numeric(true_coef,
    finite = TRUE, any.missing = FALSE,
    len = ncol(data) + 1
  )
  checkmate::assert_numeric(df, lower = 0)
}
