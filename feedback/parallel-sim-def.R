# simulate and estimate linear model with t-distributed errors for a fixed design
# inputs: reps: how many replications?
#         seed: RNG seed
#         data: data.frame containing all and only numeric covariates
#         true_coefs: coefficient vector to use for simulating new responses
#         df: degrees of freedom of the residual error t-distribution
#         cluster: (optional) a SOCKET Cluster for parallel::parLapply
#         cores: (optional) the number of processes to use for parallelization
# output: a matrix of coefficient vectors (each column is one replicate), with
#         attribute "seed" = RNG seed of the generating call for reproducibility.
simulate_parallel <- function(reps, seed, data, true_coef = 0:ncol(data),
                              df = 4, cluster = NULL, cores = NULL) {
  library(parallel)
  check_simulate_inputs(reps, seed, data, true_coef, df)
  set.seed(seed)
  # make sure RNG is set to parallel mode:
  RNGkind("L'Ecuyer-CMRG")
  if (is.null(cores)) {
    # use one less process than maximum so computer is not completely busy
    cores <- detectCores() - 1L
  }
  checkmate::assert_integerish(cores, lower = 1)

  design <- model.matrix(~., data = data)
  expected <- design %*% true_coef

  # parallel apply-function is platform-dependent, compare boot::boot()
  # for similar code...
  if (.Platform$OS.type != "windows") {
    # simply use forking (i.e.: mclapply) for parallelization for unix, mac:
    papply <- function(X, FUN, ...) {
      mclapply(X = X, FUN = FUN, mc.cores = cores, mc.set.seed = TRUE, ...)
    }
  } else {
    # use socket clusters on windows
    if (is.null(cluster)) {
      cluster <- makePSOCKcluster(cores)
      # ... and do clean up after yourself:
      on.exit(stopCluster(cl = cluster))
    }
    checkmate::assert_class(cluster, "cluster")
    # make sure RNG is set to parallel mode:
    clusterSetRNGStream(cluster)
    # define parallel-apply function for windoze:
    papply <- function(X, FUN, ...) {
      parLapply(cl = cluster, X = X, fun = FUN, ...)
    }
    # load needed objects onto cluster:
    clusterExport(
      cl = cluster,
      varlist = c(
        "simulate_once_faster", "simulate_response_faster",
        "estimate_coef_faster", "expected", "design", "df"
      ),
      # variables are found in the execution environment of the call of this
      # function and in its parents. The code below returns the environment that
      # the call is being evaluated in at runtime, see ?sys.nframe
      # the terrible(!) default for clusterExport is to use ".GlobalEnv"
      # to look for the variables in varlist, so ALWAYS set this explicitly:
      envir = sys.frame(sys.nframe())
    )
  }
  # perform parallelization by using parallel-apply on vector [1, 2, ...., reps]
  #  calling the simulate_once function (with identical arguments) each time.
  coefs <- papply(
    X = seq_len(reps),
    FUN = function(i, expected, design, df) {
      unname(simulate_once_faster(expected, design, df))
    },
    expected = expected, design = design, df = df
  )
  coefs <- do.call(cbind, coefs)
  return(structure(coefs, seed = seed))
}

# das Selbe in Grün mit foreach:
simulate_foreach <- function(reps, seed, data, true_coef = 0:ncol(data),
                             df = 4, cluster = NULL, cores = NULL) {
  library(foreach)
  library(doParallel)
  library(doRNG) #
  check_simulate_inputs(reps, seed, data, true_coef, df)
  set.seed(seed)
  if (is.null(cores)) {
    # use one less process than maximum so computer is not completely busy
    cores <- detectCores() - 1L
  }
  checkmate::assert_integerish(cores, lower = 1)

  design <- model.matrix(~., data = data)
  expected <- design %*% true_coef

  # initialize parallel backend:
  registerDoParallel(cores = cores)
  coefs <- matrix(0, nrow = length(true_coef), ncol = reps)
  coefs <- foreach(rep = seq_len(reps),
                   .export = c("simulate_once_faster", "simulate_response_faster",
                               "estimate_coef_faster", "expected", "design", "df"),
                   #.export only necessary under windows....
                   .combine = cbind) %dorng% simulate_once_faster(expected, design, df)
  # delete parallel seeds:
  attr(coefs, "rng") <- NULL
  return(structure(unname(coefs), seed = seed))
}

# das Selbe in Grün mit future:
simulate_future <- function(reps, seed, data, true_coef = 0:ncol(data),
                            df = 4, cluster = NULL, cores = NULL) {

  check_simulate_inputs(reps, seed, data, true_coef, df)
  set.seed(seed)
  if (is.null(cores)) {
    # use one less process than maximum so computer is not completely busy
    cores <- detectCores() - 1L
  }
  checkmate::assert_integerish(cores, lower = 1)

  design <- model.matrix(~., data = data)
  expected <- design %*% true_coef

  # initialize parallel backend:
  future::plan("multiprocess", workers = cores)

  coefs <- future.apply::future_replicate(n = reps,
                                          expr = simulate_once_faster(expected, design, df))

  return(structure(unname(coefs), seed = seed))
}
