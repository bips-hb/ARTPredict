#' Apply the ARTP to Binary Data
#'
#' @description Applies the ARTP to binary data given the observed covariates and the
#' response. The number of observations is n. The number of covariates
#' is m.
#'
#' @param X Binary matrix of size n x m, adjustment variables could be
#'          continuous (e.g., age)
#' @param y Binary response vector of length n
#' @param groups List of groups. Each item is a vector with the
#'               indexes of the covariates that belong to that group
#' @param adjust_vars Indexes of adjustment variables for the regression (Default = NULL)
#' @param trunc.point The truncation point used (Default = 5)
#' @param n.permutations Number of permutations (Default = 50)
#' @param verbose If TRUE, shows progress bar (Default = FALSE),
#'                FALSE in case of parallel <- TRUE
#' @param single_covariates If TRUE, covariates that do not belong to a group, get
#'                          their own individual groups (Default = TRUE)
#' @param parallel Boolean, whether to use future::plan("multiprocess")
#'                 (Default = FALSE)
#' @param nc Number of cores/workers to use for future::plan("multiprocess")
#'           (Default = 3)
#'
#' @return A list of input parameters for artp.predict:
#'         p.values.group: data frame of the p value for each group
#'         groups: the groups themselves
#'         X: the raw covariates data
#'         y: the original outcome # too large output
#'         adjust_vars: indexes of adjustment variables
#'         p.values.combined: permutation results
#'         parallel: whether fit and predict can be run using multiprocessor
#'         nc: n workers
#'         n.permutations: the number of permutations
#'
#' @examples
#' m <- 100
#' n <- 2000
#' X <- matrix(rbinom(m * n, 1, .05), ncol = m)
#' X_train <- X[1:1000, ]
#' X_test <- X[1001:2000, ]
#'
#' groups <- list(1:10, 30:40, 80:100, c(3, 5, 30, 33, 39))
#'
#' y <- sapply(1:n, function(i) {
#'   x <- X[i, ]
#'   lg <- -4 + 4 * sum(x[c(3, 5)]) + 4 * sum(x[c(30, 33, 39)])
#'   py <- 1 / (1 + exp(-lg))
#'   rbinom(1, 1, py)
#' })
#'
#' y_train <- y[1:1000]
#' y_test <- y[1001:2000]
#'
#' res <- artp.fit(X = X_train, y = y_train, groups = groups, verbose = TRUE)
#' res$p.values.group
#'
#' @export
artp.fit <- function(X, y, groups, adjust_vars,
                     trunc.point = 5, n.permutations = 50,
                     verbose = FALSE, single_covariates = TRUE,
                     parallel = FALSE, nc = 3) {
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("var_", 1:ncol(X))
  }

  if (!missing(adjust_vars)) {
    # create an iterator
    cov_ind <- setdiff(1:ncol(X), adjust_vars)
    n.cov <- length(cov_ind)
  } else {
    adjust_vars <- NULL
    n.cov <- ncol(X)
    cov_ind <- 1:n.cov
  }

  if (isTRUE(single_covariates)) {
    # any variable not in a group gets assigned their own group
    if (!is.null(names(groups))) {
      gnames <- names(groups)
    } else {
      gnames <- paste0("grp_", 1:length(groups))
    }

    variables_in_groups <- do.call(c, groups)
    variables_not_in_groups <- as.list(setdiff(cov_ind, variables_in_groups))
    groups <- c(groups, variables_not_in_groups)

    # name the new groups
    names(groups) <- c(gnames, colnames(X)[unlist(variables_not_in_groups)])
  }

  p.values.combined <- get.p.value(
    X = X, y = y, n.permutations = n.permutations, n.cov = n.cov,
    cov_ind = cov_ind, adjust_vars = adjust_vars,
    parallel = parallel, verbose = verbose, nc = nc
  )

  # map p-value indexes to group indexes
  colnames.X <- colnames(X)
  colnames.pvalue <- dimnames(p.values.combined)[[2]]
  groups.pvalue <- lapply(groups, function(group) {
    which(colnames.pvalue %in% colnames.X[group])
  })

  # go over each group
  artp.output <- sapply(groups.pvalue, function(group) {
    # get the p values obtained while permuting for the current group
    data <- t(as.matrix(p.values.combined[, group]))
    # apply the ARTP to get the p-values for the group
    ARTP(data, J = min(trunc.point, nrow(data)), n.permutations)
  })

  # create data frame with the results
  res <- data.frame(
    id = 1:length(groups),
    group.size = sapply(groups, function(group) length(group)), # the group sizes
    truncation.point = unlist(artp.output[1, ]), # the truncation points for each group
    p = unlist(artp.output[2, ]) # the p-values associated with each group
  )

  # return the results
  list(
    p.values.group = res, # p-values for each group
    groups = groups, # the groups themselves
    X = X, # the raw covariates data # too large output
    y = y, # the original outcome # too large output
    adjust_vars = adjust_vars, # indexes of adjustment variables,
    p.values.combined = p.values.combined, # permutation results
    parallel = parallel, # if fit and predict can be run using multiprocessor
    nc = nc,
    n.permutations = n.permutations # the number of permutations
  )
}

#' Get p-value for each covariate considering original and permuted outcome
#'
#' @inherit artp.fit description
#' @seealso artp.fit
#' @importFrom speedglm speedglm
get.p.value <- function(X, y,
                        n.permutations,
                        n.cov, cov_ind,
                        adjust_vars,
                        parallel = FALSE,
                        nc = 3,
                        verbose = FALSE) {
  n.permutations.done <- 0

  if (verbose) {
    pb <- utils::txtProgressBar(
      min = 0, max = n.cov * n.permutations,
      title = "no. of permutations", style = 3
    )
  }

  # compute the p-values
  p.values <- sapply(cov_ind, function(i) {
    model <- speedglm::speedglm(y ~ X[, c(adjust_vars, i), drop = FALSE],
      family = stats::binomial(link = "logit")
    )

    # update progress bar
    if (verbose) {
      n.permutations.done <<- n.permutations.done + 1
      utils::setTxtProgressBar(pb, n.permutations.done)
    }

    # get the lowest p-value of the covariates in the model
    as.numeric(stats::coefficients(summary(model))[(length(adjust_vars) + 2), 4])
  })

  # permutate the outcome (y) n.permutations time, create model
  # and obtain the p-value for each covariate

  if (isTRUE(parallel)) {
    future::plan("multiprocess", workers = nc, gc = TRUE)
  } else {
    future::plan(future::sequential)
  }

  p.values.permutations <- future.apply::future_mapply(FUN = function(k) {
    # create permutation
    set.seed(k)
    permutation <- sample(y)

    # leave one covariate out each time and store the p-value
    p.values <- sapply(cov_ind, function(i) {

      # fit the model without the covariate i
      model <- speedglm::speedglm(permutation ~ X[, c(adjust_vars, i), drop = FALSE],
        family = stats::binomial(link = "logit")
      )

      # update progress bar
      if (verbose) {
        n.permutations.done <<- n.permutations.done + 1
        utils::setTxtProgressBar(pb, n.permutations.done)
      }

      # get the lowest p-value of the covariates in the model
      as.numeric(stats::coefficients(summary(model))[(length(adjust_vars) + 2), 4])
    })
  }, k = 1:n.permutations, SIMPLIFY = TRUE, future.packages = c("speedglm"), future.seed = NULL)

  if (verbose) {
    close(pb)
  }

  # combine the actual p-values with the permutated p-values
  p.values.combined <- t(cbind(p.values, p.values.permutations))
  dimnames(p.values.combined) <- list(NULL, colnames(X)[cov_ind])

  p.values.combined
}
