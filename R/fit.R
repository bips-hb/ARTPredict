#' Apply the ARTP to Binary Data
#'
#' Applies the ARTP to binary data given the observed covariates and the
#' response. The number of observations is n. The number of covariates
#' is m.
#'
#' @param X Binary matrix of size n x m, adjustment variables are could be continous (e.g., age)
#' @param y Binary response vector of length n
#' @param groups List of groups. Each item is a vector with the
#'               indices of the covariates that belong to that group
#' @param adjust_vars Indices of adjustment variables for the regression (Default = NULL)
#' @param trunc.point The truncation point used (Default = 5)
#' @param n.permutations Number of permutations (Default = 50)
#' @param verbose If TRUE, shows progress bar (Default = FALSE)
#' @param single_covariates If TRUE, covariates that do not belong to a group, get
#'                          their own individual groups (Default = TRUE)
#' @param parallel Boolean, whether to use parallel:mcapply (Default = FALSE)
#' @param nc Number of cores to use for parallel:mcmapply (Default = 3)
#'
#' @return A data frame p.values.groups with the p value for each
#'         group
#'
#' @example
#' m = 100
#' n = 2000
#' X <- matrix(rbinom(m * n, 1, .05), ncol = m)
#' X_train <- X[1:1000,]
#' X_test <- X[1001:2000,]
#'
#' groups = list(1:10, 30:40, 80:100, c(3,5,30,33,39))
#'
#' y <- sapply(1:n, function(i) {
#'   x <- X[i, ]
#'   lg <- -4 + 4 * sum(x[c(3, 5)]) + 4 * sum(x[c(30, 33, 39)])
#'   py <- 1 / (1 + exp(-lg))
#'   rbinom(1,1,py)
#'   })
#'
#' y_train <- y[1:1000]
#' y_test  <- y[1001:2000]
#'
#' res <- artp.fit(X_train, y_train, groups = groups, verbose = T)
#' res$p.values.group
#'
#' @export
artp.fit <- function(X, y, groups, adjust_vars = NULL,
                     trunc.point = 5, n.permutations = 50,
                     verbose = FALSE, single_covariates = TRUE,
                     parallel = FALSE, nc = 3) {

  if (is.null(colnames(X))) {
    colnames(X) <- paste0('var_', 1:ncol(X))
  }

  if (!missing(adjust_vars)) {
    # recode X and n.cov iterator in case of adjust_vars
    # X.old <- X
    # X.adjust_vars <- X.old[, adjust_vars]
    # X.new <- subset(X, select = -adjust_vars, drop = TRUE)
    # n.cov <- ncol(X.new)

    # create an iterator
    cov_ind <- setdiff(1:ncol(X), adjust_vars)

  } else {
    adjust_vars <- NULL
    n.cov <- ncol(X)
    cov_ind <- 1:n.cov
  }

  if (single_covariates) {
    # Any variable not in a group gets assigned their own group
    variables_in_groups <- do.call(c, groups)
    variables_not_in_groups <- as.list(setdiff(cov_ind, variables_in_groups))
    groups <- c(groups, variables_not_in_groups)
    # uncomment to name the new groups
    # names(groups) <- c(gnames, paste0("grp", glen + 1:length(variables_not_in_groups)))
  }

  n.permutations.done <- 0

  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = n.cov * n.permutations, title = "no. of permutations", style = 3)
  }

  # compute the p-values
  p.values <- sapply(cov_ind, function(i) {
    model <- stats::glm(y ~ X[, c(adjust_vars, i)], family = stats::binomial(link = "logit"))

    # update progress bar
    if (verbose) {
      n.permutations.done <<- n.permutations.done + 1
      utils::setTxtProgressBar(pb, n.permutations.done)
    }

    # get the lowest p-value of the covariates in the model
    # stats::coefficients(summary(model))[2, 4]
    stats::coefficients(summary(model))[(length(adjust_vars) + 2), 4]
  })

  # permutate the output (y) n.permutations time, fit a logistic
  # regression and obtain the p-value for each covariate

  if (isTRUE(parallel)) {
    p.values.permutations <- parallel::mcmapply(function(k) {

      # create permutation
      set.seed(k)
      permutation <- sample(y)

      # leave one covariate out each time and store the p-value
      p.values <- sapply(cov_ind, function(i) {

        # fit the model without the covariate i
        model <- stats::glm(y ~ X[, c(adjust_vars, i)], family = stats::binomial(link = "logit"))

        # update progress bar
        if (verbose) {
          n.permutations.done <<- n.permutations.done + 1
          utils::setTxtProgressBar(pb, n.permutations.done)
        }

        # get the lowest p-value of the covariates in the model
        stats::coefficients(summary(model))[(length(adjust_vars) + 2), 4]
      })

    }, mc.cores = nc, k = 1:n.permutations, SIMPLIFY = TRUE, mc.preschedule = TRUE)
  } else {
    p.values.permutations <- sapply(1:n.permutations, function(k) {

      # create permutation
      set.seed(k)
      permutation <- sample(y)

      # leave one covariate out each time and store the p-value
      p.values <- sapply(cov_ind, function(i) {

        # fit the model without the covariate i
        model <- stats::glm(y ~ X[, c(adjust_vars, i)], family = stats::binomial(link = "logit"))

        # update progress bar
        if (verbose) {
          n.permutations.done <<- n.permutations.done + 1
          utils::setTxtProgressBar(pb, n.permutations.done)
        }

        # get the lowest p-value of the covariates in the model
        stats::coefficients(summary(model))[(length(adjust_vars) + 2), 4]
      })
    })
  }

  if (verbose) {
    close(pb)
  }

  # turn p.values.permutations into a matrix (redundant)
  # p.values.permutations <- matrix(p.values.permutations, ncol = n.permutations)

  # combine the actual p-values with the permutated p-values
  p.values.combined <- t(cbind(p.values, p.values.permutations))
  dimnames(p.values.combined) <- list(NULL, colnames(X)[cov_ind])

  # map pvalue indices to group indices
  colnames.X <- colnames(X)
  colnames.pvalue <- dimnames(p.values.combined)[[2]]

  # groups.pvalue <- lapply(groups, function(group) {
  #   which(colnames.pvalue %in% colnames.X[group])
  # })
  groups.pvalue <- groups

  # go over each group
  artp.output <- sapply(groups.pvalue, function(group) {
    # get the p values obtained while permutating for the current group
    data <- t(as.matrix(p.values.combined[, group]))
    # apply the ARTP to get the p-values for the group
    ARTP(data, J = min(trunc.point, nrow(data)), n.permutations)
  })

  # create data frame with the results
  res <- data.frame(
    id = 1:length(groups.pvalue),
    group.size = sapply(groups.pvalue, function(group) length(group)), # the group sizes
    truncation.point = unlist(artp.output[1, ]), # the truncation points for each group
    p = unlist(artp.output[2, ]) # the p-values associated with each group
  )

  # return the results
  list(
    p.values.group = res, # p-values for each group
    groups = groups, # the groups themselves
    # groups.new = groups.pvalue, # the recoded group indices in case of adjust_vars
    X = X, # the raw covariate data # too large output
    y = y, # the original output # too large output
    adjust_vars = adjust_vars, # indices of adjustment variables
    parallel = parallel, # if fit and predict can be run using mcapply
    nc = nc,
    n.permutations = n.permutations # the number of permutations
  )
}
