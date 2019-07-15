#' Apply the ARTP to Binary Data
#'
#' Applies the ARTP to binary data given the observed covariates and the
#' response. The number of observations is n. The number of covariates
#' is m.
#'
#' @param X Binary matrix of size n x m
#' @param y Binary response vector of length n
#' @param groups List of groups. Each item is a vector with the
#'               indices of the covariates that belong to that group
#' @param trunc.point The truncation point used (Default = 5)
#' @param n.permutations Number of permutations (Default = 50)
#' @param verbose If TRUE, shows progress bar (Default = FALSE)
#' @param single_covariates If TRUE, covariates that do not belong to a group, get
#'                          their own individual groups (Default = TRUE)
#'
#' @return A data frame p.values.groups with the p value for each
#'         group
#'
#' @export
artp.fit <- function(X, y, groups, trunc.point = 5, n.permutations = 50, verbose = FALSE, single_covariates = TRUE) {

  n.cov <- ncol(X)

  n.permutations.done <- 0

  if (single_covariates) {
    # Any variable not in a group gets assigned their own group
    variables_in_groups <- do.call(c, groups)
    variables_not_in_groups <- as.list(setdiff(1:n.cov, variables_in_groups))
    groups <- c(groups, variables_not_in_groups)
  }

  if (verbose) {
    pb <- txtProgressBar(min = 0, max = n.cov*n.permutations, title = "no. of permutations", style = 3)
  }

  # compute the p-values
  p.values <- sapply(1:n.cov, function(i) {
    model <- glm(y ~ X[,i], family = binomial(link = "logit"))

    # update progress bar
    if (verbose) {
      n.permutations.done <<- n.permutations.done + 1
      setTxtProgressBar(pb, n.permutations.done)
    }

    # get the lowest p-value of the covariates in the model
    coefficients(summary(model))[2,4]
    })

  # permutate the output (y) n.permutations time, fit a logistic
  # regression and obtain the p-value for each covariate
  p.values.permutations <- sapply(1:n.permutations, function(k) {

    # create permutation
    permutation <- sample(y)

    # leave one covariate out each time and store the p-value
    p.values <- sapply(1:n.cov, function(i) {

      # fit the model without the covariate i
      model <- glm(permutation ~ X[,i], family = binomial(link = "logit"))

      # update progress bar
      if (verbose) {
        n.permutations.done <<- n.permutations.done + 1
        setTxtProgressBar(pb, n.permutations.done)
      }

      # get the lowest p-value of the covariates in the model
      coefficients(summary(model))[2,4]
    })

  })

  if (verbose) {
    close(pb)
  }

  # turn p.values.permutations into a matrix
  p.values.permutations <- matrix(p.values.permutations, ncol = n.permutations)

  # combine the actual p-values with the permutated p-values
  p.values.combined <- t(cbind(p.values, p.values.permutations))
  colnames(p.values.combined) <- NULL

  # go over each group
  artp.output <- sapply(groups, function(group) {
    # get the p values obtained while permutating for the current group
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
    X = X, # the raw covariate data
    y = y, # the original output
    n.permutations = n.permutations # the number of permutations
  )
}
