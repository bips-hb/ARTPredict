#' Predict the outcome with the ARTP method
#'
#' @description Predicts the outcome based on the ARTP fit results, of binary data
#' given the observed covariates.The function uses `speedglm`, and in case of
#' non-applicability due to separation issues (in case of rare events), the function
#' uses `glm`.
#'
#' @param fit The output of fit.artp
#' @param X.new The new observation matrix (n x m)
#' @param alpha The significance level, no default
#' @param res.per.group The output of predict.per.group, no default
#' @param predict.all To predict outcome based on all groups (Default = FALSE)
#'
#' @return A list including predictions for the outcome y:
#'         selected.groups: data frame of the p value for each selected group based on alpha
#'         y.prob.per.group: logical matrix of predicted outcome per group using `predict.per.group`
#'         y.group: row sum of prediction per selected group based on alpha
#'         y.hat: the predicted response
#'
#' @seealso artp.fit
#'
#' @examples
#' set.seed(235478965)
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
#' res <- artp.fit(X_train, y_train, groups = groups, verbose = TRUE, trunc.point = 3)
#'
#' pred <- artp.predict(fit = res, X.new = X_test, alpha = .2)
#' table(pred$y.hat, y_test)
#'
#' glm.out <- glm(y_train ~ X_train, family = binomial(link = "logit"))
#' probabilities <- predict(glm.out, as.data.frame(X_test), type = "response")
#' predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
#'
#' table(predicted.classes, y_test)
#'
#'
#' @export
artp.predict <- function(fit, X.new, alpha, res.per.group, predict.all = FALSE) {
  if (missing(res.per.group)) {
    # get predictions for all groups, or in case of alpha, those 'significant'
    res.per.group <- predict.per.group(fit, X.new, alpha, predict.all)
  }

  if (isTRUE(predict.all)) {
    # select the groups that are 'significant' using subset data frame
    selected.groups <- fit$p.values.group[fit$p.values.group$p <= alpha, ]

    # no significant groups at all
    if (nrow(selected.groups) == 0) {
      stop("no significant groups found")
    }

    y.prob.per.group <- res.per.group$y.prob.per.group[, selected.groups$id]
  } else {
    selected.groups <- res.per.group$selected.groups
    y.prob.per.group <- res.per.group$y.prob.per.group
  }

  # probabilities of y being one
  y.group <- rowSums(y.prob.per.group)

  # decide whether TRUE or FALSE
  y.hat <- as.numeric((y.group > 1))

  out <- list(
    selected.groups = selected.groups,
    y.prob.per.group = res.per.group$y.prob.per.group,
    y.group = y.group,
    y.hat = y.hat
  )
  out
}

#' Predict outcome per group of covariates
#'
#' @inherit artp.predict params
#' @seealso artp.predict
#' @importFrom speedglm speedglm
predict.per.group <- function(fit, X.new, alpha, predict.all) {
  # get the old y values. Used to fit the model for each group
  X.old <- fit$X
  y.old <- fit$y
  groups <- fit$groups # the grouping with single_covariates appended
  adjust_vars <- fit$adjust_vars
  parallel <- fit$parallel
  nc <- fit$nc

  # add column names because test and train data need to have same colnames
  if (is.null(dimnames(X.new)[[2]])) {
    dimnames(X.new)[[2]] <- dimnames(X.old)[[2]]
  }

  if (isTRUE(predict.all)) {
    selected.groups <- fit$p.values.group
  } else {
    # select the groups that are 'significant' using dplyr or subset data frame
    # selected.groups <- fit$p.values.group %>% dplyr::filter(p <= alpha)
    selected.groups <- fit$p.values.group[fit$p.values.group$p <= alpha, ]

    # no significant groups at all
    if (nrow(selected.groups) == 0) {
      stop("no significant groups found")
    }
  }

  if (isTRUE(parallel)) {
    # predict for each observation the value of y
    future::plan("multiprocess", workers = nc, gc = TRUE)
  } else {
    future::plan(future::sequential)
  }

  y.prob.per.group <- future.apply::future_mapply(FUN = function(group.id) {
    # get the old data for the current group
    group <- groups[[group.id]]
    data.old <- data.frame(cbind(X.old[, c(adjust_vars, group), drop = FALSE], "y" = y.old))

    # get the new data of the current group
    data.new <- data.frame(X.new[, c(adjust_vars, group), drop = FALSE])

    # fit the model on the old data for that group alone
    model <- tryCatch({
        setTimeLimit(cpu = Inf, elapsed = Inf)
        expr <- speedglm::speedglm(y ~ ., data.old, family = stats::binomial("logit"))
      },
      error = function(e) {
        stats::glm(y ~ ., data.old, family = stats::binomial("logit"))
      }
    )

    # predict whether or not the outcome is 1 given the new data of the current group
    y.new <- stats::predict(model, newdata = data.new, type = "response")

    # decide whether TRUE or FALSE
    y.new <- (y.new > .5)
  }, group.id = selected.groups$id, SIMPLIFY = TRUE, future.packages = c("speedglm"), future.seed = NULL)

  dimnames(y.prob.per.group)[[2]] <- rownames(selected.groups)

  res.per.group <- list(
    y.prob.per.group = y.prob.per.group,
    selected.groups = selected.groups
  )

  res.per.group
}
