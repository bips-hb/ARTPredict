#' Predict the Outcome with the ARTP
#'
#' @param fit The output of fit.artp
#' @param X.new The new observation matrix
#' @param alpha The significance level
#'
#' @return The predictions for the outcome y
#'
# @import dplyr
#'
#' @export
artp.predict <- function(fit, X.new, alpha = .1) {

  # get the old y values. Used to fit the model for each group
  y.old <- fit$y
  X.old <- fit$X
  groups <- fit$groups # the grouping with single_covariates appended
  adjust_vars <- fit$adjust_vars
  parallel <- fit$parallel

  # select the groups that are 'significant'
  # selected.groups <- fit$p.values.group %>% dplyr::filter(p <= alpha)
  selected.groups <- fit$p.values.group[fit$p.values.group$p <= alpha, ]

  # no significant groups at all
  if (nrow(selected.groups) == 0) {
    stop("no significant groups found")
  }

  if (isTRUE(parallel)) {
    # predict for each observation the value of y
    y.prob.per.group <- parallel::mcmapply(function(group.id) {

      # get the old data for the current group
      group <- groups[[group.id]]
      data.old <- X.old[, c(adjust_vars, group)]

      # get the new data of the current group
      data.new <- X.new[, c(adjust_vars, group)]

      # fit the model on the old data for that group alone
      model <- stats::glm(y.old ~ data.old, family = stats::binomial("logit"))

      # predict whether or not the outcome is 1 given the new data of the current group
      y.new <- stats::predict(model, newdata = data.frame(data.new), type = "response")

      # decide whether TRUE or FALSE
      y.new <- (y.new > .5)
    }, mc.cores = fit$nc, group.id = selected.groups$id, SIMPLIFY = TRUE, mc.preschedule = TRUE)

  } else {
    # predict for each observation the value of y
    y.prob.per.group <- sapply(selected.groups$id, function(group.id) {

      # get the old data for the current group
      group <- groups[[group.id]]
      data.old <- X.old[, c(adjust_vars, group)]

      # get the new data of the current group
      data.new <- X.new[, c(adjust_vars, group)]

      # fit the model on the old data for that group alone
      model <- stats::glm(y.old ~ data.old, family = stats::binomial("logit"))

      # predict whether or not the outcome is 1 given the new data of the current group
      y.new <- stats::predict(model, newdata = data.frame(data.new), type = "response")

      # decide whether TRUE or FALSE
      y.new <- (y.new > .5)
    })
  }

  # probabilities of y being one
  y.group <- rowSums(y.prob.per.group)

  # decide whether TRUE or FALSE
  y.hat <- as.numeric((y.group > 1))

  out <- list(selected.groups = selected.groups, y.group = y.group, y.hat = y.hat)
  out
}
