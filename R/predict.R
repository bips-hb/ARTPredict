#' Predict the Outcome with the ARTP
#'
#' @param fit The output of fit.artp
#' @param X.new The new observation matrix
#' @param alpha The significance level
#'
#' @return The predictions for the outcome y
#' @export
artp.predict <- function(fit, X.new, alpha = .1) {

  # get the old y values. Used to fit the model for each group
  y.old <- fit$y
  X.old <- fit$X
  groups <- fit$groups # the grouping

  # number of observations
  n.obs <- nrow(X)

  # select the groups that are 'significant'
  selected.groups <- fit$p.values.group %>% filter(p <= alpha)

  # predict for each observation the value of y
  y.prob.per.group <- sapply(selected.groups$id, function(group.id) {

      # get the old data for the current group
      group <- groups[[group.id]]
      data.old <- X.old[, group]

      # get the new data of the current group
      data.new <- X.new[, group]

      # fit the model on the old data for that group alone
      model <- glm(y.old ~ data.old, family = binomial("logit"))

      # predict whether or not the outcome is 1 given the new data of the current group
      y.new <- predict(model, newdata = data.frame(data.new), type = "response")

      # decide whether TRUE or FALSE
      y.new <- (y.new > .5)
    })

  # probabilities of y being one
  y.prob <- rowMeans(y.prob.per.group)

  # decide whether TRUE or FALSE
  y.hat = as.numeric((y.prob > .5))
}
