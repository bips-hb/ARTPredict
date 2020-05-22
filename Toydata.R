# ARTPpredict
library(glmnet)
library(devtools)
load_all("C:/Users/foraita/Documents/BIPS/software/ARTPredict")

set.seed(34564)

### Creating a simple toy data set
m = 100
n = 1000
X <- matrix(rbinom(m * n, 1, .10), ncol = m)

# group 1, 2, and 4 influence y
groups  <- list(1:10, 30:40, 80:100, c(3,5,30,33,39))
groups2 <- list(1:10, 30:40, 80:100, c(3,5,30,33,39), c(20:29, 41:79))

y <- sapply(1:n, function(i) {
  x <- X[i, ]
  lg <- -4 + 4*sum(x[c(3, 5)]) + 4*sum(x[c(30, 33, 39)])
  py <- 1 / (1 + exp(-lg))
  rbinom(1,1,py)
})

# lg <- py <- y <- vector()
# for(i in 1:n){
#   x <- X[i, ]
#   lg[i] <- -4 + 4*sum(x[c(3, 5)]) + 4*sum(x[c(30, 33, 39)])
#   py[i] <- 1 / (1 + exp(-lg[i]))
#   y[i] <- rbinom(1,1,py[i])
#   out <- list(y=y, lg=lg, py=py)
# }
# table(out$y)
# y <- out$y

X.new <- matrix(rbinom(m * n, 1, .10), ncol = m)
y.new <- sapply(1:n, function(i) {
  x <- X[i, ]
  lg <- -4 + 4*sum(x[c(3, 5)]) + 4*sum(x[c(30, 33, 39)])
  py <- 1 / (1 + exp(-lg))
  rbinom(1,1,py)
}) 


res_F <- artp.fit(X, y, groups = groups, verbose = TRUE, single_covariates = FALSE)
res_T <- artp.fit(X, y, groups = groups, verbose = TRUE, single_covariates = TRUE)
res_FG <- artp.fit(X, y, groups = groups2, verbose = TRUE, single_covariates = FALSE)
#res$p.values.group


# prediction <- artp.predict(res, X.new, alpha = .05)
# y.hat <- as.numeric((prediction > .5))
# y.diff <- factor(2*y.new - y.hat, labels = c("FP", "TN", "TP", "FN"))
# table(y.diff)


### Function to measure precision / recall and specifity
myfunc <- function(y.hat, y.new){
  y.diff <- table(factor(2*y.new - y.hat, levels = c(-1,0,1,2),
                                          labels = c("FP", "TN", "TP", "FN")))
  
  precision <- y.diff["TP"]/(y.diff["TP"] + y.diff["FP"])
    recall <- y.diff["TP"]/(y.diff["TP"] + y.diff["FN"])
  specifity <- y.diff["TN"]/(y.diff["TN"] + y.diff["FP"])
  bAcc <- (recall + specifity) / 2
  F1score <- 2 * precision * recall /(precision + recall)
  return(list(precision = precision, recall = recall, specifity = specifity, 
              bAcc = bAcc, f1 = F1score) )
  }


### My function for simulation over alpha
mysim <- function(alpha.lauf, model, newdata, newy){
  prediction <- artp.predict(fit = model, X.new = newdata, alpha = alpha.lauf)
  y.hat <- as.numeric((prediction > .5))
  return(myfunc(y.hat, y.new = newy))
}


a <- seq(from = 0.02, to = 0.5, by = 0.01)
pred_T <- sapply(1:length(a), function(i){
          try(mysim(alpha.lauf = a[i], model = res_T, newdata = X.new, newy = y.new))
         })
df_T <- data.frame(t(pred_T), alpha = a)


pred_F <- sapply(1:length(a), function(i){
          try(mysim(alpha.lauf = a[i], model = res_F, newdata = X.new, newy = y.new))
          })
df_F <- data.frame(t(pred_F), alpha = a)


pred_FG <- sapply(1:length(a), function(i){
  try(mysim(alpha.lauf = a[i], model = res_FG, newdata = X.new, newy = y.new))
})
df_FG <- data.frame(t(pred_F), alpha = a)




attach(df_T)
plot(alpha, precision, type = "l", ylim = range(0,1), main = "with singeltons")
lines(alpha, recall, col = "green")
lines(alpha, specifity, col = "blue")
lines(alpha, bAcc, col = "orange")
lines(alpha, f1, col = "red")
detach(df_T)


attach(df_F)
plot(alpha, precision, type = "l", ylim = range(0,1), main = "without singletons")
lines(alpha, recall, col = "green")
lines(alpha, specifity, col = "blue")
lines(alpha, bAcc, col = "orange")
lines(alpha, f1, col = "red")
detach(df_F)


attach(df_FG)
plot(alpha, precision, type = "l", ylim = range(0,1), main = "singeltons in one group")
lines(alpha, recall, col = "green")
lines(alpha, specifity, col = "blue")
lines(alpha, bAcc, col = "orange")
lines(alpha, f1, col = "red")
detach(df_FG)
