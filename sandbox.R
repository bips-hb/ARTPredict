library(DDTA)
library(readr)
library(ARTPredict)


### Creating a simple toy data set
m = 100
n = 1000
X <- matrix(rbinom(m * n, 1, .05), ncol = m)

groups = list(1:10, 30:40, 80:100, c(3,5,30,33,39))

y <- sapply(1:n, function(i) {
    x <- X[i, ]
    lg <- -4 + 1.5*sum(x[c(3, 5)]) + 1.01*sum(x[c(30, 33, 39)])
    py <- 1 / (1 + exp(-lg))
    rbinom(1,1,py)
  })

res <- artp.fit(X, y, groups = groups, verbose = T)
res$p.values.group

library(ARTP)

load("artp_test_MR.RData")
fit <- artp.fit(X = X_train, y = y_train, groups = pathways, trunc.point = 5, verbose = T)

pval.gruppe

group <- pathways[[66]]
X <- X_train
y <- y_train


n.obs <- 50
n.cov <- 100
n.groups <- 25
n.permutations <- 5

y <- rbinom(n.obs, 1, prob = .5)

X <- matrix(rbinom(n.obs * n.cov, 1, prob = .5), ncol = n.cov)

g <- g_uniform(n.groups, n.cov)
s <- g()
groups <- DDTA::AssignItemsToSubsets(n.cov, n.groups, s)

X.new <- matrix(rbinom(n.obs * n.cov, 1, prob = .5), ncol = n.cov)

fit <- artp.fit(X, y, groups, n.permutations = 10, verbose = T)
prediction <- artp.predict(fit, X.new, alpha = .5)




n.cov <- ncol(X)

n.permutations.done <- 0

if (verbose) {
  pb <- txtProgressBar(min = 0, max = n.cov*n.permutations, title = "no. of permutations", style = 3)
}

# permutate the output (y) n.permutations time, fit a logistic
# regression and obtain the p-value for each covariate
p.values.permutations <- sapply(1:n.permutations, function(k) {

  # create permutation
  permutation <- sample(y)

  # leave one covariate out each time and store the p-value
  p.values <- sapply(1:n.cov, function(i) {

    # fit the model without the covariate i
    model <- glm(permutation ~ X[, -i], family = binomial(link = "logit"))

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
p.values.permutations <- matrix(p.values.permutations, nrow = n.permutations)

### get the p-values for each group by applying the ARTP function written by Malte

# go over each group
artp.output <- sapply(groups, function(group) {
  # get the p values obtained while permutating for the current group
  data <- t(as.matrix(p.values.permutations[, group]))
  # apply the ARTP to get the p-values for the group
  ARTP(data, J = min(trunc.point, nrow(data)), n.permutations - 1)
})

# create data frame with the results
res <- data.frame(
  id = 1:length(groups),
  group.size = sapply(groups, function(group) length(group[[1]])), # the group sizes
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

