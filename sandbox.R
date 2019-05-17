library(DDTA)
library(readr)
library(ARTPredict)

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

fit <- artp.fit(X, y, groups, n.permutations = 10)
prediction <- artp.predict(fit, X.new, alpha = .5)

