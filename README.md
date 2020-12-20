## ARTPredict: Prediction for Binary Outcomes with the ARTP
--------------------------------------------------------

`ARTPredict` is an `R` package for predicting binary outcomes using the ARTP method.

### Usage

Let `X` be the `n x m` binary matrix of observations and let `y` be the binary response vector of length `n`. The `groups` are given in a list where each entry is a vector that represent a single group, e.g., `groups = list(c(1,3,5), c(2,4), c(3,4))` represents three groups. Note that the last group overlaps with the first and second group.

First, one runs `artp.fit`. The resulting fit can be used for prediction with `artp.predict`, e.g.,

```R
devtools::load_all()
# or
library(ARTPredict)

### Set seed
set.seed(43)

### Create a simple dataset
m <- 100
n <- 10000
X <- matrix(rbinom(m * n, 1, .05), ncol = m)

# Create 4 groups
groups <- list(1:10, 30:40, 80:100, c(3,5,30,33,39))

# Generate binary response given the groups
y <- sapply(1:n, function(i) {
    x <- X[i, ]
    lg <- -4 + 4*sum(x[c(3, 5)]) + 4*sum(x[c(30, 33, 39)])
    py <- 1 / (1 + exp(-lg))
    rbinom(1,1,py)
  })

### Divide the simulated data in a train and test dataset
train_ind <- sort(sample(seq_len(n), size = floor(n/2)))

X_train = X[train_ind, ]
X_test = X[-train_ind, ]

y_train = y[train_ind]
y_test = y[-train_ind]

### Fit the data
fit <- artp.fit(X_train, y_train, groups = groups, verbose = TRUE)

### Predict
prediction <- artp.predict(fit, X_test, alpha = 0.2)
table(prediction$y.hat, y_test)

### Handle adjustment variables
adjust_vars <- c(22, 23, 50)
fit2 <- artp.fit(X_train, y_train, adjust_vars = adjust_vars, groups = groups, verbose = TRUE)
prediction2 <- artp.predict(fit2, X_test, alpha = 0.2)
table(prediction2$y.hat, y_test)

### Parallel, system.time and memory profiling
fit <- artp.fit(X_train, y_train, groups = groups, verbose = FALSE, parallel = TRUE)
prediction <- artp.predict(fit, X_test, alpha = 0.2)
table(prediction$y.hat, y_test)

### Compare with glm
set.seed(NULL)
set.seed(235478965)
m = 100
n = 2000
X <- matrix(rbinom(m * n, 1, .05), ncol = m)
X_train <- X[1:1000,]
X_test <- X[1001:2000,]

groups = list(1:10, 30:40, 80:100, c(3,5,30,33,39))

y <- sapply(1:n, function(i) {
  x <- X[i, ]
  lg <- -4 + 4 * sum(x[c(3, 5)]) + 4 * sum(x[c(30, 33, 39)])
  py <- 1 / (1 + exp(-lg))
  rbinom(1,1,py)
  })

y_train <- y[1:1000]
y_test  <- y[1001:2000]

res <- artp.fit(X_train, y_train, groups = groups, verbose = T, trunc.point = 3)

pred <- artp.predict(res, X_test, alpha = .2)
table(pred$y.hat, y_test)

glm.out <- glm(y_train ~ X_train, family = binomial(link = "logit"))
probabilities <- predict(glm.out, as.data.frame(X_test), type = "response")
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)

table(predicted.classes, y_test)
```

See the documentation `?artp.fit` and `?artp.predict` for more info.

### Acknowledgements

We gratefully acknowledge the financial support from the innovation fund (“Innovationsfonds”) of the Federal Joint Committee in Germany (grant number: 01VSF16020).

### Contact

Louis Dijkstra\
Leibniz Institute for Prevention Research & Epidemiology - BIPS GmbH
E-mail: dijkstra (at) leibniz-bips.de
