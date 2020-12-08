## ARTPredict: Prediction for Binary Outcomes with the ARTP 
--------------------------------------------------------

`ARTPredict` is an `R` package for predicting binary outcomes using the ARTP method. 

### Usage 

Let `X` be the `n x m` binary matrix of observations and let `y` be the binary response vector of length `n`. The `groups` are given in a list where each entry is a vector that represent a single group, e.g., `groups = list(c(1,3,5), c(2,4), c(3,4))` represents three groups. Note that the last group overlaps with the first and second group. 

First, one runs `artp.fit`. The resulting fit can be used for prediction with `artp.predict`, e.g., 

```R
library(ARTPredict)

### Creating a simple toy data set
m = 100
n = 10000

X <- matrix(rbinom(m * n, 1, .05), ncol = m)

# There are 4 groups:
groups = list(1:10, 30:40, 80:100, c(3,5,30,33,39))

# The binary responses are generated given the groups
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
prediction <- artp.predict(fit, X_test, alpha = 0.1)
``` 

See the documentation `?artp.fit` and `?artp.predict` for more info. 

### Acknowledgements

We gratefully acknowledge the financial support from the innovation fund (“Innovationsfonds”) of the Federal Joint Committee in Germany (grant number: 01VSF16020).

### Contact

Louis Dijkstra\
Leibniz Institute for Prevention Research & Epidemiology  
E-mail: dijkstra (at) leibniz-bips.de