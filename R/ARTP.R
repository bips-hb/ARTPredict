# ------------------------------------------------------------------------------
# -----------   ARTP_gene_simple   ---------------------------------------------
#' Main function, based on Yu et al.
#' Initial programming: Malte Braitmaier
#'
#' @param daten A matrix, rows contain variables,
#'              columns contain the original p-value followed
#               by p-values from permuted data
#' @param J Truncation point
#' @param n_perm Number of permutations
#'
#' @return A data frame p.values.groups with the p value for each
#'         group and truncation
#'
#' @export
ARTP <- function(daten, J, n_perm) {

  if (J > nrow(daten)) {
    stop("J must NOT be greater than nrow(daten)")
  }

  min1 <- matrix(NA, nrow = J, ncol = n_perm + 1)
  rownames(min1) <- 1:J
  colnames(min1) <- c("p_emp", paste("perm", 1:n_perm, sep = ""))

  ### sort p-values in each row
  daten.sort <- apply(daten, 2, sort)

  if (is.matrix(daten.sort)) {
    ### W_j^(b): RTP statistic
    for (i in 1:J) {
      # min1[i, 1] <- prod(daten.sort[1:i,1])
      # for (j in 1:n_perm)
      # {
      #   min1[i, j + 1] <- prod(daten.sort[1:i, j + 1])
      # }
      for (j in 0:n_perm) {
        min1[i, j + 1] <- prod(daten.sort[1:i, j + 1])
      }
    }
  } else {
    min1[1, ] <- daten.sort
  }
  order.prods <- matrix(NA, nrow = J, ncol = n_perm + 1)
  rownames(order.prods) <- 1:J

  for (i in 1:J) {
    order.prods[i, ] <- rank(min1[i, 1:(n_perm + 1)], na.last = "keep")
  }

  ### s_j^(b): p-value for W_j^(b)
  s_j <- order.prods / (n_perm + 1)

  ### which J has the minimal rank for each data set?
  minima <- rep(NA, n_perm + 1)

  for (k in 1:ncol(order.prods)) {
    minima[k] <- min(which(rank(order.prods[, k], na.last = TRUE) ==
      min(rank(order.prods[, k], na.last = TRUE))))
  }

  ### MinP^(b) = min_j s_j^(b)
  MinP <- data.frame(matrix(NA, ncol = n_perm + 2, nrow = 1))
  names(MinP) <- c("j", paste("MinP", 0:n_perm))

  for (b in 2:ncol(MinP)) {
    MinP[1, b] <- s_j[minima[b - 1], b - 1]
  }
  MinP[1, 1] <- minima[1]

  out <- data.frame(matrix(NA, nrow = 1, ncol = 2))
  out[1, 1] <- MinP[1, 1]
  out[1, 2] <- rank(as.numeric(MinP[1, 2:ncol(MinP)]))[1] / (n_perm + 1)
  names(out) <- c("j", "ARTP_p")
  out
}
