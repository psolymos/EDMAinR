# EDMAinR - Euclidean Distance Matrix Analysis in R

> A coordinate‐free approach for comparing biological shapes using landmark data

## Install

```R
devtools::install_github("psolymos/EDMAinR")
```

## Usage

```R
K <- 3 # number of landmarks
D <- 2 # dimension, 2 or 3
SigmaK <- 0.01*diag(1, K, K)
M <- matrix(c(0,1,0,0,0,1), 3, 2)
M[,1] <- M[,1] - mean(M[,1])
M[,2] <- M[,2] - mean(M[,2])

n <- 1000
Z <- matrix(nrow = n * K, ncol = D)
for (i in 1:n) {
    Z[((i - 1) * K + 1):(i * K), ] <- matrix(rnorm(K * D), nrow = K,
        ncol = D)
}
C <- chol(SigmaK)
X <- matrix(nrow = n * K, ncol = D)
for (i in 1:n) {
    X[((i - 1) * K + 1):(i * K), ] <- crossprod(C, Z[((i - 1) * K + 1):(i *
        K), ]) + M
}

(v <- vartest(X, n, K, D))
```

## References

Lele, S. R., and Richtsmeier, J. T., 1991.
Euclidean distance matrix analysis: A coordinate‐free approach for 
comparing biological shapes using landmark data.
American Journal of Physical Anthropology 86(3):415--27.
DOI: [10.1002/ajpa.1330860307](https://doi.org/10.1002/ajpa.1330860307).

Hu, L., 2007. Euclidean Distance Matrix Analysis of Landmarks Data:
Estimation of Variance. Thesis, Master of Science in Statistics,
Department of Mathematical and Statistical Sciences, 
University of Alberta, Edmonton, Alberta, Canada. Pp. 49.
