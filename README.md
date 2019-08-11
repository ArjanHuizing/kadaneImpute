# kadaneImpute
Based on the 2001 paper 'Statistical matching: A paradigm for assesing the uncertainty in the procedure' by Moriarity and Scheuren. kadaneImpute imputes missing values through linear regression, with the added ability to specify a correlation between outcomes. This can be useful when outcomes are never jointly observed, for example with potential outcomes.

The procedure can be run inside the R package `mice`, and should be specified for a block of atleast 2 continous variables with missings. Additionally, the libraries `matrixcalc`, `Rcpp` and `RcppArmadillo` are required.

# Example R code
```
data <- MASS::mvrnorm(n = 100, mu = c(0, 1), Sigma = matrix(c(1, 0.8, 0.8, 1), ncol = 2))
data[1:50, 1] <- NA
data[51:100, 2] <- NA

imp <- mice(as.data.frame(data), blocks = list(c("V1", "V2")), 
                                 method = "kadane",
                                 kadane.corr = 0.8, 
                                 kadane.match = TRUE)
plot(complete(imp))
```

# To-do
- ~~Use maximum likelihood to estimate the covariance matrix to prevent serious errors when estimating more complex data.~~
- ~~Use robust residual variance estimation~~ (More testing required)
- ~~Make code more efficient~~
- More options for correlation specification
- Clean up and document

# Changelog
- Two major bugs were squashed. The kadane method now works within mice. Proof-of-concept simulation results can be found in the Simulation folder.
- The method described by Moriarity & Scheuren (2001) to estimate the residual variance can lead to negative variance estimates. The RIEPS algorithm is used instead (Kiesl & Rassler, 2009).
- C++ based functions have sped up kadaneImpute significantly. A simple benchmark showed a ~~6x increased speed~~ **29x increased speed.**
- Now estimates covariance matrix with expectation-maximization to obtain maximum likelihood estimates when data is missing.
