# R-package for implementing the $\ell$-test

Following are instructions on installing an `R` package implementing the $\ell$-test-based procedures from Sengupta and Janson (2024). These procedures are alternatives to the traditional $t$-based inference
methods in the Gaussian linear model, having exactly the same guarantees and showing significant power gains when the coefficient vector is not too dense. In fact, $t$-test, called the $\ell$-test, achieves power close to the one-sided $t$-test
for sparse coefficient vectors without any knowledge about the true sign of the alternative coefficient.

## Installation
```R
library(dev_tools)
devtools::install_github("ssouhardya/ell_test_R")
```
