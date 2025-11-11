# R-package for implementing the $\ell$-test

Following are instructions on installing an `R` package implementing the $\ell$-test-based procedures from Sengupta and Janson (2024). These procedures are alternatives to the traditional $t$-based inference
methods in the Gaussian linear model, having exactly the same guarantees and showing significant power gains when the coefficient vector is not too dense. In fact, $t$-test, called the $\ell$-test, achieves power close to the one-sided $t$-test
for sparse coefficient vectors without any knowledge about the true sign of the alternative coefficient.

## Installation
```R
library(dev_tools)
devtools::install_github("ssouhardya/ell_test_R")
```

## Some examples

```R
source('l_testing.R')
source('adjusted_l_testing.R')
set.seed(1)

n = 100
p = 50
s =  5
A = 2.3

X = matrix(rnorm(n*p), nrow = n)
X = apply(X,2,g)
#Supplying an X with normalized columns is recommended

beta = rep(0,p)
rand_ind = sample(1:p, size = s, replace = FALSE)
j = rand_ind[1]	#index to test
beta[rand_ind] = (1-2*rbinom(s,1,0.5))*A
y = as.numeric(X%*%beta + rnorm(n)) #data

pval_l = l.test(y,X,j)	#l-test for H_j:\beta_j = 0

pval_l = l.test(y-2.3*X[,j],X,j) #for testing H_j(2.3):\beta_j = 2.3

pval_l = l.test(y,X,j, lambda_cv = 0.01) #l-test for H_j:\beta_j = 0 with a supplied lambda for cross-validation

pval_l_adjusted = l.test(y,X,j, adjusted = TRUE, lambda = 0.01) #adjusted l-test for H_j:\beta_j = 0 valid conditionally on LASSO selection using penalty 0.01, and the penalty for the test statistic chosen using cross-validation	

gamma_range = seq(from = beta[j]-10, to = beta[j]+10, length.out = 100) #the grid of \gamma values to test on

ci_l = l.ci(y,X,j, gamma_range = gamma_range, coverage = 0.95) #l-CI

ci_l_adjusted = l.ci_adjusted(y,X,j, gamma_range = gamma_range, coverage = 0.95, lambda = 0.01) #post-selection l-CI for \beta_j valid conditionally on LASSO with penalty 0.01 selecting the coefficient

```

## Reference
```
@article{SS-LJ:2024,
  title={The $\ell$-test: leveraging sparsity in the Gaussian linear model for improved inference},
  author={Sengupta, Souhardya and Janson, Lucas},
  journal={arXiv preprint arXiv:2406.18390},
  year={2024}
}
```
