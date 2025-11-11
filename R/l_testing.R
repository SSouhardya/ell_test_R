#' Evaluate the ell-CDF at a point
#'
#' @param x numeric. Point at which the ell-CDF is to be evaluated.
#' @param y numeric vector. Response variable.
#' @param X matrix or data.frame. Full column-rank, `n x p` design matrix; `nrow(X)` must equal `length(y)`.
#' @param ind integer. Index of the coefficient to test; by default an intercept is **not** added to `X`.
#' @param lambda numeric. Lambda for the selection-LASSO: a positive value, or one of
#'   * `-1` to choose `lambda.min` with `cv.glmnet()` called on `(y, X[,-ind])`,
#'   * `-2` to choose `lambda.1se` with `cv.glmnet()` called on `(y, X[,-ind])`,
#'   * `-3` to choose `lambda.min` with `cv.glmnet()` called on `(y, X)`,
#'   * `-4` to choose `lambda.1se` with `cv.glmnet()` called on `(y, X)`. Default `-1`.
#' @param lambda_cv numeric. Lambda used for evaluating the ell-distribution, same coding as `lambda`. Default `-1`.
#' @param glmnet_object list or `NULL`. If a list, then either of length-2 with object returned by `cv.glmnet()` on `(y, X)` and `(y, X[,-ind])` as first and second components,
#'   or length `d+1` list where the first component is the object returned by `cv.glmnet()` called on `(y, X)` and the `i`th is when called on `(y, X[,-i])`. Default `NULL`.
#' @param glmnet_object_type integer. Set this to `1` if `glmnet_object` is in the 2-element list form, `1` for the `d+1`-element form. If `glmnet_object` is `NULL`, this is forced to `1`. Default `1`.
#' @param adjusted logical. If `TRUE`, return the CDF conditional on LASSO selecting variable `ind`. Default `FALSE`.
#' @param return_both logical. If `TRUE`, return both unconditional and conditional CDFs. Default `FALSE`.
#' @param return_dropprob logical. If `TRUE`, also return the conditional drop probability (LASSO estimate at `ind` equals 0). Default `FALSE`.
#' @return If `return_dropprob = FALSE`, returns a numeric vector `v`. If `TRUE`, returns a list with `v` as the first and the conditional drop probability as the second component.
#'   When `return_both = FALSE`, `v` is a single value (unconditional value if `adjusted = FALSE`, conditional value if `adjusted = TRUE`).
#'   When `return_both = TRUE`, `v` has two entries: `(unconditional value of CDF, conditional value of CDF)`. If LASSO does not select variable `ind`, the conditional entry is `-1`.
#' @export
l.cdf_glmnet<-function(x,y,X,ind,lambda,lambda_cv, glmnet_object=NULL, glmnet_object_type = 1, adjusted = FALSE, return_both = FALSE, return_dropprob = FALSE){
    if(return_both){
        adjusted = TRUE
    }

    if(return_dropprob & (adjusted == FALSE) ){
    	adjusted = TRUE
    	warning('adjutsed set to TRUE since return_dropprob is TRUE')
    }
    
    n = nrow(X)
    p = ncol(X)
    
    y.std = std(y)
    Z = cbind(rep(1,n),X[,-ind])
    proj = Z%*%solve(t(Z)%*%Z)%*%t(Z)
    y.std.hat = proj %*%y.std
    
    sigma.hat.std = sqrt(sum((y.std - y.std.hat)^2))
    
    second_denom_term = sqrt(sum(((diag(n) - proj)%*%X[,ind])^2))
    
    if(is.null(glmnet_object)){
        glmnet_object_type = 1
        glmnet_object = list()
        glmnet_object[[1]] = glmnet::cv.glmnet(X,y.std, standardize=FALSE)
        glmnet_object[[2]] = glmnet::cv.glmnet(X[,-ind],y.std, standardize=FALSE)
    }
    
    ind_to_look = (glmnet_object_type-1)*(1+ind) + (2-glmnet_object_type)*2

    if(lambda == -1){
        lambda = glmnet_object[[ind_to_look]]$lambda.min
    } else if(lambda == -2){
        lambda = glmnet_object[[ind_to_look]]$lambda.1se
    } else if(lambda == -3){
        lambda = glmnet_object[[1]]$lambda.1se
    } else if(lambda == -4){
        lambda = glmnet_object[[1]]$lambda.1se
    }
    
    if(lambda_cv == -1){
        lambda_cv = glmnet_object[[ind_to_look]]$lambda.min
    } else if(lambda_cv == -2){
        lambda_cv = glmnet_object[[ind_to_look]]$lambda.1se
    } else if(lambda_cv == -3){
        lambda_cv = glmnet_object[[1]]$lambda.1se
    } else if(lambda_cv == -4){
        lambda_cv = glmnet_object[[1]]$lambda.1se
    }
    
    
    sgn<-function(x){
        if(x == 0){
            return(1)
        }
        return(sign(x))
    }
    beta.x = coef(glmnet::glmnet(X[,-ind], y.std - x*X[,ind], standardize = FALSE), s = lambda_cv, exact = TRUE, x = X[,-ind], y = y.std- x*X[,ind])
    v.x = ( -sum(X[,ind]*(y.std.hat - x*X[,ind] - Z%*%beta.x)) +n*lambda_cv*sgn(x) )/(sigma.hat.std*second_denom_term)
    
    beta.0 = coef(glmnet_object[[ind_to_look]]$glmnet.fit, s = lambda, exact = TRUE, x = X[,-ind], y = y.std)
    v1 = ( -sum(X[,ind]*(y.std.hat - Z%*%beta.0)) -n*lambda )/(sigma.hat.std*second_denom_term)
    v2 = ( -sum(X[,ind]*(y.std.hat - Z%*%beta.0)) +n*lambda )/(sigma.hat.std*second_denom_term)

    denom = 1 - (qhaar(v2,n-ncol(X), lower.tail = TRUE) - qhaar(v1,n-ncol(X), lower.tail = TRUE) )	# i am calulating this outside as this give sthe dropping probability

    uncond_prob = qhaar(v.x,n-ncol(X),lower.tail = TRUE)
    
    if(!adjusted){
        return(uncond_prob)
    }
    
    beta.hat = coef(glmnet_object[[1]]$glmnet.fit, s=lambda, exact = TRUE, x = X, y = y.std)
    

    if(beta.hat[1+ind] == 0){
        cond_prob = -1
    } else{
        numer_1 = qhaar(v.x, n-ncol(X), lower.tail = TRUE)
        numer_2 = 0
        if(v.x>v1){
            numer_2 = qhaar(min(c(v.x,v2)), n- ncol(X), lower.tail = TRUE) - qhaar(v1, n- ncol(X), lower.tail = TRUE)
        }
        cond_prob = (numer_1 - numer_2)/denom
    }

    
    if(!return_both){
        toret = cond_prob
    } else{
        toret = c(uncond_prob, cond_prob)
    }
    if(return_dropprob){
    	toret = list(toret, 1-denom)
    }
    
    return(toret)
}


#' Perform the ell-test
#'
#' Computes the ell-test p-value for a specified coefficient in a linear model.
#'
#' @param y Numeric vector. Response variable.
#' @param X Matrix or data frame. Full column-rank `n x p` design matrix; `nrow(X)` must match `length(y)`. 
#'   The intercept is included.
#' @param ind Integer. Index of the coefficient to test.
#' @param lambda Numeric. Regularization parameter for the selection LASSO. 
#'   Either a positive value or one of:
#'   * `-1` to choose `lambda.min` with `cv.glmnet()` called on `(y, X[,-ind])`,
#'   * `-2` to choose `lambda.1se` with `cv.glmnet()` called on `(y, X[,-ind])`. Default `-1`.
#' @param lambda_cv Numeric. Regularization parameter for evaluating the ell-distribution. 
#'   Same conventions as `lambda`. Default is `-1`.
#' @param glmnet_object List or `NULL`. If a list, then either:
#'   * Length-2 list with object returned by `cv.glmnet()` on `(y, X)` and `(y, X[,-ind])` as first and second components, or
#'   * Length `d+1` list where the first component is the object returned by `cv.glmnet()` called on `(y, X)` and the `i`th is when called on `(y, X[,-i])`.  
#'   Default is `NULL`, in which case the function constructs the first type.
#' @param glmnet_object_type Integer. Either `1` (2-element list form) or `2` (`d+1`-element list form). 
#'   Overridden to `1` if `glmnet_object = NULL`. Default is `1`.
#' @param adjusted Logical. If `TRUE`, returns a p-value adjusted for post-LASSO selection of the coefficient at `ind`. 
#'   Default is `FALSE`.
#' @param smoothed Logical. If `TRUE`, smooths the unconditional p-value when the LASSO estimate at `ind` is zero. 
#'   Default is `TRUE`.
#' @param return_both Logical. If `TRUE`, returns both the unconditional and adjusted ell-test p-values. 
#'   Default is `FALSE`.
#'
#' @return If `return_both = FALSE`, a single numeric value (the ell-test p-value) is returned depending on whether `adjusted` is `TRUE` or `FALSE`. 
#'   If `return_both = TRUE`, a numeric vector of length 2: `(unconditional p-value, adjusted p-value)`. 
#'   If the LASSO does not select variable `ind`, the adjusted p-value is `-1`.
#'
#' @examples
#' set.seed(1)
#'
#' n <- 100
#' p <- 50
#' s <- 5
#' A <- 2.3
#'
#' X <- matrix(rnorm(n * p), nrow = n)
#' X <- apply(X, 2, g)  # normalize columns; recommended
#'
#' beta <- rep(0, p)
#' rand_ind <- sample(1:p, size = s, replace = FALSE)
#' j <- rand_ind[1]     # index to test
#' beta[rand_ind] <- (1 - 2 * rbinom(s, 1, 0.5)) * A
#'
#' y <- as.numeric(X %*% beta + rnorm(n))
#'
#' # l-test for H0: beta_j = 0
#' pval_l <- l.test(y, X, j)
#'
#' # l-test for H0: beta_j = 2.3
#' pval_l <- l.test(y - 2.3 * X[, j], X, j)
#'
#' # l-test for H0: beta_j = 0 with supplied lambda for CV
#' pval_l <- l.test(y, X, j, lambda_cv = 0.01)
#' pval_l_adjusted = l.test(y,X,j, adjusted = TRUE, lambda = 0.01) 
#' #adjusted l-test for H_j:\beta_j = 0 valid conditionally on LASSO selection using penalty 0.01, 
#' #and the penalty for the test statistic chosen using cross-validation
#' @export
#' 
l.test<-function(y,X,ind, lambda=-1, lambda_cv=-1, glmnet_object=NULL, glmnet_object_type = 1, adjusted = FALSE, smoothed = TRUE, return_both = FALSE){
	if(return_both){
		adjusted = FALSE
	}
	n = nrow(X)
  p = ncol(X)
  if(p<=2){
  	stop('The dimension needs to be at least 3')
  }
    
  y.std = std(y)


  Z = cbind(rep(1,n),X[,-ind])
	proj = Z%*%solve(t(Z)%*%Z)%*%t(Z)
	y.std.hat = proj %*%y.std
	sigma.hat.std = sqrt(sum((y.std - y.std.hat)^2))
	V = qr.Q(qr(diag(n)-proj))[,1:(n-ncol(Z))]
	u = rnorm(n-ncol(Z))
	u = u/sqrt(sum(u^2))
	y_temp = as.vector(y.std.hat + sigma.hat.std*V%*%u)
	y_temp.std = std(y_temp)

    if(is.null(glmnet_object)){
        glmnet_object_type = 1
        glmnet_object = list()
        glmnet_object[[1]] = glmnet::cv.glmnet(X,y.std, standardize=FALSE)
        glmnet_object[[2]] = glmnet::cv.glmnet(X[,-ind],y_temp.std, standardize=FALSE)
    }
    
    ind_to_look = (glmnet_object_type-1)*(1+ind) + (2-glmnet_object_type)*2

    if(lambda == -1){
        lambda = glmnet_object[[ind_to_look]]$lambda.min
    } else if(lambda == -2){
        lambda = glmnet_object[[ind_to_look]]$lambda.1se
    } else if(lambda == -3){
        lambda = glmnet_object[[1]]$lambda.1se
    } else if(lambda == -4){
        lambda = glmnet_object[[1]]$lambda.1se
    }
    
    if(lambda_cv == -1){
        lambda_cv = glmnet_object[[ind_to_look]]$lambda.min
    } else if(lambda_cv == -2){
        lambda_cv = glmnet_object[[ind_to_look]]$lambda.1se
    } else if(lambda_cv == -3){
        lambda_cv = glmnet_object[[1]]$lambda.1se
    } else if(lambda_cv == -4){
        lambda_cv = glmnet_object[[1]]$lambda.1se
    }

    x = abs(as.numeric(coef(glmnet_object[[1]]$glmnet.fit, s = lambda_cv, exact = TRUE, x = X, y = y.std)[1+ind]) )

    if(x == 0){
    	uncond_pval = 1
    	cond_pval = 1

    	if(!adjusted){
    		if(smoothed){
    			beta.hat.null = as.vector(coef(glmnet_object[[ind_to_look]]$glmnet.fit, s = lambda_cv, exact = TRUE, y = y_temp.std, x = X[,-ind]))
    			second_denom_term = sqrt(sum((X[,ind] - proj%*%X[,ind])^2))
    			mid = -as.numeric(X[,ind]%*%(y.std.hat - Z%*%beta.hat.null))/(sigma.hat.std * second_denom_term)
    			v1 = mid - n*lambda_cv/(sigma.hat.std * second_denom_term)
    			v2 = mid + n*lambda_cv/(sigma.hat.std * second_denom_term)
    			V = qr.Q(qr(diag(n)-proj))[,1:(n-ncol(Z))]
    			u1 = sum(X[,ind]*(y.std - y.std.hat))/(sigma.hat.std*second_denom_term)
    			dist_from_mid = abs(mid - u1)
    			w1 = mid - dist_from_mid
    			w2 = mid + dist_from_mid
    			left_prob = qhaar(w1, n-ncol(Z), lower.tail = TRUE)
    			right_prob = 1 - qhaar(w2, n-ncol(Z), lower.tail = TRUE)
    			uncond_pval = left_prob + right_prob
    		}
    		if(!return_both){
    			return(uncond_pval)
    		}
    	}
    	beta.selection = abs(as.numeric(coef(glmnet_object[[1]]$glmnet.fit, s = lambda, exact = TRUE, x = X, y = y.std)[1+ind]) )
    	if(beta.selection == 0){
    		cond_pval = -1
    	}
    	if(!return_both){
    		return(cond_pval)
    	} else{
    		return(c(uncond_pval, cond_pval))
    	}
    } else{
    	pval.right = 1-l.cdf_glmnet(x,y,X,ind,lambda = lambda,lambda_cv = lambda_cv, glmnet_object = glmnet_object, glmnet_object_type = glmnet_object_type, adjusted = adjusted, return_both = return_both)
    	pval.left = l.cdf_glmnet(-x,y,X,ind,lambda = lambda,lambda_cv = lambda_cv, glmnet_object = glmnet_object, glmnet_object_type = glmnet_object_type, adjusted = adjusted, return_both = return_both)
    	pval = pval.left + pval.right
    	if(return_both){
    		if(pval.left[2] == -1){
    			pval[2] = -1
    		}
    	}
    	return(pval)
    }
}


#' computing the (un-adjusted) ell-confidence interval
#' 
#' @param y Numeric vector. Response variable.
#' @param X Matrix or data frame. Full column-rank `n x p` design matrix; `nrow(X)` must match `length(y)`. 
#'   The intercept is included.
#' @param ind Integer. Index of the coefficient to test.
#' @param gamma_range Vector. A grid of values on which the ell-test is to be inverted.
#' @param lambda_cv Numeric. Regularization parameter for evaluating the ell-distribution. 
#'   Either a positive value or one of:
#'   * `-1` to choose `lambda.min` with `cv.glmnet()` called on `(y, X[,-ind])`,
#'   * `-2` to choose `lambda.1se` with `cv.glmnet()` called on `(y, X[,-ind])`. Default `-1`.
#' @param coverage Numeric. The coverage of the confidence interval. Default is 0.95.
#' @return returns a two-dimensional vector specifying the upper and lower limits of the ell-test confidence interval (not adjusted for any LASSO selection).
#' @examples
#' set.seed(1)
#'
#' n <- 100
#' p <- 50
#' s <- 5
#' A <- 2.3
#'
#' X <- matrix(rnorm(n * p), nrow = n)
#' X <- apply(X, 2, g)  # normalize columns; recommended
#'
#' beta <- rep(0, p)
#' rand_ind <- sample(1:p, size = s, replace = FALSE)
#' j <- rand_ind[1]     # index to test
#' beta[rand_ind] <- (1 - 2 * rbinom(s, 1, 0.5)) * A
#'
#' y <- as.numeric(X %*% beta + rnorm(n))
#' 
#' gamma_range = seq(from = beta[j]-10, to = beta[j]+10, length.out = 100) 
#' #the grid of \gamma values to test on
#'ci_l = l.ci(y,X,j, gamma_range = gamma_range, coverage = 0.95) #l-CI
#' @export
#' 
l.ci<-function(y,X,ind, gamma_range, lambda_cv=-1, coverage = 0.95,  smoothed = TRUE, outer_approx = FALSE, outer_grid.length = 10){
	gamma_range = sort(gamma_range)
	g.length = length(gamma_range)
	pvals = vector(length = g.length)

	n = nrow(X)
	p = ncol(X)

	if(p<=2){
    	stop('The dimension needs to be at least 3')
    }

	lambda = lambda_cv
	for(i in 1:g.length){
		y_test = y - gamma_range[i]*X[,ind]
		pvals[i] = l.test(y_test,X, ind, lambda, lambda_cv, adjusted = FALSE, smoothed = smoothed)
	}
	inds = which(pvals > 1-coverage)
	if(length(inds) == 0){
	    warning('None of the grid elements selected. Try a different grid.')
	    return(vector(length = 0))
	}
	if(outer_approx){
		if(min(inds)!=1){
			inds1 = min(inds)-1
			left_additional = seq(from = gamma_range[inds1], to = gamma_range[min(inds)], length.out = outer_grid.length)
			pvals.left = vector(length = outer_grid.length)
			for(j in 1:outer_grid.length){
				y_test = y - left_additional[j]*X[,ind]
				pvals.left[j] = l.test(y_test,X, ind, lambda, lambda_cv, adjusted = FALSE, smoothed = smoothed)
			}
			inds.left = which(pvals.left > 1-coverage)
			inds.left = unique(c(inds.left, length(left_additional)))
			if(min(inds.left)!=1){
				inds.left = c(min(inds.left)-1, inds.left)
			}
			ci.left = left_additional[inds.left[1]]
		} else{
			ci.left = gamma_range[inds[1]]
			warning('CI hit the lower limit. Try increasing the range?')
		}
		if(max(inds)!=g.length){
			inds2 = max(inds)+1
			right_additional = seq(from = gamma_range[max(inds)], to = gamma_range[inds2], length.out = outer_grid.length)
			pvals.right = vector(length = outer_grid.length)
			for(j in 1:outer_grid.length){
				y_test = y - right_additional[j]*X[,ind]
				pvals.right[j] = l.test(y_test,X, ind, lambda, lambda_cv, adjusted = FALSE, smoothed = smoothed)
			}
			inds.right = which(pvals.right > 1-coverage)
			inds.right = unique(c(1,inds.right))
			if(max(inds.right)!= outer_grid.length){
				inds.right = c(inds.right, max(inds.right)+1)
			}
			ci.right = right_additional[inds.right[length(inds.right)]]
		} else{
			ci.right = gamma_range[inds[length(inds)]]
			warning('CI hit the upper limit. Try increasing the range?')
		}
	} else{
		temp = range( gamma_range[inds] )
		 ci.left = temp[1]
		 ci.right = temp[2]
	}
	return(c(ci.left, ci.right))
}
