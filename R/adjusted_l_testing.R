#' @export
beta_x<-function(x,y,X,ind,lamb, solver = NA){
	n = nrow(X)
	Z = cbind(rep(1,n),X[,-ind])
	beta = CVXR::Variable(ncol(X))
	obj = sum((y-x*X[,ind]-Z%*%beta)^2)/(2*n) + lamb*CVXR::p_norm(beta[-1],1)
	prob0 = CVXR::Problem(CVXR::Minimize(obj))
	beta.x = round(CVXR::solve(prob0,feastol=1e-10,reltol=1e-10,abstol=1e-10, solver = solver )$getValue(beta), digits = 9)
	return(beta.x)
}

#' @export
beta_full<-function(y,X,ind,lamb,solver = NA){
	n = nrow(X)
	Z = cbind(rep(1,n),X)
	beta = CVXR::Variable(ncol(X)+1)
	obj = sum((y-Z%*%beta)^2)/(2*n) + lamb*CVXR::p_norm(beta[-1],1)
	prob0 = CVXR::Problem(CVXR::Minimize(obj))
	beta.full = round(CVXR::solve(prob0,feastol=1e-10,reltol=1e-10,abstol=1e-10, solver = solver )$getValue(beta), digits = 9)
	return(beta.full)
}


#' @export
l.cdf_adjusted<-function(x,y,X,ind,gamma, lambda_cv , lambda , solver = NA, tail = 'left', smoothed = FALSE){
	n = nrow(X)
	p = ncol(X)

	Z = cbind(rep(1,n),X[,-ind])
	proj = Z%*%solve(t(Z)%*%Z)%*%t(Z)
	y.hat = proj %*%y
	sigma.hat = sqrt(sum(((diag(n) - proj)%*%(y - X[,ind]*gamma))^2))
	second_denom_term = sqrt(sum(((diag(n) - proj)%*%X[,ind])^2))

    sgn<-function(x){
        if(x == 0){
            return(1)
        }
        return(sign(x))
    }

    if(smoothed & (x==0)){
    	beta.x = beta_x(0,y - gamma*X[,ind],X,ind,lambda_cv, solver = solver)
    	mid = -sum(X[,ind]*(y.hat - gamma*(proj)%*%X[,ind] - Z%*%beta.x))/(sigma.hat * second_denom_term)
    	u1 = sum(X[,ind]*(y - y.hat - gamma*(diag(n) - proj)%*%X[,ind] ))/(sigma.hat*second_denom_term)
    	dist_from_mid = abs(mid - u1)
    	if(tail == 'left'){
    		v.x = mid - dist_from_mid
    	} else{
    		v.x = mid + dist_from_mid
    	}
    } else{
		beta.x = beta_x(x,y - gamma*X[,ind],X,ind,lambda_cv, solver = solver)
		v.x = ( -sum(X[,ind]*(y.hat - gamma*(proj)%*%X[,ind]- x*X[,ind] - Z%*%beta.x)) + n*lambda_cv*sgn(x) )/(sigma.hat*second_denom_term)
	}
	beta.0 = beta_x(0,y,X,ind,lambda, solver = solver)

	v1 = ( -sum(X[,ind]*(y.hat + gamma*(diag(n) - proj)%*%X[,ind] - Z%*%beta.0)) - n*lambda )/(sigma.hat*second_denom_term)
    v2 = ( -sum(X[,ind]*(y.hat + gamma*(diag(n) - proj)%*%X[,ind] - Z%*%beta.0)) + n*lambda )/(sigma.hat*second_denom_term)

    denom = 1 - (qhaar(v2,n-ncol(X), lower.tail = TRUE) - qhaar(v1,n-ncol(X), lower.tail = TRUE))
    numer_1 = qhaar(v.x, n-ncol(X), lower.tail = TRUE)
    numer_2 = 0
    if(v.x>v1){
    	numer_2 = qhaar(min(c(v.x,v2)), n- ncol(X), lower.tail = TRUE) - qhaar(v1, n- ncol(X), lower.tail = TRUE)
    }
    cond_prob = (numer_1 - numer_2)/denom

    if(tail == 'right'){
    	cond_prob = 1-cond_prob
    }
    return(cond_prob)
}


#' computing the post-LASSO-selection-adjusted ell-confidence interval
#' 
#' @param y Numeric vector. Response variable.
#' @param X Matrix or data frame. Full column-rank `n x p` design matrix; `nrow(X)` must match `length(y)`. 
#'   The intercept is included.
#' @param ind Integer. Index of the coefficient to test.
#' @param lambda Numeric. Regularization parameter for the selection LASSO. 
#' @param lambda_cv Numeric. Regularization parameter for evaluating the ell-distribution. 
#'   Either a positive value or one of:
#'   * `-1` to choose `lambda.min` with `cv.glmnet()` called on `(y, X[,-ind])`,
#'   * `-2` to choose `lambda.1se` with `cv.glmnet()` called on `(y, X[,-ind])`. Default `-1`.
#' @param coverage the coverage of the confidence interval. Default is 0.95.
#' @return returns a two-dimensional vector specifying the upper and lower limits of the ell-test confidence interval (adjusted for LASSO selection). If LASSO does not select variable `ind`, returns `NA`.
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
#' ci_l_adjusted = l.ci_adjusted(y,X,j, gamma_range = gamma_range, coverage = 0.95, lambda = 0.01) 
#' #post-selection l-CI for \beta_j valid conditionally on LASSO with penalty 0.01 
#' #selecting the coefficient
#' @export
l.ci_adjusted<-function(y,X,ind, gamma_range, lambda,  lambda_cv=-1, coverage = 0.95, display = FALSE, smoothed = FALSE, outer_approx = FALSE, outer_grid.length = 10){
	g.length = length(gamma_range)
	lambda_cv_flag = lambda_cv

	cond_pval = vector(length = length(gamma_range))
	right.pval = cond_pval
	left.pval = cond_pval
	lambda_vec = left.pval
	lambda_cv_vec = lambda_vec

	n = nrow(X)
	p = ncol(X)

	Z = cbind(rep(1,n),X[,-ind])
	proj = Z%*%solve(t(Z)%*%Z)%*%t(Z)
	y.hat = proj %*%y
	V = qr.Q(qr(diag(n)-proj))[,1:(n-ncol(Z))]
	second_denom_term = sqrt(sum(((diag(n) - proj)%*%X[,ind])^2))

	beta = CVXR::Variable(ncol(X)+1)
	obj = sum(( y - cbind(rep(1,n), X)%*%beta)^2)/(2*n) + lambda*CVXR::p_norm(beta[-1],1)
	prob = CVXR::Problem(CVXR::Minimize(obj))
	beta.hat = round(CVXR::solve(prob,feastol=1e-10,reltol=1e-10,abstol=1e-10,)$getValue(beta), digits = 9)

	if(beta.hat[1+ind] == 0){
		return(NA)	#no confidence interval if there is no selection
	}

	for(i in 1:length(gamma_range)){
		tic1 = Sys.time()
		gamma = gamma_range[i]
		if(lambda_cv_flag <0){
			sigma.hat.gamma = sqrt(sum(((diag(n) - proj)%*%(y - X[,ind]*gamma))^2))
			u = rnorm(n-ncol(Z))
			u = u/sqrt(sum(u^2))
			y_temp = as.vector(y.hat + gamma*(diag(n)-proj)%*%X[,ind] + sigma.hat.gamma*V%*%u)
			glmnet_object_cv = glmnet::cv.glmnet(X[,-ind], y_temp-gamma*X[,ind], standardize=FALSE)
			if(lambda_cv_flag == -1){
				lambda_cv = glmnet_object_cv$lambda.min
			} else if(lambda_cv_flag == -2){
				lambda_cv = glmnet_object_cv$lambda.1se
			}
		}
		lambda_vec[i] = lambda
		lambda_cv_vec[i] = lambda_cv

		beta = CVXR::Variable(ncol(X)+1)
		obj = sum(( (y - gamma*X[,ind]) -cbind(rep(1,n), X)%*%beta)^2)/(2*n) + lambda_cv*CVXR::p_norm(beta[-1],1)
		prob = CVXR::Problem(CVXR::Minimize(obj))
		beta.hat.stat = round(CVXR::solve(prob,feastol=1e-10,reltol=1e-10,abstol=1e-10)$getValue(beta), digits = 9)



		stat = abs(beta.hat.stat[1+ind])
		#print(stat)
		if(beta.hat[1+ind] == 0){
			cond_pval[i] = -1
		} else{
			right.pval[i] = l.cdf_adjusted(stat, y, X, ind, gamma, lambda_cv, lambda, smoothed = smoothed, tail = 'right')
			left.pval[i] = l.cdf_adjusted(-stat, y, X, ind, gamma, lambda_cv, lambda, smoothed = smoothed, tail = 'left')
			cond_pval[i] = right.pval[i] + left.pval[i]
		}
		toc1 = Sys.time()
		if(display){
			print(paste('For gamma value ',i))
			print(toc1 - tic1)
		}
	}

	inds = which(cond_pval>1-coverage)
	if(length(inds) == 0){
	    warning('None of the elements in the grid selected. Try increase the resolution?')
	    return(vector(length = 0))
	}

	if(outer_approx){
		if(min(inds)!=1){
			inds1 = min(inds)-1
			left_additional = seq(from = gamma_range[inds1], to = gamma_range[min(inds)], length.out = outer_grid.length)

			pvals.left = vector(length = outer_grid.length)
			for(j in 1:outer_grid.length){
				gamma = left_additional[j]
				if(lambda_cv_flag <0){
					sigma.hat.gamma = sqrt(sum(((diag(n) - proj)%*%(y - X[,ind]*gamma))^2))
					u = rnorm(n-ncol(Z))
					u = u/sqrt(sum(u^2))
					y_temp = as.vector(y.hat + gamma*(diag(n)-proj)%*%X[,ind] + sigma.hat.gamma*V%*%u)
					glmnet_object_cv = glmnet::cv.glmnet(X[,-ind], y_temp-gamma*X[,ind], standardize=FALSE)
					if(lambda_cv_flag == -1){
						lambda_cv = glmnet_object_cv$lambda.min
					} else if(lambda_cv_flag == -2){
						lambda_cv = glmnet_object_cv$lambda.1se
					}
				}			
				beta = CVXR::Variable(ncol(X)+1)
				obj = sum(( (y - gamma*X[,ind]) -cbind(rep(1,n), X)%*%beta)^2)/(2*n) + lambda_cv*CVXR::p_norm(beta[-1],1)
				prob = CVXR::Problem(CVXR::Minimize(obj))
				beta.hat.stat = round(CVXR::solve(prob,feastol=1e-10,reltol=1e-10,abstol=1e-10)$getValue(beta), digits = 9)
				stat = abs(beta.hat.stat[1+ind])

				right_tail = l.cdf_adjusted(stat, y, X, ind, gamma, lambda_cv, lambda, smoothed = smoothed, tail = 'right')
				left_tail = l.cdf_adjusted(-stat, y, X, ind, gamma, lambda_cv, lambda, smoothed = smoothed, tail = 'left')
				pvals.left[j] = left_tail + right_tail
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
				gamma = right_additional[j]
				if(lambda_cv_flag <0){
					sigma.hat.gamma = sqrt(sum(((diag(n) - proj)%*%(y - X[,ind]*gamma))^2))
					u = rnorm(n-ncol(Z))
					u = u/sqrt(sum(u^2))
					y_temp = as.vector(y.hat + gamma*(diag(n)-proj)%*%X[,ind] + sigma.hat.gamma*V%*%u)
					glmnet_object_cv = glmnet::cv.glmnet(X[,-ind], y_temp-gamma*X[,ind], standardize=FALSE)
					if(lambda_cv_flag == -1){
						lambda_cv = glmnet_object_cv$lambda.min
					} else if(lambda_cv_flag == -2){
						lambda_cv = glmnet_object_cv$lambda.1se
					}
				}			
				beta = CVXR::Variable(ncol(X)+1)
				obj = sum(( (y - gamma*X[,ind]) -cbind(rep(1,n), X)%*%beta)^2)/(2*n) + lambda_cv*CVXR::p_norm(beta[-1],1)
				prob = CVXR::Problem(CVXR::Minimize(obj))
				beta.hat.stat = round(CVXR::solve(prob,feastol=1e-10,reltol=1e-10,abstol=1e-10)$getValue(beta), digits = 9)
				stat = abs(beta.hat.stat[1+ind])

				right_tail = l.cdf_adjusted(stat, y, X, ind, gamma, lambda_cv, lambda, smoothed = smoothed, tail = 'right')
				left_tail = l.cdf_adjusted(-stat, y, X, ind, gamma, lambda_cv, lambda, smoothed = smoothed, tail = 'left')
				pvals.right[j] = left_tail + right_tail
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
