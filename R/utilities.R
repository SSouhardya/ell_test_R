#' standardize a vector
#' 
#' @param y. Numeric vector.
#' @return Returns vector of same length as `y` formed by standardizing `y`.
#' @export
std<-function(y){
	y = as.vector(y)
	m.y = mean(y)
	sd.y = sqrt(var(y)*(length(y) - 1)/length(y))
	return((y - m.y)/sd.y)
}

#' normalize a vector
#' 
#' @param v. Numeric vector.
#' @return Returns a vector of same length as `v` formed by dividing `v` with its `l2-norm`.
#' @export
g<-function(v){	
	return(v/sqrt(sum(v^2))  )
}

#' tail probabilities of the distribution of u1
#' 
#' @param q Numeric. `u1` value.
#' @param n Integer. The degrees of freedom; implies that `u` lies in an `n-1`-dimensional sphere.
#' @param lower.tail. Logical. If `TRUE`, returns the lower tail probability; if `FALSE`, returns the upper tail probability. Default is `TRUE`.
#' @param known_sigma. Logical. If `TRUE`, returns probabilities of a normal distribution with standard deviation `sigma.hat` (mentioned below) instead of the `u1`-distribution.
#'  Default is FALSE.
#' @param sigma.hat. Numeric. If `known_sigma = TRUE`, this is the standard deviation of the normal distribution. Default is 1.
#' @return Returns the percentile of `q` with respect to the `u1`-distribution or the normal distribution corresponding to whether
#'  `known_sigma` is `FALSE` or `TRUE`.
#' @export
qhaar <- function(q,n,lower.tail=TRUE,stoperr=FALSE, known_sigma = FALSE, sigma.hat = 1){
	if(known_sigma){
		return(pnorm(q*sigma.hat, lower.tail = lower.tail))
	}
	if( (abs(q)>1) & stoperr){
		p = NA
		stop("impossible haar quantile")
	}else if(q>=1){
		p = 1
	}else if(q<=-1){
		p = 0
	}else if(q==0){
		p = 0.5
	}else if(q<0){
		p = pt(sqrt((n-1)/(1/q^2-1)),n-1,lower.tail=FALSE)
	}else{
		p = 1-pt(sqrt((n-1)/(1/q^2-1)),n-1,lower.tail=FALSE)
	}
	if(!lower.tail){p = 1-p}
	p
}
