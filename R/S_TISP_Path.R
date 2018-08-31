#' S_TISP_Path Function
#'
#' 
#' This function allows you to obtain Estimation of high-dimensional 
#' multi-response regression with heterogeneous noises under 
#' eterogeneous group square-root Lasso penalty.
#'
#' @param X A block diagonal design matrix.For each block, each row represents an observation. All blocks share the same number of columns.
#' @param y response vector whose length equals to the sum of number of observations across all groups.
#' @param grps a vector to indicate which group each entry of beta belongs
#' @param k number of groups
#' @param lambdas a vector of tuning parameters of group lasso penalty 
#' @param index a vector indicates the starting point and ending point for each group. For example, if there is 100 samples in the first group and 150 samples in the second group, then it is c(1,100,101,250)
#' @export
#' @examples
#'	p <- 10
#'	n <- 20
#'	k <- 2
#'	X <- matrix(0, n*k, p*k)
#'	X[1:n, 1:p] <- rnorm(n*p)
#'	X[(n+1):(k*n), (p+1):(p*k)] <- rnorm(n*p)
#'	beta <- c(0:9, (0:9)/2)
#'	y <- X %*% beta + rnorm(n*k)*0.1
#'	grps <- rep(1:p, k)
#'	lambdas <- (1:5)/2
#'	index <- c(1, n, n+1, 2*n)
#'	betaest <- S_TISP_Path(X, y, grps, k, index, lambdas)


S_TISP_Path <- function(X, y, grps, k, index, lambdas) {
	n <- dim(X)[1]
	d <- dim(X)[2]
	k0 <- norm(X, '2') / sqrt(2)
	X <- X / k0
	y <- y / k0
	lambdas <- lambdas / k0
	tt <- index
	betaEsts <- matrix(0, d, length(lambdas))
	for(ind in 1:length(lambdas)) {
	Lambda <- lambdas[ind]
	betaInit <- rep(0, d)
	betaEst = S_TISP(X, y, Lambda, betaInit, grps,k,tt)
	betaEsts[,ind] <- betaEst
	}
return(betaEsts)
}
