S_TISP <- function(X, y, Lambda, betaInit, grps, k, tt, errBnd = 10^(-4), maxit = 4000) {
	n <- dim(X)[1]
	d <- dim(X)[2]
	p0 <- d / k
	beta_cur <- betaInit
	beta_new <- 1 + beta_cur
	iter_index <- 1
	while((iter_index <= maxit) & (max(abs(beta_new - beta_cur)) >= errBnd)) {
		if(iter_index == 1){
			beta_new <- beta_cur
		}
		iter_index <- iter_index + 1
		beta_cur <- beta_new
		resid <- y - X %*% beta_cur
		c <- rep(0, k)
		e <- rep(1, d)
		beta_new <- rep(0, d)
		temp_resid <- vector('list', k)
		for(jj in 1:k){
			temp_resid[[jj]] <- resid[tt[2*(jj-1)+1]:tt[2*jj]]
			c[jj] <- 1 / sqrt(sum(temp_resid[[jj]]^2))
			e[(p0*(jj-1)+1):(p0*jj)] <- c[jj]
		}
		A <- sum(c)
		tttt <- t(X) %*% resid
		xi_new <- e * tttt / A + beta_cur
		xi_new_norm <- xi_new
		thresh <- Lambda/A
		for(jj in 1:p0){
			xi_new_norm[which(grps == jj)] <- sqrt(sum(xi_new[which(grps == jj)] ^2))
		}
		ttInds <- which(xi_new_norm>thresh)
		beta_new[ttInds] =  xi_new[ttInds]/xi_new_norm[ttInds] * (xi_new_norm[ttInds] - thresh)
	}
return(beta_new)
}