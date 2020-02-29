element_S <- function( lower, upper, mu=0, sigma1=1)
{
  	delta1 <- (lower - mu) / sigma1 #10.2
  	delta2 <- (upper - mu) / sigma1 #Inf , upper=Inf
  	
	if(delta1 < -100) delta1 <- -100
  	if(delta2 > 100) delta2 <- 100
	dif_pnorm <- pnorm(delta2) - pnorm(delta1)
	if(dif_pnorm == 0) 
	{
		tmp1 <- 0 
		tmp2 <- 0 
		}else{
		tmp1 <- (dnorm(delta1) - dnorm(delta2)) / dif_pnorm		
		tmp2 <- (delta1*dnorm(delta1) - delta2*dnorm(delta2)) / dif_pnorm
	}

	EX <- mu + tmp1 * sigma1
  	EXX <- sigma1^2 + mu^2 + sigma1^2 * tmp2 + 2*mu*sigma1*tmp1
	
	rm(delta1, delta2, tmp1, tmp2 )
	gc()
	
	return(list(EX=EX, EXX=EXX))
}

## calculates mean and variance for conditional Gaussian
#cond_N <- function(j, Sigma, Z , Z_new, diag_element, lower, upper)
cond_N <- function(j, Sigma, Z , lower, upper)
{
	p <- ncol(Sigma)
	
	tmp <- matrix(Sigma[j, -j], 1, p-1)
	tmp1 <- solve(Sigma[-j, -j])

	#tmp1 <- glasso(Sigma[-j, -j], rho=.0001)$wi
	
	mu <- as.numeric(tmp %*% tmp1 %*% Z[-j])
	#mu <- as.vector(mu)		 
	
	#sigma1 <- as.numeric(Sigma[j, j] - tmp %*% tmp1 %*% t(tmp))
	sigma1 <- as.numeric( abs(1 - tmp %*% tmp1 %*% t(tmp)))
	sigma1 <- sqrt(sigma1) 

	obj <- element_S( lower= lower[j], upper= upper[j], mu=mu, sigma1=sigma1)      
	Z_new <- obj$EX
	diag_element <- mean(obj$EXX)
	
	rm(tmp, tmp1, mu, sigma1, obj )
	gc()
	
	return(list(Z_new=Z_new, diag_element=diag_element))
}

calculate.R.approx.internal <- function( Z, lower, upper, Sigma, ncores )
{      
	if(is.null(ncores)) ncores = detectCores() - 1
	
	n.var <- length(Z) # == 2 * p
	p <- n.var
	#nT <- nrow(y) 
	#Z_new <- rep(0, p)
	#Z_new <- matrix(0, n, p)
	#diag_element <- rep(0, p)
	
	if(ncores > 1){
	  message("Multi_core", ncores, "\n")
		cl <- makeCluster(ncores)
		#clusterEvalQ(cl, source("cond_N.R"))
		#clusterEvalQ(cl, source("element_S.R"))
		Sigma <- Sigma
		Z <- Z
		#lower_upper <- lower_upper
		cond_norm <- parLapply(cl = cl, 1:p, function(j) { 
			#cond_N(j, Sigma=Sigma, Z=Z , Z_new=Z_new, diag_element=diag_element, lower=lower, upper=upper ); 
			cond_N(j, Sigma=Sigma, Z=Z , lower=lower, upper=upper ); 
		})
		stopCluster(cl)
	}else{
		#cond_norm <- lapply(1:p, function(j){cond_N(j, Sigma, Z , Z_new, diag_element, lower=lower, upper=upper)} )
		cond_norm <- lapply(1:p, function(j){cond_N(j, Sigma, Z , lower=lower, upper=upper)} )
	}
	
	Z_new <- do.call(cbind, lapply(1:p, function(x)cond_norm[[x]]$Z_new ))
	#Z_new <-  t(sapply(1:t, function(tim){do.call(cbind, lapply(1:p, function(x)cond_norm[[tim]][[x]]$Z_new ))}))
	
	diag_element <- sapply(1:p, function(x) cond_norm[[x]]$diag_element)
	
	#ES <- t(Z_new) %*% Z_new / T
	EOmega <- t(Z_new) %*% Z_new 
	diag(EOmega) <- diag_element

	output <- list()
	output$EOmega <- EOmega
	output$Z <- Z_new
	
	rm(lower, upper, EOmega, Z_new )
	return(output)         
}