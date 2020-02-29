calculate_EM = function(time, Z, penalty, lower, upper, Omega, lam, rho, em_tol, em_iter, iter.Mstep, pen.diag.gamma, ncores= 1)
{
	p <- ncol(lower[[1]]) 
	n <- nrow(lower[[1]])
	c_em_iter = 1
	dif	<- 100
	
	while((c_em_iter <= em_iter) && (dif >= em_tol ))
	{
		if(c_em_iter > 1) Omega <- as.matrix(rbind(cbind(theta, gamma), cbind(t(gamma), theta))) #if  time_vary=none
		#if(c_em_iter > 1) Omega <- as.matrix(rbind(cbind(itheta, gamma), cbind(t(gamma), itheta))) #if  time_vary=none

		Estep <- calculate_E(t= time, Z=Z, lower=lower, upper=upper, Omega=Omega, time_vary="none", ncores=ncores)## room for improvement in terms of Z calculation.
		Z <- Estep$Z
		for(t in 2:(time-1)){ 
			Z[[t]][ , 1:(p/2)] <- Z[[t-1]][ ,((p/2) + 1):p]
		}
		if(c_em_iter == em_iter) Scc <- Estep$Scc #S_cc in the last EM will be used to calculate H function to obtain likelihood of observed data.

		Mstep <- calculate_M(n=n, t=time, lam=lam, rho=rho, penalty=penalty, Scc=Estep$Scc, Scp=Estep$Scp, Spp=Estep$Spp, iter.Mstep=iter.Mstep, pen.diag.gamma=pen.diag.gamma)	
		
		theta  <- Mstep$theta
		itheta <- Mstep$itheta
		gamma  <- Mstep$gamma
		
		if(c_em_iter < em_iter) rm(Mstep, Estep)
		c_em_iter <- c_em_iter + 1
		}
	
	SGamma <- as.matrix( Estep$Scc - (Estep$Scp %*% t(gamma)) - (gamma %*% Estep$Scp) + (gamma %*% Estep$Spp %*% t(gamma)))
	SGamma <- (SGamma + t(SGamma)) / 2
	Rgamma <- cov2cor(SGamma)	
	loglik <-  determinant(theta)$modulus[1] - sum(diag( theta %*% Rgamma))
	#loglik <- ((n * (t-1))/2) * [ determinant(theta)$modulus[1] - sum(diag( theta %*% Rgamma)) ]
	
	rm(Mstep, Estep)
	
	results <- list()
	results$theta	<- theta
	results$itheta	<- itheta
	results$gamma	<- gamma
	results$Scc		<- Scc
	##results$Omega	<- S_obj$ES
	#results$Z		<- Z
	results$loglik	<- loglik
	results$lam		<- lam
	results$rho		<- rho

	return(results)
}

