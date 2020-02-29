#NOTE for Pariya: I used matrix S3 class to generate Omega rather than Matrix S4 class. I did this because addding Matrix to the NAMESPACE
# of the package gives warnings. I may fix this later. Matrix package only has been used in initializing Omega (Not anywhere else in the package Matrix is used). So, it make sense
# to use importFrom("Matrix", "") rather than import("Matrix").


#Calculate S for jth variable
element_ini_S = function(lower, upper, mu=0, sigma=1)
{	
  	#delta1 <- (lower[j] - mu) / sigma
	#delta2 <- (upper[j] - mu) / sigma
  	delta1 <- lower #Since mu= 0 and sigma=1 then delta1= lower
  	delta2 <- upper #Since mu= 0 and sigma=1 then delta2= upper
	
	delta1[delta1 < -100] <- -100
  	delta2[delta2 > 100] <- 100
  	pnorm.upper <- apply(delta2, 1:2, function(x) pnorm(x))
  	pnorm.lower <- apply(delta1, 1:2, function(x) pnorm(x))
  	dnorm.lower <- apply(delta1, 1:2, function(x) dnorm(x))
  	dnorm.upper <- apply(delta2, 1:2, function(x) dnorm(x))
  	tmp1 <- (dnorm.lower - dnorm.upper) / (pnorm.upper - pnorm.lower)
  	EX <- mu + tmp1 * sigma
  	
  	tmp2 <- ((delta1*dnorm.lower) - (delta2 * dnorm.upper)) / (pnorm.upper - pnorm.lower)
  	EXX <- sigma^2 + mu^2 + 2 * mu * sigma * tmp1 + sigma^2 * tmp2
	rm(delta1, delta2, tmp1, tmp2, dnorm.lower, dnorm.upper, pnorm.upper,pnorm.lower)
	gc()
	
	return(list(EX=EX, EXX=EXX))
}

#Initialize Z and a *sparse* omega(two penalties) 
initialize = function( lower, upper, lam, time_vary, ncores= NULL )
{
	p <- ncol(lower[[1]]) 
	n <- nrow(lower[[1]])
	Time <- length(lower) # == ((original_t) -1)

	#########################
	##### Initialize Z ######
	#########################  
	if(ncores > 1)
	{
	  message("Multi core", ncores, "\n") 
		cl <- makeCluster(ncores)
		tmp2 <- parLapply(cl = cl, 1:Time, function(time) { 
		element_ini_S(lower[[time]], upper[[time]], mu=0, sigma=1)}); #this may crash (lapply inside multicore)
		stopCluster(cl)
	}else{
		tmp2 <- lapply(1:Time, function(time) element_ini_S(lower[[time]], upper[[time]], mu = 0, sigma= 1))#tmp2[[Time]][[p]]$EX or tmp2[[Time]][[p]]$EXX
	}
	Z <-  lapply(1:Time, function(time)tmp2[[time]]$EX) # t*p matrix
	
	############################
	##### initialize Omega #####
	############################
	# Initializing omega matrices which contain Theta and gamma
	# making 2p*2p Omega matrix for each Z where length(Z) = t-1.
	## (i): Using approximation
	diag_element <- lapply(1:Time, function(time) colMeans(tmp2[[time]]$EXX))
	EOmega <- lapply(1:Time, function(time) t(Z[[time]]) %*% Z[[time]] / n)
	for(time in 1:Time) diag(EOmega[[time]]) <- diag_element[[time]]
	#obj <- lapply(1:Time, function(time) {glasso(s=EOmega[[time]], rho=lam, maxit=1000, penalize.diagonal=FALSE)})
	#Omega <- lapply(1:Time, function(time) {Matrix((t(obj[[time]]$wi) + obj[[time]]$wi) / 2, sparse = TRUE)})
	#EOmega <- lapply(1:Time, function(time) {(t(obj[[time]]$w) + obj[[time]]$w) / 2})
	obj <- lapply(1:Time, function(time) {QUIC(S=EOmega[[time]], lam, tol=1.0e-4, msg=0, maxIter=1e4, X.init=NULL, W.init=NULL )})
	#Omega <- lapply(1:Time, function(time) {Matrix((t(obj[[time]]$X) + obj[[time]]$X) / 2, sparse = TRUE)})
	Omega <- lapply(1:Time, function(time) {(t(obj[[time]]$X) + obj[[time]]$X) / 2 })
	EOmega <- lapply(1:Time, function(time) {(t(obj[[time]]$W) + obj[[time]]$W) / 2})
	sd_marginal  <- lapply(1:Time, function(time) {sqrt(diag(EOmega[[time]]))})
	for(time in 1:Time) sd_marginal[[time]][abs(sd_marginal[[time]]) < 1e-10] <- 1e-10
	EOmega <- lapply(1:Time, function(time) {diag(1/sd_marginal[[time]]) %*% EOmega[[time]] %*% diag(1/sd_marginal[[time]])})
	Omega <- lapply(1:Time, function(time) {diag(sd_marginal[[time]]) %*% Omega[[time]] %*% diag(sd_marginal[[time]])})
	#(ii) alternative 1
	#omega = diag(p)
	#(iii) alternative 2: Ridge regression 
	
	if(time_vary == "none") 
	{	
	  message("Compute an average Gamma across all ((original_T)-1)'s \n")
		Omega <- Reduce('+', Omega)/ Time
		EOmega <- Reduce('+', EOmega)/ Time
	}else{
	  message("Each Time == ((original_T)-1) has its own Omega \n")
		Omega <- Omega
		EOmega <- EOmega
	}
	rm(tmp2, diag_element, obj, sd_marginal)
	gc()
	
	return(list(Z=Z, Omega = Omega, EOmega=EOmega))
}