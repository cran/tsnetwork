 #NOTE: we start from t=2 then all loop starts from 1:t-1 (e.g.lapply(1:(t-1), function(tim)... )

##############################
######### E_step #############
##############################
calculate_E <- function(t, Z, lower, upper, Omega, time_vary="none" , ncores=1)
{
	p <- ncol(lower[[1]]) 
	n <- nrow(lower[[1]])
	Time  <- length(lower)
	if(! Time == (t-1)) stop("check time-points \n")
	
	#******when Theta and Gamma vary across time. So we have a Gamma matrix for each individual.*****
	#if(time_vary == "both") approx <- lapply(1:(t-1) , function(i) lapply(1:n, function(i) calculate.R.approx.internal(Z=Z[[Time]][i,], lower=lower[[Time]][i,], upper=upper[[Time]][i,], Sigma= as.matrix(Omega[[Time]]), ncores=ncores))) #approximation for all Time points
	#*****when Theta and Gamma fixed across time. So we have one average Gamma matrix for all Time's.*****
	if(time_vary == "none") approx <- lapply(1:(t-1) , function(Time) lapply(1:n, function(i) calculate.R.approx.internal(Z=Z[[Time]][i,], lower=lower[[Time]][i,], upper=upper[[Time]][i,], Sigma= as.matrix(Omega), ncores=ncores))) #approximation for all Time points
	if(time_vary == "none")
	{		
		Z.new <- lapply(1:(t-1), function(Time){t(sapply( 1:n, function(i) {cbind(approx[[Time]][[i]]$Z)} ))})

		all.EOmega <- lapply(1:(t-1), function(Time){lapply(1:n, function(i){approx[[Time]][[i]]$EOmega} )})
		all.EOmega <- unlist(all.EOmega, recursive = FALSE)
	
		K <- Reduce('+', all.EOmega ) /(n*(t-1)) #if (start_t == 1) K <- Reduce('+', all.EOmega ) /(n*t)
		Scc <-  K[((p/2)+1):p, ((p/2)+1):p]
		Spp	<-  K[1:(p/2), 1:(p/2)]
		Scp	<-  K[((p/2)+1):p, 1:(p/2)]

		return(list(Z = Z.new, Scc=Scc, Spp=Spp, Scp=Scp))
		rm(approx)
	}
}	
