##############################
######### M_step #############
##############################
#lam=penalty Theta
#rho = penalty Gamma
#calculate_M = function(n, t, lam, rho, SGamma, Scc, Scp, Spp, iter.Mstep=3)
calculate_M = function(n, t, lam, rho, penalty=c("scad","lasso"), Scc, Scp, Spp, iter.Mstep, pen.diag.gamma )
{
	iter <- 1
	p <- ncol(Scc)
	
	nlam	<- (n * (t-1)) * lam ## CHECK #nlam	<- (n * (t-1)) * rho
	old.gamma = qr.solve(Scc + nlam*diag(p), Scp) #gamma 10 * 10
	if(!is.numeric(old.gamma))  old.gamma <- matrix(0,p,p)
    mab =  sum(sum(abs(old.gamma))) #initial value

	while( iter <= iter.Mstep)
	{
		SGamma <- 1/(n*(t-1))*(as.matrix( Scc - (Scp %*% t(old.gamma)) - (old.gamma %*% Scp) + (old.gamma %*% Spp %*% t(old.gamma))))
		SGamma <- (SGamma + t(SGamma)) / 2
		#diag(SGamma) <- abs(diag(SGamma))
		Rgamma <- cov2cor(SGamma) 	
		############################
		##### estimate theta########
		############################		
		if(penalty == "scad")
		{
			# Calculating SCAD penalty
			if(iter ==1 ) old.theta	<-  QUIC(S=Rgamma, lam, tol=1.0e-4, msg=0, maxIter=1e4, X.init=NULL, W.init=NULL )$X
			wt <- matrix(NA, nrow=p,ncol=p)
			a <- 3.7                                       
			for(i in 1:p)
			{
			  for(j in 1:p)
				{
					if(abs(old.theta[i,j]) <= lam ) wt[i,j] <- lam
					else { wt[i,j] <- ((a*lam-abs(old.theta[i,j]))/((a-1)))}
					#{ if((lam < abs(old.theta[i,j])) & (abs(old.theta[i,j]) <=  a*lam )) {
					#wt[i,j] <- ((a*lam-abs(old.theta[i,j]))/((a-1))) }
					#  else wt[i,j] <-  0 }
				}
			}
			diag(wt) <- 0
			quic 	<-  QUIC(S=Rgamma, wt, tol=1.0e-4, msg=0, maxIter=1e4, X.init=NULL, W.init=NULL )
		}
		if(penalty == "lasso") quic <-  QUIC(S=Rgamma, lam, tol=1.0e-4, msg=0, maxIter=1e4, X.init=NULL, W.init=NULL )
		theta <- quic$X
		Sigma <- quic$W
		
		############################
		#####estimate Gamma	########
		############################
		xtyom	<-   n*(t-1)*(Scp)  %*% theta
		warmstart = 0
		##Estimating the auotregressive coeficients based on SCAD penalty
		if(penalty == "scad")
		{
			wt1 <- matrix(NA, nrow=p,ncol=p)
			a <- 3.7
			for(i in 1:p)
			{             
			  for(j in 1:p)
				{
				  if(abs(old.gamma[i,j]) <= rho ) wt1[i,j] <- rho
				  else{if((rho < abs(old.gamma[i,j])) & (abs(old.gamma[i,j]) <  a*rho )) {
				  wt1[i,j] <- ((a*rho-abs(old.gamma[i,j]))/((a-1))) }
				   else wt1[i,j] <-0}
				}
			}
		}
		if(penalty == "lasso") wt1 <- matrix(rho, nrow= p ,ncol= p)
		rho2 <- wt1*( n * (t-1) )
		if(pen.diag.gamma == FALSE)diag(rho2) <- 0 #make it optional later. 
		total_Scc <- Scc * (t-1)*n
			
		gamma =rblasso(s=total_Scc, m=as.matrix(xtyom), om=theta, nlam=rho2, tol=1e-5, sbols=mab, maxit=100, warm=warmstart, B0=old.gamma)
		
		#Convergence criteria
		#GammaDif = sum(sum(abs(gamma - old.gamma)))
		#if(!is.numeric(gamma)) gamma <- matrix(0,p,q) 
		old.gamma = gamma
		if(penalty == "scad") old.theta = theta

		iter <- iter + 1
	}
	
	#return(list(theta = theta, gamma=gamma, Rgamma= Rgamma, loglik=loglik))
	return(list(theta = theta, itheta=Sigma,  gamma=gamma))
}


#tmp1 <- matrix(rep(res$lam1_seq,4), ncol=5, nrow=4, byrow=TRUE)
#L1 <- t(apply(res$theta_path, c(3,4), function(x) ( sum(abs((diag(x))))) ))
#pen <- (25 *tmp1)+ L1