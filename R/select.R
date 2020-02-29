#ic = c("NULL", "bic","Ebic","bic_mod","aic","gic")
select <- function(model, ic, LL.obs )
{
	#model$theta.path
	#model$gamma.path
	#model$loglik.path
	
	p <- ncol(model$theta.path[ ,,1,1])
	n <- model$n.obs
	t <- model$time
	
	#non_pen_LL_psath <- - model$loglik.path #multiply by negative because QUIC gives negative log-likelihood.
	#non_pen_LL <- matrix(model$loglik.path, byrow=TRUE, nrow= dim(non_pen_LL_path)[3], ncol=dim(non_pen_LL_path)[2] ) 
	non_pen_LL <- matrix(model$loglik.path, byrow=TRUE, nrow= dim(model$loglik.path)[3], ncol=dim(model$loglik.path)[2]) # from Q function
	
	if(LL.obs == TRUE) 
	{
		#Hfun <- t(apply( model$theta.path, c(3,4), function(x) {(n*(t-1)/2)*(determinant(x, logarithm = TRUE) - sum(diag(model$Scc %*% x)))} ))
		Hfun <- t(apply( model$theta.path, c(3,4), function(x) {(n*(t-1)/2)*(determinant(matrix(x, ncol=p, nrow=p), logarithm = TRUE)$modulus - sum(diag(model$Scc %*% matrix(x, ncol=p, nrow=p))))} ))
		LL_Y <- non_pen_LL - Hfun		
	}else{
		LL_Y <- non_pen_LL
		rm(non_pen_LL)
	}
	
	###Number of nonzero entries
	df_theta <- t(apply( model$theta.path, c(3,4), function(x) (sum(x != 0) -p ))) 
    df_gamma = t(apply( model$gamma.path, c(3,4), function(x) sum(x != 0)))
	
   	if(ic == "bic") BIC <- - (n * (t-1)) * LL_Y + (log(n*(t-1)))*(df_theta/2 + p + df_gamma) 
    if(ic == "bic_mod") bic_mod <- - (n * (t-1)) * LL_Y + log(n*(t-1))*(df_theta/2 + p +df_gamma) *log(log(p + p))
	if(ic == "Ebic") Ebic <- - (n * (t-1)) * LL_Y + (log(n*(t-1)))*(df_theta/2 + p +df_gamma) + (df_theta/2 + p + df_gamma)*4*0.5*log(p + p)    
	if(ic == "aic") aic <- - (n * (t-1)) * LL_Y + 2*(df_theta/2 + p +df_gamma)
    if(ic == "gic") gic <- - (n * (t-1)) * LL_Y + log(log(n*(t-1)))*(df_theta/2 + p +df_gamma) *log(p + p)

	if(ic == "bic") index <- which(BIC== min(BIC), arr.ind=TRUE) 
	if(ic == "bic_mod") index <- which(bic_mod== min(bic_mod), arr.ind=TRUE) 
	if(ic == "Ebic")index <- which(Ebic== min(Ebic), arr.ind=TRUE) 
	if(ic == "aic")index <- which(aic== min(aic), arr.ind=TRUE)
	if(ic == "gic")index <- which(gic== min(gic), arr.ind=TRUE)
	
	opt.theta<- model$theta.path[ , ,index[2] ,index[1]]
	opt.gamma<- model$gamma.path[ , ,index[2] ,index[1]]	
	opt.LL   <- model$loglik.path[ , ,index[1]][index[2]]	

	return(list(opt.theta = opt.theta, opt.gamma = opt.gamma, opt.LL = opt.LL, lam1.opt= model$lam1.seq[index[1]], lam2.opt=model$lam2.seq[index[2]], theta.path=model$theta.path, gamma.path=model$gamma.path, LL_path=model$loglik.path, lam1.seq= model$lam1.seq , lam2.seq= model$lam2.seq))
}


#model <- list(theta.path=res$theta_path, gamma.path=res$gamma_path, loglik.path=res$LL_path, lam1.seq= res$lam1_seq, lam2.seq= res$lam2_seq)







#Model selection  (theta, gamma)
#df.theta <- apply(theta.path, c(3,4), function(x) (sum(x!= 0)-p) / 2 ) 
#df.gamma <- apply(gamma.path, c(3,4), function(x) sum(x != 0) )
#L1		 <- apply(theta.path, c(3,4), function(x) sum(abs(x)) )
#loglik	 <- matrix(loglik.path, ncol=length(rho), nrow=length(lam) )
#lam.mat	 <- matrix(lam, nrow=length(lam), ncol=length(rho))
#BIC		 <-  - n * t *( loglik + (lam.mat * L1)) + log(n * t) * ( df.theta/2 + df.gamma + p) 
##pen.BIC	 <-  - n * t * loglik  + log(n * t) * ( df.theta/2 + df.gamma + p) 
#index    <- which(BIC== min(BIC), arr.ind=TRUE)
#opt.theta<- theta.path[ , ,index[1] ,index[2]]
#opt.gamma<- gamma.path[ , ,index[1] ,index[2]]
