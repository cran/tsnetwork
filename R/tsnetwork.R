tsnetwork <- function( dat, lower, upper, penalty= "scad", n_lam1= NULL, lam1_ratio= NULL , n_lam2= NULL, lam2_ratio= NULL, em_tol= NULL, em_iter= NULL, iter_Mstep= NULL, pen_diag_gamma= FALSE, ncores = 1)
{
 
  if(is.null(lam1_ratio))  lam1_ratio <- 0.2 ############################################## 
  if(is.null(lam2_ratio))  lam2_ratio <- 0.2 ############################################## 
  if(is.null(em_tol)) em_tol <- 0.01  ############################################ 
  if(is.null(em_iter)) em_iter <- 2   ############################################ 
  if(is.null(iter_Mstep)) iter_Mstep <- 3 #########################################
  
  
  if(is.longitudinal(dat) == TRUE)
  {
    n = get.time.repeats(dat)$repeats[[1]]
    t = get.time.repeats(dat)$time
 
    p = dim(dat)[2]
  }else{ message("Data format is not longitudinal. \n") }

  ts.dat.format <- reformat(dat= dat, t= max(t), n= n, p= p )
  lower <- ts.dat.format$lower
  upper <- ts.dat.format$upper

  p <- ncol(lower[[1]])/2 
	n <- nrow(lower[[1]])
	t <- (length(lower) + 1) #time.points from original data
	
	# penalize theta
	n_lam1 = n_lam1
	lam1_ratio = lam1_ratio
	#S <- cor(dat, use="pairwise.complete.obs", method="kendall" )
	#lam1_max = max(max(S - diag(ncol(S))), -min(S) )
	if(p == 10) lam1_max <- 0.15 #0.2
	if(p == 50) lam1_max <- 0.1 #0.25
	if(p < 10 ) lam1_max <- 0.15 #0.2
	lam1_min = lam1_ratio * lam1_max
	lam1 = exp(seq(log(lam1_max), log(lam1_min), length = n_lam1)) 

	# penalize Gamma
	n_lam2 = n_lam2
	lam2_ratio = lam2_ratio
	#lam2_max <- lam1_max - 0.05
	lam2_max <- 0.23 #0.23
	#if(p == 10) lam2_max <- 0.25
	#if(p == 50) lam2_max <- 0.5
	lam2_min = lam2_ratio * lam2_max
	lam2 = exp(seq(log(lam2_max), log(lam2_min), length = n_lam2)) ;lam2

	loglik.path <- array( NA, c(1, length(lam1), length(lam2)) )	
	theta.path <- array( NA, c(p, p, length(lam1),length(lam2)) )	
	gamma.path <- array( NA, c(p, p, length(lam1), length(lam2)) )
	for( r in 1: length(lam2) )
	{
		for( k  in 1:length(lam1))
		{
			#initials only depends on lam1 so this double for loop can be optimize.
		  message("ini is starting for lam1= ", k, " and lam2= ", r, " \n")
			ini = initialize(lower=lower, upper=upper, lam=lam1[k], time_vary="none", ncores=ncores )
			Z	= ini$Z; Z #A list which has t elements. Each element of this list contains a n*p matrix. ( Note: p is sum of the variables in time t and time (t-1) )
			Omega	= ini$Omega #A list which has n elements. Each elements is a p*p matrix (precision + gamma).
			EOmega	= ini$EOmega # do we need this?
#---------------------HERE ERROR ------------------------------------			
			message("EM is starting for lam1= ", k, " and  lam2= ", r, " ... \n")		
			EM <- calculate_EM( time = t, Z = ini$Z, penalty=penalty , lower=lower, upper=upper, Omega=ini$Omega, lam = lam1[k], rho=lam2[r], em_tol= em_tol, em_iter=em_iter, iter.Mstep=iter_Mstep, pen.diag.gamma=pen_diag_gamma, ncores= ncores)
			theta.path[ , ,k, r] <- EM$theta
			gamma.path[ , ,k, r] <- EM$gamma
			loglik.path[ ,k, r] <- EM$loglik
			Scc		<- EM$Scc #S_cc will be used to calculate H function to obtain likelihood of observed data.
			rm(EM)
			message("EM is done for lam1= ", k, " and lam2= ", r, " \n")
		}
	}
	return(list(theta.path= theta.path, gamma.path=gamma.path, loglik.path=loglik.path, lam1.seq=lam1, lam2.seq=lam2 , n.obs =n , time= t, Scc=Scc))
}