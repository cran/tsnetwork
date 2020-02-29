compare = function ( par=c("theta","gamma"), true_G, est_G) 
{
	true_G <- as.matrix(true_G)
	est_G <- as.matrix(est_G)
		
	p        <- ncol(true_G)
	true_G   <-  abs(sign(true_G))  - diag(rep(1,p))
	est_G	 <-  abs(sign(est_G)) - diag(rep(1,p)) 
	if( par == "theta" )
	{
		true_G[lower.tri(true_G, diag = TRUE)] <- 0
		est_G[lower.tri(est_G, diag = TRUE)]  <- 0
	}
	tp.all   <- (true_G != 0) * (est_G != 0) 
	fp.all   <- (true_G == 0) * (est_G != 0) 
	fn.all   <- (true_G != 0) * (est_G == 0) 
	tn.all   <- (true_G == 0) * (est_G == 0)
	
	tp       <- sum(tp.all)
	fp       <- sum(fp.all)
	fn       <- sum(fn.all)
	#tn       <- sum(tn.all)
	tn       <- sum(tn.all[upper.tri(tn.all == 1)])
	
	F1score <- (2 * tp) / (2 * tp + fp + fn)
	sen <-  tp /(tp + fn)
	spe <- tn/(tn + fp)
	mcc <-  (tp *tn-fp*fn)/ sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
	
	return( list(F1score=F1score, sen=sen, spe=spe , mcc=mcc ))
}

