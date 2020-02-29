#reformat the time-series discrete data to our format
reformat <- function(dat, t, n , p )
{
	Y_t <- Y_np(dat, time=t, p=p, n=n); Y_t  
	if(! dim(Y_t)[3] == t-1) stop("t starts from 2 until t")
	lower_upper <- apply(Y_t, 3 ,function(x) lower.upper(x) ); lower_upper #A list which has t elements.
	lower <- lapply(1:(dim(Y_t)[3]), function(Time) lower_upper[[Time]]$lower)
	upper <- lapply(1:(dim(Y_t)[3]), function(Time) lower_upper[[Time]]$upper)

	return(list(Y_t = Y_t, lower=lower, upper=upper))
}

#Reformat the longitudinal data for our analysis
Y_np <- function(Y, time, p, n){
	Y <- as.matrix(Y)
	Yt <-array(NA, c(n, 2*p, (time-1)))
	for(t in 1:(time-1))
	{
		Yt[ ,,t] <- cbind( Y[(((t-1)*n)+1):(t*n), ]  , Y[((t*n)+1):((t+1)*n), ] )
	}
	return(Yt)
}