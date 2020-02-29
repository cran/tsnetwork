## Simulate discrete time-series data from continues time-series data
#simulate Y's and Z's
gen.sim <- function(t =NULL, n= NULL, p =NULL, k=NULL, network=NULL, prob0=NULL)
{
  if(is.null(t)) t <- 10
  if(is.null(n)) n <- 100   
  if(is.null(p)) p <- 10
  if(is.null(k)) k <- 5
  if(is.null(network)) network <- "random"
  if(is.null(prob0)) prob0 <- 0.25
      
	Z <- sim(model="ar1", time=t, n.obs=n, n.var=p, prob0=prob0, network=network) 
	Y <- simDat(Z=Z$data1, k=k)
	class(Y) = "longitudinal"
	#true_theta <- Z$theta 
	#true_gamma <- Z$gamma
	#true_sigma <- Z$sigma
	#reformat the time-series discrete data to our format
	#Y <- reformat(dat= Y, t=t, n=n, p=p )
	#simulation <- list(dat= as.longitudinal(Y, repeats= n), true_theta= Z$theta, true_gamma = Z$gamma, true_sigma= Z$sigma)
	simulation <- list(dat= Y, true_theta= Z$theta, true_gamma = Z$gamma, true_sigma= Z$sigma)

	return(simulation)
}

#simulate continues time series 
sim <- function(model=c("ar1","ar2"),time, n.obs, n.var, seed=NULL, prob0=NULL,network=c("random","scale-free","hub","user_defined"),prec=NULL,gamma1=NULL,gamma2=NULL)
	{
	  model = match.arg(model)
	  network = match.arg(network)
	  t = time
	  n = n.obs
	  p = n.var
	  
	  if(is.numeric(seed)) r=0
	  else {seed = 123
		   r= round(runif(1),4)*10000}
	  if(model=="ar1") {
		if(network=="random") {
		  L = sugm.generator(n= n, d= p, graph="random", prob=prob0, seed=seed+1234+r, vis = FALSE)
		  LL = sugm.generator(n= n, d=p, graph="random", prob=prob0, seed=seed+4567+r, vis = FALSE)
		  LLL = sugm.generator(n= n, d=p, graph="random", prob=prob0, seed=seed+1564+r, vis = FALSE)

		}
		else if(network=="scale-free") {
		  L = sugm.generator(n=n, d=p, graph="scale-free", prob=prob0, seed=seed+1234+r, vis = FALSE)
		  LL = sugm.generator(n=n, d=p, graph="scale-free", prob=prob0, seed=seed+4567+r, vis = FALSE)
		  LLL = sugm.generator(n=n, d=p, graph="scale-free", prob=prob0, seed=seed+1564+r, vis = FALSE)

		}
		else if(network=="hub") {
		 L = sugm.generator(n=n,d=p,graph="hub", prob=prob0, seed=seed+1234+r, vis = FALSE)
		 LL = sugm.generator(n=n,d=p,graph="hub", prob=prob0, seed=seed+4567+r, vis = FALSE)
		 LLL = sugm.generator(n=n,d=p,graph="hub", prob=prob0, seed=seed+1564+r, vis = FALSE)

	   }
	   if(network=="user_defined"){
		   mu <- rep(0,p)
		   true_theta <- prec
		   sigma1 <- solve(prec)
		   true_gamma <- gamma1
		   B11 <- true_gamma
	   }
	  else{
	  true_theta = as.matrix(L$omega*L$theta)
	  diag(true_theta)=1     ##theta is the precision matrix
	  sigma1 <- L$sigma
	  mu <- rep(0,p)
	  thetaL = as.matrix(LL$omega*LL$theta)
	  thetaLL = as.matrix(LLL$omega*LLL$theta)
	  lwt = thetaL*(1*lower.tri(thetaL, diag =FALSE))
	  upt = thetaLL*(1*upper.tri(thetaLL, diag = FALSE))
	  B=upt+lwt
	  ua=rbinom(p,1,0.3)
	  uu=runif(p,0,1)
	  uau = ua*uu
	  diag(B)= uau
	   B11 = B
	 for(i in 1:p){
	   for(j in 1:p){
		if(B[i,j] != 0)
		   { m = rbinom(1,1,0.6)
			 if(m ==0) B11[i,j] = -B[i,j]
		   }

	  }}
	  true_gamma=B11

	  }        #Gamma is the autoregressive coefficient matrix
	  ##data generation
	  xtn <-array(NA,c(t,p,n))
	  xtt <- array(NA,c(t,p,1))
	  for(i in 1:n){
		#x0 <- rmvnorm(1, mu, sigma1, method="svd") from mvnorm package
		x0 <- rmvnorm(1, mu, sigma1) #Our own function
		for(j in 1:t){
		  #et <- rmvnorm(1, mu, sigma1, method="svd")
		  et <- rmvnorm(1, mu, sigma1)
		  xt <-  x0 %*% B11 + et
		  xtt[j,,] <- xt
		  x0 <- xt
	   }
	   xtn[,,i] <- round(xtt,3)
	 }
	 xy=matrix(aperm(xtn, c(3,1,2)), ncol=p)
	 data1 <- as.longitudinal(xy, repeats=n)
	 return(list(data1=data1,theta=true_theta, gamma=true_gamma,sigma=sigma1))
	 }
	 if(model=="ar2") {
	  if(network=="random") {
		L = sugm.generator(n=n,d=p,graph="random", prob=prob0, seed=seed+12346+r, vis = FALSE)
		LL = sugm.generator(n=n,d=p,graph="random", prob=prob0, seed=seed+45678+r, vis = FALSE)
		LLL = sugm.generator(n=n,d=p,graph="random", prob=prob0, seed=seed+43219+r, vis = FALSE)
		LL1 = sugm.generator(n=n,d=p,graph="random", prob=prob0, seed=seed+14578+r, vis = FALSE)
		LLL1 = sugm.generator(n=n,d=p,graph="random", prob=prob0, seed=seed+96879+r, vis = FALSE)
	  }
	  else if(network=="scale-free") {
	   L = sugm.generator(n=n,d=p,graph="scale-free", prob=prob0, seed=seed+12346+r, vis = FALSE)
	   LL = sugm.generator(n=n,d=p,graph="scale-free", prob=prob0, seed=seed+45678+r, vis = FALSE)
	   LLL = sugm.generator(n=n,d=p,graph="scale-free", prob=prob0, seed=seed+43219+r, vis = FALSE)
	   LL1 = sugm.generator(n=n,d=p,graph="scale-free", prob=prob0, seed=seed+14578+r, vis = FALSE)
	   LLL1 = sugm.generator(n=n,d=p,graph="scale-free", prob=prob0, seed=seed+96879+r, vis = FALSE)
	  }
	 else if(network=="hub") {
	  L = sugm.generator(n=n,d=p,graph="hub", prob=prob0, seed=seed+12346+r, vis = FALSE)
	  LL = sugm.generator(n=n,d=p,graph="hub", prob=prob0, seed=seed+45678+r, vis = FALSE)
	  LLL = sugm.generator(n=n,d=p,graph="hub", prob=prob0, seed=seed+43219+r, vis = FALSE)
	  LL1 = sugm.generator(n=n,d=p,graph="hub", prob=prob0, seed=seed+14578+r, vis = FALSE)
	  LLL1 = sugm.generator(n=n,d=p,graph="hub", prob=prob0, seed=seed+96879+r, vis = FALSE)
	 }
	 if(network=="user_defined"){
			mu <- rep(0,p)
		   true_theta <- prec
		   sigma1 <- solve(prec)
		   B1 <- gamma2
		   B2 <- gamma1
		  }
	  else{
	  true_theta = as.matrix(L$omega*L$theta)
	  diag(true_theta)=1     ##theta is the precision matrix
	  sigma1 <- L$sigma
	  mu <- rep(0,p)
	  thetaL = as.matrix(LL$omega*LL$theta)
	  thetaLL = as.matrix(LLL$omega*LLL$theta)
	  lwt = thetaL*(1*lower.tri(thetaL, diag =FALSE))
	  upt = thetaLL*(1*upper.tri(thetaLL, diag = FALSE))
	  B1c=upt+lwt
	  ua=rbinom(p,1,0.02)
	  uu=runif(p,0,1)
	  uau = ua*uu
	  diag(B1c)= uau
	   B11 = B1c
	 for(i in 1:p){
	   for(j in 1:p){
		if(B1c[i,j] != 0)
		   { m = rbinom(1,1,0.4)
			 if(m ==0) B11[i,j] = -B1c[i,j]
		   }

	  }}
	  B1=B11
	  thetaL1 = as.matrix(LL1$omega*LL1$theta)
	  thetaLL1 = as.matrix(LLL1$omega*LLL1$theta)
	  lwt = thetaL1*(1*lower.tri(thetaL1, diag =FALSE))
	  upt = thetaLL1*(1*upper.tri(thetaLL1, diag = FALSE))
	  B2c=upt+lwt
	  ua=rbinom(p,1,0.02)
	  uu=runif(p,0,1)
	  uau = ua*uu
	  diag(B2c)= uau
	   B22 = B2c
	 for(i in 1:p){
	   for(j in 1:p){
		if(B2c[i,j] != 0)
		   { m = rbinom(1,1,0.4)
			 if(m ==0) B22[i,j] = -B2c[i,j]
		   }

	  }}
	  B2=B22
	 
	 }
	#B1and B2 are the autoregressive coefficient matrices
	##data generation
	 xtn <-array(NA,c(t,p,n))
	 xtt <- array(NA,c(t,p,1))
	 for(i in 1:n){
	  #x0 <- rmvnorm(1, mu, sigma1, method="svd")
	  x0 <- rmvnorm(1, mu, sigma1)
	  x1 <- rmvnorm(1, mu, sigma1)
	  #x1 <- rmvnorm(1, mu, sigma1, method="svd")
	  for(j in 1:t){
		et <- rmvnorm(1, mu, sigma1)#, method="svd")
		xt <-   x1 %*% B2 +  x0 %*% B1 + et
		xtt[j,,] <- xt
		x0 <- x1
		x1 <- xt
	  }
	   xtn[,,i] <- round(xtt,3)
	  }
	  true_gamma <- rbind(B2, B1)
	  xy=matrix(aperm(xtn, c(3,1,2)), ncol=p)
	  data1 <- as.longitudinal(xy, repeats=n)

	 return(list(data1=data1,theta=true_theta, gamma=true_gamma,sigma=sigma1))
	}
}

#discretize the continues data to discrete time-series
simDat =  function( Z, k ){
	#3.simulate latent data=Z
	Z	<-	as.matrix(Z)
	p <- ncol(Z)
	#4.##########Determine marginals########### 
	#4.1) ***Arbitrary marginals 
	prob	  <- t(apply(matrix(runif(k*p),nrow=p, ncol=k), 1, function(x) x/sum(x)))#sum of marginal should be one for each j
	prob	  <- split(prob, row(prob))
	marginals <- lapply(prob, function(x) {qnorm(cumsum(x))[-length(x)]} )
	for(j in 1:p){
		breaks <- unlist( unique(c(min(Z[,j])-1, marginals[j],max(Z[,j])+1)) ) 
		Z[,j] <- as.integer(cut(Z[,j], breaks= breaks, right=FALSE))
	}

	#write.csv(Z, file=paste("data/allSimulations/effMrkrFam.csv", sep=""), row.names=FALSE )
	return(Z)
}

#Generate random mvnorm from Cholosky method.
#(Another possibility is with "svd" method.)
rmvnorm = function(n, mu, sigma )
{
	if(is.null(mu)) mu = rep(0, n)
	p = nrow(sigma)
	y  <- matrix(rnorm(n*p),nrow = p, ncol= n)
	chl<- chol(sigma)
	z <- t( t(chl) %*% y + mu )  ## N( mu,sigma)
return(z)
}

#make data as longitudinal 
#as.longitudinal = function(x, repeats=1, time)
#{ 
#   if (!is.matrix(x) )
#    stop("only matrices can be coerced to longitudinal")
# 
#  
#  # create repeats vector
#  if (length(repeats) == 1) 
#  {
#    if (dim(x)[1] %% repeats != 0)
#      stop("number of repeats incompatible with number of rows of data matrix")
#    repeats = rep(repeats, dim(x)[1] %/% repeats) 
#  } 
#  if (sum(repeats) != dim(x)[1])
#    stop("sum of repeats must equal the number of rows in data matrix")
#
#  # create time vector  
#  if (missing(time)) time = 1:length(repeats)
#  if (any(duplicated(time))) 
#    stop("duplicated entries in time vector") 
#  if (length(time) != length(repeats))
#    stop("length of time vector must be equal to the length of repeats vector")
#  if (any( diff(time) <= 0 ))  
#    stop("entries time vector must be monotonically increasing")
#  
# 
#  # construct longitudinal object
  #attr(x, "class") = "longitudinal"
  #attr(x, "time") = time
  #attr(x, "repeats") = repeats
  
#  rn = NULL
#  for (i in 1:length(repeats))
#  {
#    rn = c(rn, paste( time[i], seq(1, repeats[i]), sep="-") )
#  }
#  rownames(x) = rn
   
#  return(x)
#}

#from Flare package
sugm.generator <- function(n = 200, d = 50, graph = "random", v = NULL, u = NULL, g = NULL, 
                           prob = NULL, seed = NULL, vis = FALSE, verbose = TRUE){	
  gcinfo(FALSE)
  if(verbose) cat("Generating data from the multivariate normal distribution with the", graph,"graph structure...\n")
  
  if(graph!="random" && graph!="hub" && graph!="cluster" && graph!="band" && graph!="scale-free"){
    message("\"graph\" must be one of \"random\", \"hub\", \"cluster\", \"band\" and \"scale-free\" \n")
    message("More on help(sugm.generator) \n")
    return(NULL)
  }
  if(graph=="hub"||graph=="cluster"){
    if(d<4){
      message("d is too small, d>=4 required for",graph,"\n")
      message("More on help(sugm.generator) \n")
      return(NULL)
    }
  }
  if(graph=="random"||graph=="band"||graph=="scale-free"){
    if(d<3){
      message("d is too small, d>=3 required for",graph,"\n")
      message("More on help(sugm.generator) \n")
      return(NULL)
    }
  }
  if(is.null(seed)) seed = 1
  set.seed(seed)
  if(is.null(g)){
    g = 1
    if(graph == "hub" || graph == "cluster"){
      if(d > 40)	g = ceiling(d/20)
      if(d <= 40) g = 2
    }
  }
  
  if(graph == "random"){
    if(is.null(prob))	prob = min(1, 3/d)
    prob = sqrt(prob/2)*(prob<0.5)+(1-sqrt(0.5-0.5*prob))*(prob>=0.5)
  }
  
  if(graph == "cluster"){
    if(is.null(prob)){
      if(d/g > 30)	prob = 0.3
      if(d/g <= 30)	prob = min(1,6*g/d)
    }
    prob = sqrt(prob/2)*(prob<0.5)+(1-sqrt(0.5-0.5*prob))*(prob>=0.5)
  }  
  
  
  # parition variables into groups
  g.large = d%%g
  g.small = g - g.large
  n.small = floor(d/g)
  n.large = n.small+1
  g.list = c(rep(n.small,g.small),rep(n.large,g.large))
  g.ind = rep(c(1:g),g.list)
  rm(g.large,g.small,n.small,n.large,g.list)
  gc()
  
  # build the graph structure
  theta = matrix(0,d,d);
  if(graph == "band"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    for(i in 1:g){
      diag(theta[1:(d-i),(1+i):d]) = 1
      diag(theta[(1+i):d,1:(d-1)]) = 1
    }	
  }
  if(graph == "cluster"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    for(i in 1:g){
      tmp = which(g.ind==i)
      tmp2 = matrix(runif(length(tmp)^2,0,0.5),length(tmp),length(tmp))
      tmp2 = tmp2 + t(tmp2)		 	
      theta[tmp,tmp][tmp2<prob] = 1
      rm(tmp,tmp2)
      gc()
    }
  }
  if(graph == "hub"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    for(i in 1:g){
      tmp = which(g.ind==i)
      theta[tmp[1],tmp] = 1
      theta[tmp,tmp[1]] = 1
      rm(tmp)
      gc()
    }
  }
  if(graph == "random"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    
    tmp = matrix(runif(d^2,0,0.5),d,d)
    tmp = tmp + t(tmp)
    theta[tmp < prob] = 1
    #theta[tmp >= tprob] = 0
    rm(tmp)
    gc()
  }
  
  if(graph == "scale-free"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    out = .C("SFGen",dd0=as.integer(2),dd=as.integer(d),G=as.integer(theta),seed=as.integer(seed) ) #,PACKAGE="flare")
    theta = matrix(as.numeric(out$G),d,d)
  }
  if(graph=="band"||graph=="cluster"||graph=="hub"||graph=="random"||graph=="scale-free") {
    diag(theta) = 0
    omega = theta*v
    
    # make omega positive definite and standardized
    diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
    sigma = cov2cor(solve(omega))
    omega = solve(sigma)
  }
  
  # generate multivariate normal data
  x = rmvnorm(n,rep(0,d),sigma)
  sigmahat = cor(x)
  
  # graph and covariance visulization
  if(verbose) cat("done.\n")
  
  sim = list(data = x, sigma = sigma, sigmahat = sigmahat, omega = omega, 
             theta = theta, sparsity= sum(theta)/(d*(d-1)), graph.type=graph, prob = prob)
  class(sim) = "sim" 
  return(sim)
}
