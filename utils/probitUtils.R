library(mvtnorm)

source("utils/create_X.R")

getUtilities = function(X, b, n){
	matrix( X %*% b, nrow = n, ncol = nrow(X) / n, byrow = TRUE)
	}

genChoiceProbit = function(m, b, D, errsDiff = NULL, SDiff = NULL ){
	nalt = length(m$specific)
	n = nrow(D)
	
	X = getDiscreteMatrix(m, D)
	U = getUtilities(X, b, n)
	# Y1 = apply(U,1,which.max)
	# substract the first column (first alternative) to all other columns
	Udiff = U[,-1] - U[,1]
	
	if(is.null(errsDiff)){
		errsDiff = t(chol(SDiff)) %*% matrix(rnorm((nalt-1) * n), nalt -1, n)
		errsDiff = t(errsDiff)
		}
		
	#errsDiff = errsDiff / 0.01 
	Udiff = Udiff + errsDiff
	Udiff = cbind(rep(0,n), Udiff)
	Y2 = apply(Udiff,1,which.max)
	# mean(Y1 == Y2)
	Y2
	}

dc4Params = function(x,nparPr,nparReg){
	betaPr = x[1:nparPr]
	betaReg = x[(nparPr + 1):(nparPr + nparReg)]
	L = vec2mat(c(sqrt(2),x[(nparPr + nparReg + 1):length(x)]))
	return(list(betaPr = betaPr, betaReg = betaReg, L = L))
	}

# turns a lower triangular matrix into a vector of its elements
mat2vec = function(x){
	v = c()
	for(i in 1:nrow(x))
		v = c(v,x[i,1:i])
	v	
	}

# turns a vector that contains the elements of a lower triangular
# matrix into the matrix itself
vec2mat = function(x){
	n = round((-1 + sqrt(1 + 8*length(x))) / 2 + 0.1)
	L = matrix(0,n,n)
	index = 1
	for(i in 1:n){
		for(j in 1:i){
			L[i,j] = x[index]
			index = index + 1
			}
		}
	L
	}

# calculates the conditional covariance matrix of 
# normal vector conditional on observing observation
# INDEX posObs. Please note that the conditional
# covariance does not depend on the observation itself
condCov = function(sigma, posObs){
	sigma11 = sigma[-posObs,-posObs]
	sigma22 = sigma[posObs,posObs]
	sigma12 = sigma[-posObs,posObs]
	sigma21 = sigma[posObs,-posObs]
	
	sigma11 - sigma12 %*% solve(sigma22) %*% sigma21
	}

# calculate the conditional mean of a random vector
# conditional on observing obs at index posObs
condMean = function(sigma, mu, posObs, obs){
	sigma12 = sigma[-posObs,posObs]
	sigma22 = sigma[posObs,posObs]
	mu1 = mu[-posObs]
	mu2 = mu[posObs]
	
	mu1 + sigma12 %*% solve(sigma22) %*% (obs - mu2)
	}

# calculate the covariance of the differences between 
# a random vector and one of its component (indexed by
# baseUt)
cov2diff = function(n, baseUt){
	ndiff = n - 1
	
	M = matrix(0, nrow = ndiff, ncol = n)
	M[,baseUt] = -1
	M[,-baseUt] = diag(ndiff)
		
	M
	}

# generates a matrix M that transforms differences with
# respect to error term index currentBase to differences
# with respect to index newBase.
# returns identity of the two bases are the same	
reparamErrors = function(n,currentBase, newBase){
	if(currentBase == newBase)
		return(diag(n))
	
	minus1_col = newBase - (newBase > currentBase)
	zero_row = currentBase - (currentBase > newBase)
	
	M = matrix(0,n,n)
	M[,minus1_col] = -1
	M[-zero_row,-minus1_col] = diag(n-1)
	
	M
	}


simProbitLBACKUP = function(U, choice, SCond, muCond, Z, nalt){
	
	M = reparamErrors(nalt - 1, currentBase = 1, newBase = choice)
	muReparam = M %*% muCond
	SigmaReparam = M %*% SCond %*% t(M)
	LReparam = t(chol(SigmaReparam))
	
	nsim = ncol(Z)

	# P(choice) = P(err < UDiff)
	# we simulate this
	err = LReparam %*% Z
	UDiff = U[choice] - U[-choice]
	p = sum(apply(err, 2, function(e) (sum( e <  UDiff) == (nalt-1))))
	if(0 == p){
		p = 0.1
		}
	p / nsim
	}

simProbitL = function(args){
	M = reparamErrors(args$nalt - 1, currentBase = 1, newBase = args$choice)
	muReparam = M %*% args$mu
	SReparam = M %*% args$S %*% t(M)
	LReparam = t(chol(SReparam))
	
	# P(choice) = P(err < UDiff)
	err = LReparam %*% args$Z
	err = err + matrix(muReparam, nrow = nrow(err), ncol = ncol(err), byrow = FALSE)
	UDiff = args$U[ args$choice] - args$U[- args$choice]
	nsucc = max(0.1, sum(apply(err, 2, function(e) 
			(sum( e <  UDiff) == (args$nalt-1)))))
	nsucc / ncol(args$Z)
	}

pmvProbitL = function(args){
	M = reparamErrors(args$nalt - 1, currentBase = 1, newBase = args$choice)
	muReparam = as.vector(M %*% args$mu)
	SReparam = M %*% args$S %*% t(M)
	
	#UDiff = args$U[- args$choice] - args$U[args$choice]
	UDiff = args$U[ args$choice ] - args$U[ - args$choice ]
	pmvnorm(upper = UDiff,mean = muReparam, sigma = SReparam)
	}

probitL = function(args){	
	p = 0
	#if(args[["method"]] == "pmvnorm")
	#	p = pmvProbitL(args)
	#if(args[["method"]] == "simulation")
	#	p = simProbitL(args)
	#p
	pmvProbitL(args)
	}
