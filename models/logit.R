source("utils/models.R")
source("utils/create_X.R")
source("utils/R_jm_utils.R")

logit = list(
	# delta is MANDATORY
	LLVec = function(b, args){
		n = length(args$Y)
		nalt = max(args$Y) - min(args$Y) + 1	

		expU = exp(args$X %*% b)

		num = expU[seq(from = 1, length = n, by = nalt) + args$Y]
		denum = rowSums(matrix(expU, ncol = nalt, byrow = T))

		num / denum
		},

	computeArgs = function(spec, D){
		X = getDiscreteMatrix(spec, D)
		Y = D[[spec$Y]]
		Y = Y - min(Y)
		delta = getElem(spec, name = "delta", default = 0.001) 
		list(X = X, Y = Y, delta = delta)
		},

	computeStart = function(spec, D){
		rep(0, getNumVar(spec, D))
		},

	computeOther = function(spec, D){
		list(names = getNames(spec, D))
		}
	) # logit function list
	
logitApply = function(b, X, nalt, n){
	expU = matrix(exp( X %*% b), n, nalt, byrow = TRUE)
	# divides each column by the sum of each rows
	
	apply(expU, 2, function(col) col / rowSums(expU))
	} 
	
#genChoice = function(specLogit, cx, D){
#	n = nrow(D)
#	nalt = length(specLogit$specific)
#	X = getDiscreteMatrix(specLogit, D)
#	U = matrix(X %*% cx + rgumbel(n*nalt), nrow = n, ncol = nalt, byrow = T)
#	apply(U,1,which.max)	
#	}
