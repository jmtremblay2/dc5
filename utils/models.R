source("utils/R_jm_utils.R")
source("utils/print.R")
source("utils/LLfunctions.R")
source("utils/SD.R")	

modelFit = function(modelFns, args, start, other.params, spec){
	# check for bounds
	big_num = 9999999
	lb = rep(-big_num, length(start))
	ub = rep(big_num, length(start))
	if("lb" %in% other.params)
		lb = other.params[["lb"]]
	if("ub" %in% other.params)
		ub = other.params[["ub"]]

	O = optim(fn = LLWrapper, gr = LLGradWrapper, method = "BFGS"
			, hessian = F, par = start, control = list(fnscale = -1
			, reltol = spec$reltol), args = args, f = modelFns$LLVec)
	beta_hat = O$par
	
	if(spec[["SD"]] == "none")
		SD = NoneSD(beta_hat)
	if(spec[["SD"]] == "hessian")
		SD = HessianSD(beta_hat, modelFns, args)
	if(spec[["SD"]] == "bootstrap")
		SD = BootstrapSD(beta_hat, args, modelFns, spec, start, lb, ub)

	t = beta_hat / SD
	p = 2*pnorm(-abs(t))
	results = data.frame(name = other.params[["names"]]
			, beta_hat = beta_hat, SD = SD, t = t, p = p)
	maxLL = O$value
	list(results = results, maxLL = maxLL)				
	}
	
model = function(modelFns, spec, D){
	args = modelFns$computeArgs(spec, D)
	# if the user did not specify a value for delta we 
	# need to make one up
	if( ! "delta" %in% names(args))
		args[["delta"]] = 0.01
	
	# relative tolerance to stop the solver
	if( ! "reltol" %in% names(spec))
		spec[["reltol"]] = 1e-4

	# print details or not
	if( ! "verbose" %in% names(spec))
		spec[["verbose"]] = FALSE
	args[["verbose"]] = spec[["verbose"]]
	
	# get user specified starting value or default one
	start = NULL
	if("start" %in% names(spec))
		start = scan(spec[["start"]])
	if( ! "start" %in% names(spec))
		start = modelFns$computeStart(spec, D)
	
	other.params = modelFns$computeOther(spec, D)
	
	# if SD is not specified, calculate hessian based SDs
	if( ! "SD" %in% names(spec))
		spec[["SD"]] = "hessian"

	M = modelFit(modelFns, args, start, other.params, spec)
	
	# if output is requested, we plot the LL
	if("output" %in% names(spec))
		plotLL(spec$plotrange, spec$numpoints, M$results$beta_hat,
				spec$output, f, args, M$results$name)

	class(M) = "dcModel"
	return(M)
	}
