# Dicrete-Continuous model version 4
# fits a probit and a regression with all error terms
# correlated
#
# Jean-Michel Tremblay
# University of Maryland
# jeanmi.tremblay@gmail.com

library(mvtnorm)

source("utils/create_X.R")
source("utils/probitUtils.R")

source("models/logit.R")
#source("models/dc4.R")

dc4 = list(
	LLVec = function(b, args){
		params = dc4Params(b, nparPr = ncol(args$XPr)
				, nparReg = ncol(args$XReg))
		L = params$L
		YReg = args$YReg
		YPr = args$YPr
		S = L %*% t(L)
		
		# LL for the regression
		epsilon = YReg - (args$XReg %*% as.matrix(params$betaReg))
		sigma2 = S[nrow(S), ncol(S)]
		LReg = dnorm(epsilon, 0, sigma2)
		
		# utilities
		U = matrix(args$XPr %*% params$betaPr, nrow = args$n
				, ncol = args$nalt, byrow = TRUE)

		# we simulate for the LL of the probit, 
		# we need a loop because at every line we need to reparametrize
		LPr = rep(0, length(args$YPr))
		for(i in 1:length(args$YPr)){
			SigmaDiffCond = condCov(S, args$nalt)
			muDiffCond = condMean(S, mu = rep(0, args$nalt)
					, posObs = args$nalt, obs = epsilon[i])
			p = 0
			#S = SigmaDiffCond
			#muCond = muDiffCond
			#Uline = U[i,]
			probitArgs = list(nalt = args$nalt, choice = YPr[i], mu = muDiffCond
					, S = SigmaDiffCond, U = U[i,]) 
			#print(args)
			if(args$method == "simulation")
				p = simProbitL(c(list(Z = args$Z[[i]]), probitArgs))
			if(args$method != "simulation")
				p = probitL(probitArgs)
			LPr[i] = p
			}
		LReg * LPr
		},
	
	computeArgs = function(spec, D){
		XReg = as.matrix(D[, spec$reg])
		YReg = D[, spec$YReg ]
		XPr = getDiscreteMatrix(spec$probit, D)
		YPr = D[, spec$YPr]
		YPr = YPr - min(YPr) + 1	
		nalt = length(spec$probit$specific)
		n = nrow(D)
		args = list(XReg = XReg, YReg = YReg, XPr = XPr, YPr = YPr
				, nalt = nalt, n = n, method = spec$method)

		# if we simulate the probas, generate the error terms here
		if(spec$method == "simulation"){
			args[["Z"]] = list()
			for(i in 1:n)
				args[["Z"]][[i]] = matrix(rnorm(spec$nsim * (nalt-1))
						, (nalt -1), spec$nsim)
			}
		args
		},
		
	computeStart = function(spec, D){
		nalt = length(spec$probit$specific)
		L = model(logit, c(list("Y" = spec$YPr), spec$probit), D)
		formReg = as.formula(paste(spec$YPr, paste(spec$reg,collapse="+")
				, sep = " ~ 0 +"))
		R = lm(formReg, D)
		betaPr = L$results$beta_hat / (pi / sqrt(6))
		betaReg = R$coefficients
		sigma2 = mean(R$residuals ** 2)
		
		# if errors of utilities ~ N(0,I) then differences ~ N(0, I + 1)
		S = diag(nalt)
		S[nalt, nalt] = sigma2
		S[-nalt, -nalt] = S[-nalt, -nalt] + 1
		L = t(chol(S))
		# we fix the first element of L so we remove it
		c(betaPr, betaReg, mat2vec(L)[-1])
		},
		
	computeOther = function(spec, D)	{
		nalt = length(spec$probit$specific)
	
		# record the names of coefficients
		namesPr = getNames(spec$probit, D)
		namesReg = spec$reg
		namesSigma = c()
		for(i in 2:nalt)
			for(j in 1:i)
				namesSigma = c(namesSigma,paste("L_",i,j,sep=""))
		names = c(namesPr, namesReg, namesSigma)
		list(names = names)
		},
	
	findStartL = function(m){
		i = 1
		while(m$results$name[i] != "L_21")
			i = i + 1
		i
		},	
		
	reparam = function(m){
		startL = dc4$findStartL(m)
		endL = length(m$results$beta_hat)
		L = vec2mat(c(sqrt(2), m$results$beta_hat[startL:endL]))
		S = L %*% t(L)
		m2 = m
		m2$results = m2$results[1:(startL -1),]
		m2$SigmaDiffFirst = S
		class(m2) = "dc4"
		m2
		},
		
	comparePmvSim = function(m, spec, D, output){
		# calculate LL for Pmv
		specPmv = spec
		specPmv[["method"]] = "pmvnorm"
		
		specSim = spec
		specSim[["method"]] = "simulation"
		if(! "nsim" %in% names(specSim))
			specSim[["nsim"]] = 500
			
		argsPmv = dc4$computeArgs(specPmv, D)
		argsSim = dc4$computeArgs(specSim, D)
		
		pSim = dc4$LLVec(m$results$beta_hat, argsSim)
		pPmv = dc4$LLVec(m$results$beta_hat, argsPmv)
		
		pdf(output)
		plot(pSim, pPmv, xlab = "probas computed with simulation"
				, ylab = "probas computed with pmv"
				, main = "comparision of invididual probabilities"
				)
		abline(lm(pPmv ~ 0 + pSim))
		dev.off()
		},
		
	outputUtilities = function(m, spec, D, output){
		#b = as.vector(m$results$beta_hat)
		args = dc4$computeArgs(spec,D)
		#print(args$X)
		nparPr = length(spec$probit$common)
		for(i in spec$probit$specific)
			nparPr = nparPr + length(i)
			
		betaPr = m$results$beta_hat[1:nparPr]
		U = getUtilities(args$XPr, betaPr, nrow(D))
		U = round(1000*U) / 1000
		p = dc4$LLVec(m$results$beta_hat, args)
		p = round(1000*p) / 1000
		choice = D[,spec[["YPr"]]]
		
		out = data.frame(U = U, choice = choice, p = p)
		write.csv(out, output, quote = FALSE, row.names = FALSE)
		},
		
	comparePis = function(b1, b2, spec, D, output){
		args = dc4$computeArgs(spec, D)
		p1 = dc4$LLVec(b1, args)
		p2 = dc4$LLVec(b2, args)
		pdf(output)
		plot(p1,p2,xlab="p_i with first set of coefficients",
				ylab="p_i with second set", main="comparision of 2 cx sets")
		abline(lm(p2 ~ p1))
		dev.off()
		},
	
	apply = function(b, spec, D){
		XPr = create_X(spec$probit$common, spec$probit$specific, D)	
		XReg = as.matrix(D[,reg])
		params = dc4Params(b,ncol(XPr),ncol(XReg))	
		nalt = length(spec$probit$specific)
		nerr = nalt -1 +1
		n = nrow(D)

		# error terms
		err = params$L %*% matrix(rnorm(n*nerr),nrow = nerr, ncol = n)
		errsDiff = t(err[1:(nalt-1),])
		errsReg = err[nalt,]

		# calculate dependant varialbe
		b = params$betaPr
		YPr = genChoiceProbit(spec$probit, params$betaPr, D, errsDiff)
		YReg = XReg %*% params$betaReg + errsReg
		list(Ydisc = YPr, Yreg = YReg)
		}
	
	)
