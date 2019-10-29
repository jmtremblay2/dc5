applyDCBatch = function(b, modelFns, spec, vars_to_mod, consts, B,D){
	nalt = -1
	if("Yop" %in% names(spec))
		nalt = length(unique(D[,spec$Yop]))
	if("YPr" %in% names(spec))
		nalt = length(spec$probit$specific)

	# data frame to store resultsi of apply
	applySummary = data.frame(matrix(0
			, nrow = length(vars_to_mod) * length(consts)
			, ncol = 2 + nalt +1))
	names(applySummary) = c("var","const",paste("alt",1:nalt,sep="_"),"reg")
	
	# index of applySummary row to update
	rowIndex = 1
	for(v in vars_to_mod){
		for(k in consts){
			# change the data
			D2 = D
			D2[,v] = D[,v] * k

			# will sum the results of the apply here, then compute mean
			applyDisc = rep(0,nalt)
			applyReg = 0
				# do the apply B time
			for(i in 1:B){
				a = modelFns$apply(b,spec,D2)
				applyDisc = applyDisc + as.vector(table(a$Ydisc))
				applyReg = applyReg + mean(a$Yreg)
				}	
			# divide by number of applies to get average results
			applyDisc = applyDisc / B
			applyReg = applyReg / B
			applyReg = round(100000 * applyReg) / 100000
			# write results into data frame and increment row to update
			applySummary[rowIndex,] = c(v,k,applyDisc,applyReg)
			rowIndex = rowIndex + 1
			}
		}
	applySummary
	}
