create_X = function(common,specific,D){
	# Creates a data.frame to estimate utility based models
	# 
	# Args:
	#		common: 
	#			list of vectors containing the name
	#			of variables that share the same coefficient
	#			accross alternative
	#			the vectors must have the same number of components
	#		specific: 
	#			list of vectors of variables whose
	#			coefficients are specific to the alternative
	#		D: 	
	#			the data.frame containing the variables
	#
	#	Returns:
	#		a data.frame containing n*nalt row for which
	#		the utility of alternative i equals line_i * beta
	nalt = length(specific)
	n = nrow(D)
	nrow = nalt * n
	npar = length(common)
	for(i in specific)
		npar = npar + length(i)
	
	X = matrix(0,nrow = nrow, ncol = npar)
	col_index = 1
	
	# add columns for variables with common cx
	for(com in common){	
		X[,col_index] = as.vector(t(as.matrix(D[,com])))
		col_index = col_index + 1
		}
		
	# add columns for variables with specific cx
	line_index = 1
	for(spec in specific){
		
		for(v in spec){
			X[seq(from = line_index, by = nalt, length = n),col_index] = D[,v]
			col_index = col_index + 1
			}
		line_index = line_index + 1
		}

	labels = c()
	for(var_group in common)
		labels = c(labels,paste(var_group,collapse="_"))
	for(var_group in specific)
		labels = c(labels,var_group)
	colnames(X) = labels
	
	return(X)
	}
	
create_X_diff = function(common,specific,D){
	nalt = length(specific)
	n = nrow(D)
	
	X = create_X(common,specific,D)

	# substract the first line
	for(i in 2:nalt){
		altIndex = seq(from = i, by = nalt, length = n)
		firstAlt =  seq(from = 1, by = nalt, length = n)
		X[altIndex,] = X[altIndex,] - X[firstAlt,]
		}
	X[-seq(from = 1, by = nalt, length = n),]
}

getDiscreteMatrix = function(m, D){
	create_X(m$common, m$specific, D)
	}

getNumVar = function(m, D){
	ncol(create_X(m$common, m$specific, D))
	}
	
getNames = function(m, D){
	colnames(create_X(m$common, m$specific, D))
	}
