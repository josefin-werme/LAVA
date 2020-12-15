
### Joint logistic model for input processing of binary phenotypes ###

fit.logistic = function(G, X, R, N.orig, case.prop, snp.stat, snp.N=NULL, n.iter=25, svar.thresh=1.75, phen.id, loc.id, dropped=c()) {
	#x n.iter=25; svar.thresh=1.75; N.orig=loc$N[i]; case.prop = input$info$prop_cases[i]; snp.stat=loc.sum[[i]]$STAT; snp.N=loc.sum[[i]]$N; phen.id=input$info$phenotype[i]; loc.id=loc$id; dropped=dropped
	K = dim(G)[2]; K.block = dim(X)[2]; N.ref = dim(G)[1];
	N1 = case.prop * N.orig; N.ratio = N.ref/N.orig
	
	if (is.null(snp.N)) snp.N = rep(N.orig, K.block)
	beta.marg = matrix(NA, nrow=K.block, ncol=2) # third column stores a boolean indicator about the inferred scaling of the SNP (standardized vs unstandardized), fourth column stores error value
	for (i in 1:K.block) beta.marg[i,] = find.beta(X[,i], case.prop*snp.N[i], snp.N[i], snp.stat[i]) #reconstruct marginal model b0 and b1 for each SNP
	wty0 = N1*N.ratio; wty1 = get.wty(X, R, beta.marg) #adding the first for the model intercept, wty0 = t(1) %*% Y = sum(Y) = N1
	
	comp.use = which(!(1:K %in% dropped))
	while (length(comp.use) > 1) { # shouldn't get anywhere near 1 (unless locus is tiny to begin with I suppose), but just to guarantee the loop ends; K is rechecked against min.K in the calling locus.info function
		W = cbind(1, G[,comp.use])
		wty = c(wty0, wty1[comp.use])
		
		beta = c(log(N1/(N.orig-N1)), rep(0, length(comp.use))) # set starting values to null model for now
		for (i in 1:n.iter) {
			mu = 1 / (1+exp(-W%*%beta))
			s = as.numeric(mu * (1-mu))
			wsw.inv = try(solve(t(W) %*% diag(s) %*% W),silent=T)	# going to just assume this is invertible, but since cor(W) = I it would be very strange if it weren't
			if (class(wsw.inv)[1]=="try-error") { print(paste0("DEV.ERROR: wsw inversion problem for phenotype ",phen.id,", locus ",loc.id)); return(NA) } # TODO: REMOVE
			beta = beta + wsw.inv %*% (wty - t(W) %*% mu)
		}
		V = diag(wsw.inv)[-1]					# vector of sampling variances (excluding intercept)
		svar.ratio = max(V) / quantile(V, 0.5)	# using ratio of max to median; will normally be very close to 1
		if (svar.ratio < svar.thresh) break 	# if above threshold, drop the PC with highest sampling variance and rerun
		comp.use = comp.use[-which.max(V)]
	}
	
	# h2
	var.Y = case.prop*(1-case.prop) / (N.orig / (N.orig - 1))
	var.Wb = sum(wty1[comp.use]^2) / ((N.orig-1)*N.ratio)^2
	
	return(
		list(
			beta=beta[-1,],
			var=mean(V) * N.ratio,
			dropped=which(!(1:K %in% comp.use)),
			h2.obs = 1 - (1 - var.Wb / var.Y) * ((N.orig-1) / (N.orig - length(comp.use) - 1))
		)
	)
}

# helper function
# beta is a K.block x 2 matrix containing the b0 and b1 for the logistic regression model of each SNP
get.wty = function(X, R, beta) {   
	K = dim(beta)[1]; xty = rep(NA, K)
	for (i in 1:K) {
		mu = 1 / (1+exp(-(beta[i,1] + beta[i,2]*X[,i]))) 
		xty[i] = sum(X[,i]*mu)
	}
	
	return(t(R) %*% xty)
}


# This function reconstructs the b0 and b1 for the logistic regression of a single SNP from its test statistic
find.beta = function(x, N1, N.orig, stat, tolerance=1e-5, reduction=0.25) {
	tbl.x = table(x); tbl.x = tbl.x / sum(tbl.x) * N.orig; 
	val.x = as.numeric(names(tbl.x))
	n.x = as.numeric(tbl.x)
	x = cbind(1, val.x)
	
	make.b1 = function(curr, value, step, N1, stat, x, n.x, reduction, tolerance) {
		curr$b1 = value
		curr$b0 = find.b0(curr$b0, curr$b1, N1, n.x, x[,2], step=abs(step), reduction=reduction, tolerance=max(abs(step), tolerance))
		curr$err = stat.error(curr$b0, curr$b1, stat, x, n.x)    
		return(curr)
	}
	
	update.b1 = function(curr, step, N1, stat, x, n.x, reduction, tolerance) {
		return(make.b1(curr, curr$b1 + step, step, N1, stat, x, n.x, reduction, tolerance))
	}
	
	step = 0.01 * sign(stat); try.reverse = F
	curr = make.b1(list(b0 = log(N1/(N.orig-N1))), 0, step, N1, stat, x, n.x, reduction, tolerance)
	while (abs(step) > tolerance/10 && curr$err > tolerance) {
		prop = update.b1(curr, step, N1, stat, x, n.x, reduction, tolerance)
		if (prop$err < curr$err) {
			curr = prop
			try.reverse = F
		} else {
			if (try.reverse) {
				step = -step
				try.reverse = F
			} else {
				step = step * reduction
				try.reverse = (curr$b1 != 0) 
			}
		}
	}
	
	return(c(curr$b0, curr$b1))
}

# finds the b0 corresponding to the given b0, minimizing the difference between sum of model-fitted probabilities and sum of sample probabilities (ie. sum.error()) 
find.b0 = function(b0, b1, N1, n.x, val.x, step=0.01, reduction=0.25, tolerance=0.01) {
	err = sum.error(b0, b1, N1, n.x, val.x)
	while (step != 0 && abs(err) > tolerance) {
		prop.b0 = b0 + -step*sign(err)
		prop.err = sum.error(prop.b0, b1, N1, n.x, val.x)
		if (abs(prop.err) < abs(err)) {
			b0 = prop.b0
			err = prop.err
		} else {
			step = step * reduction
		}
	}
	return(b0)
}

stat.error = function(b0, b1, stat, x, n.x) {
	mu = 1/(1+exp(-(b0 + b1*x[,2])))      
	xsx = t(x) %*% diag(mu*(1-mu)*n.x) %*% x
	var = xsx[1] / (xsx[1]*xsx[4] - xsx[2]^2)    
	return(ifelse(var > 0, abs(b1/sqrt(var)-stat)/abs(stat), Inf))
}   

sum.error = function(b0, b1, N1, n.x, val.x) {
	mu = 1/(1+exp(-(b0 + b1*val.x)))      
	return((sum(mu*n.x) - N1)/N1)
}