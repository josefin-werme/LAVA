

process.binary = function(stat, N, freq, case.prop) {
	freq[freq > 0.5] = 1 - freq[freq > 0.5]; N.case = N * case.prop; no.snps = length(stat)
	
	corrs = c()
	for (i in 1:no.snps) {
		count.x = N[i] * c((1 - freq[i])^2, 2 * freq[i] * (1 - freq[i]), freq[i]^2); x = c(0,1,2)
		beta = find.beta(count.x, N.case[i], N[i], stat[i]) 
		
		mu = 1 / (1+exp(-(beta[1] + beta[2] * x))) 
		xty = sum(x * count.x * mu)
		if (is.finite(xty)) {
			sx = sum(count.x * x); sx2 = sum(count.x * x^2); sy = N.case[i]; 
			var.x = (sx2 - sx^2 / N[i]) / (N[i] - 1)
			var.y = (sy - sy^2 / N[i])  / (N[i] - 1)
			cov.xy = (xty - sx*sy/N[i]) / (N[i] - 1)
			if (var.x > 0 & var.y > 0) corrs[i] = cov.xy / sqrt(var.x*var.y)	
		}
	}
	corrs[!is.finite(corrs)] = NA
	return(corrs)
}


# This function reconstructs the b0 and b1 for the logistic regression of a single SNP from its test statistic
find.beta = function(count.x, N1, N.orig, stat, tolerance=1e-5, reduction=0.25) {
	x = cbind(1, c(0,1,2))
	
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
	curr = make.b1(list(b0 = log(N1/(N.orig-N1))), 0, step, N1, stat, x, count.x, reduction, tolerance)
	while (abs(step) > tolerance/10 && curr$err > tolerance) {
		prop = update.b1(curr, step, N1, stat, x, count.x, reduction, tolerance)
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