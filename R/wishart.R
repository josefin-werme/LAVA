### All functions used for p-value computation ###

# the integral.func argument should just be one of bivariate.integral, multivariate.integral or pcov.integral
# for pcov.integral, omega and sigma should be formatted such that the first two pheno are X and Y, and the remainder is Z
integral.p = function(integral.func, K, omega, sigma, min.iter=10000, adap.thresh=c(1e-4, 1e-6)) {
	tot.iter = min.iter * 10^(0:length(adap.thresh))
	adap.thresh = c(adap.thresh,0)  # adding dummy 0 at the end to simplify loop code
	
	p = 1; curr.iter = 0
	for (i in 1:length(tot.iter)) {
		add.iter = tot.iter[i] - curr.iter
		add.p = integral.func(K, omega, sigma, n.iter=add.iter)
		p = (curr.iter*p + add.iter*add.p) / tot.iter[i]
		curr.iter = tot.iter[i]
		if (all(is.na(p)) || all(p[!is.na(p)] >= adap.thresh[i])) break
	}
	return(p)
}


# tests gamma_j = 0 for each element j of gamma
# input omega and sigma should be the whole matrix, including Y
# computing observed and null gammas internally, to save on number of required input arguments
# not doing checks on eg. invertability here, is assumed to have been checked and addressed before
multivariate.integral = function(K, omega, sigma, n.iter=1000) {
	P = dim(sigma)[1]; Px = P - 1
	omega.x = omega[-P,-P]; omega.xy = omega[-P,P]
	sig.xys = solve(sigma[1:Px,1:Px]) %*% sigma[1:Px,P] / K
	var.y = as.numeric(sigma[P,P] - sigma[P,1:Px] %*% solve(sigma[1:Px,1:Px]) %*% sigma[1:Px,P]) / K^2
	
	sigma.use = matrix(0,P+Px,P+Px); sigma.use[1:Px,1:Px] = sigma[1:Px,1:Px]
	theta = matrix(0,P+Px,P+Px); theta[1:Px+Px,1:Px+Px] = K*omega.x; theta[P+Px,P+Px] = K
	
	draws = matrixsampling::rwishart(n.iter, K, Sigma=sigma.use, Theta=theta)
	param = apply(draws, 3, multi.cond.stats, K, sigma, sig.xys, var.y)
	C1 = param[1:(Px^2),]
	C2 = param[(1:Px)+(Px^2),]
	C3 = param[(1:Px)+(Px^2+Px),]
	sds = param[(1:Px)+(Px^2+2*Px),]
	
	gamma.ss.obs = diag(sqrt(diag(omega.x))) %*% solve(omega.x) %*% omega.xy
	p.out = rep(NA, Px)
	for (index in 1:Px) {
		gamma.null = rep(0,Px); gamma.null[-index] = solve(omega.x[-index,-index]) %*% omega.xy[-index]
		tau.null = as.numeric(omega[P,P] - t(omega.xy[-index]) %*% solve(omega.x[-index,-index]) %*% omega.xy[-index])
		tau.null.sqrt = ifelse(tau.null > 0, sqrt(tau.null), 0) # setting negative tau to 0 here
		
		M = gamma.null %*% C1[1:Px+(index-1)*Px,] + tau.null.sqrt * C2[index,] + C3[index,]
		p.out[index] = conditional.norm(gamma.ss.obs[index], M, sds[index,])
	}
	return(p.out)
}


multi.cond.stats = function(draw, K, sigma, sig.xys, var.y) {
	Px = dim(sigma)[1]-1; i.eps = 1:Px; i.delta = i.eps + Px; i.y = 2*Px+1
	dtd.x = matrix(draw[i.eps,i.eps] + draw[i.eps,i.delta] + draw[i.delta,i.eps] + draw[i.delta,i.delta], ncol=Px)
	
	omega.x = matrix(dtd.x/K - sigma[1:Px,1:Px], ncol=Px)
	omega.x.inv = tryCatch(solve(omega.x),error=function(x){omega.x*NA}) # silently put to NA if not invertible
	O.x = suppressWarnings(diag(sqrt(diag(omega.x)), ncol=Px) %*% omega.x.inv)
	
	sds = suppressWarnings(sqrt(diag(var.y * O.x %*% dtd.x %*% t(O.x))))
	
	C1 = (draw[i.delta,i.delta] + draw[i.delta,i.eps]) %*% t(O.x) / K
	C2 = O.x %*% (draw[i.delta,i.y] + draw[i.eps,i.y])/K
	C3 = O.x %*% ((draw[i.delta,i.eps] + draw[i.eps,i.eps]) %*% sig.xys - sigma[1:Px,Px+1])
	
	return(c(C1,C2,C3,sds))
}


bivariate.integral = function(K, omega, sigma, n.iter=1000, add.reverse=T) {
	if (!add.reverse) {
		omega.null = diag(diag(omega))
		sig.use = matrix(0,3,3); sig.use[1,1] = sigma[1,1]
		theta = matrix(0,3,3); theta[-1,-1] = omega.null*K
		
		sig.xy = sigma[1,2]
		sig.xys = sig.xy/sigma[1,1]
		var.y = sigma[2,2] - (sigma[1,2]^2)/sigma[1,1]
		
		params = apply(matrixsampling::rwishart(n.iter, K, Sigma=sig.use, Theta=theta), 3, bivar.cond.stats, K=K, sig.xy, sig.xys, var.y)   # first row is means, second is SDs
		return(conditional.norm(omega[1,2], params[1,], params[2,]))
	} else {
		p1 = bivariate.integral(K, omega, sigma, n.iter/2, add.reverse=F)
		p2 = bivariate.integral(K, omega[2:1,2:1], sigma[2:1,2:1], n.iter/2, add.reverse=F)    
		return((p1+p2)/2)
	}
}

# this is an internal function for the apply in integral.p(), defined here for clarity
# draw will be the 3x3 matrix drawn from the wishart
bivar.cond.stats = function(draw, K, sig.xy, sig.xys, var.y) {
	m = draw[2,3] + draw[1,3] + sig.xys*(draw[1,2] + draw[1,1])
	m = m/K - sig.xy
	
	v = var.y * (draw[2,2] + 2*draw[1,2] + draw[1,1])    
	v = v / K^2
	v = ifelse(v <= 0, NA, sqrt(v))
	return(c(m,v))
}

conditional.norm = function(obs, means, sds) {
	obs = abs(obs)
	prob = suppressWarnings(pnorm(obs, mean=means, sd=sds, lower.tail=F))
	prob = prob + suppressWarnings(pnorm(-obs, mean=means, sd=sds, lower.tail=T))
	return(mean(prob, na.rm=T))
}

pcov.integral = function(K, omega, sigma, n.iter=1000, add.reverse=T, xy.index=NULL, z.index=NULL) {
	P = dim(omega)[1] 
	if (P <= 2) return(NA)                                      # fail quietly if nothing to condition on
	if (is.null(xy.index)) xy.index = 1:2                       # if not specified, assume first two pheno are X and Y
	if (is.null(z.index)) z.index = which(!(1:P %in% xy.index)) # if not specified, assume it's all the other pheno in Omega
	
	index = check.index(xy.index, z.index, dim(omega)[1])
	if (is.null(index)) return(NA)                              # failing quietly if index is null
	
	if (!add.reverse) {
		omega = omega[index,index]; sigma = sigma[index,index]    # put y last (index == P), x second to last (index == Pw)
		Pw = P - 1; Pz = Pw - 1;
		
		fit.x = omega[1:Pz,Pw] %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,Pw]
		fit.y = omega[1:Pz,P] %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,P]
		if (omega[Pw,Pw] <= fit.x) omega[Pw,Pw] = fit.x / 0.99999               # scale var up to have R2 slightly below 1
		fit.xy = omega[1:Pw,P] %*% solve(omega[1:Pw,1:Pw]) %*% omega[1:Pw,P]
		if (omega[P,P] <= fit.xy) omega[P,P] = fit.xy / 0.99999                 # scale var up to have R2 slightly below 1
		
		var.y = as.numeric(sigma[P,P] - sigma[P,1:Pw] %*% solve(sigma[1:Pw,1:Pw]) %*% sigma[1:Pw,P]) / K^2
		sig.xys = solve(sigma[1:Pw,1:Pw]) %*% sigma[1:Pw,P] / K
		gamma.parts = list(
			x = solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,Pw],
			y = solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,P],
			fit.x = fit.x*K,
			dw.dz.gamma = c(omega[1:Pz,P], omega[1:Pz,Pw] %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,P])*K
		)
		
		sigma.use = matrix(0,P+Pw,P+Pw); sigma.use[1:Pw,1:Pw] = sigma[1:Pw,1:Pw]
		theta = matrix(0,P+Pw,P+Pw); theta[1:Pz+Pw,1:Pz+Pw] = K*omega[1:Pz,1:Pz];
		diag(theta)[Pw+Pw:P] = K*c(omega[Pw,Pw] - fit.x, omega[P,P] - fit.y)
		
		draws = matrixsampling::rwishart(n.iter, K, Sigma=sigma.use, Theta=theta)
		params = apply(draws, 3, pcov.cond.stats, K, sigma, gamma.parts, sig.xys, var.y)
		
		pcov.obs = omega[Pw,P] - t(omega[1:Pz,Pw]) %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,P]
		p = conditional.norm(pcov.obs, params[1,], params[2,])
		return(p)
	} else {
		p1 = pcov.integral(K, omega, sigma, n.iter/2, add.reverse=F, xy.index=xy.index, z.index=z.index)      # moved index parameters to the back
		p2 = pcov.integral(K, omega, sigma, n.iter/2, add.reverse=F, xy.index=rev(xy.index), z.index=z.index) # moved index parameters to the back
		return((p1+p2)/2)
	}
}






##########################
#  comparison functions  #
##########################
# TODO: remove all these?


#just hardcoding four pheno
wish.sampler.4var = function(K, omega, sigma, index, n.iter=1000) {
	params = generate.params(omega, sigma, index)
	obs = params$gamma.ss[index]
	perm = apply(matrixsampling::rwishart(n.iter, K, Sigma=sigma, Theta=K*params$omega.null), 3, est.gamma.ss, K, sigma, index)
	p = mean(abs(perm) > abs(obs))
	return(p)
}

mvn.sampler.4var = function(K, omega, sigma, index, n.iter=1000) {
	params = generate.params(omega, sigma, index)
	delta = make.delta.4var(K, params$omega.null)
	obs = params$gamma.ss[index]
	perm = rep(NA, n.iter)
	for (i in 1:n.iter) {
		delta.hat = delta + rmvnorm(K, sigma=sigma) 
		perm[i] = est.gamma.ss(t(delta.hat) %*% delta.hat, K, sigma, index)
	}
	p = mean(abs(perm) > abs(obs), na.rm=T)
	return(p)
}

est.gamma.ss = function(dtd, K, sigma, index) {
	omega.hat = dtd/K - sigma
	gamma.hat = solve(omega.hat[1:3,1:3]) %*% omega.hat[1:3,4]
	gamma.ss.hat = diag(sqrt(diag(omega.hat)[-4])) %*% gamma.hat
	return(gamma.ss.hat[index])
}

generate.params = function(omega, sigma, null.index) {
	omega.x = omega[1:3,1:3]
	gamma = solve(omega.x) %*% omega[1:3,4]
	gamma.ss = diag(sqrt(diag(omega.x))) %*% gamma
	gamma.null = rep(0,3)
	gamma.null[-null.index] = solve(omega.x[-null.index,-null.index]) %*% omega[1:3,4][-null.index]
	omega.null = omega
	omega.null[1:3,4] = omega.null[4,1:3] = omega.x %*% gamma.null
	return(list(gamma.ss=gamma.ss, omega.null=omega.null))
}

# helper function for null delta for mvn bivariate
make.delta.4var = function(K, omega) {
	P = 4
	delta = scale(rmvnorm(K, sigma=diag(P)), scale=F)
	eig = eigen(t(delta) %*% delta/K)
	delta = delta %*% eig$vectors %*% diag(1/sqrt(eig$values))

	eig = eigen(omega)    
	delta = delta %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
	return(delta)  
}


wish.sampler.bivar = function(K, omega, sigma, n.iter=1000) {
	obs = omega[1,2]
	omega.null = diag(diag(omega))
	
	prod.xy = apply(matrixsampling::rwishart(n.iter, K, Sigma=sigma, Theta=K*omega.null), 3, function(x) {x[1,2]})
	omega.xy = prod.xy/K - sigma[1,2]
	p = mean(abs(omega.xy) > abs(obs))
	return(p)
}


mvn.sampler.bivar = function(K, omega, sigma, n.iter=1000) {
	delta = make.delta.bivar(K, diag(omega))
	
	obs = omega[1,2]
	prod.xy = rep(NA, n.iter)
	for (i in 1:n.iter) {
		delta.hat = delta + rmvnorm(K, sigma=sigma)    
		prod.xy[i] = sum(delta.hat[,1]*delta.hat[,2])
	}
	omega.xy = prod.xy/K - sigma[1,2]
	p = mean(abs(omega.xy) > abs(obs))
	return(p)
}

# helper function for null delta for mvn bivariate
make.delta.bivar = function(K, vars) {
	P = 2
	delta = scale(rmvnorm(K, sigma=diag(P)), scale=F)
	eig = eigen(t(delta) %*% delta/K)
	delta = sweep(delta %*% eig$vectors %*% diag(1/sqrt(eig$values)), 2, sqrt(vars), FUN="*")
	return(delta) 
}
