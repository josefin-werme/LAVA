# The functions do some bounds checking on correlations truncating them to just below the bound to stop the matrixsampling::rwishart function from aborting
# In general, its assumed that the input has been validated in advance (ie. omega.x is invertible, etc.)

### BIVARIATE ###
ci.bivariate = function(K, omega, sigma, n.iter=10000) {  
	S = diag(sqrt(diag(omega)))
	corrs = solve(S) %*% omega %*% solve(S)
	corrs[corrs >= 1] = 0.99999; corrs[corrs <= -1] = -0.99999; diag(corrs) = 1
	omega = S %*% corrs %*% S
	
	P = dim(omega)[1]; tri = lower.tri(corrs)
	out = data.frame(pheno1=col(corrs)[tri], pheno2=row(corrs)[tri], r=corrs[tri], rho.lower=NA, rho.upper=NA, r2.lower=NA, r2.upper=NA)
	draws = tryCatch(matrixsampling::rwishart(n.iter, K, Sigma=sigma/K, Theta=omega), error=function(e){return(NULL)})
	if (!is.null(draws)) {
		if (P == 2) {
			func = function(draw, sigma) {o = draw-sigma; o[1,2]/sqrt(o[1,1]*o[2,2])}
			r = matrix(suppressWarnings(apply(draws, 3, func, sigma)), nrow=1)
		} else {
			func = function(draw, sigma) {cov2cor(draw-sigma)[lower.tri(sigma)]}
			r = suppressWarnings(apply(draws, 3, func, sigma))
		}
		# quantiles rho
		qq = round(apply(r, 1, quantile, c(0.025, 0.975), na.rm=T),5)
		qq[qq < -1] = -1; qq[qq > 1] = 1
		out$rho.lower = qq[1,]; out$rho.upper = qq[2,]
		
		# quantiles r2
		qq = round(apply(r^2, 1, quantile, c(0.025, 0.975), na.rm=T),5)
		qq[qq < 0] = 0; qq[qq > 1] = 1
		out$r2.lower = qq[1,]; out$r2.upper = qq[2,]
		if (sign(out$rho.lower)!=sign(out$rho.upper)) out$r2.lower = 0  # set r2 lower to 0 if rho CI spans 0
	}
	
	return(out)
}


### MULTILPE REG ###
# expects omega.x to be invertible
ci.multivariate = function(K, omega, sigma, n.iter=10000) {
	P = dim(omega)[1]
	S = diag(sqrt(diag(omega)))
	corrs = solve(S) %*% omega %*% solve(S)
	corrs[corrs >= 1] = 0.99999; corrs[corrs <= -1] = -0.99999; diag(corrs) = 1
	omega = S %*% corrs %*% S
	fit = omega[-P,P] %*% solve(omega[-P,-P]) %*% omega[-P,P]
	if (fit >= omega[P,P]) omega[P,P] = fit/0.99999
	# increasing omega_Y to fit if r2 > 1; setting r2 slightly below 1 in that case, since otherwise the matrixsampling::rwishart function will fail 
	
	gamma.ss = solve(corrs[-P,-P]) %*% corrs[-P,P]
	r2 = max(0, fit/omega[P,P])
	
	draws = tryCatch(matrixsampling::rwishart(n.iter, K, Sigma=sigma/K, Theta=omega), error=function(e){return(NULL)}) 
	if (!is.null(draws)) {
		est = apply(draws, 3, estimate.std, sigma)  
		qq = apply(est, 1, quantile, c(0.025, 0.975), na.rm=T)
	} else {
		qq = matrix(NA, nrow=2, ncol=P)
	}
	qq.r2 = qq[,P]; qq.r2[qq.r2 < 0] = 0; qq.r2[qq.r2 > 1] = 1;
	
	ci = list(
		gamma = round(data.frame(est=gamma.ss, ci.low=qq[1,-P], ci.high=qq[2,-P]),5), 
		r2 = round(data.frame(est=r2, ci.low=qq.r2[1], ci.high=qq.r2[2]),5)
	)
	return(ci)
}


estimate.std = function(draw, sigma) {
	P = dim(sigma)[1]
	o = cov2cor(draw-sigma)
	g = solve(o[-P,-P]) %*% o[-P,P]
	r2 = o[-P,P] %*% g
	return(c(g,r2))
}



### PARTIAL COR ###
# xy.index must be two unique index values, z.index any number (>0) of values not in xy.index
# expects omega.z to be invertible
# returns vector with three values: estimate (for reference), ci.low, ci.high
ci.pcor = function(K, xy.index, z.index, omega, sigma, n.iter=10000) {
	index = check.index(xy.index, z.index, dim(omega)[1]); out = rep(NA,3)
	if (is.null(index)) return(out) # just failing quietly for now
	
	omega = omega[index,index]; sigma = sigma[index,index]  
	P = dim(omega)[1]; Pw = P - 1; Pz = Pw - 1;
	S = diag(sqrt(diag(omega)))
	corrs = solve(S) %*% omega %*% solve(S)
	corrs[corrs >= 1] = 0.99999; corrs[corrs <= -1] = -0.99999; diag(corrs) = 1
	omega = S %*% corrs %*% S
	
	fit.x = omega[1:Pz,Pw] %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,Pw]
	fit.y = omega[1:Pz,P] %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,P]  
	if (omega[Pw,Pw] <= fit.x) omega[Pw,Pw] = fit.x / 0.99999 #scale var up to have R2 slightly below 1
	fit.xy = omega[1:Pw,P] %*% solve(omega[1:Pw,1:Pw]) %*% omega[1:Pw,P]
	if (omega[P,P] <= fit.xy) omega[P,P] = fit.xy / 0.99999 #scale var up to have R2 slightly below 1
	
	p.cov = omega[Pw,P] - t(omega[1:Pz,Pw]) %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,P]
	p.vars = c(omega[Pw,Pw] - fit.x, omega[P,P] - fit.y)
	out[1] = ifelse(all(p.vars > 0), p.cov/sqrt(prod(p.vars)), NA)
	
	draws = tryCatch(matrixsampling::rwishart(n.iter, K, Sigma=sigma/K, Theta=omega), error=function(e){return(NULL)}) 
	if (!is.null(draws)) {
		est = apply(draws, 3, estimate.pcor, sigma)  
		out[2:3] = quantile(est, c(0.025, 0.975), na.rm=T)    
		out[out < -1] = -1; out[out > 1] = 1
	}
	return(round(out,5))
}                    
