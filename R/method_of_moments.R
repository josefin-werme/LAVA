#### METHOD OF MOMENTS 
estimate.moments = function(delta, sigma, K, null.idx=0, null.omega=NULL, only.omega=F) {
	#x delta=locus$delta; sigma=locus$sigma; K=locus$K; null.idx=0; only.omega=F
	delta = as.matrix(delta)
	P = ncol(delta)
	
	if (P!=1) {; Px = P-1; x = 1:Px; y = P	}
	
	# Omega
	if (null.idx==0) { omega = t(delta)%*%delta / K - sigma } else { omega = null.omega }
	if (only.omega) return(omega)
	
	# Gamma / tau
	gamma = t(omega[y,x] %*% solve(omega[x,x]))
	tau = omega[y,y] - omega[y,x] %*% solve(omega[x,x]) %*% omega[x,y]
	
	if (null.idx) {
		gamma.null = rep(0,P)
		gamma.null[-null.idx] = gamma
		return(list(gamma=gamma.null, tau=tau))
	} else {
		# standardise
		gamma.std = rep(NA,Px); for (p in x) { gamma.std[p] = gamma[p]*sqrt(omega[p,p]/omega[y,y]) }
		tau.std = tau / omega[y,y]
		r2 = 1 - tau.std
		return(list(omega=omega, omega.xy=omega[P,-P], gamma=c(gamma), gamma.std=c(gamma.std), tau=tau, tau.std=tau.std, r2=r2))
	}
}


# separate this into a estimate.omega() function, and a estimate.params() function


