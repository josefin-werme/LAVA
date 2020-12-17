
pcov.cond.stats = function(draw, K, sigma, gamma, sig.xys, var.y) {
	Pw = dim(sigma)[1]-1; Pz = Pw - 1
	i.eps = 1:Pw; i.eps.z = 1:Pz; i.eps.x = Pw
	i.delta = 1:Pw + Pw; i.delta.z = 1:Pz + Pw; i.x = 2*Pw; i.y = i.x+1
	
	dtd.w = matrix(draw[i.eps,i.eps] + draw[i.eps,i.delta] + draw[i.delta,i.eps] + draw[i.delta,i.delta], ncol=Pw)
	dtd.w[Pw,Pw] = dtd.w[Pw,Pw] + 2*t(gamma$x) %*% draw[i.delta.z,i.eps.x] + gamma$fit.x
	dtd.w[-Pw,Pw] = dtd.w[-Pw,Pw] + (draw[i.delta.z,i.delta.z] + draw[i.eps.z,i.delta.z]) %*% gamma$x
	dtd.w[Pw,-Pw] = dtd.w[-Pw,Pw]
	
	omega.w = matrix(dtd.w/K - sigma[1:Pw,1:Pw], ncol=Pw)
	omega.z.inv = tryCatch(solve(omega.w[-Pw,-Pw]),error=function(x){omega.w[-Pw,-Pw]*NA})  # silently put to NA if not invertible
	b = matrix(c(-(omega.z.inv %*% omega.w[1:Pz,Pw]), 1), ncol=1)
	
	dhw.dy = gamma$dw.dz.gamma + draw[i.eps,i.delta.z] %*% gamma$y + draw[i.eps,i.y]
	dw.ew = rbind(draw[i.delta.z,i.eps], t(gamma$x) %*% draw[i.delta.z,i.eps] + draw[i.x,i.eps])
	
	M = dhw.dy/K + (dw.ew + draw[i.eps,i.eps]) %*% sig.xys - sigma[1:Pw,Pw+1]
	V = t(b) %*% dtd.w %*% b * var.y  
	
	return(c(t(b) %*% M, ifelse(V >= 0, sqrt(V), NA)))
}

# checks indices, puts xy at end
check.index = function(xy.index, z.index, P) {
	xy.index = unique(round(xy.index))
	z.index = unique(round(z.index))
	index = as.numeric(c(z.index,xy.index))
	if (any(is.na(index) | index > P | index <= 0) || length(xy.index) != 2 || length(z.index) == 0 || any(z.index %in% xy.index)) index = NULL
	return(index) 
}

estimate.pcor = function(draw, sigma) {
	i.y = dim(sigma)[1]; i.x = i.y - 1; i.z = 1:(i.x-1)
	o = draw - sigma
	
	cov = o[i.x,i.y] - t(o[i.z,i.x]) %*% solve(o[i.z,i.z]) %*% o[i.z,i.y]   
	vars = c(
		o[i.x,i.x] - t(o[i.z,i.x]) %*% solve(o[i.z,i.z]) %*% o[i.z,i.x], 
		o[i.y,i.y] - t(o[i.z,i.y]) %*% solve(o[i.z,i.z]) %*% o[i.z,i.y]
	) 
	r = suppressWarnings(cov / sqrt(prod(vars)))
	return(r)
}

conditional.norm = function(obs, means, sds) {
	obs = abs(obs)
	prob = suppressWarnings(pnorm(obs, mean=means, sd=sds, lower.tail=F))
	prob = prob + suppressWarnings(pnorm(-obs, mean=means, sd=sds, lower.tail=T))
	return(mean(prob, na.rm=T))
}

# single indices x and y, vector of indices z
partial.cor = function(omega, x, y, z) {
	p.cov = partial.cov(omega, x, y, z)
	return(p.cov/sqrt(partial.var(omega, x, z) * partial.var(omega, y, z)))
}

# single indices x and y, vector of indices z
partial.cov = function(omega, x, y, z) {omega[x,y] - t(omega[z,x]) %*% solve(omega[z,z]) %*% omega[z,y]}

partial.var = function(omega, x, z) {omega[x,x] - t(omega[z,x]) %*% solve(omega[z,z]) %*% omega[z,x]}