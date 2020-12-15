
### 20-09-17: updated pcov.integral so that it can be passed to the new integral.p function
pcov.integral = function(K, omega, sigma, n.iter=1000, add.reverse=T, xy.index=NULL, z.index=NULL) { #C17: moved xy.index and z.index to the back and set default values, so input format is identical to multivariate.integral() and bivariate.integral()
  P = dim(omega)[1]   #C17: moved this line up, was before "Pw = ..." before
  if (P <= 2) return(NA) #C17: fail quietly if nothing to condition on
  if (is.null(xy.index)) xy.index = 1:2 #C17: if not specified, assume first two pheno are X and Y
  if (is.null(z.index)) z.index = which(!(1:P %in% xy.index)) #C17: if not specified, assume it's all the other pheno in Omega
  
  index = check.index(xy.index, z.index, dim(omega)[1])
  if (is.null(index)) return(NA)#just failing quietly for now
  
  if (!add.reverse) {
    omega = omega[index,index]; sigma = sigma[index,index] #put y last (index == P), x second to last (index == Pw)
    Pw = P - 1; Pz = Pw - 1;
    
    fit.x = omega[1:Pz,Pw] %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,Pw]
    fit.y = omega[1:Pz,P] %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,P]
    if (omega[Pw,Pw] <= fit.x) omega[Pw,Pw] = fit.x / 0.99999 #scale var up to have R2 slightly below 1
    fit.xy = omega[1:Pw,P] %*% solve(omega[1:Pw,1:Pw]) %*% omega[1:Pw,P]
    if (omega[P,P] <= fit.xy) omega[P,P] = fit.xy / 0.99999 #scale var up to have R2 slightly below 1
    
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
    p1 = pcov.integral(K, omega, sigma, n.iter/2, add.reverse=F, xy.index=xy.index, z.index=z.index) #C17: moved index parameters to the back
    p2 = pcov.integral(K, omega, sigma, n.iter/2, add.reverse=F, xy.index=rev(xy.index), z.index=z.index) #C17: moved index parameters to the back
    return((p1+p2)/2)
  }
}



pcov.cond.stats = function(draw, K, sigma, gamma, sig.xys, var.y) {
  Pw = dim(sigma)[1]-1; Pz = Pw - 1
  i.eps = 1:Pw; i.eps.z = 1:Pz; i.eps.x = Pw
  i.delta = 1:Pw + Pw; i.delta.z = 1:Pz + Pw; i.x = 2*Pw; i.y = i.x+1
  
  dtd.w = matrix(draw[i.eps,i.eps] + draw[i.eps,i.delta] + draw[i.delta,i.eps] + draw[i.delta,i.delta], ncol=Pw)
  dtd.w[Pw,Pw] = dtd.w[Pw,Pw] + 2*t(gamma$x) %*% draw[i.delta.z,i.eps.x] + gamma$fit.x
  dtd.w[-Pw,Pw] = dtd.w[-Pw,Pw] + (draw[i.delta.z,i.delta.z] + draw[i.eps.z,i.delta.z]) %*% gamma$x
  dtd.w[Pw,-Pw] = dtd.w[-Pw,Pw]
   
  omega.w = matrix(dtd.w/K - sigma[1:Pw,1:Pw], ncol=Pw)
  omega.z.inv = tryCatch(solve(omega.w[-Pw,-Pw]),error=function(x){omega.w[-Pw,-Pw]*NA}) #silently put to NA if not invertible
  b = matrix(c(-(omega.z.inv %*% omega.w[1:Pz,Pw]), 1), ncol=1)
 
  dhw.dy = gamma$dw.dz.gamma + draw[i.eps,i.delta.z] %*% gamma$y + draw[i.eps,i.y]
  dw.ew = rbind(draw[i.delta.z,i.eps], t(gamma$x) %*% draw[i.delta.z,i.eps] + draw[i.x,i.eps])

  M = dhw.dy/K + (dw.ew + draw[i.eps,i.eps]) %*% sig.xys - sigma[1:Pw,Pw+1]
  V = t(b) %*% dtd.w %*% b * var.y  
  
  return(c(t(b) %*% M, ifelse(V >= 0, sqrt(V), NA)))
}

#checks indices, puts xy at end
check.index = function(xy.index, z.index, P) {
  xy.index = unique(round(xy.index))
  z.index = unique(round(z.index))
  index = as.numeric(c(z.index,xy.index))
  if (any(is.na(index) | index > P | index <= 0) || length(xy.index) != 2 || length(z.index) == 0 || any(z.index %in% xy.index)) index = NULL  #jw
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

# this is identical to the one in the previous script files, just included here for completeness
conditional.norm = function(obs, means, sds) {
  obs = abs(obs)
  prob = suppressWarnings(pnorm(obs, mean=means, sd=sds, lower.tail=F))
  prob = prob + suppressWarnings(pnorm(-obs, mean=means, sd=sds, lower.tail=T))
  return(mean(prob, na.rm=T))
}

