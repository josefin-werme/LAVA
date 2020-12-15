### 20/09/17: updated integral.p that requires integral.func to be specified (rather than autodetected)

### 20/09/17: updated integral.p that requires integral.func to be specified (rather than autodetect)
#C17: the integral.func argument should just be one of bivariate.integral, multivariate.integral or pcov.integral
#C17: for pcov.integral, omega and sigma should be formatted such that the first two pheno are X and Y, and the remainder is Z
integral.p = function(integral.func, K, omega, sigma, min.iter=10000, adap.thresh=c(1e-4, 1e-6)) { #C17: moved integral.func from a variable defined inside the function to a required function argument
  #C17: these first two lines can go away, left them in now for reference of what changed
  #  P = as.numeric(ncol(omega)) #C17: deleted
  #  integral.func = ifelse(P == 2, bivariate.integral, multivariate.integral) #C17: deleted
  tot.iter = min.iter * 10^(0:length(adap.thresh))
  adap.thresh = c(adap.thresh,0) # adding dummy 0 at the end to simplify loop code
  
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


#tests gamma_j = 0 for each element j of gamma
#input omega and sigma should be the whole matrix, including Y
#computing observed and null gammas internally, to save on number of required input arguments
# -> not doing checks on eg. invertability here, is assumed to have been checked and addressed before
#NB: I haven't really tested whether this works if Px = 1 (though not sure why it would ever be used for that)
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
    tau.null.sqrt = ifelse(tau.null > 0, sqrt(tau.null), 0) #setting negative tau to 0 here

    M = gamma.null %*% C1[1:Px+(index-1)*Px,] + tau.null.sqrt * C2[index,] + C3[index,]
    p.out[index] = conditional.norm(gamma.ss.obs[index], M, sds[index,])
  }
  return(p.out)
}

multi.cond.stats = function(draw, K, sigma, sig.xys, var.y) {
  Px = dim(sigma)[1]-1; i.eps = 1:Px; i.delta = i.eps + Px; i.y = 2*Px+1
  dtd.x = matrix(draw[i.eps,i.eps] + draw[i.eps,i.delta] + draw[i.delta,i.eps] + draw[i.delta,i.delta], ncol=Px)

  omega.x = matrix(dtd.x/K - sigma[1:Px,1:Px], ncol=Px)
  omega.x.inv = tryCatch(solve(omega.x),error=function(x){omega.x*NA}) #silently put to NA if not invertible
  O.x = suppressWarnings(diag(sqrt(diag(omega.x)), ncol=Px) %*% omega.x.inv)
  
  sds = suppressWarnings(sqrt(diag(var.y * O.x %*% dtd.x %*% t(O.x))))

  C1 = (draw[i.delta,i.delta] + draw[i.delta,i.eps]) %*% t(O.x) / K
  C2 = O.x %*% (draw[i.delta,i.y] + draw[i.eps,i.y])/K
  C3 = O.x %*% ((draw[i.delta,i.eps] + draw[i.eps,i.eps]) %*% sig.xys - sigma[1:Px,Px+1])

  return(c(C1,C2,C3,sds))
}

#this was previously the integral.p() function
#moved the functionality to swap order for half the permutations into this function
bivariate.integral = function(K, omega, sigma, n.iter=1000, add.reverse=T) {
  if (!add.reverse) {
    omega.null = diag(diag(omega))
    sig.use = matrix(0,3,3); sig.use[1,1] = sigma[1,1]
    theta = matrix(0,3,3); theta[-1,-1] = omega.null*K
    
    sig.xy = sigma[1,2]
    sig.xys = sig.xy/sigma[1,1]
    var.y = sigma[2,2] - (sigma[1,2]^2)/sigma[1,1]
    
    params = apply(matrixsampling::rwishart(n.iter, K, Sigma=sig.use, Theta=theta), 3, bivar.cond.stats, K=K, sig.xy, sig.xys, var.y) #first row is means, second is SDs
    return(conditional.norm(omega[1,2], params[1,], params[2,]))
  } else {
    p1 = bivariate.integral(K, omega, sigma, n.iter/2, add.reverse=F)
    p2 = bivariate.integral(K, omega[2:1,2:1], sigma[2:1,2:1], n.iter/2, add.reverse=F)    
    return((p1+p2)/2)
  }
}

#this is an internal function for the apply in integral.p(), defined here for clarity
#draw will be the 3x3 matrix drawn from the wishart, see Word doc for details
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

###################################


##########################
#  comparison functions  #
##########################

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

#just helper function for null delta for mvn bivariate
make.delta.4var = function(K, omega) {P = 4
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

#just helper function for null delta for mvn bivariate
make.delta.bivar = function(K, vars) {P = 2
  delta = scale(rmvnorm(K, sigma=diag(P)), scale=F)
  eig = eigen(t(delta) %*% delta/K)
  delta = sweep(delta %*% eig$vectors %*% diag(1/sqrt(eig$values)), 2, sqrt(vars), FUN="*")
  return(delta)  
}
