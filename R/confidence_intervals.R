#The functions do some bounds checking on correlations truncating them to just below the bound to stop the matrixsampling::rwishart function from aborting
#(I will pull apart the matrixsampling::rwishart function at some point to improve efficiency further, we can fix this then too)
#In general though, assumes these to have been validated in advance (ie. omega.x is invertible, etc.)
#Output includes the observed value used for reference (note that for corrs/R2 of 1 (or -1) this will be slightly below 1 because of the bounds correcting,
# so I'd generally use original estimates here and not directly output the estimate value from these functions)


#CIs are for correlations in omega
#It works for any size of omega, if more than 2x2 it computes CI's for each correlation
#It's more efficient to do that in one go rather than do all the pairwise correlations separately
#!!! HOWEVER: it seems that due to a check inside the matrixsampling::rwishart function it won't work if omega as a whole is not invertible (it will complain about Theta not being positive)
#We can fix that later once we deconstruct the matrixsampling::rwishart function itself, but for now might be best only to use this for individual correlations,
#just iterating over all pairs instead
ci.bivariate = function(K, omega, sigma, n.iter=10000) {  # K = locus$K; omega = locus$omega; sigma = locus$sigma
  S = diag(sqrt(diag(omega)))
  corrs = solve(S) %*% omega %*% solve(S)
  corrs[corrs >= 1] = 0.99999; corrs[corrs <= -1] = -0.99999; diag(corrs) = 1
  omega = S %*% corrs %*% S

  P = dim(omega)[1]; tri = lower.tri(corrs)
  out = data.frame(pheno1=col(corrs)[tri], pheno2=row(corrs)[tri], r=corrs[tri])
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
    qq = apply(r, 1, quantile, c(0.025, 0.975), na.rm=T)
    qq[qq < -1] = -1; qq[qq > 1] = 1
    out$ci.rho.low = qq[1,]
    out$ci.rho.high = qq[2,]
    
    # quantiles r2
    qq = apply(r^2, 1, quantile, c(0.025, 0.975), na.rm=T)  #JW
    qq[qq < -1] = -1; qq[qq > 1] = 1
    out$ci.r2.low = qq[1,]
    out$ci.r2.high = qq[2,]
    if (sign(out$ci.rho.low)!=sign(out$ci.rho.high)) out$ci.r2.low = 0  # set r2 CI to 0 if rho CI spans 1
  }
  return(round(out,5))  #jw (round)
}

#expects omega.x to be invertible
# K=locus$K; omega=locus$omega; sigma=locus$sigma; n.iter=10000
ci.multivariate = function(K, omega, sigma, n.iter=10000) {
  P = dim(omega)[1]
  S = diag(sqrt(diag(omega)))
  corrs = solve(S) %*% omega %*% solve(S)
  corrs[corrs >= 1] = 0.99999; corrs[corrs <= -1] = -0.99999; diag(corrs) = 1
  omega = S %*% corrs %*% S
  fit = omega[-P,P] %*% solve(omega[-P,-P]) %*% omega[-P,P]
  if (fit >= omega[P,P]) omega[P,P] = fit/0.99999
  #increasing omega_Y to fit if r2 > 1; setting r2 slightly below 1 in that case, since otherwise the matrixsampling::rwishart function will fail 
  
  gamma.ss = solve(corrs[-P,-P]) %*% corrs[-P,P]
  r2 = max(0, fit/omega[P,P])

  draws = tryCatch(matrixsampling::rwishart(n.iter, K, Sigma=sigma/K, Theta=omega), error=function(e){return(NULL)}) 
  if (!is.null(draws)) {
    est = apply(draws, 3, estimate.std, sigma)  
    qq = apply(est, 1, quantile, c(0.025, 0.975), na.rm=T)
    #qq[qq < -1] = -1; qq[qq > 1] = 1
  } else {
    qq = matrix(NA, nrow=2, ncol=P)
  }
  qq.r2 = qq[,P]; qq.r2[qq.r2 < 0] = 0; qq.r2[qq.r2 > 1] = 1;

  ci = list(
    gamma = round(data.frame(est=gamma.ss, ci.low=qq[1,-P], ci.high=qq[2,-P]),5),  #jw (round)
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
