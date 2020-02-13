"heidel.diag" <- function (x, eps = 0.1, pvalue=0.05) 
{
  if (is.mcmc.list(x)) 
    return(lapply(x, heidel.diag, eps))
  x <- as.mcmc(as.matrix(x))
  HW.mat0 <- matrix(0, ncol = 6, nrow = nvar(x))
  dimnames(HW.mat0) <- list(varnames(x),
                            c("stest", "start", "pvalue", "htest",
                              "mean", "halfwidth"))
  HW.mat <- HW.mat0
  for (j in 1:nvar(x)) {
    start.vec <- seq(from=start(x), to = end(x)/2, by=niter(x)/10)
    Y <- x[, j, drop = TRUE]    
    n1 <- length(Y)
    ## Schruben's test for convergence, applied sequentially
    ##
    S0 <- spectrum0.ar(window(Y, start=end(Y)/2))$spec
    converged <- FALSE
    for (i in seq(along = start.vec)) {
      Y <- window(Y, start = start.vec[i])
      n <- niter(Y)
      ybar <- mean(Y)
      B <- cumsum(Y) - ybar * (1:n)
      Bsq <- (B * B)/(n * S0)
      I <- sum(Bsq)/n
      if(converged <- !is.na(I) && pcramer(I) < 1 - pvalue)
        break
    }
    ## Recalculate S0 using section of chain that passed convergence test
    S0ci <- spectrum0.ar(Y)$spec
    halfwidth <- 1.96 * sqrt(S0ci/n)
    passed.hw <- !is.na(halfwidth) & (abs(halfwidth/ybar) <= eps)
    if (!converged || is.na(I) || is.na(halfwidth)) {
      nstart <- NA
      passed.hw <- NA
      halfwidth <- NA
      ybar <- NA
    }
    else {
      nstart <- start(Y)
    }
    HW.mat[j, ] <- c(converged, nstart, 1 - pcramer(I), 
                     passed.hw, ybar, halfwidth)
  }
  class(HW.mat) <- "heidel.diag"
  return(HW.mat)
}

"print.heidel.diag" <-
  function (x, digits = 3, ...) 
{
  HW.title <- matrix(c("Stationarity", "test", "start", "iteration",
                       "p-value", "", 
                       "Halfwidth", "test", "Mean", "", "Halfwidth", ""),
                     nrow = 2)
  y <- matrix("", nrow = nrow(x), ncol = 6)
  for (j in 1:ncol(y)) {
    y[, j] <- format(x[, j], digits = digits)
  }
  y[, c(1, 4)] <- ifelse(x[, c(1, 4)], "passed", "failed")
  y <- rbind(HW.title, y)
  vnames <- if (is.null(rownames(x))) 
    paste("[,", 1:nrow(x), "]", sep = "")
  else rownames(x)
  dimnames(y) <- list(c("", "", vnames), rep("", 6))
  print.default(y[, 1:3], quote = FALSE, ...)
  print.default(y[, 4:6], quote = FALSE, ...)
  invisible(x)
}

"spectrum0.ar" <- function(x)
{
  x <- as.matrix(x)
  v0 <- order <- numeric(ncol(x))
  names(v0) <- names(order) <- colnames(x)
  z <- 1:nrow(x)
  for (i in 1:ncol(x))
  {
      lm.out <- lm(x[,i] ~ z)
      if (identical(all.equal(sd(residuals(lm.out)), 0), TRUE)) {
          v0[i] <- 0
          order[i] <- 0
      }
      else {
          ar.out <- ar(x[,i], aic=TRUE)
          v0[i] <- ar.out$var.pred/(1 - sum(ar.out$ar))^2
          order[i] <- ar.out$order
      }
  }
  return(list(spec=v0, order=order))
}

effectiveSize <- function(x, crosschain = FALSE)
{
  if (is.mcmc.list(x))
    {
      if (crosschain) {
        ans <- effectiveSize_crosschain(x)
      }else{
        ##RGA changed to sum across all chains
        ess <- do.call("rbind",lapply(x,effectiveSize))
        ans <- apply(ess,2,sum)
      }
    }
  else
    if (crosschain) message('Option crosschain only works with multiple chains; running default ESS calculations.')
    {
      x <- as.mcmc(x)
      x <- as.matrix(x)
      spec <- spectrum0.ar(x)$spec
      ans <- ifelse(spec==0, 0, nrow(x) * apply(x, 2, var)/spec)
    }
  return(ans)
}

"spectrum0" <- function(x, max.freq=0.5, order=1, max.length=200)
{
  x <- as.matrix(x)
  if (!is.null(max.length) && nrow(x) > max.length) {
    batch.size <- ceiling(nrow(x)/max.length)
    if (is.R()) {
      x <- aggregate(ts(x, frequency=batch.size), nfreq = 1, FUN=mean)
    }
    else {
      x <- aggregate(ts(x, frequency=batch.size), nf = 1, fun=mean)
    }
  }
  else {
    batch.size <- 1
  }
  
  out <- do.spectrum0(x, max.freq=max.freq, order=order)
  out$spec <- out$spec * batch.size
  return(out)
}

"do.spectrum0" <- function(x, max.freq=0.5, order=1)
{
  ## Estimate spectral density of time series x at frequency 0.
  ## spectrum0(x)/length(x) estimates the variance of mean(x)
  ##
  ## NB We do NOT use the same definition of spectral density
  ## as in spec.pgram.
  ##
  fmla <- switch(order+1,
                 spec ~ one,
                 spec ~ f1,
                 spec ~ f1 + f2)
  if(is.null(fmla))
    stop("invalid order")

  N <- nrow(x)
  Nfreq <- floor(N/2)
  freq <- seq(from = 1/N, by = 1/N, length = Nfreq)
  f1 <- sqrt(3) * (4 * freq - 1)
  f2 <- sqrt(5) * (24 * freq^2 - 12 * freq + 1)
  v0 <- numeric(ncol(x))
  for(i in 1:ncol(x)) {
    y <- x[,i]
    if (var(y) == 0) {
      v0[i] <- 0
    }
    else {
      yfft <- fft(y)
      spec <- Re(yfft * Conj(yfft))/ N
      spec.data <- data.frame(one = rep(1, Nfreq), f1=f1, f2=f2,
                              spec = spec[1 + (1:Nfreq)],
                              inset = I(freq<=max.freq))
      
      glm.out <- glm(fmla, family=Gamma(link="log"), data=spec.data)
      v0[i] <- predict(glm.out, type="response",
                       newdata=data.frame(spec=0,one=1,f1=-sqrt(3),f2=sqrt(5)))
    }
  }
  return(list(spec=v0))
}

"pcramer" <- function (q, eps=1.0e-5)
{
  ## Distribution function of the Cramer-von Mises statistic
  ##
  log.eps <- log(eps)
  y <- matrix(0, nrow=4, ncol=length(q))
  for(k in 0:3) {
    z <- gamma(k + 0.5) * sqrt(4*k + 1)/(gamma(k+1) * pi^(3/2) * sqrt(q))
    u <- (4*k + 1)^2/(16*q)
    y[k+1,] <- ifelse(u > -log.eps, 0, z * exp(-u) * besselK(x = u, nu=1/4))
  }
  return(apply(y,2,sum))
}

effectiveSize_crosschain <- function(x) {
  if (!is.mcmc.list(x)) stop('x must be an mcmc.list')
  chains <- nchain(x)
  iters <- niter(x)
  vars <- nvar(x)
  varnames <- colnames(x[[1]]) 
  
  M_fail <- iters > 2
  M <- iters
  while (M_fail) {
    M_test <- M
    while (M_test%%2 == 0) M_test <- M_test/2
    while (M_test%%3 == 0) M_test <- M_test/3
    while (M_test%%5 == 0) M_test <- M_test/5
    if (M_test<=1) M_fail <- FALSE
    else M <- M + 1
  }
  
  acs <- lapply(x,function(c){
    apply(c,2,function(y) {
      padded <- c(y - mean(y), rep.int(0, M*2 - iters))
      xf <- fft(padded)
      ac <- fft(Conj(xf) * xf, inverse = TRUE)
      return(Re(ac)[1:iters]/(iters^2 * 2))
    })
  })
  acs <- array(unlist(acs), dim=c(iters, vars, chains), dimnames = list(NULL, varnames, NULL))
  
  chain_means <- lapply(x, apply, 2, mean)
  mean_vars <- apply(acs[1,,],1,mean) * iters / (iters-1)
  between_vars <- apply(do.call(rbind,chain_means),2,var)
  var_hats <- between_vars + mean_vars * (iters-1) / iters
  
  rho_hats_t <- 1 - sweep(-sweep(apply(acs, c(1,2), mean), 2, mean_vars),2,var_hats,'/')
  rho_hats_t[1,] <- rep(1,vars)
  
  rho_hats_even <- rho_hats_t[1:iters %% 2 == 1,]
  rho_hats_odd <- rho_hats_t[1:iters %% 2 == 0,]
  if (iters %% 2 == 1) rho_hats_odd <- rbind(rho_hats_odd, rep(0,vars))
  rho_hats_pair <- rho_hats_even + rho_hats_odd
  
  prior_pairs <- rho_hats_pair[1,]
  for (i in 2:(ceiling(iters/2)-2)) {
    rho_hats_pair[i,] <- pmin(prior_pairs,rho_hats_pair[i,])
    rho_hats_pair[i,] <- pmax(rho_hats_pair[i,], rep(0, vars))
    prior_pairs <- rho_hats_pair[i,]
  }
  rho_hats_pair[ceiling(iters/2)-1,] <- rep(0,vars)
  rho_hats_pair[ceiling(iters/2),] <- rep(0,vars)
  
  for (j in 1:vars) {
    for (i in 2:ceiling(iters/2)-1) {
      if (rho_hats_pair[i,j]==0) {
        rho_hats_pair[i,j] <- ifelse( rho_hats_even[i,j] > 0, rho_hats_even[i,j]/2, 0)
        break
      }
    }
  }
  
  tau_hats <- -1 + 2*apply(rho_hats_pair, 2, sum)
  tau_hats <- pmax(tau_hats, 1/log10(iters*chains))
  neffs <- iters * chains / tau_hats
  neffs
}



















