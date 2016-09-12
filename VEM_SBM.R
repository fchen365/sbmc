safe.tau = function(beta, alpha,
                    f1 = NA, P1 = NA,
                    is.cov, cov.type) {
  if (is.cov) {
    if (cov.type == 'norm')
      tem1 = t(t(beta * f1) * alpha)
    else if (cov.type == 'bin')
      tem1 = t(t(beta * P1) * alpha)
  }
  else
    tem1 = t(t(beta) * alpha)
  tem2 = apply(tem1, 1, sum)
  wrong.list = which(tem2 == 0)
  for (i in wrong.list) {
    if (is.cov) {
      if (cov.type == 'norm') {
        tem1[i,] = f1[i,] * alpha
        tem2[i] = f1[i,] %*% alpha
      } else if (cov.type == 'bin') {
        tem1[i,] = P1[i,] * alpha
        tem2[i] = P1[i,] %*% alpha
      }
    } else {
      tem1[i,] = rep(1 / length(tem1[i,]), length(tem1[i,])) + runif(n = length(tem1[i,]), min = -0.1, max = 0.1)
      tem2[i] = sum(tem1[i,])
    }
  }
  
  return(tem1 / tem2)
  
}

# ----------------------------------------------
# Variational EM regarding or disregarding covariates
# ----------------------------------------------

VEM_SBM = function (As, K, iter = 75,
                    X = NA, is.cov = FALSE, cov.type = 'norm',
                    init.type = 'scp') {

  # - As        (sparse) adjacency matrix
  # - K         number of communities
  # - X         covariate data
  # - iter      number of iteration
  # - is.cov    logical, indicating whether to take covariates into consideration
  # - cov.type  type of covariates data, 'norm' for normal dist., 'bin' for bernoulli
  # - inittype  initial label method, possible value: 'scp', 'casc'...
  
  if (is.cov) {
    if (is.na(X)[1])
      stop('Covariates data can not be omitted.')
    if (cov.type == 'bin' && all(which(X == 1 | X == 0)) == FALSE)
      stop('Covariates data must be either 0 or 1 in case of binary type.')
  }
  
  #   sparse matrix case
  #   library(Matrix)
  #   As = Matrix(As, sparse = T)
  
  # abtain initial label
  setwd("~/R")
  # source('spectralClusteringMethods.R')
  source('initLabel.r')
  init = initLabel(
    As = As, K = K, X = X, cov.type = cov.type, type = init.type
  )
  Z = init$label
  
  time1 = Sys.time()
  
  n = dim(As)[1]
  R = dim(X)[2]
  epsilon = 1e-2
  
  # initial tau according to initial label
  tau = matrix(0.2 / K, n, K)
  for (i in 1:n)
    tau[i,Z[i]] = tau[i,Z[i]] + 0.8
  
  for (iterate in 1:iter) {
    # M-Step
    
    # pi-update
    pi_1 = matrix(0, K, K)
    for (k in 1:K)
      for (l in 1:K) {
        pi_1[k,l] = as.numeric(t(tau[,k]) %*% As %*% tau[,l] / sum(tau[,k] %o% tau[,l]))
      }
    
    # alpha-update
    alpha = apply(tau, 2, sum) / sum(tau);
    
    if (is.cov)
      if (cov.type == 'norm') {
        # mu-update
        mu = t(tau) %*% X / apply(tau, 2, sum)
        
        # Sigma-update
        w = t(t(tau) / apply(tau, 2, sum))
        Sigma = list()
        for (k in 1:K) {
          Sigma[[k]] = matrix(0, R, R)
          for (i in 1:n)
            Sigma[[k]] = Sigma[[k]] + w[i,k] * (X[i,] - mu[k,]) %*% t(X[i,] - mu[k,])
        }
      } else if (cov.type == 'bin') {
        # mu-update
        mu = t(tau) %*% X / (t(tau) %*% X + apply(tau, 2, sum))
      }
    
    # E-step
    
    # beta-update
    A = As - diag(diag(As))
    eta = log(pi_1 / (1 - pi_1))
    rho = log(1 - pi_1)
    beta = exp(t(t(A %*% tau %*% t(eta)) + as.numeric(rep(1,n) %*% tau %*% t(rho))))
    # beta = exp( A %*% tau %*% t(eta) + rep(1,n) %*% tau %*% t(rho) )
    
    # density function f1-update
    if (is.cov)
      if (cov.type == 'norm') {
        library(mnormt)
        f1 = matrix(0, n, k)
        for (k in 1:K)
          f1[,k] = dmnorm(x = X, mean = mu[k,], varcov = Sigma[[k]])
        
      } else if (cov.type == 'bin') {
        P1 = matrix(0, n, k)
        for (k in 1:K)
          P1[,k] = apply(
            X = X %*% diag(plogis(mu[k,]), nrow = R)
            + (1 - X) %*% diag(1 - plogis(mu[k,]), nrow = R),
            MARGIN = 1, FUN = prod
          )
        
      }
    
    # tau-update
    tau1 = tau
    tau = safe.tau(
      beta = beta, alpha = alpha,
      f1 = f1, P1 = P1,
      is.cov = is.cov, cov.type = cov.type
    )
    
    if (max(abs(tau1 - tau)) <= epsilon)
      break
  }
  
  time2 = Sys.time()
  
  fun = function(x) {
    return(which(x == max(x))[1])
  }
  result = list(
    runtime = as.numeric(time2 - time1), #+ init$runtime,
    itertimes = iterate,
    groups = apply(tau, 1, fun), initLabel = Z,
    tau = tau, alpha = alpha, pi_1 = pi_1
  )
  if (is.cov)
    if (cov.type == 'norm') {
      result$f1 = f1
      result$mu = mu
      result$Sigma = Sigma
    } else if (cov.type == 'bin') {
      result$P1 = P1
    }
  
  result
}

# M = gen_demo(1500, K = 3, R = 3, cov.type = 'bin')
# res1 = VEM_SBM(As = M$adj, K = M$blocks, is.cov = FALSE, init = 'scp', iter = 15)
# res1$groups
# res2 = VEM_SBM(As = M$adj, K = M$blocks, X = M$Cov, is.cov = TRUE, init = 'casc', iter = 10)
# res2$groups
# res3 = VEM_SBM(As = M$adj, K = M$blocks, X = M$Cov, is.cov = TRUE, cov.type = 'bin', init = 'scp')
# res3$groups