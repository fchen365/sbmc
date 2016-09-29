# ----------------------------------------------
# Helper function: Generating a random positive-definite matrix with user-specified positive eigenvalues
# ----------------------------------------------
PosDef = function (n, ev = runif(n, 0, n / 10)) {
  # If eigenvalues are not specified, they are generated from a uniform distribution
  # n - dimension of matrix
  # ev - eigenvectors
  
  Z = matrix(rnorm(n ^ 2), n, n)
  decomp = qr(Z)
  Q = qr.Q(decomp)
  R = qr.R(decomp)
  d = diag(R)
  ph = d / abs(d)
  O = Q %*% diag(ph)
  Z = t(O) %*% diag(ev) %*% O
  return(Z)
}

# pdmat = PosDef(n=5, ev=1:5)
# eigen(pdmat)$val

# ----------------------------------------------
# Helper function: Generating n random orthogonal vectors
# ----------------------------------------------

Orth = function (n) {
  Z = matrix(rnorm(n ^ 2), n, n)
  decomp = qr(Z)
  Q = qr.Q(decomp)
  return(Q)
}

# ----------------------------------------------
# Helper function: Mean vectors for continuous covariates (normal) for each blocks,
# or location parameters of binary data for each blocks
# ----------------------------------------------

getMV = function(K, R, angle = pi / 4) {
  
  # K - number of verter to generate
  # R - dimension of vertor
  
  if (K <= R) {
    if (K == 1)
      mu = matrix(
        data = c(1, rep(0,R - 1)),
        nrow = K, ncol = R, byrow = T
      ) %*% Orth(R)
    else if (K == 2)
      mu = matrix(
        data = c(1, rep(0,R - 1), cos(angle), sin(angle), rep(0, R - 2)),
        nrow = K, ncol = R, byrow = T
      ) %*% Orth(R)
    else if (K == 3)
      mu = matrix(
        data = c(
          1, rep(0, R - 1),
          cos(angle), sin(angle), rep(0, R - 2),
          cos(angle), sin(angle) * cos(angle) / (1 + cos(angle)), sqrt((cos(angle) - cos(2 * angle)) / (1 + cos(angle))),
          rep(0,R - 3)
        ),
        nrow = K, ncol = R,byrow = T
      ) %*% Orth(R)
    else
      mu = Orth(R)[1:K,]
  } else {
    mu = matrix(0,K,R)
    for (k in 1:K)
      mu[k,] = runif(R, 2 * k - 2, 2 * k - 1)
  }
  
  mu
}

# ----------------------------------------------
# Helper function: Generation of one SBM network (and Covariate data)
# ----------------------------------------------

gen_demo = function (n, K, size = NA, R = NA, cov.type = 'norm',
                     sym = TRUE, # maxiter = 100,
                     ratio = 0.25, dif = 0, ExpDeg = 7, angle = pi / 4) {

  # - n       number of nodes
  # - K       number of communities
  # - size    vector of sizes of each block, the default is equal size for each one
  # - sym     type of adjacency matrix, 'FALSE' (default) refers to normal adjacency matrix
  # - R       numerical integer indicates number of Covariates dimension
  # - ratio   ratio of connection probability, i.e, q over p
  # - matiter max iteration times of generation of demo data
  # - dif     proportion of different covariates in different blocks
  # - angle   the angle between each mean vector
  
  if (is.na(size)) {
    if (n %% K != 0) 
      stop('Wrong size!')
    size = rep(n / K, K)
  }
  else if (n != sum(size) || length(size) != K) 
      stop('Wrong size!')
  
  # generate blocks assignment matrix Z
  Z = diag(K) %x% rep(c(1:K),size)
  for (k in 1:K)
    Z = Z[-which(Z[,k] != k & Z[,k] != 0),]
  Z[which(Z != 0)] = 1
  
  # adjacency matrix
  library(Matrix)
  p = ExpDeg / (n * (1 - ratio) / K + n * ratio) # diagonal entries
  P = matrix(ratio * p, K, K)
  diag(P) = rep(p,K)
  M = 1 * (matrix(runif(n * n), n, n) < Z %*% P %*% t(Z))
  if (sym)
    M[lower.tri(M)] = t(M)[lower.tri(M)]
  
  #   # Test isolated nodes
  #     for (i in 1:maxiter) {
  #       M = 1 * (matrix(runif(n * n), n, n) < Z %*% P %*% t(Z))
  #       if (sym)
  #         M[lower.tri(M)] = t(M)[lower.tri(M)]
  #       if (length(which(apply(M,1,sum) * apply(M,2,sum) == 0)) == 0)
  #         break;
  #     }
  #     if (i == maxiter) {
  #       print('failed.')
  #       return()
  #     }
  
#   # Plot adjacency matrix
#   library(igraph)
#   plot(graph.adjacency(M),
#        vertex.color=rep(c(1:K), size), vertex.size=3, vertex.label=NA,
#        edge.width=0.01, edge.arrow.size=0.05)
#   title(paste('in-out ratio = ', ratio))
  
  # Covariates data
  if (!is.na(R))
    if (cov.type == 'norm') {
      Sigma = PosDef(R) # fix the same covariance matrix for each block
      mu = getMV(K = K, R = R, angle = angle) # Mean vector
      library(mnormt)
      # set.seed(123)
      X = matrix(ncol = R)
      for (k in 1:K) {
        X = rbind(X, rmnorm(size[k] - round(size[k] * dif), mu[k,], Sigma))
        X = rbind(X, rmnorm(round(size[k] * dif), mu[k %% K + 1,], Sigma))
      }
      X = X[-1,]
    } else if (cov.type == 'bin') {
      range = R
      mu = range * getMV(K = K, R = R, angle = angle)
      X = matrix(ncol = R)
      for (k in 1:K) {
        X = rbind(X,
                  matrix(
                    data = as.numeric(runif((
                      size[k] - round(size[k] * dif)
                    ) * R) < plogis(mu[k,])),
                    nrow = size[k] - round(size[k] * dif), ncol = R, byrow = T
                  ))
        X = rbind(X,
                  matrix(
                    data = as.numeric(runif(round(size[k] * dif) * R) < plogis(mu[k %% K + 1,])),
                    nrow = round(size[k] * dif), ncol = R, byrow = T
                  ))
      }
      X = X[-1,]
    }
  
  result = list(
    adj = M,
    nodes = n,
    size = size,
    blocks = K,
    cov.type = cov.type
  )
  if (!is.na(R))
    result$Cov = X
  
  return(result)
}

# A = gen_demo(n = 600, K = 3, R = 3, sym = T)
# M = gen_demo(n = 600, K = 3, R = 3, ratio = 0.5, sym = T)
# M = gen_demo(n = 600, K = 3, R = 3, dif = 0.2, sym = T)
# M = gen_demo(n = 600, K = 3, R = 3, ExpDeg = 40, sym = T)
# M = gen_demo(n = 600, K = 3, R = 3, angle = pi / 4, sym = T)
# M = gen_demo(n = 600, K = 3, R = 3, ratio = 0.5, dif = 0, ExpDeg = 7, angle = pi / 4, sym = T)
# M = gen_demo(n = 600, K = 3, R = 3, cov.type = 'bin', sym = T)