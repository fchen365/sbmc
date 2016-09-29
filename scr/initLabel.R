# ---------------------------------------------------------------------
# Generate initial labeling
# ---------------------------------------------------------------------

initLabel = function (As, K, 
                      X = NA, cov.type = 'norm', 
                      type = 'scp') {

  # - As    (sparse) adjacency matrix
  # - K     number of communities
  # - type  'scp' spectral clustering with perturbations, 
  #         'casc' refers to Covariates assisted spectral clustering
  
  if (type == 'casc') {  
    
    if (is.na(X)[1])
      stop('Covariates data can not be omitted.')
    if (cov.type == 'bin' && all(which(X == 1 | X == 0)) == FALSE)
      stop('Covariates data must be either 0 or 1 in case of binary type.')

    setwd("~/R")
    
    time1 = Sys.time()
    
#     source('spectralClusteringMethods.R')
#     if (cov.type == 'norm')
#       method = 'adjacency'
#     else if (cov.type == 'bin') 
#       method = 'regLaplacian'
#     label = getCascClusters(
#       adjacencyMat = As, nBlocks = K, method = 'adjacency',
#       covariates = X, hTuningParam = TRUE
#     )
    
    source('casc.R')
    label = casc(adjMat = As, covMat = X, nBlocks = K, nIter = 10)$cluster
    
    time2 = Sys.time()
    
    result = list(runtime = as.numeric(time2 - time1), label = label)
    return(result)
  }
  else if (type == 'scp') {
    rhoPert = 0.25
    degPert = 0.01
    
    # sparse matrix case
    library(Matrix)
    degh = Matrix()
    As = Matrix(As, sparse = T);
    
    # As = as.matrix(As)
    n = dim(As)[1]
    avgDeg = mean(apply(As, 1, sum))
    
    time1 = Sys.time()
    alpha0 = rhoPert * avgDeg;
    
    # degh = Matrix(sum(As, 2) + alpha0, sparse = T)
    degh = apply(As, 1, sum) + alpha0
    Gh = diag(degh ^ (-0.5))
    Gh[which(is.infinite(Gh))] = 0
    
    # bh = Matrix(Gh %*% rep(1, n), sparse = T)
    bh = Gh %*% rep(1, n)
    bhn = (alpha0 / n) * bh
    Lh = Gh %*% As %*% Gh
    
    if (dim(As)[1] == K)
      kmIDX = c(1:K)
    else if (dim(As)[1] == K + 1) {
      U = eigen(x = Lh, symmetric = T)$vectors
      kmIDX = kmeans(x = U[, 2:K], centers = K, nstart = max(K, 10))$cluster
    } else {
      func = function(x, extra = NULL) {
        as.numeric(Lh %*% x + as.numeric(t(bh) %*% x) * bhn)
      }
      
      library(igraph)
      U = arpack(
        func = func, sym = TRUE, options = list(
          n = n, nev = K, ncv = min(K + 2, n - 1), which = "LM", maxiter = 1000000
        )
      )$vectors
      
#       chooseCRANmirror(ind = 95)
#       install.packages('rARPACK')
#       library(rARPACK)
#       U = eigs(A = func, k = K, which = 'LM', n = n)$vectors
      
      kmIDX = kmeans(x = U[, 2:K], centers = K, nstart = max(K, 10))$cluster
    }
    
    time2 = Sys.time()
    
    result = list(runtime = as.numeric(time2 - time1), label = kmIDX)
    return(result)
  }
}

# A = gen_demo(600, 3, sym = T)
# init = initLabel(As = M$adj, K = M$blocks, type = 'scp')
# init = initLabel(As = M$adj, K = M$blocks, X = M$Cov, type = 'casc')
