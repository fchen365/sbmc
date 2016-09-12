# ----------------------------------------------
# Function to return blockmodel results
# ----------------------------------------------
getRes = function(Ego, CG, K, is.plot = FALSE) {
  # As - adjacency matrix
  # X - covariates data
  # K - number of blocks
  
  setwd('D:/University/Scholarship/CSST/R')
  source('initLabel.r')
  source('VEM_SBM.R')
  source('MR_cal.R')
  
  nMet = 5
  nVer = dim(Ego$adj)[1]
  nCir = length(CG$cir)
  
  # Standard communities
  std = rep(0, nVer)
  CirLab = initLabel(As = CG$adj, K = K, type = 'scp')$label
  for (i in 1:nCir)
    std[CG$cir[[i]]] = CirLab[i]
  
  # Plot original ego network
  if (is.plot) {
    library(igraph)
    mygraph = graph.adjacency(adjmatrix = Ego$adj, mode = 'undirected', weighted = TRUE)
    plot.igraph(
      x = mygraph,
      vertex.label = NA, vertex.size = 4, vertex.color = std,
      edge.arrow.size = 0.5 # edge.width = sqrt(E(mygraph)$weight),
    )
    title('True Community Structure')
  }
  
  res = matrix(0, nMet, nVer)
  res[1,] = initLabel(As = Ego$adj, K = K,
                      type = 'scp')$label
  res[2,] = VEM_SBM(
    As = Ego$adj, K = K,
    is.cov = FALSE, init.type = 'scp'
  )$groups
  res[3,] = VEM_SBM(
    As = Ego$adj, K = K,
    X = Ego$Cov, is.cov = T,
    init.type = 'scp'
  )$groups
  res[4,] = initLabel(
    As = Ego$adj, K = K,
    X = Ego$Cov, type = 'casc'
  )$label
  res[5,] = VEM_SBM(
    As = Ego$adj, K = K,
    X = Ego$Cov, is.cov = TRUE,
    init.type = 'casc'
  )$groups
#   res[6,] = VEM_SBM(
#     As = Ego$adj, K = K,
#     X = Ego$Cov, is.cov = T, cov.type = 'bin',
#     init.type = 'scp'
#   )$groups
  
  MR = rep(0, nMet)
  for (i in 1:nMet)
    MR[i] = MR_cal(C = res[i,], std = std, K = K)
  MR
}