# ----------------------------------------------
# Function to read ego network data
# ----------------------------------------------

EgoReader = function(egoID, threshold = 3) {
  # egoID - ID of ego center node
  # threshold - treshold of biggest size of circle to be emit (unused)
  
  library(igraph)
  
  # filename
  filename = as.character(egoID)
  fn.e = paste(filename, '.edges', sep = '')
  fn.c = paste(filename, '.circles', sep = '')
  fn.ef = paste(filename, '.egofeat', sep = '')
  fn.f = paste(filename, '.feat', sep = '')
  
  # vertices data
  ego.vertices = as.matrix(read.table(fn.e)) + 1
  ego.circles = strsplit(x = readLines(fn.c), split = '\t')
  for (j in 1:length(ego.circles))
    ego.circles[[j]] = as.numeric(ego.circles[[j]][-1]) + 1
  ego.circles = unlist(ego.circles)
  ego.vertices = sort(unique(c(
    as.numeric(ego.vertices), egoID + 1, ego.circles
  )))
  miss = setdiff(c(1:max(ego.vertices)), ego.vertices)
  
  # edges data and adjacency matrix
  ego.edges = read.graph(fn.e, format = 'edgelist', directed = F)
#   for (j in ego.vertices)
#     if (j != egoID + 1)
#       ego.edges = ego.edges + edge(egoID + 1, j)
  adj = as.matrix(get.adjacency(ego.edges))
  adj = adj[-miss,-miss]
  adj[which(adj != 0)] = 1
  
  # circles data
  ego.circles = strsplit(x = readLines(fn.c), split = '\t')
  for (j in 1:length(ego.circles)) {
    ego.circles[[j]] = as.numeric(ego.circles[[j]][-1]) + 1
    ego.circles[[j]] = match(ego.circles[[j]], ego.vertices)
  }
  
  # feature data
  feat = rbind(cbind(egoID, as.matrix(read.table(fn.ef))),
               as.matrix(read.table(fn.f)))
  feat = feat[as.numeric(na.omit(match(ego.vertices, feat[,1]+1))),]
  feat = feat[order(feat[,1]), 2:dim(feat)[2]]
  
  # generate covariates data using PCA
  # library(psych)
  # nCov = fa.parallel(feat, n.iter = 5)$ncomp # parallel analysis
  # pr = principal(feat, nfactors = nCov)
  pr = princomp(feat)
  nCov = length(which(pr$sdev > 3.5 * mean(pr$sdev)))
  Cov = feat %*% pr$loadings[, 1:nCov]
  
  result = list(adj = adj,
                Cov = Cov,
                cir = ego.circles)
  result
}
