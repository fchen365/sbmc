library(igraph)

setwd('D:/University/Scholarship/CSST/R')
source('EgoReader.R')
source('CirGraph.R')
source('getRes.R')

egolist = c(0, 107, 348, 414, 686, 698, 1684, 1912, 3437, 3980)
nEgo = length(egolist) # number of egos
nMet = 5

MR = matrix(0, nEgo, nMet)
nVer = nBlk = nCov = nCir = ratio = ExpDeg = rep(0, nEgo)
nBlk = c(3, 8, 3, 3, 3, 2, 3, 6, 3, 2)
# nBlk = rep(4, 10)

times = 5

for (e in 1:nEgo) {
  setwd("~/R")
  Ego = EgoReader(egoID = egolist[e])
  CG = CirGraph(adj = Ego$adj, cir = Ego$cir)
  nVer[e] = dim(Ego$adj)[1]
  nCov[e] = dim(Ego$Cov)[2]
  nCir[e] = length(CG$cir)
  # nBlk[e] = sqrt(nCir[e]) %/% 1
  ratio[e] = CG$ratio
  ExpDeg[e] = mean(apply(Ego$adj, 1, sum))

  for (t in 1:times) {
    res = getRes(Ego = Ego, CG = CG, K = nBlk[e], is.plot = FALSE)
    MR[e,] = MR[e,] + res
  }
  MR[e,] = MR[e,] / times
}

print(cbind(
  nVer, nCov, nCir, nBlk,
  ratio = round(ratio, digits = 3), ExpDeg = round(ExpDeg, digits = 3),
  round(MR, digits = 3)
))
print(round(apply(MR, 2, mean), digits = 3))
names = switch(
  as.character(nMet),
  '5' = c('SCP', 'VEM_scp', 'VEM_Cov_scp', 'CASC', 'VEM_Cov_casc'),
  '6' = c(
    'SCP', 'VEM_scp', 'VEM_Cov_scp', 'CASC', 'VEM_Cov_casc', 'VEM_(bin)Cov_scp'
  )
)
boxplot(MR, col = 'bisque', names = names, ylab = 'Misclassified Rate')
title('Comparison of five methods')

# # plot the original graph
# mygraph = graph.adjacency(adjmatrix = Ego$adj, mode = 'undirected', weighted = TRUE)
# plot.igraph(
#   x = mygraph,
#   vertex.label = NA, vertex.size = 5, vertex.color = res,
#   edge.arrow.size = 0.5 # edge.width = sqrt(E(mygraph)$weight),
# )
# print('Please assign the number of blocks: ')
# nBlk[e] = as.numeric(readline())


# plot the weighted circle graph
# mygraph = graph.adjacency(adjmatrix = CG$adj, mode = 'undirected', weighted = TRUE)
# nodeSize = sapply(CG$cir, length)
# plot.igraph(
#   x = mygraph,
#   vertex.label = NA, vertex.color = CirLab, vertex.size = sqrt(nodeSize),
#   edge.arrow.size = 0.5 # edge.width = sqrt(E(mygraph)$weight),
# )
# print('Please assign the number of blocks: ')
# nBlk[e] = as.numeric(readline())
