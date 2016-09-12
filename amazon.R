setwd('~/R')
item = read.table('group.txt', header = T)
edge = read.table('similar.txt', header = T)
revw = read.table('review.txt', header = T)

nodes = intersect(item$PID, revw$PID)
nNodes = length(nodes)
from = item$PID[match(edge$from, item$ASIN)]
from = match(from, nodes)
to = item$PID[match(edge$to, item$ASIN)]
to = match(to, nodes)
edges = na.omit(cbind(from, to))

library(igraph)
edges = graph(edges = t(edges))
adj = get.adjacency(edges)
cov = cbind(as.numeric(as.factor(item$GROUP))[nodes], revw$TOTAL, revw$AVG)

setwd('~/R')
source('VEM_SBM.R')
nBlocks = 10
VEM_SBM(As = adj, K = nBlocks, iter = 50, X = cov, is.cov = T, cov.type = 'norm', init.type = 'casc')