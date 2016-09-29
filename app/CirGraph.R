# ----------------------------------------------
# Function to generate circle graph. V1: Split the overlaping circles
# ----------------------------------------------

CirGraph = function(adj, cir) {
  # edges - adjacency matrix
  # cir - list of circles
  # threshold - upper threshold size of circles to be omitted
  
  # Split overlaping circles
  flag = TRUE
  nCir = length(cir)
  while (flag) {
    cnt = 0
    for (i in 1:(nCir - 1))
      for (j in (i + 1):nCir)
        if (length(intersect(cir[[i]], cir[[j]])) > 0)
          if (all(cir[[i]] %in% cir[[j]])) {
            cir[[j]] = setdiff(cir[[j]], cir[[i]])
          } else if (all(cir[[j]] %in% cir[[i]])) {
            cir[[i]] = setdiff(cir[[i]], cir[[j]])
          } else {
            cnt = cnt + 1
            cir[[nCir + cnt]] = intersect(cir[[i]], cir[[j]])
            cir[[i]] = setdiff(cir[[i]], cir[[nCir + cnt]])
            cir[[j]] = setdiff(cir[[j]], cir[[nCir + cnt]])
          }
    nCir = nCir + cnt
    if (cnt == 0)
      flag = FALSE
  }
  cir = cir[which(sapply(cir, length) > 0)]
  nCir = length(cir)
  
  # circle graph adjacency, emit isolated circles and then iGraph
  CG.adj = matrix(0, nCir, nCir)
  for (i in 1:(nCir - 1))
    for (j in (i + 1):nCir)
      CG.adj[i,j] = CG.adj[j,i] = sum(adj[cir[[i]], cir[[j]]])
  
  result = list(cir = cir,
                adj = CG.adj)
  result
}

# ----------------------------------------------
# Function to generate circle graph. V2: Directly remove the overlaping node
# ----------------------------------------------

CirGraph = function(adj, cir, threshold = 3) {
  # edges - adjacency matrix
  # cir - list of circles
  # threshold - upper threshold size of circles to be omitted
  
  # Split overlaping circles
  flag = TRUE
  nCir = length(cir)
  while (flag) {
    cnt = 0
    for (i in 1:(nCir - 1))
      for (j in (i + 1):nCir)
        if (length(intersect(cir[[i]], cir[[j]])) > 0)
          if (all(cir[[i]] %in% cir[[j]])) {
            cir[[j]] = setdiff(cir[[j]], cir[[i]])
          } else if (all(cir[[j]] %in% cir[[i]])) {
            cir[[i]] = setdiff(cir[[i]], cir[[j]])
          } else {
            cnt = cnt + 1
            intersection = intersect(cir[[i]], cir[[j]])
            if (sample(c(0,1),1) == 0)
              cir[[i]] = setdiff(cir[[i]], intersection)
            else
              cir[[j]] = setdiff(cir[[j]], intersection)
          }
    if (cnt == 0)
      flag = FALSE
  }
  cir = cir[which(sapply(cir, length) > 0)]
  nCir = length(cir)
  
  # Omit tiny circles
  cir = cir[which(sapply(cir, length) > threshold)]
  nCir = length(cir)
  
  # circle graph adjacency, emit isolated circles and then iGraph
  CG.adj = matrix(0, nCir, nCir)
  for (i in 1:nCir)
    for (j in i:nCir)
      CG.adj[i,j] = CG.adj[j,i] = sum(adj[cir[[i]], cir[[j]]])
  
  cirSize = sapply(cir, length)
  nVer = length(unique(unlist(cir)))
  r.in = mean(diag(CG.adj) / cirSize)
  r.out = mean((apply(CG.adj, 1, sum) - diag(CG.adj)) / (nVer - cirSize))
  ratio = r.out / r.in
  diag(CG.adj) = 0
  
  result = list(cir = cir,
                adj = CG.adj, 
                ratio = ratio)
  result
}

# # Detect containing relationship
# omit = c()
# for (i in 1:nCir)
#   for (j in 1:nCir)
#     if (i != j && all(cir[[i]] %in% cir[[j]]))
#       omit = c(omit, i)
# cir = cir[unique(omit)]
# nCir = length(cir)
#

# # Omit isolated circles
# cirSize = sapply(cir, length)
# iso = which(apply(CG.adj, 1, sum) == 0 & circsize < threshold2)
# cir = cir[-iso]
# nCir = length(cir)
# cirSize = cirSize[-iso]
# CG.adj = CG.adj[-iso,-iso]
# nCir = length(cir)
# CG = graph.adjacency(adjmatrix = CG.adj, mode = 'undirected', weighted = TRUE)
