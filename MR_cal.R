# ----------------------------------------------
# Helper function: Function to minimize misclassified cost
# ----------------------------------------------
dfs = function (op, used, cost, n, K) {
  min = Inf
  for (i in 1:K)
    if (used[i] == FALSE) {
      op[n] = i
      used[i] = TRUE
      if (n == K) {
        total = 0
        for (j in 1:K)
          total = total + cost[j,op[j]]
        min = min(min,total)
      }else
        min = min(min, dfs(
          op = op, used = used, cost = cost, n = n + 1, K = K
        ))
      op[n] == 0;
      used[i] = FALSE
    }
  return(min)
}

# ----------------------------------------------
# Function to Calculate the misclassified rate
# ----------------------------------------------
MR_cal = function (C, size = rep(length(C) / K, K), N = sum(size),
                   K = length(size), std = rep(c(1:K), size)) {
  # - C     vector of grouping result
  # - std   standard communities
  # - size  vector of size of blocks
  # - N     number of nodes
  # - K     number of blocks
  
  # cnt[k,l] - cost of assigning k to position of l
  cnt = matrix(0,K,K)
  for (k in 1:K)
    for (l in 1:K)
      cnt[k,l] = length(which(C[which(std == l)] != k &
                                C[which(std == l)] != 0))
    
#     dfs(
#       op = rep(0,K), used = logical(K), cost = cnt, n = 1, K = K
#     ) / N
    
    result = 0
    for (k in 1:K) {
      result = result + min(cnt)
      Index = match(min(cnt), cnt)
      cnt[ ,Index %/% K + (Index %% K > 0) ] = rep(Inf, K)
      cnt[ Index %% K + K * (Index %% K == 0), ] = rep(Inf, K)
    }
    result / N
    
}

# MR_cal = function (C, size, N = sum(size), K = length(size)) {
#   # Calculate the misclassified rate
#   # - C     vector of grouping result
#   # - size  vector of size of blocks
#   # - N     number of nodes
#   # - K     number of blocks
#
#   cnt = matrix(0,K,K)
#   for (k in 1:K) {
#     for (l in 1:K) {
#       if (l == 1)
#         head = 1
#       else
#         head = sum(size[1:(l - 1)]) + 1
#       tail = sum(size[1:l])
#       cnt[k,l] = length(which(C[head:tail] != k))
#     }
#   }
#
#   return(dfs(
#     op = rep(0,K), used = logical(K), cost = cnt, n = 1, K = K
#   ) / N)
# }
#
# MR_cal_sp = function(c, st, K = max(c), times = max(10, K)) {
#   min = Inf
#   tran = rep(0, K)
#   for (i in 1:times) {
#     perm = sample(c(1:K), K)
#     for (j in 1:K)
#       tran[which(c == j)] = perm[j]
#     min = min(min, length(which(tran != st & st != 0)) / length(st))
#   }
#   min
# }


# MR_cal(C = c(1,2,3,2,2,3,1,2,2,1,3,3), size = rep(4,3))
