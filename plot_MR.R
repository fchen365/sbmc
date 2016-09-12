# ----------------------------------------------
# Plot the misclassified rate vs. number of blocks
# ----------------------------------------------
plot.MR = function (N, K, R = 3, cov.type = 'norm', type, times = 25, start = 0, pace, step) {

  #   - N      number of nodes
  #   - K      number of blocks
  #   - cov.type  type of covariates data
  #   - type   'ratio'   indicates comparison of in-out ratio,
  #            'dif'     means proportion of different covariates in different blocks,
  #            'ExpDeg'  indicates camparison of expected degree
  #            'nCov'    refers to numbers of covariates, must be larger than K
  #            'angle'   means angle of every two mean vector should be varied as comparison
  #            'nNode'   number of nodes will be varied and runtime will be compared
  #   - times  times to simulation
  #   - step   step of covariates # addition, required
  #   - pace   pace of covariates # addition, required
  
  setwd("~/R")
  source('gen_demo.r')
  source('initLabel.r')
  source('VEM_SBM.R')
  source('MR_cal.R')
  
  nMet = 5
  MR = rt = matrix(0, nMet, step) # misclassified rate
  
  # default parameters
  ratio = 0.3
  dif = 0
  ExpDeg = 7
  nCov = 3
  angle = pi / 4
  
  for (s in 1:step) {
    for (t in 1:times) {
      # generate simulation data
      if (type == 'ratio')
        M = gen_demo(
          n = N, K = K, sym = T, R = nCov, cov.type = cov.type,
          ratio = start + s * pace, dif = dif, ExpDeg = ExpDeg, angle = angle
        )
      else if (type == 'dif')
        M = gen_demo(
          n = N, K = K, sym = T, R = nCov, cov.type = cov.type,
          ratio = ratio, dif = start + s * pace, ExpDeg = ExpDeg, angle = angle
        )
      else if (type == 'ExpDeg')
        M = gen_demo(
          n = N, K = K, sym = T, R = nCov, cov.type = cov.type,
          ratio = ratio, dif = dif, ExpDeg = start + s * pace, angle = angle
        )
      else if (type == 'nCov')
        M = gen_demo(
          n = N, K = K, sym = T, R = max(start + s * pace, 3), cov.type = cov.type,
          ratio = ratio, dif = dif, ExpDeg = ExpDeg, angle = angle
        )
      else if (type == 'angle')
        M = gen_demo(
          n = N, K = K, sym = T, R = nCov, cov.type = cov.type,
          ratio = ratio, dif = dif, ExpDeg = ExpDeg, angle = start + s * pace
        )
      else if (type == 'nNode')
        M = gen_demo(
          n = 10 ^ (start + s * pace) - 10 ^ (start + s * pace) %% K, 
          K= K, sym = T, R = nCov, cov.type = cov.type,
          ratio = ratio, dif = dif, ExpDeg = ExpDeg, angle = angle
        )
      else
        stop('Wrong type!')
      
      # result comparison
      
      # initial lable: spectral clustering with perturbation
      res1 = initLabel(As = M$adj, K = M$blocks, type = 'scp')
      MR[1,s] = MR[1,s] + MR_cal(res1$label, size = M$size)
      rt[1,s] = rt[1,s] + res1$runtime
      
      # VEM disregarding covariates
      res2 = VEM_SBM(
        As = M$adj, K = M$blocks,
        X = M$Cov, is.cov = FALSE,
        init.type = 'scp'
      )
      MR[2,s] = MR[2,s] + MR_cal(C = res2$groups, size = M$size)
      rt[2,s] = rt[2,s] + res2$runtime
      
      # VEM with covariates over scp initial labels
      res3 = VEM_SBM(
        As = M$adj, K = M$blocks,
        X = M$Cov, is.cov = TRUE, cov.type = cov.type,
        init.type = 'scp'
      )
      MR[3,s] = MR[3,s] + MR_cal(C = res3$groups, size = M$size)
      rt[3,s] = rt[3,s] + res3$runtime
      
      # initial lable: covariates assisted spectral clustering
      res4 = initLabel(
        As = M$adj, K = M$blocks,
        X = M$Cov, type = 'casc', cov.type = cov.type
      )
      MR[4,s] = MR[4,s] + MR_cal(C = res4$label, size = M$size)
      rt[4,s] = rt[4,s] + res4$runtime
      
      # VEM with covariates over casc initial label
      res5 = VEM_SBM(
        As = M$adj, K = M$blocks,
        X = M$Cov, is.cov = TRUE, cov.type = cov.type,
        init.type = 'casc'
      )
      MR[5,s] = MR[5,s] + MR_cal(C = res5$groups, size = M$size)
      rt[5,s] = rt[5,s] + res5$runtime
    }
    MR[,s] = MR[,s] / times
    rt[,s] = rt[,s] / times
    print(s)
  }
  
  if (type == 'nNode')
    data = log(rt,10)
  else
    data = MR
  # Plot Misclassified Rate
  setEPS()
  xlab = type
  if (type == 'nNode')
    ylab = 'RT'
  else
    ylab = 'MR'
  postscript(
    paste(
      ylab, '-', type,
      ' (ratio=', ratio,
      ', times=', times,
      ', type=', cov.type,
      ').eps', sep = ''
    )
  )
  plot(
    c(start + pace, start + step * pace), c(min(data), max(data)), type = 'n',
    xlab = xlab, ylab = ylab, main = paste(ylab, '-', type)
  )
  for (i in 1:nMet)
    lines(seq(
      from = start + pace, to = start + step * pace, by = pace
    ), data[i,], col = i + 1)
  
  if (type == 'ExpDeg' || type == 'angle')
    location = 'topright'
  else if (type == 'ratio' || type == 'nCov')
    location = 'topleft'
  else if (type == 'dif' || type == 'nNode')
    location = 'bottomright'
  legend(
    location, title = "Line Type",c('SCP', 'VEM', 'VEM_Cov_SCP', 'CASC', 'VEM_Cov_CASC'),
    lty = rep(1, 5), col = c(2:6), cex = 0.75
  )
  dev.off()
  
  MR
  rt
  
}

# plot.MR(N = 600, K = 3, R = 3, type = 'ratio', times = 25, pace = 0.05, step = 15)
# plot.MR(N = 1500, K = 3, R = 3, type = 'ratio', times = 50, start = -0.03, pace = 0.04, step = 20)
# plot.MR(N = 600, K = 3, R = 3, type = 'dif', times = 25, pace = 0.1, step = 6)
# plot.MR(N = 600, K = 3, R = 3, type = 'ExpDeg', times = 25, pace = 5, step = 5)
# plot.MR(N = 600, K = 3, R = 3, type = 'angle', times = 25, pace = pi/12, step = 6)
# plot.MR(N = 600, K = 3, R = 3, type = 'nCov', times = 25, pace = 2, step = 5)
# plot.MR(N = 1500, K = 3, R = 3, type = 'nNode', cov.type = 'norm', times = 1, start = 2, pace = 0.5, step = 6)