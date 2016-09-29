source('~/R/plot_MR.R')

# plot.MR(N = 1500, K = 3, R = 3, type = 'ratio', times = 25, start = -0.03, pace = 0.04, step = 18)
# plot.MR(N = 1500, K = 3, R = 3, type = 'dif', times = 35, pace = 0.1, step = 8)
# plot.MR(N = 1500, K = 3, R = 3, type = 'ExpDeg', times = 35, start = 0, pace = 3, step = 10)
# plot.MR(N = 1500, K = 3, R = 3, type = 'angle', times = 35, start = -pi/12, pace = pi/12, step = 7)
# plot.MR(N = 1500, K = 3, R = 3, type = 'nCov', times = 35, start = 0, pace = 2, step = 8)
# plot.MR(N = 1500, K = 3, R = 3, type = 'ratio', times = 25, start = -0.03, pace = 0.04, step = 18, cov.type = 'bin')
plot.MR(N = 1500, K = 3, R = 3, type = 'nNode', cov.type = 'norm', times = 10, start = 2, pace = 0.5, step = 5)
