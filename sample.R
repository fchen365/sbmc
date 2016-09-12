A = diag(c(1:5))
print(A)
diag(A)

source('~/R/VEM_SBM.R')

setEPS()
postscript('sample.eps')
plot(diag(A), main = 'sample: 1->5')
dev.off()
