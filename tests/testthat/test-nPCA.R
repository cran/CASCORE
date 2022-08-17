context("nPCA")

test_that("nPCA stops when it should", {
  expect_error( runningmean(0, c(0,0)) )
})

set.seed(2022)
n = 10; K = 2
P  = matrix(c(1/2, 1/4, 1/4, 1/2), byrow = TRUE, nrow = K)
distribution = c(1, 2)
l = sample(distribution, n, replace=TRUE, prob = c(1/2, 1/2))
Pi = matrix(0, n, 2)
for (i in 1:n){
  Pi[i, l[i]] = 1
  }
Omega = Pi %*% P %*% t(Pi)
Adj = matrix(runif(n*n, 0, 1), nrow = n)
Adj = Adj - Omega
temp = Adj
Adj[which(temp > 0)] = 0
Adj[which(temp <= 0)] = 1
diag(Adj) = 0
Adj[lower.tri(Adj)] = t(Adj)[lower.tri(Adj)]

test_that("This function returns a list of predicted membership of all nodes", {
  expect_length(nPCA(Adj, K), n)
  expect_length(unique(nPCA(Adj, K)), K)
})
