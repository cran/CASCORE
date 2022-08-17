#' @title Ordinary Principle Component Analysis.
#' @description \emph{Ordinary Principle Component Analysis (oPCA)}, also known as spectral clustering
#'   on the adjacency matrix is a classical spectral clustering method that applies \code{k-means} on
#'   the first \eqn{K} leading (unit-norm) eigenvectors of the adjacency matrix of a graph.
#' @param Adj A 0/1 adjacency matrix.
#' @param K A positive integer, indicating the number of underlying communities in
#'   graph \code{Adj}.
#' @param itermax \code{k-means} parameter, indicating the maximum number of
#'   iterations allowed. The default value is 100.
#' @param startn \code{k-means} parameter. If centers is a number, how many
#'   random sets should be chosen? The default value is 10.
#' @return A label vector.
#'
#' @importFrom stats kmeans runif
#' @references Chung, F. R., & Graham, F. C. (1997). \emph{Spectral graph theory (Vol. 92)}
#'   \emph{American Mathematical Soc..}
#' @examples
#'
#' # Simulate the Network
#' n = 10; K = 2;
#' theta = 0.4 + (0.45-0.05)*(seq(1:n)/n)^2; Theta = diag(theta);
#' P  = matrix(c(0.8, 0.2, 0.2, 0.8), byrow = TRUE, nrow = K)
#' set.seed(2022)
#' l = sample(1:K, n, replace=TRUE); # node labels
#' Pi = matrix(0, n, K) # label matrix
#' for (k in 1:K){
#'   Pi[l == k, k] = 1
#' }
#' Omega = Theta %*% Pi %*% P %*% t(Pi) %*% Theta;
#' Adj = matrix(runif(n*n, 0, 1), nrow = n);
#' Adj = Omega - Adj;
#' Adj = 1*(Adj >= 0)
#' diag(Adj) = 0
#' Adj[lower.tri(Adj)] = t(Adj)[lower.tri(Adj)]
#' oPCA(Adj, 2)
#' @export

oPCA = function(Adj, K, itermax = 100, startn = 10){
  # Inputs:
  # 1) Adj: an n by n symmetric adjacency matrix whose diagonals = 0 and positive entries = 1.
  # 2) K: a positive integer which is no larger than n. This is the predefined number of communities.

  # Optional Arguments for Kmeans:
  # 1) itermax: the maximum number of iterations allowed.
  # 2) nstart: R will try startn different random starting assignments and then select the one with the lowest within cluster variation.

  # Outputs:
  # 1) a factor indicating nodes' labels. Items sharing the same label are in the same community.

  if(!isSymmetric(Adj)) stop("Error! Adj is not symmetric!")
  if(K > dim(Adj)[1]) stop("Error! More communities than nodes!")
  if(K %% 1 != 0) stop("Error! K is not an integer!")
  if(K <= 0) stop("Error! Nonpositive K!")

  o.eigen = eigen(Adj)
  O = o.eigen$vectors
  O = O[, 1:K]
  result = tryCatch({kmeans(O, K, iter.max = 100, nstart = 10)}, error = function(x)
  {kmeans(O, K, iter.max = 1000, nstart = 10, algorithm = "Lloyd")})
  est = as.factor(result$cluster)
  return(est)
}
