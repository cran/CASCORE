#' @title Network-based Clustering.
#' @description \emph{Network-based Clustering} is a spectral clustering method that focuses
#'   solely on the topological structure of a network. It employs \code{k-means} on the first
#'   \eqn{K} leading eigenvectors of the weighted adjacency matrix of a graph, with each
#'   eigenvector normalized to have unit magnitude.
#' @param Adj A 0/1 adjacency matrix.
#' @param tau An optional tuning parameter, the default value is the mean of adajacency matrix.
#' @param K A positive integer, indicating the number of underlying communities in
#'   graph \code{Adj}.
#' @param itermax \code{k-means} parameter, indicating the maximum number of
#'   iterations allowed. The default value is 100.
#' @param startn \code{k-means} parameter. If centers is a number, how many
#'   random sets should be chosen? The default value is 10.
#' @return A label vector.
#'
#' @importFrom stats kmeans runif
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
#' Net_based(Adj, 2)
#' @export

Net_based <- function(Adj, K, tau = NULL, itermax = NULL, startn = NULL){
  if(!isSymmetric(Adj)) stop("Error! Adj is not symmetric!")
  if(K > dim(Adj)[1]) stop("Error! More communities than nodes!")
  if(K %% 1 != 0) stop("Error! K is not an integer!")
  if(K <= 0) stop("Error! Nonpositive K!")

  if(is.null(tau)) tau = mean(Adj);

  n <- dim(Adj)[1]
  A_tau = Adj + tau * matrix(1, n, n)/n
  s = rowSums(A_tau)
  s =  s^(-1/2)
  S = diag(s)
  Z = S %*% A_tau %*% S
  #D <- diag(rowSums(A_tau))
  #L_tau <- (D + tau*J/n)^{-1/2} %*% Adj %*% (D + tau*J/n)^{-1/2}
  #L <- (D + tau*diag(n))^{-1/2} %*% Adj %*% (D + tau*diag(n))^{-1/2}
  g.eigen <-  eigen(Z)
  R = g.eigen$vectors
  R = R[, 1: K]
  R <- t(apply(R, 1, function(x) x/sqrt(sum(x^2))))

  # apply Kmeans to assign nodes into communities
  if(!is.null(itermax) & !is.null(startn)){
    result = kmeans(R, K, iter.max = itermax, nstart = startn) #apply kmeans on ratio matrix
  }
  if(!is.null(itermax) & is.null(startn)){
    result = kmeans(R, K, iter.max = itermax, nstart = 10) #apply kmeans on ratio matrix
  }
  if(is.null(itermax) & !is.null(startn)){
    result = kmeans(R, K, iter.max = 100, nstart = startn) #apply kmeans on ratio matrix
  }
  else{
    result = kmeans(R, K, iter.max = 100, nstart = 10) #apply kmeans on ratio matrix
  }

  est = as.factor(result$cluster)
  return(est)
}
