#' @title Covariate Assisted Spectral Clustering.
#' @description \emph{CASC} clusters graph nodes by applying spectral clustering with the assistance from
#' node covariates.
#' @details \emph{CASC} is a community detection algorithm for networks with node covariates, proposed
#'   in \emph{Covariate-assisted spectral clustering} of Binkiewicz, et al. (2017). \emph{CASC} applies
#'   \code{k-means} on the first \code{K} leading eigenvectors of the balanced matrix between the Laplacian
#'   matrix and the covariate matrix.
#' @param Adj A 0/1 adjacency matrix.
#' @param Covariate A covariate matrix. The rows correspond to nodes and the columns correspond to covariates.
#' @param K A positive integer, indicating the number of underlying communities in graph \code{Adj}.
#' @param alphan A tuning parameter to balance between the contributions of the graph and the covariates.
#' @param itermax \code{k-means} parameter, indicating the maximum number of iterations allowed.
#'   The default value is 100.
#' @param startn \code{k-means} parameter. If centers is a number, how many random sets should
#'   be chosen? The default value is 10.
#'
#' @return  \item{estall}{A lavel vector.}
#'
#' @importFrom pracma eig Norm
#' @importFrom stats kmeans runif
#'
#' @references Binkiewicz, N., Vogelstein, J. T., & Rohe, K. (2017). \emph{Covariate-assisted spectral clustering}.
#'   \emph{Biometrika, 104(2), 361-377.}\cr\doi{10.1093/biomet/asx008}\cr
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
#' caseno = 4; Nrange = 10; Nmin = 10; prob1 = 0.9; p = n*4;
#' Q = matrix(runif(p*K, 0, 1), nrow = p, ncol = K)
#' Q = sweep(Q,2,colSums(Q),`/`)
#' W = matrix(0, nrow = n, ncol = K);
#' for(jj in 1:n) {
#'   if(runif(1) <= prob1) {W[jj, 1:K] = Pi[jj, ];}
#'   else W[jj, sample(K, 1)] = 1;
#' }
#' W = t(W)
#' D0 = Q %*% W
#' X = matrix(0, n, p)
#' N = switch(caseno, rep(100, n), rep(100, n), round(runif(n)*Nrange+ Nmin),
#'   round(runif(n)* Nrange+Nmin))
#' for (i in 1: ncol(D0)){
#'   X[i, ] = rmultinom(1, N[i], D0[, i])
#' }
#' CASC(Adj, X, 2)
#' @export

CASC = function(Adj, Covariate, K, alphan = 5, itermax = 100, startn = 10){
  s = rowSums(Adj)
  s = s + mean(s)
  s =  s^(-1/2)
  S = diag(s)
  Z = S %*% Adj %*% S
  net.eigen = eigen(Z%*%Z)
  ca = Covariate %*% t(Covariate);
  ca.eigen = eigen(ca);
  alphalower = (net.eigen$values[K] - net.eigen$values[K+1])/ca.eigen$values[1];
  alphaupper = net.eigen$values[1]/(ca.eigen$values[K] - ca.eigen$values[K+1]);
  d = rep(0, alphan);
  alpha = seq(alphalower, alphaupper, length.out = alphan);
  est = matrix(0, alphan, dim(Adj)[1])

  for(ii in 1:alphan){
    casc.eigen = eigen(Z%*%Z + alpha[ii]*ca);
    U = casc.eigen$vectors[,1:K];
    Unorm = apply(U, 1, Norm);
    indu = which(Unorm > 0);
    U = U[indu, ]/Unorm[indu]
    result = kmeans(U, K, iter.max = itermax, nstart = startn);
    d[ii] = result$tot.withinss;
    est[ii, indu] = as.factor(result$cluster)
  }
  result = est[which.min(d), ]
  return(result)
}
