#' @title Covariate Assisted Spectral Clustering on Ratios of Eigenvectors.
#' @description Using ratios-of-eigenvectors to detect underlying communities in networks with node covariates.
#' @details \emph{CASCORE} is fully established in \emph{Covariate-Assisted Community Detection
#'   on Sparse Networks} of Hu & Wang (2022). \emph{CASCORE} detects the latent community structure
#'   under the covariate assisted degree corrected stochastic block model (CADCSBM), and it allows
#'   the disagreement between the community structures indicated in the graph and the covariates,
#'   respectively. \code{K-means} is applied on the entry-wise ratios between first leading eigenvector
#'   and each of the other \eqn{K} leading eigenvectors of the combined matrix of the adjacency matrix and the
#'   covariate matrix, to reveal the underlying memberships.
#' @param Adj A 0/1 adjacency matrix.
#' @param Covariate A covariate matrix. The rows correspond to nodes and the columns correspond to covariates.
#' @param K A positive integer, indicating the number of underlying communities in graph \code{Adj}.
#' @param alpha A numeric vector, each element of which is a tuning parameter to weigh the covariate matrix.
#' @param alphan The number of candidates \eqn{\alpha}. The default number is 5.
#' @param itermax \code{k-means} parameter, indicating the maximum number of
#'   iterations allowed. The default value is 100.
#' @param startn \code{k-means} parameter. If centers is a number, how many
#'   random sets should be chosen? The default value is 10.
#' @return \item{estall}{A lavel vector.}
#'
#' @importFrom pracma eig Norm
#' @importFrom stats kmeans runif median
#'
#' @references Hu, Y., & Wang, W. (2022). \emph{Covariate-Assisted Community Detection on Sparse Networks}.
#'   \cr\doi{10.48550/arXiv.2208.00257}\cr
#'
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
#' CASCORE(Adj, X, 2)

#' @export

CASCORE = function(Adj, Covariate, K, alpha = NULL, alphan = 5, itermax = 100, startn = 10){
  # Inputs:
  # 1) Adj: an n by n symmetric adjacency matrix whose diagonals = 0 and positive entries = 1.
  # 2) Covariate: an n by p covariate matrix
  # 3) K: a positive integer which is no larger than n. This is the predefined number of communities.
  # 3) alpha: a positive number to tune the weight of covariate matrix

  # Optional Arguments for Kmeans:
  # 1) itermax: the maximum number of iterations allowed.
  # 2) nstart: R will try startn different random starting assignments and then select the one with the lowest within cluster variation.

  # Outputs:
  # 1) a factor indicating nodes' labels. Items sharing the same label are in the same community.

  if(!isSymmetric(Adj)) stop("Error! Adj is not symmetric!")
  if(K > dim(Adj)[1]) stop("Error! More communities than nodes!")
  if(K %% 1 != 0) stop("Error! K is not an integer!")
  if(K < 2) stop("Error: There must be at least 2 communities!")
  if(dim(Adj)[1] != dim(Covariate)[1]) stop("Error! Incompatible!")
  #  if(alpha < 0) stop("Negative Alpha")

  #Regularity check
  estall = rep(NA, dim(Adj)[1]);
  netrow = rowSums(Adj);
  covrow = rowSums(abs(Covariate));
  ind_reg = which(netrow != 0 | covrow != 0)
  Adj = Adj[ind_reg, ind_reg];
  Covariate = Covariate[ind_reg, ];

  ##Algorithm
  n = dim(Adj)[1]
  d = rowSums(Adj);
  X = Adj %*% Covariate

  diagcomp = cbind(1, median(d)/(d+0.5));
  lambda = apply(diagcomp, 1, min)

  if(is.null(alpha)){
    d.net = sort(abs(eig(Adj)), decreasing = TRUE);
    alphaupper = d.net[1]*log(n)/mean(d);
    alphalower = max(0.05, d.net[K]/4);

    alpha = seq(alphalower, alphaupper, length.out = alphan);
    #    print(alpha)
    est1 = matrix(0, alphan, n); est2 = est1;
    prop1 = rep(0, alphan); prop2 = rep(0, alphan)

    for(i in 1:alphan){
      Newmat = X + alpha[i]*diag(lambda)%*%Covariate;
      zz = Newmat%*%t(Newmat)
      c = eigen(zz)
      vec = c$vectors
      vecn = vec[,1:K]/apply(vec[,1:K], 1, Norm);
      result = kmeans(vecn, K, iter.max = itermax, nstart = startn);
      if (result$ifault==4) { result = kmeans(X, K,  iter.max = itermax, nstart = startn, algorithm="Lloyd"); }
      prop2[i] = result$tot.withinss;
      est2[i,] = as.factor(result$cluster);
    }
    est = est2[which.min(prop2), ];
    #    print(prop2)
  }
  else{
    Newmat = X + alpha*diag(lambda)%*%Covariate;
    c.svd = svd(Newmat, nu = K, nv = K)
    vec = c.svd$u
    vecn = vec/apply(vec, 1, Norm);
    result = kmeans(vecn, K, iter.max = itermax, nstart = startn);
    est = as.factor(result$cluster);
  }
  estall[ind_reg] = est;
  return(estall)
}
