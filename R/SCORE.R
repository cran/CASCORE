#' @title Spectral Clustering On Ratios-of-Eigenvectors.
#' @description Using ratios-of-eigenvectors to detect underlying communities.
#' @details \emph{SCORE} is fully established in \emph{Fast community detection by
#'   SCORE} of Jin (2015). \emph{SCORE} uses the entry-wise ratios between the
#'   first leading eigenvector and each of the other \eqn{K-1} leading eigenvectors for
#'   clustering. It is noteworthy that SCORE only works on connected graphs,
#'   in other words, it does not allow for isolated vertices.
#' @param G A 0/1 adjacency matrix of a connected graph.
#' @param K A positive integer, indicating the number of underlying communities in graph \code{G}.
#' @param itermax \code{k-means} parameter, indicating the maximum number of
#'   iterations allowed. The default value is 100.
#' @param startn \code{k-means} parameter. If centers is a number, how many
#'   random sets should be chosen? The default value is 10.
#' @return \item{estall}{A lavel vector.}
#'
#' @importFrom stats kmeans runif
#'
#' @references Jin, J. (2015). \emph{Fast community detection by score}.
#'   \emph{The Annals of Statistics 43 (1), 57–89.}\cr\doi{10.1214/14-AOS1265}\cr
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
#' library(igraph)
#' is.igraph(Adj) # [1] FALSE
#' ix = components(graph.adjacency(Adj))
#' componentLabel = ix$membership
#' giantLabel = which(componentLabel == which.max(ix$csize))
#' Giant = Adj[giantLabel, giantLabel]
#' SCORE(Giant, 2)
#'
#' @export


####################################################
######## Spectral Clustering Method: SCORE #########
####################################################

# Assume there are n nodes and K communities
# Before applying SCORE, need to:
# 1) transform the network graph into an n by n adjacency matrix. It has following properties:
#    i)   symmetrix
#    ii)  diagonals = 0
#    iii) positive entries = 1.


##### SCORE #####
# spectral clustering On ratios-of-eigenvectors
SCORE = function(G, K, itermax = NULL, startn = NULL){
  # Inputs:
  # 1) G: an n by n symmetric adjacency matrix whose diagonals = 0 and positive entries = 1.
  # 2) K: a positive integer which is no larger than n. This is the predefined number of communities.

  # Optional Arguments for Kmeans:
  # 1) itermax: the maximum number of iterations allowed.
  # 2) nstart: R will try startn different random starting assignments and then select the one with the lowest within cluster variation.

  # Outputs:
  # 1) a factor indicating nodes' labels. Items sharing the same label are in the same community.

  # Remark:
  # SCORE only works on connected graphs, i.e., no isolated node is allowed.

  # exclude all wrong possibilities:
  if(!isSymmetric(G)) stop("Error! G is not symmetric!")
  if(K > dim(G)[1]) stop("Error! More communities than nodes!")
  if(K %% 1 != 0) stop("Error! K is not an integer!")
  if(K  <= 0) stop("Error! Nonpositive K!")

  g.eigen = eigen(G)
  if(sum(g.eigen$vectors[, 1]==0) > 0) stop("Error! Zeroes in the first column")
  R = g.eigen$vectors[, -1]
  R = R[, 1: (K-1)]
  R = R / g.eigen$vectors[, 1]

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
