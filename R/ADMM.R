#' @title Penalized Optimization Framework for Community Detection in Networks with Covariates.
#' @description Semidefinite programming for optimizing the inner product between combined network and the
#'   solution matrix.
#' @details \emph{ADMM} is proposed in \emph{Covariate Regularized Community Detection in Sparse Graphs}
#'   of Yan & Sarkar (2021). \emph{ADMM} relies on semidefinite programming (SDP) relaxations for detecting
#'   the community structure in sparse networks with covariates.
#' @param Adj A 0/1 adjacency matrix.
#' @param Covariate A covariate matrix. The rows correspond to nodes and the columns correspond to covariates.
#' @param lambda A tuning parameter to weigh the covariate matrix.
#' @param K A positive integer, indicating the number of underlying communities in graph \code{Adj}.
#' @param alpha A number. The elementwise upper bound in the SDP.
#' @param rho The learning rate of ADMM.
#' @param TT The maximum of iteration.
#' @param tol The tolerance for stopping criterion.
#' @param quiet An optional inoput. Whether to print result at each step.
#' @param report_interval An optional inoput. The frequency to print intermediate result.
#' @param r An optional inoput. The expected rank of the solution, leave NULL if no constraint is required.
#' @return \item{estall}{A lavel vector.}
#'
#' @importFrom pracma eig Norm
#' @importFrom stats kmeans runif rnorm
#'
#' @references Yan, B., & Sarkar, P. (2021). \emph{Covariate Regularized Community Detection in Sparse Graphs}.
#'   \emph{Journal of the American Statistical Association, 116(534), 734-745}.
#'   \cr\doi{10.1080/01621459.2019.1706541}\cr
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
#' ADMM(Adj, X, lambda = 0.2, K = K, alpha = 0.5, rho = 2, TT = 100, tol = 5)

#' @export

ADMM = function(Adj, Covariate, lambda, K, alpha, rho, TT, tol, quiet = NULL, report_interval = NULL, r = NULL){

  # Inputs:     Adj: adjacency matrix
  #             Covariate: n x p covaraite matrix
  #             lambda: tuning parameter between graph and covariates
  #             K: number of clusters
  #             alpha: elementwise upper bound in the SDP
  #             rho: learning rate of ADMM
  #             TT:   max iteration
  #             tol: tolerance for stopping criterion
  #             quiet: whether to print result at each step
  #             report_interval: frequency to print intermediate result
  #             r: expected rank of the solution, leave blank if no constraint is required.
  # Outputs:    estall: the label vector

  As = Adj + lambda* Covariate %*% t(Covariate)

  n = dim(As)[1]
  U <- V <- matrix(0, n, n)

  #Initialization - spectral with perturbation
  v = eigen(As)$vectors[, 1: K]
  e = diag(eigen(As)$values[1: K])
  X = v %*% t(v) + 0.1*matrix(rnorm(n*n), nrow = n)
  Y = v %*% t(v) + 0.1*matrix(rnorm(n*n), nrow = n)
  Z = v %*% t(v) + 0.1*matrix(rnorm(n*n), nrow = n)

  As_rescaled = (1/rho) * As;

  if(is.null(report_interval)) report_interval = 1
  if(is.null(quiet)) quiet = FALSE
  if(is.null(r)) r = Inf

  if (is.infinite(TT)) {
    delta = matrix(0, 1000, 1)
    infeas = matrix(0, 1000, 1)
  }
  else {
    delta = matrix(0, TT, 1)
    infeas = matrix(0, TT, 1)
  }

  dt = matrix(0, 1, 3);
  t = 1;
  CONVERGED = FALSE;

  while(CONVERGED == FALSE & t<= TT){
    Xold = X
    X = projAXb(0.5*(Z - U + Y - V + As_rescaled), K, n);
    Z = X + U
    Z[Z < 0] = 0
    Z[Z > alpha] = alpha

    Y = projSp(X + V);
    U = U + X - Z;
    V = V + X - Y;
    delta[t] = norm(X - Xold) / norm(Xold);
    infeas[t] = (sqrt(sum(diag(t(X - Y) * (X - Y)))) + sqrt(sum(diag(t(X - Z) * (X - Z))))) / sqrt(sum(diag(t(X)*X)));
    CONVERGED = max(delta[t], infeas[t]) < tol;
    t = t + 1;
  }

  T_term = t - 1
  estall = rsc(X, K, 'adj')
  return(estall)
}
