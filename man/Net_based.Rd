% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Net_based.R
\name{Net_based}
\alias{Net_based}
\title{Network-based Clustering.}
\usage{
Net_based(Adj, K, tau = NULL, itermax = NULL, startn = NULL)
}
\arguments{
\item{Adj}{A 0/1 adjacency matrix.}

\item{K}{A positive integer, indicating the number of underlying communities in
graph \code{Adj}.}

\item{tau}{An optional tuning parameter, the default value is the mean of adajacency matrix.}

\item{itermax}{\code{k-means} parameter, indicating the maximum number of
iterations allowed. The default value is 100.}

\item{startn}{\code{k-means} parameter. If centers is a number, how many
random sets should be chosen? The default value is 10.}
}
\value{
A label vector.
}
\description{
\emph{Network-based Clustering} is a spectral clustering method that focuses
  solely on the topological structure of a network. It employs \code{k-means} on the first
  \eqn{K} leading eigenvectors of the weighted adjacency matrix of a graph, with each
  eigenvector normalized to have unit magnitude.
}
\examples{

# Simulate the Network
n = 10; K = 2;
theta = 0.4 + (0.45-0.05)*(seq(1:n)/n)^2; Theta = diag(theta);
P  = matrix(c(0.8, 0.2, 0.2, 0.8), byrow = TRUE, nrow = K)
set.seed(2022)
l = sample(1:K, n, replace=TRUE); # node labels
Pi = matrix(0, n, K) # label matrix
for (k in 1:K){
  Pi[l == k, k] = 1
}
Omega = Theta \%*\% Pi \%*\% P \%*\% t(Pi) \%*\% Theta;
Adj = matrix(runif(n*n, 0, 1), nrow = n);
Adj = Omega - Adj;
Adj = 1*(Adj >= 0)
diag(Adj) = 0
Adj[lower.tri(Adj)] = t(Adj)[lower.tri(Adj)]
Net_based(Adj, 2)
}
