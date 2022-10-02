#' @title normalizeSym
#' @param A A 0/1 adjacency matrix.
#' @param L An optional input.
#' @return \item{An}{The result matrix}.
#' @noRd
#' @keywords internal

normalizeSym = function(A, L = NULL){
  if (is.null(L)) L = 0
  n = dim(A)[1]
  d = colSums(A)
  d[which(d == 0)] = 1e10
  if (L == 0) {
    d = 1/sqrt(d)
    D = diag(d)
    An = D %*% A %*% D
    An = 1/2*(An + t(An))
  }
  else{
    d = 1 / d
    D=diag(d)
    An=D %*% A;
  }
  return(An)
}
