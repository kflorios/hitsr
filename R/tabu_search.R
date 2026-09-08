#' Run Discrete Tabu Search
#'
#' @param X Numeric matrix of dimensions (T, p_in).
#' @param y Integer vector of length T containing binary target indicators.
#' @param b0 Double, intercept parameter.
#' @param d Integer parameter.
#' @param iSeed Integer, random number generator seed.
#' @return A list containing best attributes, coefficients, and best score.
#' @useDynLib hitsr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
run_tabu_search <- function(X, y, b0 = -1.0, d = 10, iSeed = 123456) {
  X_mat <- as.matrix(X)
  storage.mode(X_mat) <- "double"
  
  y_vec <- as.integer(y)
  
  res <- cpp_run_tabu_search(X_mat, y_vec, b0, as.integer(d), as.integer(iSeed))
  return(res)
}
