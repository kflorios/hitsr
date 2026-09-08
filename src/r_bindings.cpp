#include <Rcpp.h>
#include "tabu_search.hpp"

extern "C" int run_tabu_search(
    double* X,
    int T,
    int p_in,
    int* y_in,
    double b0,
    int d,
    int iSeed,
    int* out_attributes,
    double* out_coeffs,
    double* out_score
);

// [[Rcpp::export]]
Rcpp::List cpp_run_tabu_search(Rcpp::NumericMatrix X_in, 
                              Rcpp::IntegerVector y_in, 
                              double b0, 
                              int d, 
                              int iSeed) 
{
    // X_in expected shape in R: (T, p_in)
    // R column-major (T, p_in) matches C++ row-major (p_in, T) in flat memory layout.
    int T = X_in.nrow();
    int p_in = X_in.ncol();

    if (y_in.size() != T) {
        Rcpp::stop("Length of y must match number of rows in X");
    }

    int p = p_in - 1;
    if (p <= 0 || p > pmax) {
        Rcpp::stop("p_in - 1 must be between 1 and pmax");
    }

    std::vector<int> attrs(p);
    std::vector<double> coeffs(p);
    double score = 0.0;

    // Direct pointer cast to memory buffer
    double* X_ptr = const_cast<double*>(&(X_in[0]));
    int* y_ptr = const_cast<int*>(&(y_in[0]));

    int ret = run_tabu_search(
        X_ptr, T, p_in, y_ptr, b0, d, iSeed,
        attrs.data(), coeffs.data(), &score
    );

    if (ret != 0) {
        Rcpp::stop("run_tabu_search failed with internal error code: %d", ret);
    }

    // Convert 0-based C++ attribute indices to 1-based R indices
    Rcpp::IntegerVector r_attrs(p);
    Rcpp::NumericVector r_coeffs(p);
    for (int i = 0; i < p; ++i) {
        r_attrs[i] = attrs[i] + 1;
        r_coeffs[i] = coeffs[i];
    }

    return Rcpp::List::create(
        Rcpp::Named("attributes") = r_attrs,
        Rcpp::Named("coeffs")     = r_coeffs,
        Rcpp::Named("score")      = score
    );
}