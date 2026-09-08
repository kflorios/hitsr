#ifndef TABU_SEARCH_HPP
#define TABU_SEARCH_HPP

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <omp.h>

// -----------------------------------------------------------------------------
// Global constants (Fortran PARAMETER)
// -----------------------------------------------------------------------------
static const int Tmax = 50000;
static const int pmax = 20;
static const int Itmax = 1000;

// -----------------------------------------------------------------------------
// Function prototypes
// -----------------------------------------------------------------------------

// bsort: sort nat[0..np-1] into snat[0..np-1]
void bsort(int np, const int nat[], int snat[]);

// Gaussian elimination (Fortran gelim)
void gelim(int np,
           double sk[pmax][pmax],
           double r1[pmax],
           double u[pmax]);

// EvaluateScore1: compute score h
void EvaluateScore1(
    const double w[pmax],
    int T,
    int p,
    double &h,
    const double Z[pmax][Tmax],
    double b0,
    const int data_y[Tmax]
);

// evaluate_neighbourhood: OpenMP-parallelized neighborhood evaluation
void evaluate_neighbourhood(
    int i_idx,
    int p,
    int n,
    const int natsol[pmax],
    const int Mat[pmax][Itmax],
    const double Z[pmax][Tmax],
    const double y[Tmax],
    const int id[Tmax],
    const double xlbound[pmax],
    const double xubound[pmax],
    double b0,
    double MaxScore,
    int ntenure,
    int csubstitutionstotal,
    int T,
    double &dz_i,
    int natcur_i[pmax],
    double &MaxScoreEver_i,
    int natever_i[pmax],
    int &row_i,
    int &cbd_i,
    int &cinfsd_i,
    const int data_y[Tmax]
);

#endif // TABU_SEARCH_HPP
