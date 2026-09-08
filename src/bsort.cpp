#include "tabu_search.hpp"

// -----------------------------------------------------------------------------
// bsort: sort nat[0..np-1] into snat[0..np-1] in increasing order
// Fortran original:
//   do i=1,np
//   do i=1,np-1
//      do j=i+1,np
// -----------------------------------------------------------------------------
void bsort(int np, const int nat[], int snat[])
{
    // copy input
    for (int i = 0; i < np; ++i)
        snat[i] = nat[i];

    // simple bubble/selection hybrid as in Fortran
    for (int i = 0; i < np - 1; ++i) {
        for (int j = i + 1; j < np; ++j) {
            if (snat[i] > snat[j]) {
                int tmp = snat[i];
                snat[i] = snat[j];
                snat[j] = tmp;
            }
        }
    }
}

