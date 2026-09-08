#include "tabu_search.hpp"

// -----------------------------------------------------------------------------
// gelim: Gaussian elimination (Fortran version preserved exactly)
// 0-based C++ translation of the Fortran routine:
//
//   do i=1,np-1
//      do j=i+1,np
//         p = sk(j,i)/sk(i,i)
//         r1(j) = r1(j) - p*r1(i)
//         do k=1,np
//             sk(j,k) = sk(j,k) - p*sk(i,k)
//         enddo
//      enddo
//   enddo
//
//   Back substitution, then restore original sk and r1
// -----------------------------------------------------------------------------
void gelim(int np,
           double sk[pmax][pmax],
           double r1[pmax],
           double u[pmax])
{
    double bsk[pmax][pmax];
    double br1[pmax];

    // Save original sk and r1 (Fortran behavior)
    for (int i = 0; i < np; ++i) {
        br1[i] = r1[i];
        for (int j = 0; j < np; ++j)
            bsk[i][j] = sk[i][j];
    }

    // Forward elimination
    for (int i = 0; i < np - 1; ++i) {
        for (int j = i + 1; j < np; ++j) {
            double p = sk[j][i] / sk[i][i];
            r1[j] -= p * r1[i];
            for (int k = 0; k < np; ++k)
                sk[j][k] -= p * sk[i][k];
        }
    }

    // Back substitution
    for (int l = np - 1; l >= 0; --l) {
        u[l] = r1[l] / sk[l][l];
        for (int j = l - 1; j >= 0; --j) {
            double p1 = sk[j][l] / sk[l][l];
            sk[j][l] = 0.0;
            r1[j] -= p1 * r1[l];
        }
    }

    // Restore original sk and r1 (Fortran semantics)
    for (int i = 0; i < np; ++i) {
        r1[i] = br1[i];
        for (int j = 0; j < np; ++j)
            sk[i][j] = bsk[i][j];
    }
}
