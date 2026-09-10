#include "tabu_search.hpp"

// Portable RNG
static inline double get_random() {
    return double(rand()) / double(RAND_MAX);
}

// -----------------------------------------------------------------------------
// Public API entry point callable from Python (pybind11)
// -----------------------------------------------------------------------------
extern "C" int run_tabu_search(
    double* X,
    int T,
    int p_in,
    int* y_in,
    double b0,
    double d,
    int iSeed,                 // <-- NEW ARGUMENT
    int* out_attributes,
    double* out_coeffs,
    double* out_score
)
{
	srand(iSeed);
	//srand(12345);   // or use time(NULL) if you want nondeterministic runs
	//srand((unsigned)time(NULL)); //If you want randomness each run:	

    int p = p_in - 1;   // number of attributes (excluding intercept)

    static double Z[pmax][Tmax];   // now pmax and Tmax are known
    static double y[Tmax];
    static int data_y[Tmax];
    static int id[Tmax];

    if (p <= 0 || p > pmax)
        return -10;
    if (T <= 0 || T > Tmax)
        return -11;

    // Core data
    /*double Z[pmax][Tmax];
    double y[Tmax];
    int data_y[Tmax];
    int id[Tmax];*/

    // Working arrays
    double A[pmax][pmax], b[pmax], bb[pmax];
    double mA[pmax][pmax], mb[pmax], bmb[pmax];
    double x[pmax], score1, MaxScore0, MaxScore, MaxScoreEver;
    double Af[pmax][pmax], bf[pmax], bbf[pmax];
    double ArgMax[pmax];

    int nat[pmax], natsol[pmax], natcur[pmax], natever[pmax];
    int Mat[pmax][Itmax], snat[pmax], snatsol[pmax];

    double xlbound[pmax], xubound[pmax];

    int ntenure, niterover, niternotimproved;
    int csubstitutionstotal, cbdsubseq, cbd, cinfsd;
    double dzmax;
    int row;

    // -------------------------------------------------------------------------
    // Copy inputs into local arrays
    // -------------------------------------------------------------------------
    for (int i = 0; i < T; ++i) {
        id[i] = i;
        data_y[i] = y_in[i];
        y[i] = double(y_in[i]);
    }

    for (int j = 0; j < p_in; ++j) {
        for (int i = 0; i < T; ++i) {
            Z[j][i] = X[j*T + i];
        }
    }

    // -------------------------------------------------------------------------
    // Bounds (you can later tie these to d if desired)
    // -------------------------------------------------------------------------
    for (int j = 0; j < p; ++j) {
        //xlbound[j] = -1.0e300;
        //xubound[j] =  1.0e300;
        xlbound[j] = -d;
        xubound[j] =  d;
    }

    // -------------------------------------------------------------------------
    // Initial nat (first p indices)
    // -------------------------------------------------------------------------
    /*for (int j = 0; j < p; ++j)
        nat[j] = j;*/
	for (int j = 0; j < p; ++j) {
    double rr = get_random();
    nat[j] = int(std::floor(rr * T));
}
    printf("Initial nat: ");
	for (int j = 0; j < p; ++j) printf("%d ", nat[j]);
	printf("\n");


    bsort(p, nat, snat);
    for (int j = 0; j < p; ++j)
        natsol[j] = snat[j];

    // -------------------------------------------------------------------------
    // Build mb
    // -------------------------------------------------------------------------
    for (int j = 0; j < p; ++j)
        mb[j] = (y[id[natsol[j]]] >= 0.0 ? 1.0 : 0.0);

    // -------------------------------------------------------------------------
    // Build mA and bmb
    // -------------------------------------------------------------------------
    for (int j = 0; j < p; ++j)
        for (int i = 0; i < p; ++i)
            mA[i][j] = (1.0 - 2.0 * mb[i]) * Z[j + 1][ id[natsol[i]] ];

    for (int i = 0; i < p; ++i)
        bmb[i] = b0 * (-1.0 + 2.0 * mb[i]) * Z[0][ id[natsol[i]] ];

    // -------------------------------------------------------------------------
    // Solve initial system
    // -------------------------------------------------------------------------
    gelim(p, mA, bmb, x);

    for (int i = 0; i < p; ++i) {
        if (x[i] < xlbound[i] || x[i] > xubound[i]) {
            return -2; // infeasible initial solution
        }
    }

    // -------------------------------------------------------------------------
    // Evaluate initial score
    // -------------------------------------------------------------------------
    EvaluateScore1(
        x, T, p, score1,
        Z, b0, data_y
    );
    MaxScore0 = score1;
    MaxScore = MaxScore0;
    MaxScoreEver = MaxScore0;

    // -------------------------------------------------------------------------
    // Tabu parameters
    // -------------------------------------------------------------------------
    ntenure = 200;
    niterover = 200;
    niternotimproved = 200;

    for (int j = 0; j < p; ++j)
        for (int i = 0; i < ntenure; ++i)
            Mat[j][i] = -1;

    csubstitutionstotal = 0;
    cbdsubseq = 0;

    // -------------------------------------------------------------------------
    // Main tabu loop
    // -------------------------------------------------------------------------
    while (true) {

        cbd = 0;
        cinfsd = 0;
        dzmax = -1.0;
        row = 0;

        bsort(p, natsol, snatsol);
        for (int i = 0; i < p; ++i)
            natsol[i] = snatsol[i];

        for (int i = 0; i < p; ++i)
            b[i] = (y[id[natsol[i]]] >= 0.0 ? 1.0 : 0.0);

        for (int j = 0; j < p; ++j)
            for (int i = 0; i < p; ++i)
                A[i][j] = (1.0 - 2.0 * b[i]) * Z[j + 1][ id[natsol[i]] ];

        for (int i = 0; i < p; ++i)
            bb[i] = b0 * (-1.0 + 2.0 * b[i]) * Z[0][ id[natsol[i]] ];

        int row_i_arr[pmax];
        int cbd_i_arr[pmax];
        int cinfsd_i_arr[pmax];
        double dz_i_arr[pmax];
        double MaxScoreEver_i_arr[pmax];
        int natcur_i_arr[pmax][pmax];
        int natever_i_arr[pmax][pmax];

        #pragma omp parallel for default(shared)
        for (int i = 0; i < p; ++i) {
            evaluate_neighbourhood(
                i, p, T,
                natsol, Mat,
                Z, y, id,
                xlbound, xubound,
                b0, MaxScore,
                ntenure, csubstitutionstotal, T,
                dz_i_arr[i], natcur_i_arr[i],
                MaxScoreEver_i_arr[i], natever_i_arr[i],
                row_i_arr[i], cbd_i_arr[i], cinfsd_i_arr[i],
                data_y
            );
        }

        for (int i = 0; i < p; ++i) {
            row += row_i_arr[i];
            cbd += cbd_i_arr[i];
            cinfsd += cinfsd_i_arr[i];

            if (MaxScoreEver_i_arr[i] > MaxScoreEver) {
                MaxScoreEver = MaxScoreEver_i_arr[i];
                for (int j = 0; j < p; ++j)
                    natever[j] = natever_i_arr[i][j];
            }

            if (dz_i_arr[i] > dzmax) {
                dzmax = dz_i_arr[i];
                for (int j = 0; j < p; ++j)
                    natcur[j] = natcur_i_arr[i][j];
            }
        }

        csubstitutionstotal++;

        if (cbd == 0)
            cbdsubseq++;

        if (csubstitutionstotal == niterover ||
            cbdsubseq == niternotimproved)
            break;

        for (int j = 0; j < p; ++j)
            natsol[j] = natcur[j];

        MaxScore += dzmax;

        bsort(p, natsol, snatsol);
        for (int j = 0; j < p; ++j)
            Mat[j][ntenure + csubstitutionstotal] = snatsol[j];

		printf("iteration...%d...:",csubstitutionstotal);
		printf("current natsol: ");
		for (int j = 0; j < p; ++j) printf("%d ", natsol[j]);
		printf("\n");
		
    }

    // -------------------------------------------------------------------------
    // Final evaluation
    // -------------------------------------------------------------------------
    for (int j = 0; j < p; ++j)
        bf[j] = (y[id[natever[j]]] >= 0.0 ? 1.0 : 0.0);

    for (int j = 0; j < p; ++j)
        for (int l = 0; l < p; ++l)
            Af[l][j] = (1.0 - 2.0 * bf[l]) * Z[j + 1][ id[natever[l]] ];

    for (int j = 0; j < p; ++j)
        bbf[j] = b0 * (-1.0 + 2.0 * bf[j]) * Z[0][ id[natever[j]] ];

    gelim(p, Af, bbf, x);
    EvaluateScore1(
        x, T, p, score1,
        Z, b0, data_y
    );

    for (int j = 0; j < p; ++j) {
        out_attributes[j] = natever[j];
        out_coeffs[j] = x[j];
    }

    *out_score = score1;

    return 0;
}
