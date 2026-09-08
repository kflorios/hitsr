#include "tabu_search.hpp"

// -----------------------------------------------------------------------------
// evaluate_neighbourhood
// 0-based C++ translation of the Fortran routine.
//
// This routine evaluates all possible substitutions of attribute i_idx
// with every candidate id[j], j = 0..n-1, except those already in natsol.
//
// It performs:
//   - triviality check
//   - tabu check
//   - local system solve (gelim)
//   - bounds check
//   - score computation (EvaluateScore1)
//   - tracking of best dz_i and MaxScoreEver_i
//
// This version is called inside an OpenMP parallel loop in main.cpp.
// -----------------------------------------------------------------------------
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
)
{
	row_i = 0;
	cbd_i = 0;
	cinfsd_i = 0;
	dz_i = -1.0;
	MaxScoreEver_i = -1.0e99;

	int nat_local[pmax];
	int snat_local[pmax];
	double b_local[pmax];
	double A_local[pmax][pmax];
	double bb_local[pmax];
	double x_local[pmax];

	for (int j = 0; j < n; ++j)
	{
		bool trivial = false;
		for (int iii = 0; iii < p; ++iii)
		{
			if (natsol[iii] == id[j])
			{
				trivial = true;
				break;
			}
		}
		if (trivial)
			continue;

		row_i++;

		for (int iii = 0; iii < p; ++iii)
			nat_local[iii] = (iii == i_idx ? id[j] : natsol[iii]);

		bsort(p, nat_local, snat_local);

		bool tabu = false;
		for (int k = 0; k < ntenure; ++k)
		{
			int ct = 0;
			for (int iii = 0; iii < p; ++iii)
			{
				if (snat_local[iii] == Mat[iii][csubstitutionstotal + k])
					ct++;
			}
			if (ct == p)
			{
				tabu = true;
				break;
			}
		}
		if (tabu)
			continue;

		for (int iii = 0; iii < p; ++iii)
			b_local[iii] = (y[nat_local[iii]] >= 0.0 ? 1.0 : 0.0);

		for (int iii = 0; iii < p; ++iii)
		{
			for (int jjj = 0; jjj < p; ++jjj)
			{
				A_local[jjj][iii] =
					(1.0 - 2.0 * b_local[jjj]) *
					Z[iii + 1][id[nat_local[jjj]]];
			}
		}

		for (int iii = 0; iii < p; ++iii)
		{
			bb_local[iii] =
				b0 * (-1.0 + 2.0 * b_local[iii]) *
				Z[0][id[nat_local[iii]]];
		}

		gelim(p, A_local, bb_local, x_local);

		bool infeasible = false;
		for (int iii = 0; iii < p; ++iii)
		{
			if (x_local[iii] < xlbound[iii] ||
				x_local[iii] > xubound[iii])
			{
				cinfsd_i++;
				infeasible = true;
				break;
			}
		}
		if (infeasible)
			continue;

		double score1;
		EvaluateScore1(
			x_local, T, p, score1,
			Z, b0, data_y
		);

		double dz = score1 - MaxScore;

		if (dz > 0.0)
			cbd_i++;

		if (score1 > MaxScoreEver_i)
		{
			MaxScoreEver_i = score1;
			for (int iii = 0; iii < p; ++iii)
				natever_i[iii] = nat_local[iii];
		}

		if (dz > dz_i)
		{
			dz_i = dz;
			for (int iii = 0; iii < p; ++iii)
				natcur_i[iii] = nat_local[iii];
		}
	}
}
