#include "tabu_search.hpp"

// -----------------------------------------------------------------------------
// EvaluateScore1
//   sump(i) = sum_j beta(j)*Z(j+1,i) + b0*Z(1,i)
//   firep(i) = sign(sump(i))
//   sfire = sum_i data_y(i)*firep(i)
//   h = (1 + sfire/T) / 2
// -----------------------------------------------------------------------------
void EvaluateScore1(
	const double w[pmax],
	int T,
	int p,
	double &h,
	const double Z[pmax][Tmax],
	double b0,
	const int data_y[Tmax]
)
{
	double beta[pmax];
	int firep[Tmax];

	for (int j = 0; j < p; ++j)
		beta[j] = w[j];

	int sfire = 0;

	for (int i = 0; i < T; ++i)
	{
		double s = 0.0;

		for (int j = 0; j < p; ++j)
			s += beta[j] * Z[j + 1][i];

		s += b0 * Z[0][i];

		firep[i] = (s >= 0.0 ? 1 : -1);

		sfire += data_y[i] * firep[i];
	}

	h = (1.0 + double(sfire) / double(T)) * 0.5;
}
