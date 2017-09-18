#include "mex.h"
#include "MatlabMatrix.h"
#include "vpchip.h"
#include "CInterp.h"
#include "nr3_opt.h"

#ifdef OMP
#include "omp.h"
#endif

using namespace CInterp;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Parameters
	GET_DBL(gamma);
	GET_DBL(kMin);
	GET_DBL(TOL_OPT);

	// Get variables from matlab workspace
	GET_INT(ePts);
	GET_INT(kPts);
	GET_DV_VIEW(kGrid);
	GET_DV_VIEW(eGrid);

	// Prices
	GET_DBL(w);
	GET_DBL(r);
	
	// Get vFuture
	GET_DV_VIEW(vFutureCoefs);

	// Output
	GET_DM_VIEW(kpPolicy, 2);
	GET_DM_VIEW(v_new, 2);

	// openmp
#ifdef OMP
	GET_INT(NUM_THREADS);
	omp_set_num_threads(NUM_THREADS);
#endif

	// Loop over e and k
#pragma omp parallel for collapse (2)
	for (int i_e = 1; i_e <= ePts ; i_e++)
	{
		for (int i_k = 1; i_k <= kPts ; i_k++)
		{
			double e = eGrid(i_e);
			double k = kGrid(i_k);

			// Compute budget
			double budget = k*(1 + r) + w*e;

			// "Nested" function to evaluate value
			auto compute_value_of_kp = [&](double kp)
			{
				double c = budget - kp;
				double v = log(c) + search_eval4_1d(_kGrid, kPts, _vFutureCoefs, kp, i_e-1);

				return -v;
			};

			// Call the minimization routine
			v_new(i_e, i_k) = -find_local_min(kMin, budget - 1e-12, 1e-15, TOL_OPT, compute_value_of_kp, &kpPolicy(i_e, i_k));

		}
	}
}
