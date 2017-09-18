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
	GET_INT(ePts);
	GET_INT(betaPts);
	GET_INT(kPts);

	GET_DV_VIEW(eGrid);
	GET_DV_VIEW(betaGrid);
	GET_DV_VIEW(kGrid);
	GET_DBL(kMin);
	GET_DBL(TOL_OPT);

	// tensor
	GET_DM_VIEW(budgetMesh, 3);

	// Coefs
	GET_DV_VIEW(valueKSplineCoefs);

	// Output
	GET_DM_VIEW(value_sub, 3);
	GET_DM_VIEW(kpPolicy_sub, 3);

	// openmp
#ifdef OMP
	GET_INT(NUM_THREADS);
	omp_set_num_threads(NUM_THREADS);
#endif

	// Loop over e and k
#pragma omp parallel for collapse (3)
	for (int i_e = 1; i_e <= ePts; i_e++)
	{
		for (int i_beta = 1; i_beta <= betaPts; i_beta++)
		{
			for (int i_k = 1; i_k <= kPts; i_k++)
			{
				double budget = budgetMesh(i_e, i_beta, i_k);

				// "Nested" function to evaluate value
				auto compute_value_of_kp = [&](double kp)
				{
					double c = budget - kp;
					// e goes first in memory
					double v = log(c) + search_eval4_1d(_kGrid, kPts, _valueKSplineCoefs, kp, (i_beta - 1)*ePts + (i_e - 1));

					return -v;
				};

				// Call the minimization routine
				value_sub(i_e, i_beta, i_k) = -find_local_min(kMin, budget - 1e-12, 1e-15, TOL_OPT, compute_value_of_kp, &kpPolicy_sub(i_e, i_beta, i_k));
			}
		}
	}
}
