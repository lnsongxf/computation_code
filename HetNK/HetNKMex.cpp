// The Consumer Credit Channel
// Author: Wenlan Luo (luowenlan at gmail.com), Georgetown University
// qqsMex.cpp
// Last updated: 2016-2-18

#include <signal.h>

#include "mex.h"
#include "math.h"
#include "CInterp.h"
#include "MatlabMatrix.h"
#include "nr3_opt.h"
#include <string.h>
#include <algorithm>
#include <limits>

#ifdef USE_OMP
#include <omp.h>
#endif

#define MAX(a,b)(((a)>(b))?(a):(b))
#define MIN(a,b)(((a)>(b))?(b):(a))

void my_function_to_handle_aborts(int signal_number)
{
	/*
	printf("Break here\n");
	exit(-1);
	*/
	char ErrMsg[200];
	sprintf(ErrMsg, "Abort from CMEX.\nLINE: %d\nFILE: %s\n", __LINE__, __FILE__);
	mexErrMsgTxt(ErrMsg);
}

#define CRRA(c) ( pow((c),1-Sigma)/(1-Sigma) )

template<typename T>
T find_max(int nData, T* data)
{
	T currentMax = -1e20;
	for (int i = 0; i < nData; i++)
	{
		currentMax = MAX(data[i], currentMax);
	}
	return currentMax;
}

void SOLVE_N_BP0();


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Handle errors
	signal(SIGABRT, &my_function_to_handle_aborts);

	// OPENMP
	GET_INT(NumThreads);
#ifdef USE_OMP
	omp_set_num_threads(NumThreads);
#endif

	// Get task number
	GET_INT(MEX_TASK);
	GET_INT(MEX_SOLVE_N_BP0);

	if (MEX_TASK == MEX_SOLVE_N_BP0)
	{
		SOLVE_N_BP0();
		return;
	}
	mexErrMsgTxt("No task executed");
}

void SOLVE_N_BP0()
{
	GET_DV_VIEW(bPolicy);
	GET_DV_VIEW(EMesh);
	GET_DV_VIEW(BudgetWithoutLabor);
	GET_DV_VIEW(nZero);
	GET_INT(NumProblems);
	GET_DBL(Gamma);
	GET_DBL(Chi);
	GET_DBL(Nu);
	GET_DBL(wage);

#pragma omp parallel for schedule(dynamic)
	for (int i = 1; i <= NumProblems ; i++)
	{
		if (bPolicy(i) <= 0) {
			auto focForLabor = [&](double n) {
				double c = BudgetWithoutLabor(i) + wage*EMesh(i)*n;
				double foc = wage*EMesh(i)*pow(c, -Gamma) - Chi*pow(n, 1 / Nu);

				return foc;
			};

			int exitFlag;
			int iter;
			nZero(i) = zbrent(focForLabor, 0.0, 10.0, 1e-10, 100, &exitFlag, &iter);
		}
	}
}
