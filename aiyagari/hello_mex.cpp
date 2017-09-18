#include "mex.h"
#include "MatlabMatrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Get variables from matlab workspace
	GET_INT(inputSize);
	GET_DV_VIEW(input);
	GET_DV_VIEW(output);

	printf("hello, mex!\n");

	for (int i = 1; i <= inputSize ; i++)
	{
		output(i) = input(i) * 2;
	}
}
