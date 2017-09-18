// Construct coefficients for one-dimension pchip.
// Constructions follow matlab. Coefficients store like Intel MKL
// Author: Wenlan Luo, Tsinghua University, luowenlan@gmail.com

#include <string.h>

#define SIGN(x) (((x)>0)?1:(((x)<0)?-1:0))
#define MAX(x,y) (((x)>=(y))?(x):(y))
#define MIN(x,y) (((x)<=(y))?(x):(y))
#define ABS(x) (((x)>=0)?(x):(-x))

inline
void pchipslopes(double* x, double* y, double* del, int nx, double* h, double* d) {

    memset(d, 0, nx*sizeof(double));

    // 0..nx-1

    // interior
    // 1..nx-2
#pragma simd
    for (int i = 0;
        i < nx - 2;
        ++i)
    {
        if (SIGN(del[i]) * SIGN(del[i + 1]) > 0) {
            double hs = h[i] + h[i + 1];
            double w1 = (h[i] + hs) / (3 * hs);
            double w2 = (hs + h[i + 1]) / (3 * hs);
            double dmax = MAX(ABS(del[i]), ABS(del[i + 1]));
            double dmin = MIN(ABS(del[i]), ABS(del[i + 1]));
            d[i + 1] = dmin/(w1*del[i]/dmax + w2*(del[i + 1]/dmax));
        }
    }


    // slopes at end point
    d[0] = ((2 * h[0] + h[1])*del[0] - h[0] * del[1]) / (h[0] + h[1]);
    if (SIGN(d[0]) != SIGN(del[0]))
        d[0] = 0;
    else if ((SIGN(del[0]) != SIGN(del[1])) && (ABS(d[0]) > ABS(3 * del[0])))
        d[0] = 3 * del[0];

    d[nx - 1] = ((2 * h[nx - 2] + h[nx - 3])*del[nx - 2] - h[nx - 2] * del[nx - 3]) / (h[nx - 2] + h[nx - 3]);
    if (SIGN(d[nx - 1]) != SIGN(del[nx - 2]))
        d[nx - 1] = 0;
    else if (SIGN(del[nx - 2]) != SIGN(del[nx - 3]) && (ABS(d[nx - 1]) > ABS(3 * del[nx - 2])))
        d[nx - 1] = 3 * del[nx - 2];

}

/*
/// @param[in] x: grid
/// @param[in] nx: length of grid
/// @param[in] y: vector value, of size [ny,nx], row major form
/// @param[in] del: vector to store del, of size [ny,nx]
/// @param[in] h: vector to store h, of size [1,nx]
/// @param[in] slopes: vector to store slopes, of size [ny,nx]
/// @param[out] coefs: returned coefs
void pchip(double* x, int nx, double* yvec, int ny, double* del, double* h, double* slopes, double* coefs)
{
#pragma simd
	for (int ix = 0; ix < nx - 1; ix++)
	{
		h[ix] = x[ix + 1] - x[ix];
	}

	for (int iy = 0; iy < ny; iy++)
	{
		double* y = yvec + iy*nx;
#pragma simd
		for (int ix = 0; ix < nx - 1; ix++)
		{
			del[ix] = (y[ix + 1] - y[ix]) / h[ix];
		}
		pchipslopes(x, y, del, nx, h, slopes);
#pragma simd
		for (int ix = 0; ix < nx - 1; ix++)
		{
			double dzzdx = (del[ix] - slopes[ix]) / h[ix];
			double dzdxdx = (slopes[ix + 1] - del[ix]) / h[ix];
			// Intel coefs is stored [ny,nx-1,4], row major form
			double* pcoefs = coefs + iy*(nx-1) * 4 + ix * 4;
			pcoefs[0] = y[ix];
			pcoefs[1] = slopes[ix];
			pcoefs[2] = 2 * dzzdx - dzdxdx;
			pcoefs[3] = (dzdxdx - dzzdx) / h[ix];
		}
	}
}
*/

#define MAX_GRID_SIZE 1000

/// @param[in] x: grid
/// @param[in] nx: length of grid
/// @param[in] y: vector value, of size [ny,nx], row major form
/// @param[in] del: vector to store del, of size [ny,nx]
/// @param[in] h: vector to store h, of size [1,nx]
/// @param[in] slopes: vector to store slopes, of size [ny,nx]
/// @param[out] coefs: returned coefs
void pchip_omp(double* x, int nx, double* yvec, int ny, double* coefs)
{
	double h[MAX_GRID_SIZE];

#pragma simd
	for (int ix = 0; ix < nx - 1; ix++)
	{
		h[ix] = x[ix + 1] - x[ix];
	}

#pragma omp parallel for
	for (int iy = 0; iy < ny; iy++)
	{
		double del[MAX_GRID_SIZE];
		double slopes[MAX_GRID_SIZE];
		double* y = yvec + iy*nx;
#pragma simd
		for (int ix = 0; ix < nx - 1; ix++)
		{
			del[ix] = (y[ix + 1] - y[ix]) / h[ix];
		}
		pchipslopes(x, y, del, nx, h, slopes);
#pragma simd
		for (int ix = 0; ix < nx - 1; ix++)
		{
			double dzzdx = (del[ix] - slopes[ix]) / h[ix];
			double dzdxdx = (slopes[ix + 1] - del[ix]) / h[ix];
			// Intel coefs is stored [ny,nx-1,4], row major form
			double* pcoefs = coefs + iy*(nx-1) * 4 + ix * 4;
			pcoefs[0] = y[ix];
			pcoefs[1] = slopes[ix];
			pcoefs[2] = 2 * dzzdx - dzdxdx;
			pcoefs[3] = (dzdxdx - dzzdx) / h[ix];
		}
	}
}

void pchip(double* x, int nx, double* yvec, int ny, double* coefs)
{
	double h[MAX_GRID_SIZE];

#pragma simd
	for (int ix = 0; ix < nx - 1; ix++)
	{
		h[ix] = x[ix + 1] - x[ix];
	}

	for (int iy = 0; iy < ny; iy++)
	{
		double del[MAX_GRID_SIZE];
		double slopes[MAX_GRID_SIZE];
		double* y = yvec + iy*nx;
#pragma simd
		for (int ix = 0; ix < nx - 1; ix++)
		{
			del[ix] = (y[ix + 1] - y[ix]) / h[ix];
		}
		pchipslopes(x, y, del, nx, h, slopes);
#pragma simd
		for (int ix = 0; ix < nx - 1; ix++)
		{
			double dzzdx = (del[ix] - slopes[ix]) / h[ix];
			double dzdxdx = (slopes[ix + 1] - del[ix]) / h[ix];
			// Intel coefs is stored [ny,nx-1,4], row major form
			double* pcoefs = coefs + iy*(nx-1) * 4 + ix * 4;
			pcoefs[0] = y[ix];
			pcoefs[1] = slopes[ix];
			pcoefs[2] = 2 * dzzdx - dzdxdx;
			pcoefs[3] = (dzdxdx - dzzdx) / h[ix];
		}
	}
}


void construct_linear_interp(double* x, int nx, double* yvec, int ny, double* coefs)
{
	double h[MAX_GRID_SIZE];

#pragma simd
	for (int ix = 0; ix < nx - 1; ix++)
	{
		h[ix] = x[ix + 1] - x[ix];
	}

	for (int iy = 0; iy < ny; iy++)
	{
		double* y = yvec + iy*nx;
#pragma simd
		for (int ix = 0; ix < nx - 1; ix++)
		{
			// Intel coefs is stored [ny,nx-1,2], row major form
			double* pcoefs = coefs + iy*(nx - 1) * 2 + ix * 2;
			pcoefs[0] = y[ix];
			pcoefs[1] = (y[ix + 1] - y[ix]) / h[ix];
		}
	}
}

void construct_linear_interp_omp(double* x, int nx, double* yvec, int ny, double* coefs)
{
	double h[MAX_GRID_SIZE];

#pragma simd
	for (int ix = 0; ix < nx - 1; ix++)
	{
		h[ix] = x[ix + 1] - x[ix];
	}

#pragma omp parallel for
	for (int iy = 0; iy < ny; iy++)
	{
		double* y = yvec + iy*nx;
#pragma simd
		for (int ix = 0; ix < nx - 1; ix++)
		{
			// Intel coefs is stored [ny,nx-1,2], row major form
			double* pcoefs = coefs + iy*(nx - 1) * 2 + ix * 2;
			pcoefs[0] = y[ix];
			pcoefs[1] = (y[ix + 1] - y[ix]) / h[ix];
		}
	}
}
