#pragma once

#include <blitz/array.h>
#include <MatlabMatrix.h>

#define MAXDIM 10

namespace CInterp {
	/*
	* recursive evaluation doesn't allow compiler to optimize
	*/
	inline double receval4_1d(double* Coefs, double* XSite, int* CellOfSite, int Shift)
	{
		double r;
		double* pCoefs = Coefs + Shift + (*CellOfSite + 1) * 4;
		r = *(--pCoefs);
		r *= *XSite;
		r += *(--pCoefs);
		r *= *XSite;
		r += *(--pCoefs);
		r *= *XSite;
		r += *(--pCoefs);
		return r;
	}

	/*
	* recursive evaluation doesn't allow compiler to optimize
	*/
	inline double receval2_1d(double* Coefs, double* XSite, int* CellOfSite, int Shift)
	{
		double r;
		double* pCoefs = Coefs + Shift + (*CellOfSite + 1) * 2;
		r = *(--pCoefs);
		r *= *XSite;
		r += *(--pCoefs);
		return r;
	}

	inline double eval4_1d(double* XGrid, double* Coefs, int* CoefsSize,
		double* XSite, int* CellOfSite, int VecIdx)
	{
		double XSiteToLeft = *XSite - XGrid[*CellOfSite];
		return receval4_1d(Coefs, &XSiteToLeft, CellOfSite, VecIdx*CoefsSize[1]);
	}

	inline double receval4_1d_011(double* Coefs, double* XSite, int* CellOfSite, int Shift, double* g)
	{
		double r;

		double* pCoefs = Coefs + Shift + (*CellOfSite + 1) * 4;

		--pCoefs;
		r = *pCoefs;
		*g = *pCoefs * 3;

		--pCoefs;
		r *= *XSite;
		r += *pCoefs;
		*g *= *XSite;
		*g += *pCoefs * 2;

		--pCoefs;
		r *= *XSite;
		r += *pCoefs;
		*g *= *XSite;
		*g += *pCoefs;

		--pCoefs;
		r *= *XSite;
		r += *pCoefs;
		return r;
	}

	inline double receval4_1d_111(double* Coefs, double* XSite, int* CellOfSite, int Shift, double* g, double* h)
	{
		double r;

		double* pCoefs = Coefs + Shift + (*CellOfSite + 1) * 4;

		--pCoefs;
		r = *pCoefs;
		*g = *pCoefs * 3;
		*h = *pCoefs * 6;

		--pCoefs;
		r *= *XSite;
		r += *pCoefs;
		*g *= *XSite;
		*g += *pCoefs * 2;
		*h *= *XSite;
		*h += *pCoefs * 2;

		--pCoefs;
		r *= *XSite;
		r += *pCoefs;
		*g *= *XSite;
		*g += *pCoefs;

		--pCoefs;
		r *= *XSite;
		r += *pCoefs;
		return r;
	}

	inline double eval4_1d_111(double* XGrid, double* Coefs, int* CoefsSize,
		double* XSite, int* CellOfSite, int VecIdx, double* g, double* h)
	{
		double XSiteToLeft = *XSite - XGrid[*CellOfSite];
		return receval4_1d_111(Coefs, &XSiteToLeft, CellOfSite, VecIdx*CoefsSize[1], g, h);
	}

	/**
	* Evaluate cubic sline at one particular vector function recursively
	* \param[in] Coefs Coefficients of interpolation, order (vector function, X1, X2, ...).
	* \param[in] CoefsSize Size of coefficients of the above order.
	* \param[in] CurrentDim How many dimension left to be interpolated.
	* \param[in] XSite Evaluation site, already substracted from the left knot points.
	* \param[in] CellOfSite Which pieces does the i-th coordinate of XSite fall in.
	*/
	inline double receval2(double* coefs, int* coefsSize, int xDim, int idim, double* xsite, int* cellOfSite, int shift)
	{
		// recursive evaluation of a piece of coefficient at certain dimension
		double r;
		if (idim == xDim - 1) {
			// last dimension
			double* pCoefs = coefs + shift + (*cellOfSite + 1) * 2;
			r = *(--pCoefs);
			r *= *xsite;
			r += *(--pCoefs);
		}
		else {
			shift += (*cellOfSite + 1) * 2;
			r = receval2(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
			r *= *xsite;
			r += receval2(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
		}
		return r;
	}

	/**
	* Evaluate cubic sline at one particular vector function recursively
	* \param[in] Coefs Coefficients of interpolation, order (vector function, X1, X2, ...).
	* \param[in] CoefsSize Size of coefficients of the above order.
	* \param[in] CurrentDim How many dimension left to be interpolated.
	* \param[in] XSite Evaluation site, already substracted from the left knot points.
	* \param[in] CellOfSite Which pieces does the i-th coordinate of XSite fall in.
	*/
	inline double receval4(double* coefs, int* coefsSize, int xDim, int idim, double* xsite, int* cellOfSite, int shift)
	{
		// recursive evaluation of a piece of coefficient at certain dimension
		double r;
		if (idim == xDim - 1) {
			// last dimension
			double* pCoefs = coefs + shift + (*cellOfSite + 1) * 4;
			r = *(--pCoefs);
			r *= *xsite;
			r += *(--pCoefs);
			r *= *xsite;
			r += *(--pCoefs);
			r *= *xsite;
			r += *(--pCoefs);
		}
		else {
			shift += (*cellOfSite + 1) * 4;
			r = receval4(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
			r *= *xsite;
			r += receval4(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
			r *= *xsite;
			r += receval4(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
			r *= *xsite;
			r += receval4(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
		}
		return r;
	}

	inline double receval4_011(double* coefs, int* coefsSize, int xDim, int idim, double* xsite, int* cellOfSite, int shift, double* grad)
	{
		// recursive evaluation of a piece of coefficient at certain dimension
		double r;
		int sOrder = 4;

		if (idim == xDim - 1) {
			// last dimension
			double* pcoeff = coefs + shift + (*cellOfSite + 1) * sOrder;
			r = *(--pcoeff);

			// gradient
			for (int j = 0; j < xDim; ++j) {
				if (j == idim)
					grad[j] = r*(sOrder - 1);
				else
					grad[j] = r;
			}

			for (int i = 0; i < sOrder - 1; ++i) {
				r *= *xsite;
				r += *(--pcoeff);
				for (int j = 0; j < xDim; ++j) {
					if (j == idim) {
						if (i < sOrder - 2) {
							grad[j] *= *xsite;
							grad[j] += *(pcoeff)* (sOrder - 2 - i);
						}
					}
					else {
						grad[j] *= *xsite;
						grad[j] += *(pcoeff);
					}
				} // j
			} // i
			return r;
		} // idim == xDim-1
		else {
			double gradPrime[MAXDIM];
			shift += (*cellOfSite + 1) * sOrder;
			r = receval4_011(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2], gradPrime);

			for (int j = 0; j < xDim; ++j) {
				if (j == idim)
					grad[j] = gradPrime[j] * (sOrder - 1);
				else
					grad[j] = gradPrime[j];
			}

			for (int i = 0; i < sOrder - 1; ++i) {
				r *= *xsite;
				r += receval4_011(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2], gradPrime);
				for (int j = 0; j < xDim; ++j) {
					if (j == idim) {
						if (i < sOrder - 2) {
							grad[j] *= *xsite;
							grad[j] += gradPrime[j] * (sOrder - 2 - i);
						}
					}
					else {
						grad[j] *= *xsite;
						grad[j] += gradPrime[j];
					}
				}
			}
			return r;
		}
	}

	inline double eval4(int xDim, double** xGrid, double* coefs, int* coefsSize,
		double* xSite, int* cellOfSite, int vecIdx)
	{
		double xSiteToLeft[MAXDIM];
		for (int j = 0; j < xDim; ++j)
		{
			xSiteToLeft[j] = xSite[j] - xGrid[j][cellOfSite[j]];
		}
		return receval4(coefs, coefsSize, xDim, 0, xSiteToLeft, cellOfSite, vecIdx*coefsSize[1]);
	}

	inline double eval4_011(int xDim, double** xGrid, double* coefs, int* coefsSize,
		double* xSite, int* cellOfSite, int vecIdx, double* g)
	{
		double xSiteToLeft[MAXDIM];
		for (int j = 0; j < xDim; ++j)
		{
			xSiteToLeft[j] = xSite[j] - xGrid[j][cellOfSite[j]];
		}
		return receval4_011(coefs, coefsSize, xDim, 0, xSiteToLeft, cellOfSite, vecIdx*coefsSize[1], g);
	}


	inline void hunt(double xx[], int n, double x, int *jlo)
	{
		int jm, jhi, inc;

		if (*jlo <= 0 || *jlo > n) {
			*jlo = 0;
			jhi = n + 1;
		}
		else {
			inc = 1;
			if (x >= xx[*jlo]) {
				if (*jlo == n) return;
				jhi = (*jlo) + 1;
				while (x >= xx[jhi]) {
					*jlo = jhi;
					inc += inc;
					jhi = (*jlo) + inc;
					if (jhi > n) {
						jhi = n + 1;
						break;
					}
				}
			}
			else {
				if (*jlo == 1) {
					*jlo = 0;
					return;
				}
				jhi = (*jlo)--;
				while (x < xx[*jlo]) {
					jhi = (*jlo);
					inc <<= 1;
					if (inc >= jhi) {
						*jlo = 0;
						break;
					}
					else *jlo = jhi - inc;
				}
			}
		}
		while (jhi - (*jlo) != 1) {
			jm = (jhi + (*jlo)) >> 1;
			if (x >= xx[jm])
				*jlo = jm;
			else
				jhi = jm;
		}
		if (x == xx[n]) *jlo = n - 1;
		if (x == xx[1]) *jlo = 1;
	}

	inline void locate(double xx[], int n, double x, int *j)
	{
		int ju, jm, jl;

		jl = 0;
		ju = n + 1;
		while (ju - jl > 1) {
			jm = (ju + jl) >> 1;
			if (x >= xx[jm])
				jl = jm;
			else
				ju = jm;
		}
		if (x == xx[1]) *j = 1;
		else if (x == xx[n]) *j = n - 1;
		else *j = jl;
	}

#define INTERP_DIM MAXDIM
	inline double search_eval4(int xDim, int* xPts, double** xGrid, double* coefs, int* coefsSize,
		double* xSite, int* cellOfSite, int vecIdx)
	{
		double xSiteToLeft[INTERP_DIM];
		for (int j = 0; j < xDim; ++j)
		{
			locate(xGrid[j], xPts[j] - 2, xSite[j], cellOfSite + j);
			xSiteToLeft[j] = xSite[j] - xGrid[j][cellOfSite[j]];
		}
		return receval4(coefs, coefsSize, xDim, 0, xSiteToLeft, cellOfSite, vecIdx*coefsSize[1]);
	}

	inline double search_eval2(int xDim, int* xPts, double** xGrid, double* coefs, int* coefsSize,
		double* xSite, int* cellOfSite, int vecIdx)
	{
		double xSiteToLeft[INTERP_DIM];
		for (int j = 0; j < xDim; ++j)
		{
			locate(xGrid[j], xPts[j] - 2, xSite[j], cellOfSite + j);
			xSiteToLeft[j] = xSite[j] - xGrid[j][cellOfSite[j]];
		}
		return receval2(coefs, coefsSize, xDim, 0, xSiteToLeft, cellOfSite, vecIdx*coefsSize[1]);
	}

	inline double hunt_eval4(int xDim, int* xPts, double** xGrid, double* coefs, int* coefsSize,
		double* xSite, int* cellOfSite, int vecIdx)
	{
		double xSiteToLeft[INTERP_DIM];
		for (int j = 0; j < xDim; ++j)
		{
			hunt(xGrid[j], xPts[j] - 2, xSite[j], cellOfSite + j);
			xSiteToLeft[j] = xSite[j] - xGrid[j][cellOfSite[j]];
		}
		return receval4(coefs, coefsSize, xDim, 0, xSiteToLeft, cellOfSite, vecIdx*coefsSize[1]);
	}

	inline double nosearch_eval2(int xDim, int* xPts, double** xGrid, double* coefs, int* coefsSize,
		double* xSiteToLeft, int* cellOfSite, int vecIdx)
	{
		return receval2(coefs, coefsSize, xDim, 0, xSiteToLeft, cellOfSite, vecIdx*coefsSize[1]);
	}

	inline double nosearch_eval4(int xDim, int* xPts, double** xGrid, double* coefs, int* coefsSize,
		double* xSiteToLeft, int* cellOfSite, int vecIdx)
	{
		return receval4(coefs, coefsSize, xDim, 0, xSiteToLeft, cellOfSite, vecIdx*coefsSize[1]);
	}

	inline double hunt_eval4_011(int xDim, int* xPts, double** xGrid, double* coefs, int* coefsSize,
		double* xSite, int* cellOfSite, int vecIdx, double* g)
	{
		double xSiteToLeft[INTERP_DIM];
		for (int j = 0; j < xDim; ++j)
		{
			hunt(xGrid[j], xPts[j] - 2, xSite[j], cellOfSite + j);
			xSiteToLeft[j] = xSite[j] - xGrid[j][cellOfSite[j]];
		}
		return receval4_011(coefs, coefsSize, xDim, 0, xSiteToLeft, cellOfSite, vecIdx*coefsSize[1], g);
	}

	template<typename T>
	inline int locate2(T* xx, int n, T x)
	{
		// adjust xx[] to 1-based, 
		xx--;
		// return j if x \in [xx[j],xx[j+1])
		// exception, return n-1 if x==xx[n]; only return n if x > xx[n]
		int j;
		int ju, jm, jl;

		jl = 0;
		ju = n + 1;
		while (ju - jl > 1) {
			jm = (ju + jl) >> 1;
			if (x >= xx[jm])
				jl = jm;
			else
				ju = jm;
		}
		if (x == xx[1]) j = 1;
		else if (x == xx[n]) j = n - 1;
		else j = jl;

		return j;
	}

	/**** Matlab-alike Interface ***********/
#define INTERP_DIM MAXDIM
	inline double search_eval4_1d(double* xGrid, int xPts, double* coefs,
		double site, int idx)
	{
		int cellOfSite;
		locate(xGrid, xPts - 2, site, &cellOfSite);
		double XSiteToLeft = site - xGrid[cellOfSite];
		return receval4_1d(coefs, &XSiteToLeft, &cellOfSite, idx*(xPts - 1) * 4);
	}

	inline double nosearch_eval4_1d(int xPts, double* coefs, double XSiteToLeft, int cellOfSite, int idx)
	{
		return receval4_1d(coefs, &XSiteToLeft, &cellOfSite, idx*(xPts - 1) * 4);
	}

	inline double hunt_eval4_1d(double* xGrid, int xPts, double* coefs,
		double site, int idx, int* cellOfSite)
	{
		locate(xGrid, xPts - 2, site, cellOfSite);
		double XSiteToLeft = site - xGrid[*cellOfSite];
		return receval4_1d(coefs, &XSiteToLeft, cellOfSite, idx*(xPts - 1) * 4);
	}

	inline double search_eval4_1d_011(double* xGrid, int xPts, double* coefs,
		double site, int idx, double* g)
	{
		int cellOfSite;
		locate(xGrid, xPts - 2, site, &cellOfSite);
		double XSiteToLeft = site - xGrid[cellOfSite];
		return receval4_1d_011(coefs, &XSiteToLeft, &cellOfSite, idx*(xPts - 1) * 4, g);
	}

	inline double hunt_eval4_1d_011(double* xGrid, int xPts, double* coefs,
		double site, int idx, int* cellOfSite, double* g)
	{
		hunt(xGrid, xPts - 2, site, cellOfSite);
		double XSiteToLeft = site - xGrid[*cellOfSite];
		return receval4_1d_011(coefs, &XSiteToLeft, cellOfSite, idx*(xPts - 1) * 4, g);
	}

	inline void search_1d(double* xGrid, int xPts, double site, double* XSiteToLeft, int* cellOfSite)
	{
		locate(xGrid, xPts - 2, site, cellOfSite);
		*XSiteToLeft = site - xGrid[*cellOfSite];
	}

	inline void hunt_1d(double* xGrid, int xPts, double site, double* XSiteToLeft, int* cellOfSite)
	{
		hunt(xGrid, xPts - 2, site, cellOfSite);
		*XSiteToLeft = site - xGrid[*cellOfSite];
	}
	
	inline double nosearch_eval2_1d(int xPts, double* coefs, double XSiteToLeft, int cellOfSite, int idx)
	{
		return receval2_1d(coefs, &XSiteToLeft, &cellOfSite, idx*(xPts - 1) * 2);
	}

	inline double search_eval2_1d(double* xGrid, int xPts, double* coefs,
		double site, int idx)
	{
		int j;
		int ju, jm, jl;
		int n = xPts - 1;

		jl = 0;
		ju = n + 1;
		while (ju - jl > 1) {
			jm = (ju + jl) >> 1;
			if (site >= xGrid[jm])
				jl = jm;
			else
				ju = jm;
		}
		if (site == xGrid[1]) j = 1;
		else if (site == xGrid[n]) j = n - 1;
		else j = jl;

		double XSiteToLeft = site - xGrid[j];

		int Shift = idx*(xPts - 1) * 2;
		double* pCoefs = coefs + Shift + j*2;
		double r = pCoefs[0] + pCoefs[1] * XSiteToLeft;
		return r;
	}

	template <typename T, int dim>
	T interp1(blitz::Array<T, dim> xGrid, blitz::Array<T, dim> val, T x)
	{
		int nGrid = xGrid.numElements();
		T* _xGrid = xGrid.data();
		int nVal = val.numElements();
		T* _val = val.data();
		assert(nGrid == nVal);
		// Search cell of site
		int cellOfSite = locate2(_xGrid, nGrid, x);

		if (cellOfSite <= 0)
			return _val[0];

		if (cellOfSite >= nGrid)
			return _val[nGrid - 1];

		T rightWeight = (x - _xGrid[cellOfSite - 1]) / (_xGrid[cellOfSite] - _xGrid[cellOfSite - 1]);
		T rslt = _val[cellOfSite - 1] * (1 - rightWeight) + _val[cellOfSite]*rightWeight;

		return rslt;
	}
}

