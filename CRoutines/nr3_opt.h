/**********************************************************
 * This is a modified interface to NR3 optimizers.
 * It's a c interface.
 * It allows to pass function parameters.
 * Function should have the prototype
 * for golden and brent: void f(double x, double* fx, void* params)
 * for dbrent: void f(double x, double* fx, double* dfx, void* params)
 *********************************************************************/
#include <cmath>
#include <cstdio>

/******************************************************************************/

/******************************************************************************/

template <typename Lambda>
inline
double find_local_min ( double a, double b, double eps, double t, 
  Lambda f, double *x )

/******************************************************************************/
/*
  Purpose:

    LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].

  Discussion:

    The method used is a combination of golden section search and
    successive parabolic interpolation.  Convergence is never much slower
    than that for a Fibonacci search.  If F has a continuous second
    derivative which is positive at the minimum (which is not at A or
    B), then convergence is superlinear, and usually of the order of
    about 1.324....

    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
    F is never evaluated at two points closer than TOL.  

    If F is a unimodal function and the computed values of F are always
    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
    LOCAL_MIN approximates the abscissa of the global minimum of F on the 
    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.  

    If F is not unimodal, then LOCAL_MIN may approximate a local, but 
    perhaps non-global, minimum to the same accuracy.

    Thanks to Jonathan Eggleston for pointing out a correction to the 
    golden section step, 01 July 2013.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 July 2013

  Author:

    Orignal FORTRAN77 version by Richard Brent.
    C version by John Burkardt.

  Reference:

    Richard Brent,
    Algorithms for Minimization Without Derivatives,
    Dover, 2002,
    ISBN: 0-486-41998-3,
    LC: QA402.5.B74.

  Parameters:

    Input, double A, B, the endpoints of the interval.

    Input, double EPS, a positive relative error tolerance.
    EPS should be no smaller than twice the relative machine precision,
    and preferably not much less than the square root of the relative
    machine precision.

    Input, double T, a positive absolute error tolerance.

    Input, double F ( double x ), a user-supplied
    function whose local minimum is being sought.

    Output, double *X, the estimated value of an abscissa
    for which F attains a local minimum value in [A,B].

    Output, double LOCAL_MIN, the value F(X).
*/
{
  double c;
  double d;
  double e;
  double fu;
  double fv;
  double fw;
  double fx;
  double m;
  double p;
  double q;
  double r;
  double sa;
  double sb;
  double t2;
  double tol;
  double u;
  double v;
  double w;
/*
  C is the square of the inverse of the golden ratio.
*/
  c = 0.5 * ( 3.0 - sqrt ( 5.0 ) );

  sa = a;
  sb = b;
  *x = sa + c * ( b - a );
  w = *x;
  v = w;
  e = 0.0;
  fx = f ( *x );
  fw = fx;
  fv = fw;

  for ( ; ; )
  { 
    m = 0.5 * ( sa + sb ) ;
    tol = eps * fabs ( *x ) + t;
    t2 = 2.0 * tol;
/*
  Check the stopping criterion.
*/
    if ( fabs ( *x - m ) <= t2 - 0.5 * ( sb - sa ) )
    {
      break;
    }
/*
  Fit a parabola.
*/
    r = 0.0;
    q = r;
    p = q;

    if ( tol < fabs ( e ) )
    {
      r = ( *x - w ) * ( fx - fv );
      q = ( *x - v ) * ( fx - fw );
      p = ( *x - v ) * q - ( *x - w ) * r;
      q = 2.0 * ( q - r );
      if ( 0.0 < q )
      {
        p = - p;
      }
      q = fabs ( q );
      r = e;
      e = d;
    }

    if ( fabs ( p ) < fabs ( 0.5 * q * r ) && 
         q * ( sa - *x ) < p && 
         p < q * ( sb - *x ) )
    {
/*
  Take the parabolic interpolation step.
*/
      d = p / q;
      u = *x + d;
/*
  F must not be evaluated too close to A or B.
*/
      if ( ( u - sa ) < t2 || ( sb - u ) < t2 )
      {
        if ( *x < m )
        {
          d = tol;
        }
        else
        {
          d = - tol;
        }
      }
    }
/*
  A golden-section step.
*/
    else
    {
      if ( *x < m )
      {
        e = sb - *x;
      }
      else
      {
        e = sa - *x;
      }
      d = c * e;
    }
/*
  F must not be evaluated too close to X.
*/
    if ( tol <= fabs ( d ) )
    {
      u = *x + d;
    }
    else if ( 0.0 < d )
    {
      u = *x + tol;
    }
    else
    {
      u = *x - tol;
    }

    fu = f ( u );
/*
  Update A, B, V, W, and X.
*/
    if ( fu <= fx )
    {
      if ( u < *x )
      {
        sb = *x;
      }
      else
      {
        sa = *x;
      }
      v = w;
      fv = fw;
      w = *x;
      fw = fx;
      *x = u;
      fx = fu;
    }
    else
    {
      if ( u < *x )
      {
        sa = u;
      }
      else
      {
        sb = u;
      }

      if ( fu <= fw || w == *x )
      {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      }
      else if ( fu <= fv || v == *x || v == w )
      {
        v = u;
        fv = fu;
      }
    }
  }
  return fx;
}

template <typename Lambda>
inline
double find_zero ( double a, double b, double machep, double t, Lambda f)

	/******************************************************************************/
	/*
	  Purpose:

		ZERO seeks the root of a function F(X) in an interval [A,B].

	  Discussion:

		The interval [A,B] must be a change of sign interval for F.
		That is, F(A) and F(B) must be of opposite signs.  Then
		assuming that F is continuous implies the existence of at least
		one value C between A and B for which F(C) = 0.

		The location of the zero is determined to within an accuracy
		of 6 * MACHEPS * abs ( C ) + 2 * T.

		Thanks to Thomas Secretin for pointing out a transcription error in the
		setting of the value of P, 11 February 2013.

	  Licensing:

		This code is distributed under the GNU LGPL license.

	  Modified:

		11 February 2013

	  Author:

		Original FORTRAN77 version by Richard Brent.
		C version by John Burkardt.

	  Reference:

		Richard Brent,
		Algorithms for Minimization Without Derivatives,
		Dover, 2002,
		ISBN: 0-486-41998-3,
		LC: QA402.5.B74.

	  Parameters:

		Input, double A, B, the endpoints of the change of sign interval.

		Input, double MACHEP, an estimate for the relative machine
		precision.

		Input, double T, a positive error tolerance.

		Input, double F ( double x ), a user-supplied function whose zero
		is being sought.

		Output, double ZERO, the estimated value of a zero of
		the function F.
	*/
{
	double c;
	double d;
	double e;
	double fa;
	double fb;
	double fc;
	double m;
	double p;
	double q;
	double r;
	double s;
	double sa;
	double sb;
	double tol;
	/*
	  Make local copies of A and B.
	*/
	sa = a;
	sb = b;
	fa = f(sa);
	fb = f(sb);

	c = sa;
	fc = fa;
	e = sb - sa;
	d = e;

	for (; ; )
	{
		if (fabs(fc) < fabs(fb))
		{
			sa = sb;
			sb = c;
			c = sa;
			fa = fb;
			fb = fc;
			fc = fa;
		}

		tol = 2.0 * machep * fabs(sb) + t;
		m = 0.5 * (c - sb);

		if (fabs(m) <= tol || fb == 0.0)
		{
			break;
		}

		if (fabs(e) < tol || fabs(fa) <= fabs(fb))
		{
			e = m;
			d = e;
		}
		else
		{
			s = fb / fa;

			if (sa == c)
			{
				p = 2.0 * m * s;
				q = 1.0 - s;
			}
			else
			{
				q = fa / fc;
				r = fb / fc;
				p = s * (2.0 * m * q * (q - r) - (sb - sa) * (r - 1.0));
				q = (q - 1.0) * (r - 1.0) * (s - 1.0);
			}

			if (0.0 < p)
			{
				q = -q;
			}
			else
			{
				p = -p;
			}

			s = e;
			e = d;

			if (2.0 * p < 3.0 * m * q - fabs(tol * q) &&
				p < fabs(0.5 * s * q))
			{
				d = p / q;
			}
			else
			{
				e = m;
				d = e;
			}
		}
		sa = sb;
		fa = fb;

		if (tol < fabs(d))
		{
			sb = sb + d;
		}
		else if (0.0 < m)
		{
			sb = sb + tol;
		}
		else
		{
			sb = sb - tol;
		}

		fb = f(sb);

		if ((0.0 < fb && 0.0 < fc) || (fb <= 0.0 && fc <= 0.0))
		{
			c = sa;
			fc = fa;
			e = sb - sa;
			d = e;
		}
	}
	return sb;
}

#define NRANSI
#define R 0.61803399
#define C (1.0-R)
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

template <typename Lambda>
double golden(double ax, double bx, double cx, Lambda f, void* params, double tol, double *xmin)
{
	double f1, f2, x0, x1, x2, x3;
	int iter = 0;

	x0 = ax;
	x3 = cx;
	if (fabs(cx - bx) > fabs(bx - ax)) {
		x1 = bx;
		x2 = bx + C*(cx - bx);
	}
	else {
		x2 = bx;
		x1 = bx - C*(bx - ax);
	}
	f1 = f(x1, params);
	f2 = f(x2, params);
	while (fabs(x3 - x0) > tol) {
		iter++;
		if (f2 < f1) {
			SHFT3(x0, x1, x2, R*x1 + C*x3)
				SHFT2(f1, f2, f(x2, params))
		}
		else {
			SHFT3(x3, x2, x1, R*x2 + C*x0)
				SHFT2(f2, f1, f(x1, params))
		}
	}
	if (f1 < f2) {
		*xmin = x1;
		return f1;
	}
	else {
		*xmin = x2;
		return f2;
	}
}

template <typename Lambda>
double brent(double ax, double bx, double cx, Lambda f, void* params, double tol, double *xmin)
{
	int iter;
	double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
	double e = 0.0;

	a = (ax < cx ? ax : cx);
	b = (ax > cx ? ax : cx);
	x = w = v = bx;
	fw = fv = fx = f(x, params);
	while (1) {
		xm = 0.5*(a + b);
		tol2 = 2.0*(tol1 = tol*fabs(x) + ZEPS);
		if (fabs(x - xm) <= (tol2 - 0.5*(b - a))) {
			*xmin = x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r = (x - w)*(fx - fv);
			q = (x - v)*(fx - fw);
			p = (x - v)*q - (x - w)*r;
			q = 2.0*(q - r);
			if (q > 0.0) p = -p;
			q = fabs(q);
			etemp = e;
			e = d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a - x) || p >= q*(b - x))
				d = CGOLD*(e = (x >= xm ? a - x : b - x));
			else {
				d = p / q;
				u = x + d;
				if (u - a < tol2 || b - u < tol2)
					d = SIGN(tol1, xm - x);
			}
		}
		else {
			d = CGOLD*(e = (x >= xm ? a - x : b - x));
		}
		u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
		fu = f(u, params);
		if (fu <= fx) {
			if (u >= x) a = x; else b = x;
			SHFT(v, w, x, u)
				SHFT(fv, fw, fx, fu)
		}
		else {
			if (u < x) a = u; else b = u;
			if (fu <= fw || w == x) {
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if (fu <= fv || v == x || v == w) {
				v = u;
				fv = fu;
			}
		}
	}
	*xmin = x;
	return fx;
}

/*
double dbrent(double ax, double bx, double cx, double (*f)(double, double*, void*), void* params, double tol, double *xmin)
{
	int iter,ok1,ok2;
	double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
	double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
	double dfx, dfu;


	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;

	fw=fv=fx=(*f)(x, &dfx, params);
	dw=dv=dx=dfx;
	for (iter=1;iter<=ITMAX;iter++) {
	xm=0.5*(a+b);
	tol1=tol*fabs(x)+ZEPS;
	tol2=2.0*tol1;
	if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
		*xmin=x;
		return fx;
	}
	if (fabs(e) > tol1) {
		d1=2.0*(b-a);
		d2=d1;
		if (dw != dx) d1=(w-x)*dx/(dx-dw);
		if (dv != dx) d2=(v-x)*dx/(dx-dv);
		u1=x+d1;
		u2=x+d2;
		ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
		ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
		olde=e;
		e=d;
		if (ok1 || ok2) {
		if (ok1 && ok2)
			d=(fabs(d1) < fabs(d2) ? d1 : d2);
		else if (ok1)
			d=d1;
		else
			d=d2;
		if (fabs(d) <= fabs(0.5*olde)) {
			u=x+d;
			if (u-a < tol2 || b-u < tol2)
			d=SIGN(tol1,xm-x);
		} else {
			d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
		} else {
		d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
	} else {
		d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	}
	if (fabs(d) >= tol1) {
		u=x+d;
		fu=(*f)(u, &dfu, params);
	} else {
		u=x+SIGN(tol1,d);
		fu=(*f)(u, &dfu, params);
		if (fu > fx) {
		*xmin=x;
		return fx;
		}
	}
	du=dfu;
	if (fu <= fx) {
		if (u >= x) a=x; else b=x;
		MOV3(v,fv,dv, w,fw,dw)
		MOV3(w,fw,dw, x,fx,dx)
		MOV3(x,fx,dx, u,fu,du)
	} else {
		if (u < x) a=u; else b=u;
		if (fu <= fw || w == x) {
		MOV3(v,fv,dv, w,fw,dw)
			MOV3(w,fw,dw, u,fu,du)
		} else if (fu < fv || v == x || v == w) {
		MOV3(v,fv,dv, u,fu,du)
		}
	}
	}
	return 0.0;
}
*/

#undef C
#undef R
#undef SHFT2
#undef SHFT3
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef MOV3
#undef NRANSI
#undef SIGN

template <typename Lambda, typename T>
inline
double rtbis(Lambda func, T x1, T x2, T xacc, int maxiter, int* exitflag, int* iter)
{
	int j;
	double dx, f, fmid, xmid, rtb;
	*exitflag = -2;

	f = func(x1);
	fmid = func(x2);
	if (f*fmid >= 0.0) {
		*exitflag = -2;
		// Wenlan: handle corner
		if (f < 0.0)
			return x1;
		else
			return x2;
	}
	rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
	for ((*iter) = 1; (*iter) <= maxiter; (*iter)++) {
		fmid = func(xmid = rtb + (dx *= 0.5));
		if (fmid <= 0.0) rtb = xmid;
		if (fabs(dx) < xacc || fmid == 0.0) {
			*exitflag = 1;
			return rtb;
		}
	}
	*exitflag = 0;
	return 0.0;
}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define EPS 3.0e-8
template <typename Lambda, typename T>
T zbrent(Lambda func, T x1, T x2, T tol, int maxiter, int* exitflag, int* iter)
{
	T a = x1, b = x2, c = x2, d, e, min1, min2;
	T fa = func(a), fb = func(b), fc, p, q, r, s, tol1, xm;
	*exitflag = 1;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && b < 0.0)) {
		*exitflag = -2;
		// Wenlan: handle corner
		if (fa < 0.0)
			return x1;
		else
			return x2;
	}
	fc = fb;
	for ((*iter) = 1; (*iter) <= maxiter; (*iter)++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c = a;
			fc = fa;
			e = d = b - a;
		}
		if (fabs(fc) < fabs(fb)) {
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		tol1 = 2.0*EPS*fabs(b) + 0.5*tol;
		xm = 0.5*(c - b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s = fb / fa;
			if (a == c) {
				p = 2.0*xm*s;
				q = 1.0 - s;
			}
			else {
				q = fa / fc;
				r = fb / fc;
				p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
				q = (q - 1.0)*(r - 1.0)*(s - 1.0);
			}
			if (p > 0.0) q = -q;
			p = fabs(p);
			min1 = 3.0*xm*q - fabs(tol1*q);
			min2 = fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e = d;
				d = p / q;
			}
			else {
				d = xm;
				e = d;
			}
		}
		else {
			d = xm;
			e = d;
		}
		a = b;
		fa = fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1, xm);
		fb = func(b);
	}
	*exitflag = 0;
	return 0.0;
}
#undef SIGN

template <typename Lambda, class T>
T rtsafe(Lambda funcd, T x0, T x1, T x2, T tol, int maxiter, int* exitflag, int* iter) {
	T xh, xl;
	*exitflag = 1;

	T fl = funcd(x1, 0);
	T fh = funcd(x2, 0);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
		*exitflag = -2;
		// Wenlan: handle corner
		if (fl < 0.0)
			return x1;
		else
			return x2;
	}
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl = x1;
		xh = x2;
	}
	else {
		xh = x1;
		xl = x2;
	}
	T rts = x0;
	T dxold = abs(x2 - x1);
	T dx = dxold;
	T df;
	T f = funcd(rts, &df);
	for ((*iter) = 1; (*iter) <= maxiter; (*iter)++) {
		if ((((rts - xh)*df - f)*((rts - xl)*df - f) > 0.0)
			|| (abs(2.0*f) > abs(dxold*df))) {
			dxold = dx;
			dx = 0.5*(xh - xl);
			rts = xl + dx;
			if (xl == rts) return rts;
		}
		else {
			dxold = dx;
			dx = f / df;
			T temp = rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (abs(dx) < tol) return rts;
		f = funcd(rts, &df);
		if (f < 0.0)
			xl = rts;
		else
			xh = rts;
	}
	*exitflag = 0;
	return 0.0;
}

