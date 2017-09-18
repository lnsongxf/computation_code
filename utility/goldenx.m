function [x1,f1,exitflag] = goldenx(f,a,b,maxiter,tol,varargin)

% GOLDENX Computes local maximum of univariate function on interval via Golden Search
%   Vectorized version of Golden
% USAGE
%   [x,fval] = goldenx(f,a,b,maxiter,tol,varargin)
% INPUTS
%   f         : name of function of form fval=f(x)
%   a,b       : left, right endpoints of interval
%   maxiter   : maximum number of iterations
%   tol       : convergence tolerance
%   P1,P2,... : optional additional arguments for f
% OUTPUTS
%   x       : local maximum of f
%   fval    : function value estimate
% exitflag  : conditions of exit (1 = converged; 0 = max iteration achieved)

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

% Revised by Jinhui Bai, January 15, 2005;Feb 14, 2005

% check bracketing point 
if isequal(size(a),size(b))==0 
    error('Bracketing vectors must be of the same length.')
end

if any(((b-a)>0)==0)
    error('Lower and upper interval points are not in the correct sequence.')
end

% main code
alpha1 = (3-sqrt(5))/2;             % golden ratio at the left
alpha2 = 1-alpha1;                  % golden ratio at the right
alpha12 = alpha2 - alpha1;          % distance between two ratio
d  = b-a;                           % length of the interval
x1 = a+alpha1*d;                    % left section point
x2 = a+alpha2*d;                    % right section point
s  = ones(size(x1));                % the sign of difference between x1 and x2
f1 = feval(f,x1,varargin{:});       % fcn value at left point
f2 = feval(f,x2,varargin{:});       % fcn value at right point

iter = 1;
while any(d>tol) && (iter <= maxiter)
  d = d*alpha2;                     % shrink interval length
  i = f2>f1;                        % indicator for f2>f1
  x1(i) = x2(i);                    % update boundary point
  f1(i) = f2(i);                    % update fcn value  
  d12 = alpha12*d;                  % distance between two section points
  x2 = x1+s.*(i-(~i)).*d12;         % assign value to x2 by distance
  f2 = feval(f,x2,varargin{:});     % evaluate x2 
  s = sign(x2-x1);                  % indicator on left or right boundary
  iter = iter + 1;
end

% Return the larger of the two choice
i = f2>f1; x1(i) = x2(i);  f1(i) = f2(i);

% check that endpoints are less than the maximum found
% take the maximum of the two boundaries
fa = feval(f,a,varargin{:});        % value at the left boundary
fb = feval(f,b,varargin{:});        % value at the right boundary
i = fb > fa; a(i) = b(i);  fa(i) = fb(i);
% take the maximum of x1 and a
i = fa > f1; x1(i) = a(i);  f1(i) = fa(i); 

if all(d <= tol)
    exitflag = 1;
else
    exitflag = 0;
end;
