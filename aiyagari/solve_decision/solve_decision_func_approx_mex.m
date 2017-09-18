function [kpPolicy,v] = solve_decision_func_approx_mex(w,r,params)
% Aiyagari (1994)
% Solve decision problem
% Wenlan Luo, luowenlan@gmail.com
MAXITER_VFI=params.MAXITER_VFI; NUM_THREADS=params.NUM_THREADS; PRINT_FREQ=params.PRINT_FREQ; TOL_EQ=params.TOL_EQ; TOL_OPT=params.TOL_OPT; TOL_VFI=params.TOL_VFI; alpha=params.alpha; beta=params.beta; delta=params.delta; eGrid=params.eGrid; ePts=params.ePts; eRange=params.eRange; eRho=params.eRho; eSigma=params.eSigma; eTrans=params.eTrans; gamma=params.gamma; kGrid=params.kGrid; kMax=params.kMax; kMin=params.kMin; kPts=params.kPts; kShift=params.kShift; 

% Initialize value function
v = zeros(ePts,kPts);

% Start iterations
metric = 1;
iter = 0;
while (metric > TOL_VFI && iter<=MAXITER_VFI)
    % Continution value
    vFuture = beta*eTrans*v;
    
    % Allocate space for value and policy
    kpPolicy = zeros(ePts,kPts);
    v_new = zeros(ePts,kPts);
    
    % Construct interp coefficients
    vFutureInterp = tensor_pchip({kGrid},vFuture);
    vFutureInterpCForm = myppual(vFutureInterp);
    vFutureCoefs = vFutureInterpCForm.coefs;
    
    % Call mex
    vfi_mex;
    
    % Update v
    metric = max(abs(v_new(:)-v(:)));
    v = v_new;
    iter = iter+1;
    
    if mod(iter,PRINT_FREQ) == 0
        fprintf('iter: %d, metric: %g\n',iter,metric);
    end
end

end