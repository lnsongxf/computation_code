function [kpPolicy,vp] = solve_decision_egm(w,r,params)
% Aiyagari (1994)
% Solve decision problem
% Wenlan Luo, luowenlan@gmail.com
MAXITER_VFI=params.MAXITER_VFI; NUM_THREADS=params.NUM_THREADS; PRINT_FREQ=params.PRINT_FREQ; TOL_EQ=params.TOL_EQ; TOL_OPT=params.TOL_OPT; TOL_VFI=params.TOL_VFI; alpha=params.alpha; beta=params.beta; delta=params.delta; eGrid=params.eGrid; ePts=params.ePts; eRange=params.eRange; eRho=params.eRho; eSigma=params.eSigma; eTrans=params.eTrans; gamma=params.gamma; kGrid=params.kGrid; kMax=params.kMax; kMin=params.kMin; kPts=params.kPts; kShift=params.kShift; 

% Initialize marginal value
budget = reshape(kGrid,[1,kPts])*(1+r) + reshape(eGrid,[ePts,1])*w;
vp = up(budget,gamma)*(1+r);

% Start iterations
metric = 1;
iter = 0;
solverOptions = optimset('TolX',TOL_OPT);
while (metric > TOL_VFI && iter<=MAXITER_VFI)
    % Continution value
    vpFuture = beta*eTrans*vp;
    
    % Endogenous grid
    kTilde = (up_inverse(vpFuture,gamma) + reshape(kGrid,[1,kPts]) - reshape(eGrid,[ePts,1])*w) / (1+r);
    
    for i_e=1:ePts
        % Interpolate kTilde over k
        kpPolicy(i_e,:) = interp1(kTilde(i_e,:),kGrid,kGrid,'pchip');
        % We know kp(i_e,k)=0 for k <= kTilde(i_e,0), and the interpolation
        % will give kp(i_e,k)<0 for k<=kTilde(i_e,0)
    end
    kpPolicy = max(kpPolicy,kMin);
    
    % update vp using envelope theorem
    c = budget - kpPolicy;
    vp_new = up(c,gamma)*(1+r);
    
    % Update vp
    metric = max(abs(vp_new(:)-vp(:)));
    vp = vp_new;
    iter = iter+1;
    
    if mod(iter,PRINT_FREQ) == 0
        fprintf('iter: %d, metric: %g\n',iter,metric);
    end
end

fprintf('iter: %d, metric: %g\n',iter,metric);

end