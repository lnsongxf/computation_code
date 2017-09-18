function [kpPolicy,v] = solve_decision_grid_search(w,r,params)
% Aiyagari (1994)
% Solve decision problem using grid search
% Wenlan Luo, luowenlan@gmail.com
NUM_THREADS=params.NUM_THREADS; PRINT_FREQ=params.PRINT_FREQ; TOL_EQ=params.TOL_EQ; TOL_OPT=params.TOL_OPT; TOL_VFI=params.TOL_VFI; alpha=params.alpha; beta=params.beta; delta=params.delta; eGrid=params.eGrid; ePts=params.ePts; eRange=params.eRange; eRho=params.eRho; eSigma=params.eSigma; eTrans=params.eTrans; gamma=params.gamma; kGrid=params.kGrid; kMax=params.kMax; kMin=params.kMin; kPts=params.kPts; kShift=params.kShift;

% Initialize value function
v = zeros(ePts,kPts);

% Compute consumption for each (kp,e,k).
% matalb supports broadcasting after R2016b
c = reshape(kGrid,[1,1,kPts])*(1+r) + reshape(eGrid,[1,ePts,1])*w - reshape(kGrid,[kPts,1,1]);

% Start iterations
metric = 1;
iter = 0;
while (metric > TOL_VFI && iter<=inf)
    % Compute value for each (kp,e,k)
    vFuture = permute(beta*eTrans*v,[2,1]);
    valueOfKp = u(c,gamma) + reshape(vFuture,[kPts,ePts,1]);
    valueOfKp(c<0) = -1e20;
    
    % Compute max w.r.t. kp
    [v_new,argmaxKpIdx] = max(valueOfKp,[],1);
    v_new = reshape(v_new,[ePts,kPts]);
    
    % Extract policy
    kpPolicy = kGrid(reshape(argmaxKpIdx,[ePts,kPts]));
    
    % Update v
    metric = max(abs(v_new(:)-v(:)));
    v = v_new;
    iter = iter+1;
    
    if mod(iter,PRINT_FREQ) == 0
        fprintf('iter: %d, metric: %g\n',iter,metric);
    end
end

end