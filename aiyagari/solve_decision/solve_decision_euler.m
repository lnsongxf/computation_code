function [kpPolicy,vp] = solve_decision_euler(w,r,params)
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
    
    % Allocate space for value and policy
    kpPolicy = zeros(ePts,kPts);
    
    % For each e and k
    for i_e = 1:ePts
        for i_k = 1:kPts
            k = kGrid(i_k);
            e = eGrid(i_e);
            kpUpper = k*(1+r) + e*w-1e-12;
            
            if euler_of_kp(kMin,i_e,i_k,w,r,vpFuture,params) > 0
                kpPolicy(i_e,i_k) = kMin;
            else
                kpPolicy(i_e,i_k) = ...
                    fzero(@(kp) euler_of_kp(kp,i_e,i_k,w,r,vpFuture,params),[kMin,kpUpper],solverOptions);
            end
        end
    end
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

end

function residual = euler_of_kp(kp,i_e,i_k,w,r,vpFuture,params)
gamma = params.gamma;
eGrid = params.eGrid;
kGrid = params.kGrid;

k = kGrid(i_k);
e = eGrid(i_e);

% Compute consumption given kp
c = k*(1+r) + e*w - kp;

residual = up(c,gamma) - interp1(kGrid,vpFuture(i_e,:),kp,'linear','extrap');
end