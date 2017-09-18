function [kpPolicy,v] = solve_decision_func_approx(w,r,params)
% Aiyagari (1994)
% Solve decision problem
% Wenlan Luo, luowenlan@gmail.com
MAXITER_VFI=params.MAXITER_VFI; NUM_THREADS=params.NUM_THREADS; PRINT_FREQ=params.PRINT_FREQ; TOL_EQ=params.TOL_EQ; TOL_OPT=params.TOL_OPT; TOL_VFI=params.TOL_VFI; alpha=params.alpha; beta=params.beta; delta=params.delta; eGrid=params.eGrid; ePts=params.ePts; eRange=params.eRange; eRho=params.eRho; eSigma=params.eSigma; eTrans=params.eTrans; gamma=params.gamma; kGrid=params.kGrid; kMax=params.kMax; kMin=params.kMin; kPts=params.kPts; kShift=params.kShift; 

% Initialize value function
v = zeros(ePts,kPts);

% Start iterations
metric = 1;
iter = 0;
solverOptions = optimset('TolX',TOL_OPT);
while (metric > TOL_VFI && iter<=MAXITER_VFI)
    % Continution value
    vFuture = beta*eTrans*v;
    
    % Allocate space for value and policy
    kpPolicy = zeros(ePts,kPts);
    v_new = zeros(ePts,kPts);
    
    % For each e and k
    for i_e = 1:ePts
        for i_k = 1:kPts
            k = kGrid(i_k);
            e = eGrid(i_e);
            kpUpper = k*(1+r) + e*w;
            [kpPolicy(i_e,i_k),v_new(i_e,i_k)] = ...
                fminbnd(@(kp) value_of_kp(kp,i_e,i_k,w,r,vFuture,params),kMin,kpUpper,solverOptions);
        end
    end
    
    % value_of_kp returns negative of v;
    v_new = -v_new;
    
    % Update v
    metric = max(abs(v_new(:)-v(:)));
    v = v_new;
    iter = iter+1;
    
    if mod(iter,PRINT_FREQ) == 0
        fprintf('iter: %d, metric: %g\n',iter,metric);
    end
end

end

function negativeV = value_of_kp(kp,i_e,i_k,w,r,vFuture,params)
gamma = params.gamma;
eGrid = params.eGrid;
kGrid = params.kGrid;

k = kGrid(i_k);
e = eGrid(i_e);

% Compute consumption given kp
c = k*(1+r) + e*w - kp;

% Negative consumption, return infeasible value
if c<0
    negativeV = 1e20;
    return;
end


v = u(c,gamma) + interp1(kGrid,vFuture(i_e,:),kp,'linear','extrap');

negativeV = -v;
end