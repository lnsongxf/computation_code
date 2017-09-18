function [kpPolicy,v] = solve_decision_func_approx_vec(w,r,params)
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
    
    kLower = kMin*ones(ePts,kPts);
    budget = reshape(kGrid,[1,kPts])*(1+r) + reshape(eGrid,[ePts,1])*w;
    
    % For each e and k
    [kpPolicy,v_new] = goldenx(@(kp)value_of_kp(kp(:),budget,vFuture,params),kLower(:),budget(:)-1e-12,1000,TOL_OPT);
    kpPolicy = reshape(kpPolicy,[ePts,kPts]);
    v_new = reshape(v_new,[ePts,kPts]);
    
    % Update v
    metric = max(abs(v_new(:)-v(:)));
    v = v_new;
    iter = iter+1;
    
    if mod(iter,PRINT_FREQ) == 0
        fprintf('iter: %d, metric: %g\n',iter,metric);
    end
end

end

function v = value_of_kp(kp,budget,vFuture,params)
gamma = params.gamma;
eGrid = params.eGrid;
kGrid = params.kGrid;
ePts = params.ePts;
kPts = params.kPts;

% Compute consumption given kp
kp = reshape(kp,[ePts,kPts]);
c = budget - kp;

% interp1 does the vectorized interpolation, return size [ePts,kpPts]
v = zeros(ePts,kPts);
for i_e = 1:ePts
    v(i_e,:) = u(c(i_e,:),gamma) + interp1(kGrid,vFuture(i_e,:),kp(i_e,:),'pchip','extrap');
end

v = v(:);
end