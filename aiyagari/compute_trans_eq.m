function trans = compute_trans_eq(ssEq,dist0,params)
MAXITER_VFI=params.MAXITER_VFI; NUM_THREADS=params.NUM_THREADS; PRINT_FREQ=params.PRINT_FREQ; T=params.T; TOL_DIST=params.TOL_DIST; TOL_EQ=params.TOL_EQ; TOL_OPT=params.TOL_OPT; TOL_VFI=params.TOL_VFI; alpha=params.alpha; beta=params.beta; delta=params.delta; eGrid=params.eGrid; ePts=params.ePts; eRange=params.eRange; eRho=params.eRho; eSigma=params.eSigma; eTrans=params.eTrans; gamma=params.gamma; kGrid=params.kGrid; kMax=params.kMax; kMin=params.kMin; kPts=params.kPts; kShift=params.kShift; numAgents=params.numAgents; transT=params.transT; 

% Extract steady state margial value and distirbution
vp0 = ssEq.vp;
kpPolicy0 = ssEq.kpPolicy;
L = ssEq.aggRslt.L;

% Allocate space for transition path
vp_t = zeros([size(vp0),transT+1]);
dist_t = zeros([size(dist0),transT+1]);
kpPolicy_t = zeros([size(kpPolicy0),transT]);
K_t = ssEq.aggRslt.K*ones(1,transT);
K_t_new = K_t;

% Initiate value and ditribution
vp_t(:,:,transT+1) = vp0;
dist_t(:,:,1) = dist0;

% Start loop by searching K_t
metric = 1;
iter = 0;
TOL_TRANS_EQ = 1e-6;

while metric > TOL_TRANS_EQ
    % Compute prices given K_t
    KL_t = K_t / L;
    r_t = alpha* KL_t.^(alpha-1) - delta;
    w_t = (1-alpha)* KL_t.^(alpha);
    
    % Backward iteration for solving the decision problem
    for t=transT:-1:1
        % Continution value
        vpFuture = beta*eTrans*vp_t(:,:,t+1);
        [kpPolicy_t(:,:,t), vp_t(:,:,t)] = solve_decision_egm_1(vpFuture,w_t(t),r_t(t),params);
    end
    
    % Forward iteration for distribution
    for t=1:1:transT
        dist_t(:,:,t+1) = simulate_histogram_1(dist_t(:,:,t),kpPolicy_t(:,:,t),w_t(t),r_t(t),params);
        aggRslt = aggregate_histogram(dist_t(:,:,t),w_t(t),r_t(t),params);
        K_t_new(t) = aggRslt.K;
    end
    
    % Compute metric;
    metric = max(abs(K_t(:)-K_t_new(:)));
    iter = iter+1;
    % Do a mannual line search
    TRANS_UPDATE_SPEED = 0.2;
    PRINT_FREQ_TRANS_EQ = 20;
    K_t = K_t_new*TRANS_UPDATE_SPEED + K_t*(1-TRANS_UPDATE_SPEED);
    
    % Print something
    if mod(iter,PRINT_FREQ_TRANS_EQ)==0
        fprintf('iter: %d, metric: %g\n',iter, metric);
    end
end
fprintf('iter: %d, metric: %g\n',iter, metric);

trans = v2struct(vp_t,dist_t,kpPolicy_t,K_t);
end