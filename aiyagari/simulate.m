function dist = simulate(kpPolicy,rn,w,r,params)
MAXITER_VFI=params.MAXITER_VFI; NUM_THREADS=params.NUM_THREADS; PRINT_FREQ=params.PRINT_FREQ; T=params.T; TOL_EQ=params.TOL_EQ; TOL_OPT=params.TOL_OPT; TOL_VFI=params.TOL_VFI; alpha=params.alpha; beta=params.beta; delta=params.delta; eGrid=params.eGrid; ePts=params.ePts; eRange=params.eRange; eRho=params.eRho; eSigma=params.eSigma; eTrans=params.eTrans; gamma=params.gamma; kGrid=params.kGrid; kMax=params.kMax; kMin=params.kMin; kPts=params.kPts; kShift=params.kShift; numAgents=params.numAgents; 

% Unpack shocks
eShockIdx = rn.eShockIdx;

% Construct time series
k_t = zeros(numAgents,T);

% Simulate forward
for t=1:T-1
    % Find samples with shock index
    for i_e=1:ePts
        idxOfShockJ = find(eShockIdx(:,t)==i_e);
        % Interpolate using corresponding policy function
        k_t(idxOfShockJ,t+1) = interp1(kGrid,kpPolicy(i_e,:),k_t(idxOfShockJ,t),'pchip');
    end
end

inc_t = k_t*r + w*eGrid(eShockIdx);

dist.k_t = k_t;
dist.inc_t = inc_t;
end