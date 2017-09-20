function [kpPolicy,vp] = solve_decision_egm_1(vpFuture,w,r,params)
MAXITER_VFI=params.MAXITER_VFI; NUM_THREADS=params.NUM_THREADS; PRINT_FREQ=params.PRINT_FREQ; T=params.T; TOL_DIST=params.TOL_DIST; TOL_EQ=params.TOL_EQ; TOL_OPT=params.TOL_OPT; TOL_VFI=params.TOL_VFI; alpha=params.alpha; beta=params.beta; delta=params.delta; eGrid=params.eGrid; ePts=params.ePts; eRange=params.eRange; eRho=params.eRho; eSigma=params.eSigma; eTrans=params.eTrans; gamma=params.gamma; kGrid=params.kGrid; kMax=params.kMax; kMin=params.kMin; kPts=params.kPts; kShift=params.kShift; numAgents=params.numAgents; transT=params.transT; 

% Endogenous grid
kTilde = (up_inverse(vpFuture,gamma) + reshape(kGrid,[1,kPts]) - reshape(eGrid,[ePts,1])*w) / (1+r);
kpPolicy = zeros(ePts,kPts);
for i_e=1:ePts
    % Interpolate kTilde over k
    kpPolicy(i_e,:) = interp1(kTilde(i_e,:),kGrid,kGrid,'pchip');
    % We know kp(i_e,k)=0 for k <= kTilde(i_e,0), and the interpolation
    % will give kp(i_e,k)<0 for k<=kTilde(i_e,0)
end
kpPolicy = max(kpPolicy,kMin);

% update vp using envelope theorem
budget = reshape(kGrid,[1,kPts])*(1+r) + reshape(eGrid,[ePts,1])*w;
c = budget - kpPolicy;
vp  = up(c,gamma)*(1+r);
end