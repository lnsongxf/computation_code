function simuRslt = simulate(vfiRslt, rn, params)
% Replication of Krusell and Smith (1998)
% Author: Wenlan Luo
% SIMULATE: simulate samples forward
NUM_THREADS=params.NUM_THREADS; PRINT_FREQ=params.PRINT_FREQ; PRINT_FREQ_SIMULATE=params.PRINT_FREQ_SIMULATE; PROFILE=params.PROFILE; TOL_OPT=params.TOL_OPT; TOL_VFI=params.TOL_VFI; ZETrans=params.ZETrans; alpha=params.alpha; betaGrid=params.betaGrid; betaPts=params.betaPts; betaTrans=params.betaTrans; delta=params.delta; eBetaEPrimeBetaPrimeZZprime=params.eBetaEPrimeBetaPrimeZZprime; eBetaTrans=params.eBetaTrans; eBetaTransConditionalOnZ=params.eBetaTransConditionalOnZ; eBetaTransInv=params.eBetaTransInv; eEPrimeZZPrime=params.eEPrimeZZPrime; eGrid=params.eGrid; ePts=params.ePts; eTransConditionalOnZ=params.eTransConditionalOnZ; fullDiscountedTrans=params.fullDiscountedTrans; fullTrans=params.fullTrans; i_beta=params.i_beta; i_betap=params.i_betap; i_e=params.i_e; i_ep=params.i_ep; i_z=params.i_z; i_zp=params.i_zp; kBarGrid=params.kBarGrid; kBarPts=params.kBarPts; kGrid=params.kGrid; kMax=params.kMax; kMin=params.kMin; kPts=params.kPts; noWorkWage=params.noWorkWage; numOfAgents=params.numOfAgents; numOfPeriods=params.numOfPeriods; numOfPeriodsRemoved=params.numOfPeriodsRemoved; params=params.params; phi=params.phi; zGrid=params.zGrid; zPts=params.zPts; zTrans=params.zTrans; 

% Unpack from structure
kpPolicy = vfiRslt.kpPolicy;
zIdx = rn.zIdx;
eBetaIdx = rn.eBetaIdx;

% Construct full coefficient w.r.t. k and kBar
% permute kpPolicy to have order (column major) [kBar,k,e,beta,z]
kpPolicyPermute = permute(kpPolicy,[5,4,2,3,1]);
kpPolicyKBarSplineCoefs = myppualMKL_CMEX(int32(-NUM_THREADS), {kGrid,kBarGrid}, kpPolicyPermute, [], int32([4,4]), int32(zPts*ePts*betaPts), [], [], [], []);
kpPolicyKBarSplineCoefsByZ = reshape(kpPolicyKBarSplineCoefs,[],zPts);

% Initiate simulation
k_t = 10*ones(numOfAgents,numOfPeriods);
kBar_t = zeros(1,numOfPeriods);
for t=1:numOfPeriods
    currentZIdx = zIdx(t);
    kBar_t(t) = mean(k_t(:,t));
    
    if t<numOfPeriods
        % Do a reduction at kBar
        kpPolicyKBarSplineCoefsAtCurrentZ = kpPolicyKBarSplineCoefsByZ(:,currentZIdx);
        kpPolicyKSplineCoefs = myppualMKL_CMEX(int32(NUM_THREADS), {kBarGrid}, kpPolicyKBarSplineCoefsAtCurrentZ, ...
            [], int32(4), int32(ePts*betaPts*(kPts-1)*4), [], kBar_t(t), [], [], []);
        
        % Do an evaluation at k
        k_t(:,t+1) = myppualMKL_CMEX(int32(NUM_THREADS), {kGrid}, kpPolicyKSplineCoefs, ...
            [], int32(4), int32(ePts*betaPts), [], k_t(:,t)', [], int32(eBetaIdx(:,t)')-1, []);
    end
    
    if mod(t,PRINT_FREQ_SIMULATE)==0
        fprintf('%-10s%-10s%-10s%-10s%-10s\n','t','kMin','kMax','kBar','Gini');
        fprintf('%-10d%-10g%-10g%-10g%-10g\n',t,min(k_t(:,t)),max(k_t(:,t)),kBar_t(t),ginicoeff(k_t(:,t)));
    end
end
simuRslt = v2struct(k_t,kBar_t);
end