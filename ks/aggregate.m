function aggRslt = aggregate(simuRslt,rn,params)
NUM_THREADS=params.NUM_THREADS; PRINT_FREQ=params.PRINT_FREQ; PRINT_FREQ_SIMULATE=params.PRINT_FREQ_SIMULATE; PROFILE=params.PROFILE; TOL_OPT=params.TOL_OPT; TOL_VFI=params.TOL_VFI; ZETrans=params.ZETrans; alpha=params.alpha; betaGrid=params.betaGrid; betaPts=params.betaPts; betaTrans=params.betaTrans; delta=params.delta; eBetaEPrimeBetaPrimeZZprime=params.eBetaEPrimeBetaPrimeZZprime; eBetaTrans=params.eBetaTrans; eBetaTransConditionalOnZ=params.eBetaTransConditionalOnZ; eBetaTransInv=params.eBetaTransInv; eEPrimeZZPrime=params.eEPrimeZZPrime; eGrid=params.eGrid; ePts=params.ePts; eTransConditionalOnZ=params.eTransConditionalOnZ; fullDiscountedTrans=params.fullDiscountedTrans; fullTrans=params.fullTrans; i_beta=params.i_beta; i_betap=params.i_betap; i_e=params.i_e; i_ep=params.i_ep; i_z=params.i_z; i_zp=params.i_zp; kBarGrid=params.kBarGrid; kBarPts=params.kBarPts; kGrid=params.kGrid; kMax=params.kMax; kMin=params.kMin; kPts=params.kPts; noWorkWage=params.noWorkWage; numOfAgents=params.numOfAgents; numOfPeriods=params.numOfPeriods; numOfPeriodsRemoved=params.numOfPeriodsRemoved; params=params.params; phi=params.phi; zGrid=params.zGrid; zPts=params.zPts; zTrans=params.zTrans;

% Unpack things
k_t = simuRslt.k_t;
kBar_t = simuRslt.kBar_t;

zIdx = rn.zIdx;

% Run regression for conditional on zIdx
% Collect samples
kBar_t(1:numOfPeriodsRemoved) = [];
zIdx(1:numOfPeriodsRemoved) = [];

r2 = zeros(1,zPts);
sse = zeros(1,zPts);
coefs = zeros(2,zPts);

for i_z=1:zPts
    sampleIdxZ = find(zIdx(1:end-1)==i_z);
    kBarZ = kBar_t(sampleIdxZ);
    kBarPrimeZ = kBar_t(sampleIdxZ+1);
    
    mdl = fitlm(log(kBarZ),log(kBarPrimeZ));
    
    coefs(:,i_z) = mdl.Coefficients.Estimate;
    r2(i_z) = mdl.Rsquared.Ordinary;
    sse(i_z) = mdl.SSE;
end

aggRslt = v2struct(coefs,r2,sse);

end