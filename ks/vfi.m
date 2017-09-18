function [vfiResult, exitFlag] = vfi(params,warmUp)
NUM_THREADS=params.NUM_THREADS; PRINT_FREQ=params.PRINT_FREQ; PRINT_FREQ_SIMULATE=params.PRINT_FREQ_SIMULATE; PROFILE=params.PROFILE; TOL_EQ=params.TOL_EQ; TOL_OPT=params.TOL_OPT; TOL_VFI=params.TOL_VFI; ZETrans=params.ZETrans; alpha=params.alpha; betaGrid=params.betaGrid; betaPts=params.betaPts; betaTrans=params.betaTrans; delta=params.delta; eBetaEPrimeBetaPrimeCondZZPrime=params.eBetaEPrimeBetaPrimeCondZZPrime; eBetaEPrimeBetaPrimeZZPrime=params.eBetaEPrimeBetaPrimeZZPrime; eGrid=params.eGrid; ePts=params.ePts; fullDiscountedTrans=params.fullDiscountedTrans; fullTrans=params.fullTrans; i_beta=params.i_beta; i_betap=params.i_betap; i_e=params.i_e; i_ep=params.i_ep; i_z=params.i_z; i_zp=params.i_zp; kBarGrid=params.kBarGrid; kBarPts=params.kBarPts; kGrid=params.kGrid; kMax=params.kMax; kMin=params.kMin; kPts=params.kPts; noWorkWage=params.noWorkWage; numOfAgents=params.numOfAgents; numOfPeriods=params.numOfPeriods; numOfPeriodsRemoved=params.numOfPeriodsRemoved; phi=params.phi; zGrid=params.zGrid; zPts=params.zPts; zTrans=params.zTrans; 

value = zeros(zPts, ePts, betaPts, kPts, kBarPts);
kpPolicy = zeros(zPts, ePts, betaPts, kPts, kBarPts);

value_new = zeros(size(value));
metric = 1;
iter = 0;

if nargin>1 && ~isempty(warmUp)
    v2struct(warmUp);
end

[eMesh,betaMesh,kMesh] = ndgrid(eGrid,betaGrid,kGrid);

while (metric > TOL_VFI)
    for i_z = 1:zPts
        % do a construction w.r.t. kBar
        valueKBarSplineCoefs = myppualMKL_CMEX(int32(-8), {kBarGrid}, reshape(value(i_z,:), [ePts*betaPts*kPts,kBarPts])', [], int32(4), int32(ePts*betaPts*kPts), [], [], [], []);
        for i_kBar = 1:kBarPts
            kBar = kBarGrid(i_kBar);
            % predict kBarPrime using guess of transition rules
            kBarPrime = exp(phi(i_z,1)*log(kBar) + phi(i_z,2));
            lBar = exp(phi(i_z,3)*log(kBar) + phi(i_z,4));
            % prices
            w = zGrid(i_z) * (1-alpha) * (kBar/lBar)^alpha;
            r = zGrid(i_z) * alpha * (kBar/lBar)^(alpha-1) - delta;
            
            % interpoate at kBarPrime
            valueAtKBarPrime = myppualMKL_CMEX(int32(8), {kBarGrid}, valueKBarSplineCoefs, [], int32(4), int32(ePts*betaPts*kPts), [], kBarPrime, [], int32([1:ePts*betaPts*kPts]')-1, []);
            
            % treat k as continuous variable, and beta and e as vector
            % functions
            valueKSplineCoefs = myppualMKL_CMEX(int32(-8), {kGrid}, reshape(valueAtKBarPrime, [], kPts)', [], int32(4), int32(ePts*betaPts), [], [], [], []);
            
            budgetMesh = eMesh*w + kMesh*(1+r);
            switch PROFILE
                case 'Beta'
                    budgetMesh = budgetMesh + (eMesh==0)*noWorkWage;
            end
            
            % Allocate space for value_sub, kpPolicy_sub given i_z and kbar
            value_sub = zeros(ePts,betaPts,kPts);
            kpPolicy_sub = zeros(ePts,betaPts,kPts);
            
            vfi_mex;
            
            value_new(i_z, :, :, :, i_kBar) = reshape(value_sub, [1 ePts betaPts kPts 1]);
            kpPolicy(i_z, :, :, :, i_kBar) = reshape(kpPolicy_sub, [1 ePts betaPts kPts 1]);
        end
    end
    
    % integrate over z e beta
    expectedValue_new = fullDiscountedTrans * reshape(value_new, zPts*ePts*betaPts, []);
    
    % update metric
    metric = max(abs(value(:) - expectedValue_new(:)));
    value = reshape(expectedValue_new, [zPts ePts betaPts kPts kBarPts]);
    
    iter = iter + 1;
    
    if mod(iter,PRINT_FREQ) == 0
        fprintf('iter: %d, metric: %g\n', iter, metric);
    end
end

fprintf('iter: %d, metric: %g\n', iter, metric);

vfiResult = v2struct(kpPolicy,value);
exitFlag = 0;
end
