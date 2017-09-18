function rn = gen_rn(params)
numOfPeriods = params.numOfPeriods;
numOfAgents = params.numOfAgents;
eBetaEPrimeBetaPrimeCondZZPrime = params.eBetaEPrimeBetaPrimeCondZZPrime ;
zTrans = params.zTrans;
% Generate zShock
zIdx = gen_discrete_markov_rn(zTrans, 1, numOfPeriods, 1);

% Generate eBetaShock0
eBetaIdx0 = ones(numOfAgents,1);

% Generate betaEShock together
eBetaIdx = gen_discrete_markov_rn_with_agg(eBetaEPrimeBetaPrimeCondZZPrime,numOfAgents,numOfPeriods,eBetaIdx0,zIdx);

rn = v2struct(zIdx,eBetaIdx);
end