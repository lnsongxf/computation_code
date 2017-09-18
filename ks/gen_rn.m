function rn = gen_rn(params)
numOfPeriods = params.numOfPeriods;
numOfAgents = params.numOfAgents;
eBetaTransConditionalOnZ = params.eBetaTransConditionalOnZ ;
eBetaTransInv = params.eBetaTransInv;
zTrans = params.zTrans;
% Generate zShock
zIdx = gen_discrete_markov_rn(zTrans, 1, numOfPeriods, 1);

% Generate eBetaShock0
un = rand(numOfAgents,1);
[~,eBetaIdx0] = histc(un,[0 cumsum(eBetaTransInv)]);

% Generate betaEShock together
eBetaIdx = gen_discrete_markov_rn_with_agg(eBetaTransConditionalOnZ,numOfAgents,numOfPeriods,eBetaIdx0,zIdx);

rn = v2struct(zIdx,eBetaIdx);
end