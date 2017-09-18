function rn = gen_rn(params)
eTrans = params.eTrans;
numAgents = params.numAgents;
T = params.T;
ePts = params.ePts;

% Compute invariant distribution of eTrans
eTransStationary = eTrans^1000;
eTransStationary = eTransStationary(1,:);

% Generate uniform random number
un = rand(numAgents,1);

[eShock0Histogram,eShockIdx0] = histc(un,[0,cumsum(eTransStationary)]);

eShockIdx = gen_discrete_markov_rn(eTrans,numAgents,T,eShockIdx0);

eShockHistogram = histc(eShockIdx,1:ePts) ./ numAgents;

rn = v2struct(eShockIdx);
end