function EqSs = EQ_SS(Params,EqWarmUp)
v2struct(Params);

% Initiate Value
wage = z0 / Mu;
G = G0;
R = R0;
Dividend = 0;
tau = 0;
EVPrime = 1./(wage*EMesh + BMesh);
Dist = zeros(sizeDist);
% Dist(:,1) = 1/EPts;
Dist(:) = 1/numel(Dist);
Metric = 1;
Iter = 0;

if nargin>1
    v2struct(EqWarmUp);
end

EVMetric = 1;
DistMetric = 1;
AggMetric = 1;

while (Metric > TolEqSs)
    EVMetric = 1;
    while (EVMetric > TolEV)
        [EVPrime_new,bPolicy,nPolicy,cPolicy] = VFI(EVPrime,Params,tau,Dividend,wage,R);
        % Compute Metric
        EVMetric = max(abs(EVPrime_new(:) - EVPrime(:)));
        % Update
        EVPrime = EVPrime_new;
    end
    
    DistMetric = 1;
    [bPolicyDist,nPolicyDist,cPolicyDist] = EXPAND_TO_DIST(bPolicy,nPolicy,wage,R,tau,Dividend,Params);
    while (DistMetric > TolEqSs)
        Dist_new = UPDATE_DIST(Dist,bPolicyDist,Params);
        DistMetric = sum(abs(Dist_new(:) - Dist(:)));
        Dist = Dist_new;
    end
    
    % Aggregate
    [N,B,B2,tau_new,C,Y,Dividend_new,DistE,hours] = AGGREGATE_SS(bPolicyDist,nPolicyDist,cPolicyDist,Dist,Params,z0,G,R,wage);
    AggMetric = max([
        abs(Dividend_new - Dividend)
        abs(tau_new - tau)
        ]);
    
    tau = tau_new*UpdateSpeed + tau*(1-UpdateSpeed);
    Dividend = Dividend_new*UpdateSpeed + Dividend*(1-UpdateSpeed);
    
    Metric = max([EVMetric,DistMetric,AggMetric]);
    Iter = Iter+1;
    
    % Print something
    if mod(Iter,PrintFreq) == 0
        fprintf('%16s%16s%16s%16s%16s\n','Iter','Metric','EVMetric','DistMetric','AggMetric');
        fprintf('%16d%16g%16g%16g%16g\n',Iter,Metric,EVMetric,DistMetric,AggMetric);
    end
end

L = N;
Pi = 1;
ii = R-1;

EqSs = v2struct(EVPrime,bPolicy,nPolicy,cPolicy,Dist,tau,Dividend,N,B,C,Y,wage,Pi,L,R,DistE,ii,hours,G);

end