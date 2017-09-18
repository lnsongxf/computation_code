function eq = compute_partial_eq_given_r(r,rn,params)
v2struct(params);
alpha = params.alpha;
delta = params.delta;

% Compute the implied KL ratio given r
KLRatio = ( (r+delta) / alpha ) ^ (1/(alpha-1));
% Compute implied wage
w = (1-alpha) * KLRatio^alpha;

% Compute policie functions
[kpPolicy,vp] = solve_decision_egm(w,r,params);
% Simulate
dist = simulate(kpPolicy,rn,w,r,params);
% Aggregate
aggRslt = aggregate(dist,rn,w,r,params);

eq.kpPolicy = kpPolicy;
eq.vp = vp;
eq.aggRslt = aggRslt;
eq.dist = dist;
eq.w = w;
eq.KLRatio = KLRatio;
eq.r = r;
end