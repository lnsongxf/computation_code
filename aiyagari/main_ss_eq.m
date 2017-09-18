% Aiyagari (1994)
% Gateway to solve steady state
function main_ss_eq
params = setup;
write_params(params);

rng(0823) % Set the seed
rn = gen_rn(params);

rMin = 0.000;
rMax = 0.01;

params.PRINT_FREQ = 1e3;
options = optimset('TolX',params.TOL_EQ);
rEq = fzero(@(r) compute_r_implied(r,rn,params)-r, [rMin,rMax], options);

ssEq = compute_partial_eq_given_r(rEq,rn,params);

save('mat/ss','ssEq','params','rn');
end

function rImplied = compute_r_implied(r,rn,params)

eq = compute_partial_eq_given_r(r,rn,params);
K = eq.aggRslt.K_t(end);
L = eq.aggRslt.L;

delta = params.delta;
alpha = params.alpha;

KLRatio = K/L;

rImplied = alpha* KLRatio^(alpha-1) - delta;

% Print some information
fprintf('r: %g, rImplied: %g\n',r,rImplied);
end