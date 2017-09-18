% Aiyagari (1994)
% Gateway to solve partial equilibrium varying interest rate

params = setup;
write_params(params);

rng(0823) % Set the seed
rn = gen_rn(params);

r = 0.005;
eq = compute_partial_eq_given_r(r,rn,params);

figure;
plot(eq.aggRslt.K_t);
xlabel('t');
title('Aggregate Capital over Time');