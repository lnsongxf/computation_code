% Aiyagari (1994)
% Gateway for solving decision problem
% Wenlan Luo, luowenlan@gmail.com

clear;

r = 0.01;
w = 1;

params = setup;
write_params(params);
params.PRINT_FREQ = 1;
params.MAXITER_VFI = 10;

[kpPolicy,v] = solve_decision_func_approx_vec(w,r,params);

figure;
plot(params.kGrid(1:150),v(:,1:150)');
title('value');

figure;
plot(params.kGrid(1:100),kpPolicy(:,1:100)');
title('kpPolicy');