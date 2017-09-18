% Aiyagari (1994)
% Gateway for solving decision problem
% Wenlan Luo, luowenlan@gmail.com

clear;
addpath('../CRoutines');

r = 0.01;
w = 1;

params = setup;
write_params(params);
params.PRINT_FREQ = 50;
params.MAXITER_VFI = inf;

tic;
[kpPolicy,v] = solve_decision_func_approx_mex(w,r,params);
toc;

figure;
plot(params.kGrid(1:150),v(:,1:150)');
title('value');

figure;
plot(params.kGrid(1:100),kpPolicy(:,1:100)');
title('kpPolicy');