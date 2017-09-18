% Aiyagari (1994)
% Gateway for solving decision problem using grid search
% Wenlan Luo, luowenlan@gmail.com

clear;

r = 0.01;
w = 1;

params = setup;

[kpPolicy,v] = solve_decision_grid_search(w,r,params);

figure;
plot(params.kGrid(1:150),v(:,1:150)');
title('value');

figure;
plot(params.kGrid(1:100),kpPolicy(:,1:100)');
title('kpPolicy');