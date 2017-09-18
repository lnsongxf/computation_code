% Aiyagari (1994)
% Gateway to inspect steady state properties

% Load eq
ss = load('mat/ss.mat');
eq = ss.ssEq;
params = ss.params;

% %{
% Policy
figure;
plot(params.kGrid(1:120),eq.kpPolicy(:,1:120)');
title('kpPolicy');

% Policy zoom in
figure;
plot(params.kGrid(1:70),eq.kpPolicy([1 7],1:70)');
title('kpPolicy zoom in');
legend('\epsilon_{low}','\epsilon_{high}');

% Gini
ginicoeff(eq.dist.inc_t(:,end))

ginicoeff(eq.dist.k_t(:,end))

% Histogram
figure;
hist3([ss.rn.eShockIdx(:,end),eq.dist.k_t(:,end)],'Nbins',[7 50]);
xlabel('\epsilon'); ylabel('k');
%}

% Non-stochastic simulation
params.PRINT_FREQ = 1000;
dist = simulate_histogram(kpPolicy,eq.w,eq.r,params);
% Plot histogram
[eMesh,kMesh] = ndgrid(params.eGrid,params.kGrid);
surf(eMesh,kMesh,dist);