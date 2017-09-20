% Aiyagari (1994)
% Gateway to solve transition path
clear;
params = setup;
write_params(params);
% Load steady state
ss = load('mat/ss.mat');
ssEq = ss.ssEq;

% Destroy everyone's capital to half
kPts = params.kPts;
ePts = params.ePts;
kGrid = params.kGrid;
kHalfGrid = repmat( reshape(kGrid/2,[],kPts), [ePts,1]);
dist0 = simulate_histogram_1(ssEq.dist,kHalfGrid,ssEq.w,ssEq.r,params);

% Inspect distribution before and after
%{
close all;
[eMesh,kMesh] = ndgrid(params.eGrid,params.kGrid);
figure;
surf(eMesh,kMesh,ssEq.dist);
figure;
surf(eMesh,kMesh,dist0);
%}

% trans = compute_trans_eq(ssEq,ssEq.dist,params);
trans = compute_trans_eq(ssEq,dist0,params);

%% Insepct transition equilibrium
figure;
plot(trans.K_t);
title('Capital');

% 
eGrid = params.eGrid;
kGrid = params.kGrid;
[eMesh,kMesh] = ndgrid(eGrid,kGrid);
transT = params.transT;
gini_t = zeros(transT,1);
for t=1:transT
    gini_t(t) = gini_histogram(reshape(trans.dist_t(:,:,t),[],1), kMesh(:), 1e5, 0823);
end
figure;
plot(gini_t);
title('Wealth Gini');
