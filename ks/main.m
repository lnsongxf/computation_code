% Replication of Krusell and Smith (1998)
% Author: Wenlan Luo
% MAIN: Gateway
% clear;

params = setup;
params = common(params);
write_params(params);

rn = gen_rn(params);

% tic;
% vfiRslt = vfi(params);
% toc;

% vfiRslt = load('mat/vfiRslt.mat');
% vfiRslt = vfiRslt.vfiRslt;

metric = 1;
iter = 0;

phi = params.phi;
phi_new = params.phi;
vfiRslt = [];
while metric > params.TOL_EQ
    params.phi = phi;
    fprintf('VFI...\n');
    tic;
    vfiRslt = vfi(params,vfiRslt);
    toc;
    fprintf('Simualte...\n');
    tic;
    simuRslt = simulate(vfiRslt, rn, params);
    toc;
    aggRslt = aggregate(simuRslt,rn,params);
    
    phi_new(1,1) = aggRslt.coefs(2,1);
    phi_new(2,1) = aggRslt.coefs(2,2);
    phi_new(1,2) = aggRslt.coefs(1,1);
    phi_new(2,2) = aggRslt.coefs(1,2);
    
    metric = max(abs(phi_new(:)-phi(:)));
    iter = iter+1;
    fprintf('iter: %d, metric: %g\n',iter,metric);
    phi = phi_new*0.1 + phi*(1-0.1);
    phi
end



