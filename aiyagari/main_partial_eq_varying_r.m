% Aiyagari (1994)
% Gateway to solve partial equilibrium varying interest rate

params = setup;
write_params(params);

rng(0823) % Set the seed
rn = gen_rn(params);

% Set an interest grid
rGrid = [0.000 0.005 0.007 0.008 0.009 0.0095 0.01];
% Solve the capital supplied by family
KGrid = zeros(1,length(rGrid));
for i_r = 1:length(rGrid)
    r = rGrid(i_r);
    eq = compute_partial_eq_given_r(r,rn,params);
    KGrid(i_r) = eq.aggRslt.K_t(end);
    
    % Labor supply does not vary
    L = eq.aggRslt.L;
end

% Compute a capital demand curve, that is KLRatio that gives interest rate
rGridDemand = [0.00:0.005:0.02];
delta = params.delta;
alpha = params.alpha;
KLRatioImplied = ( (rGridDemand+delta) / alpha ) .^ (1/(alpha-1));
KImplied = KLRatioImplied * L;

% Plot the "supply" and "demand" curve in the figure;
beta = params.beta;
figure;hold on;
plot(KGrid,rGrid,'k-');
plot(KImplied,rGridDemand,'k--');
% A flat line for 1/beta-1
plot([0,KGrid(end)],[1/beta-1 1/beta-1],'k:');
text(1,1/beta-1,'r=1/\beta-1','VerticalAlignment','bottom');
legend('Asset','Capital');
xlim([0 70]);
