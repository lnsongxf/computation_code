function params = setup
% Replication of Krusell and Smith (1998)
% Author: Wenlan Luo
% SETUP: Setup parameters


PROFILE = 'Benchmark';

addpath('../utility');

%
alpha = 0.36; % capital share in production
delta = 0.025; % depreciation rate

% aggregate shocks
zGrid = [0.99 1.01];
zTrans = [
    0.875 0.125
    0.125 0.875
    ];
zPts = length(zGrid);

% individual shocks
noWorkWage = 0.07; % unemployment insurance in stochastic beta model
eGrid = [0 0.3271];
ePts = length(eGrid);

% joint transition matrix
ZETrans = [
    0.5250 0.3500 0.0312 0.0938
    0.0389 0.8361 0.0021 0.1229
    0.0938 0.0312 0.2917 0.5833
    0.0091 0.1159 0.0243 0.8507
    ];

% conditional transition matrix

% beta states
switch PROFILE
    case 'Benchmark'
        betaGrid = 0.99;
        betaTrans = 1;
    case 'Beta'
        betaGrid = [0.9858 0.9894 0.9930];
        betaTrans = [
            0.9950 0.0050 0.0000
            0.0006 0.9988 0.0006
            0.0000 0.0050 0.9950
            ];
end
betaPts = length(betaGrid);

% grid for asset holdings
kGrid = csvread('points.kgd');
switch PROFILE
    case 'Benchmark'
        kGrid = kGrid + 2.4 + 1e-3;
end
kPts = length(kGrid);
kMin = min(kGrid);
kMax = max(kGrid);

% pseudo states for aggregate capital
kBarGrid = [11.1 11.52 11.94 12.36];
kBarPts = length(kBarGrid);

% agents
numOfAgents = 1e4;
numOfPeriods = 11000;
numOfPeriodsRemoved = 1000;

% initial guess for transition rule, Phi has the order [z order]
switch PROFILE
    case 'Benchmark'
        phi = [
            0.96479 0.0837134 0 -1.22282
            0.962935 0.0936052 0 -1.15836
            ];
    case 'Beta'
        phi = [
            0.960769 0.0944097 0 -1.22282
            0.960787 0.0990515 0 -1.15836
            ];
end

%
TOL_VFI = 1e-7;
TOL_OPT = 1e-10;
TOL_EQ = 1e-6;
PRINT_FREQ = 100;
PRINT_FREQ_SIMULATE = 1000;
NUM_THREADS = feature('numcores')*2;

params = v2struct();
end