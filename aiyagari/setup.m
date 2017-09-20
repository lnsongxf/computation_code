function params = setup
% Aiyagari (1994)
% Setting model parameters
% Wenlan Luo, luowenlan@gmail.com

addpath('./solve_decision');
addpath('../utility');

% Preferences
beta = 0.98;
gamma = 1; % log utility

% Labor shocks
eRho = 0.9;
eSigma = 0.1; % The unconditional standard deviation;
ePts = 7;
eRange = 3;
[eTrans,eGrid] = markovappr(eRho,eSigma,eRange,ePts);
eGrid = exp(eGrid);

% Production
alpha = 0.36;
delta = 0.025;

% Grid
kMin = 0.0;
kMax = 5e3;
kPts = 200;
kShift = 1e-2;
kGrid = exp(linspace( log(kMin+kShift),log(kMax+kShift),kPts ))-kShift;
kGrid(1) = kMin; %Rouding error

% For simulations
numAgents = 2e3;
T = 1.0e3;

% Transition path
transT = 200;

% Computation related
TOL_OPT = 1e-12;
TOL_VFI = 1e-12;
TOL_DIST = 1e-12;
MAXITER_VFI = inf;
TOL_EQ = 1e-12;
PRINT_FREQ = 50;
NUM_THREADS = feature('numcores')*2;

% Pack parameters
params = v2struct();
end