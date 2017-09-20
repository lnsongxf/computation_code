function NewParams = SETUP(Params)
v2struct(Params);

%% Preamble
run '../set_path.m';
addpath('generated_equation');
DEF_CONSTANTS;

%% Parameters
Beta = 0.986; % Discount factor
Gamma = 2; % CRRA
Nu = 0.5; % Frisch
Chi = 1; % Labor weight

BYRatio = 1.4;
Mu = 1.2; % Markup
Eta = Mu/(Mu-1); % Elasticity

Theta = 0.15; % Price revision rate

RhoE = 0.966;
VarE = 0.017;

%% Grid
% Bond
BGrid = textread('points.kgd');
BGrid = BGrid(:,1);
BGrid = BGrid - BGrid(1);
BPts = length(BGrid);
BMin = BGrid(1);
BMax = BGrid(end);

% Expand Bond Grid for simulation
BDistGrid = [];
for i = 0:3
    BDistGrid = [BDistGrid;
        (4-i)*1/4*BGrid(1:end-1)' + i*1/4*BGrid(2:end)'
        ];
end
% BDistGrid = BGrid;
BDistGrid = BDistGrid(:);
BDistGrid = [BDistGrid; BGrid(end)];
BDistPts = length(BDistGrid);

% Discount factor
BetaGrid = 0.986;
BetaTrans = 1;
BetaPts = 1;
if isfield(Params,'HetBeta')
    BetaGrid = [0.1 0.986];
    BetaPts = length(BetaGrid);
    BetaTrans = [
        0.98,0.02;
        0.02,0.98
        ];
end

% Income shock
[ETrans, EGrid] = rouwen(RhoE,0,sqrt(VarE / (1-RhoE^2)),3);
ETrans = ETrans';
EGrid = exp(EGrid);
EPts = length(EGrid);
TauEGrid = [0 0 1]; % Tax who?

% Aggregate
z0 = 1;
R0 = 1 + 0.02/4;
if isfield(Params,'Complete')
    Beta = 1/R0;
elseif isfield(Params,'Incomplete')
else
    error('Specify Comp/Incomp');
end
G0 = 0;
GYRatio = 0;

% Monetary policy
Phi = 1.5; % Taylor rule

%% Transition path
TransPeriods = 250;
G_t = zeros(1,TransPeriods);
z_t = ones(1,TransPeriods);
MShock_t = zeros(1,TransPeriods);
BetaShock_t = ones(1,TransPeriods);
R0_t = R0 * ones(1,TransPeriods);

%% For computation
TolEqSs = 1e-12;
TolEV = 1e-12;
TolSol = 1e-14;
TolEqTrans = 5e-6;
NumThreads = 4;
PrintFreq = 1;
UpdateSpeed = 0.5;
TransUpdateSpeed = 0.5;
EtaUpdateSpeed = 0.1;
FsolveOptions = optimoptions('fsolve','SpecifyObjectiveGradient',true,'MaxIterations',inf,'MaxFunctionEvaluations',inf,'Display','iter','Algorithm','trust-region-dogleg','FunctionTolerance',1e-6,'StepTolerance',1e-6);

clear Params;
NewParams = v2struct;
end