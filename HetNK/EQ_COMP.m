function EqTransCompMkt = EQ_NOMPEG_COMP(Params,EqSs,EqWarmUp)
v2struct(Params);

eta2 = EqSs.eta2;

%% Transition
% Unknowns
C_t = EqSs.C*ones(1,TransPeriods);
N_t = EqSs.N*ones(1,TransPeriods);
Y_t = EqSs.Y*ones(1,TransPeriods);
Dividend_t = EqSs.Dividend*ones(1,TransPeriods);
ii_t = (EqSs.II-1)*ones(1,TransPeriods);
wage_t = EqSs.wage*ones(1,TransPeriods);
R_t = R0 * ones(1,TransPeriods);
tau_t = EqSs.tau*ones(1,TransPeriods);
Pi_t = EqSs.Pi*ones(1,TransPeriods);

if nargin>2
    v2struct(EqWarmUp);
end

x = [C_t;N_t;Y_t;Dividend_t;ii_t;wage_t;Pi_t;];
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'MaxIterations',inf,'MaxFunctionEvaluations',inf,'Display','iter','Algorithm','trust-region-dogleg','FunctionTolerance',1e-8,'StepTolerance',1e-8);
xx = fsolve(@(x) eq_nompeg_rotemberg(x,Params,EqSs,BetaShock_t,MShock_t,z_t,G_t), x(:), options);
x = reshape(xx,7,[]);

iC = 1;
iN = 2;
iY = 3;
iDividend = 4;
iii = 5;
iwage = 6;
iPi = 7;

C_t = x(iC,:);
N_t = x(iN,:);
Y_t = x(iY,:);
Dividend_t = x(iDividend,:);
ii_t = x(iii,:);
wage_t = x(iwage,:);
Pi_t = x(iPi,:);

R_t = (1+ii_t) ./ [Pi_t(2:end) EqSs.Pi];
TaylorI_t = R0 + Phi*(Pi_t-1) - 1;

EqTransCompMkt = v2struct(C_t,N_t,Y_t,Dividend_t,ii_t,wage_t,Pi_t,R_t,tau_t,TaylorI_t);
end

