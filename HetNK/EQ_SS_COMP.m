function EqSsCompMkt = EQ_SS_COMP(Params,EqSs)
v2struct(Params);
% v2struct(EqSs);

% From hours to labor efficiency
EMean = EqSs.N / EqSs.hours;

%% Steady state
% Calibrate eta2 to equalize steady state output
eta2 = EMean^(-1-1/Nu);
Ss0 = [EqSs.C,EqSs.N,EqSs.Y,eta2];

options = optimoptions(@fsolve,'Display','Iter','FunctionTolerance',1e-12,'OptimalityTolerance',1e-12);
Ss = fsolve(@(x)eq_resid_ss(x,Params,EqSs), Ss0, options);

eta2 = Ss(4);

% Correct for complete market
EqSs.R = R0;
EqSs.II = R0;
EqSs.tau = EqSs.B + EqSs.G - EqSs.B/EqSs.R;
EqSs.C = Ss(1);
EqSs.N = Ss(2);
EqSs.Y = Ss(3);
EqSs.Dividend = EqSs.Y - EqSs.wage*EqSs.N;

EqSsCompMkt = EqSs;
EqSsCompMkt.eta2 = eta2;
end

function resid = eq_resid_ss(x,Params,EqSs)
v2struct(Params);

C = x(1);
N = x(2);
Y = x(3);
eta2 = x(4);

resid(1) = C^(-Gamma)*EqSs.wage - eta2*N^(1/Nu);
resid(2) = Y - C - G0;
resid(3) = Y - z0*N;
resid(4) = Y - EqSs.Y;
end