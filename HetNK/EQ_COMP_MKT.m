function [EqSsCompMkt, EqTransCompMkt] = EQ_COMP_MKT(Params,EqSs)
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
EqSs.PA = Mu*EqSs.wage/z0*EqSs.Y / (1 - (1-Theta)*Beta);
EqSs.PB = EqSs.Y / (1 - (1-Theta)*Beta);
EqSs.Dividend = EqSs.Y - EqSs.wage*EqSs.N;


%% Transition
% Unknowns
C_t = EqSs.C*ones(1,TransPeriods);
N_t = EqSs.N*ones(1,TransPeriods);
Y_t = EqSs.Y*ones(1,TransPeriods);
Dividend_t = EqSs.Dividend*ones(1,TransPeriods);
II_t = EqSs.II * ones(1,TransPeriods);
wage_t = EqSs.wage*ones(1,TransPeriods);
% r_t = r0 * ones(1,TransPeriods);
R_t = R0_t;
pRatio_t = EqSs.pRatio*ones(1,TransPeriods);
S_t = EqSs.S*ones(1,TransPeriods);
tau_t = EqSs.tau*ones(1,TransPeriods);
PA_t = EqSs.PA*ones(1,TransPeriods);
PB_t = EqSs.PB*ones(1,TransPeriods);
Pi_t = EqSs.Pi*ones(1,TransPeriods);

x = [C_t;N_t;Y_t;Dividend_t;II_t;wage_t;Pi_t;R_t;pRatio_t;S_t;tau_t;PA_t;PB_t];
MinorMetric = 1;
while (MinorMetric>1e-8)
    for t=1:TransPeriods
        % Unpack things at current period
        % var =[C,N,Y,Dividend,II,wage,Pi,r,pRatio,S,tau,PA,PB];
        C = x(1,t);
        N = x(2,t);
        Y = x(3,t);
        Dividend = x(4,t);
        II = x(5,t);
        wage = x(6,t);
        Pi = x(7,t);
        R = x(8,t);
        pRatio = x(9,t);
        S = x(10,t);
        tau = x(11,t);
        PA = x(12,t);
        PB = x(13,t);
        
        if t>1
            C_m1 = x(1,t-1);
            N_m1 = x(2,t-1);
            Y_m1 = x(3,t-1);
            Dividend_m1 = x(4,t-1);
            II_m1 = x(5,t-1);
            wage_m1 = x(6,t-1);
            Pi_m1 = x(7,t-1);
            R_m1 = x(8,t-1);
            pRatio_m1 = x(9,t-1);
            S_m1 = x(10,t-1);
            tau_m1 = x(11,t-1);
            PA_m1 = x(12,t-1);
            PB_m1 = x(13,t-1);
        else
            C_m1 = EqSs.C;
            N_m1 = EqSs.N;
            Y_m1 = EqSs.Y;
            Dividend_m1 = EqSs.Dividend;
            II_m1 = EqSs.II;
            wage_m1 = EqSs.wage;
            Pi_m1 = EqSs.Pi;
            R_m1 = EqSs.R;
            pRatio_m1 = EqSs.pRatio;
            S_m1 = EqSs.S;
            tau_m1 = EqSs.tau;
            PA_m1 = EqSs.PA;
            PB_m1 = EqSs.PB;
        end
        
        if t<TransPeriods
            C_1 = x(1,t+1);
            N_1 = x(2,t+1);
            Y_1 = x(3,t+1);
            Dividend_1 = x(4,t+1);
            II_1 = x(5,t+1);
            wage_1 = x(6,t+1);
            Pi_1 = x(7,t+1);
            R_1 = x(8,t+1);
            pRatio_1 = x(9,t+1);
            S_1 = x(10,t+1);
            tau_1 = x(11,t+1);
            PA_1 = x(12,t+1);
            PB_1 = x(13,t+1);
        else
            C_1 = EqSs.C;
            N_1 = EqSs.N;
            Y_1 = EqSs.Y;
            Dividend_1 = EqSs.Dividend;
            II_1 = EqSs.II;
            wage_1 = EqSs.wage;
            Pi_1 = EqSs.Pi;
            R_1 = EqSs.R;
            pRatio_1 = EqSs.pRatio;
            S_1 = EqSs.S;
            tau_1 = EqSs.tau;
            PA_1 = EqSs.PA;
            PB_1 = EqSs.PB;
        end
        
        eta1 = 1;
        % eta2 = eta2;
        
        Beta_C = Beta;
        z = z_t(t);
        G = G_t(t);
        B = EqSs.B;
        RFix = R0_t(t);
        DistE3 = 1;
        
        % Evaluate function
        eq_func;
        F(:,t) = TEMP;
        
        % Evaluate Jacobian
        jac_func;
        Coef(:,:,t) = TEMP;
        jacF_func;
        Coef_F(:,:,t) = TEMP;
        jacB_func;
        Coef_B(:,:,t) = TEMP;
    end
    % Compute gradient
    dX = solve_block(Coef_B,Coef,Coef_F,F);
    x = x+dX;
    
    MinorMetric = max(abs(F(:)))
end
ptr = 0;

%1 C
C_t = x(ptr+1,:);
ptr = ptr+1;

%2 N
N_t = x(ptr+1,:);
ptr = ptr+1;

%3 Y
Y_t = x(ptr+1,:);
ptr = ptr+1;

%4 D
Dividend_t = x(ptr+1,:);
ptr = ptr+1;

%5 i
II_t = x(ptr+1,:);
ptr = ptr+1;

%6 wage
wage_t = x(ptr+1,:);
ptr = ptr+1;

%7 Pi
Pi_t = x(ptr+1,:);
ptr = ptr+1;

%8 r
R_t = x(ptr+1,:);
ptr = ptr+1;

%9 pRatio
pRatio_t = x(ptr+1,:);
ptr = ptr+1;

%10 S
S_t = x(ptr+1,:);
ptr = ptr+1;

%11 tau
tau_t = x(ptr+1,:);
ptr = ptr+1;

%12 PA
PA_t = x(ptr+1,:);
ptr = ptr+1;

%13 PB
PB_t = x(ptr+1,:);
ptr = ptr+1;

EqSsCompMkt = EqSs;
EqSsCompMkt.eta2 = eta2;
EqTransCompMkt = v2struct(C_t,N_t,Y_t,Dividend_t,II_t,wage_t,Pi_t,R_t,pRatio_t,S_t,tau_t,PA_t,PB_t);
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
% resid(4) = eta2 - 1;
end