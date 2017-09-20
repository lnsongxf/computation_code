function EqTransCompMkt = EQ_COMP_MKT_ZLB(Params,EqSs)
v2struct(Params);

eta2 = EqSs.eta2;

%% Transition
% Unknowns
C_t = EqSs.C*ones(1,TransPeriods);
N_t = EqSs.N*ones(1,TransPeriods);
Y_t = EqSs.Y*ones(1,TransPeriods);
Dividend_t = EqSs.Dividend*ones(1,TransPeriods);
ii_t = (EqSs.II-1)*ones(1,TransPeriods);
ii_t(1:CommitPeriods) = 0;
% ii_t = zeros(1,TransPeriods);
wage_t = EqSs.wage*ones(1,TransPeriods);
R_t = R0 * ones(1,TransPeriods);
% r_t = r0_t;
pRatio_t = EqSs.pRatio*ones(1,TransPeriods);
S_t = EqSs.S*ones(1,TransPeriods);
tau_t = EqSs.tau*ones(1,TransPeriods);
PA_t = EqSs.PA*ones(1,TransPeriods);
PB_t = EqSs.PB*ones(1,TransPeriods);
Pi_t = EqSs.Pi*ones(1,TransPeriods);

x = [C_t;N_t;Y_t;Dividend_t;ii_t;wage_t;Pi_t;pRatio_t;S_t;PA_t;PB_t];
iC = 1;
iN = 2;
iY = 3;
iDividend = 4;
iii = 5;
iwage = 6;
iPi = 7;
ipRatio = 8;
iS = 9;
iPA = 10;
iPB = 11;
% xL = repmat([1e-3;1e-3;1e-3;1e-3;1e-3;1e-3;-inf;-inf;1e-3;1e-3;1e-3;1e-3;1e-3],[1 TransPeriods]);
MinorMetric = 1;
while (MinorMetric>1e-8)
    for t=1:TransPeriods
        % Unpack things at current period
        % var =[C,N,Y,Dividend,ii,wage,Pi,r,pRatio,S,tau,PA,PB];
        C = x(iC,t);
        N = x(iN,t);
        Y = x(iY,t);
        Dividend = x(iDividend,t);
        ii = x(iii,t);
        wage = x(iwage,t);
        Pi = x(iPi,t);
        pRatio = x(ipRatio,t);
        S = x(iS,t);
        PA = x(iPA,t);
        PB = x(iPB,t);
        
        if t>1
            C_m1 = x(iC,t-1);
            N_m1 = x(iN,t-1);
            Y_m1 = x(iY,t-1);
            Dividend_m1 = x(iDividend,t-1);
            ii_m1 = x(iii,t-1);
            wage_m1 = x(iwage,t-1);
            Pi_m1 = x(iPi,t-1);
            pRatio_m1 = x(ipRatio,t-1);
            S_m1 = x(iS,t-1);
            PA_m1 = x(iPA,t-1);
            PB_m1 = x(iPB,t-1);
        else
            C_m1 = EqSs.C;
            N_m1 = EqSs.N;
            Y_m1 = EqSs.Y;
            Dividend_m1 = EqSs.Dividend;
            ii_m1 = EqSs.II - 1;
            wage_m1 = EqSs.wage;
            Pi_m1 = EqSs.Pi;
            pRatio_m1 = EqSs.pRatio;
            S_m1 = EqSs.S;
            PA_m1 = EqSs.PA;
            PB_m1 = EqSs.PB;
        end
        
        if t<TransPeriods
            C_1 = x(iC,t+1);
            N_1 = x(iN,t+1);
            Y_1 = x(iY,t+1);
            Dividend_1 = x(iDividend,t+1);
            ii_1 = x(iii,t+1);
            wage_1 = x(iwage,t+1);
            Pi_1 = x(iPi,t+1);
            pRatio_1 = x(ipRatio,t+1);
            S_1 = x(iS,t+1);
            PA_1 = x(iPA,t+1);
            PB_1 = x(iPB,t+1);
        else
            C_1 = EqSs.C;
            N_1 = EqSs.N;
            Y_1 = EqSs.Y;
            Dividend_1 = EqSs.Dividend;
            ii_1 = EqSs.II - 1;
            wage_1 = EqSs.wage;
            Pi_1 = EqSs.Pi;
            pRatio_1 = EqSs.pRatio;
            S_1 = EqSs.S;
            PA_1 = EqSs.PA;
            PB_1 = EqSs.PB;
        end
        
        eta1 = 1;
        % eta2 = eta2;
        
        % Beta = Beta_t(t);
        BetaShock = BetaShock_t(t);
        MShock = MShock_t(t);
        z = z_t(t);
        G = G_t(t);
        B = EqSs.B;
        ii0 = 0;
        
        % Evaluate function
        eq_taylor_func;
        F(:,t) = TEMP;
        
        % Evaluate Jacobian
        jac_taylor_func;
        Coef(:,:,t) = TEMP;
        jacF_taylor_func;
        Coef_F(:,:,t) = TEMP;
        jacB_taylor_func;
        Coef_B(:,:,t) = TEMP;
        
        %{
        if t<=CommitPeriods
            % Evaluate function
            eq_zlb_func;
            F(:,t) = TEMP;
            
            % Evaluate Jacobian
            jac_zlb_func;
            Coef(:,:,t) = TEMP;
            jacF_zlb_func;
            Coef_F(:,:,t) = TEMP;
            jacB_zlb_func;
            Coef_B(:,:,t) = TEMP;
        else
            % Evaluate function
            eq_taylor_func;
            F(:,t) = TEMP;
            
            % Evaluate Jacobian
            jac_taylor_func;
            Coef(:,:,t) = TEMP;
            jacF_taylor_func;
            Coef_F(:,:,t) = TEMP;
            jacB_taylor_func;
            Coef_B(:,:,t) = TEMP;
        end
        %}
    end
    % Compute gradient
    dX = solve_block(Coef_B,Coef,Coef_F,F);
    
    step = 0.25;
    if (MinorMetric < 0.001)
        step = 1;
    end
    x = x+step*dX;
    
    MinorMetric = max(abs(F(:)));
end
ptr = 0;


C_t = x(iC,:);
N_t = x(iN,:);
Y_t = x(iY,:);
Dividend_t = x(iDividend,:);
ii_t = x(iii,:);
wage_t = x(iwage,:);
Pi_t = x(iPi,:);
pRatio_t = x(ipRatio,:);
S_t = x(iS,:);
PA_t = x(iPA,:);
PB_t = x(iPB,:);

R_t = (1+ii_t) ./ [Pi_t(2:end) EqSs.Pi];

TaylorI_t = R0 + Phi*(Pi_t-1) - 1;

EqTransCompMkt = v2struct(C_t,N_t,Y_t,Dividend_t,ii_t,wage_t,Pi_t,R_t,pRatio_t,S_t,tau_t,PA_t,PB_t,TaylorI_t);
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