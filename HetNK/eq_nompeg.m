function [F,J] = eq_nompeg(xx,Params,EqSs,BetaShock_t,MShock_t,z_t,G_t)
x = reshape(xx,11,[]);
F = zeros(size(x));

v2struct(Params);
eta2 = EqSs.eta2;

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
    BetaShock = BetaShock_t(t);
    MShock = MShock_t(t);
    z = z_t(t);
    G = G_t(t);
    B = EqSs.B;
    ii0 = EqSs.R-1;

    if t<=PegPeriods
        % Evaluate function
        eq_nom_func;
        F(:,t) = TEMP;
        
        % Evaluate Jacobian
        jac_nom_func;
        Coef(:,:,t) = TEMP;
        jacF_nom_func;
        Coef_F(:,:,t) = TEMP;
        jacB_nom_func;
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
end

J = solve_block_sparse(Coef_B,Coef,Coef_F,F);
F = F(:);
end

