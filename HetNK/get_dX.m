function [dX,F] = get_dX(x,Params,EqSs,eta1_t,eta2_t,BetaShock_t,MShock_t,z_t,G_t,tax_unit,B)
x = reshape(x,8,[]);
F = zeros(size(x));

v2struct(Params);

iC = 1;
iN = 2;
iY = 3;
iDividend = 4;
iii = 5;
iwage = 6;
iPi = 7;
itau = 8;

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
    tau = x(itau,t);
    
    if t>1
        C_m1 = x(iC,t-1);
        N_m1 = x(iN,t-1);
        Y_m1 = x(iY,t-1);
        Dividend_m1 = x(iDividend,t-1);
        ii_m1 = x(iii,t-1);
        wage_m1 = x(iwage,t-1);
        Pi_m1 = x(iPi,t-1);
        tau_m1 = x(itau,t-1);
    else
        C_m1 = EqSs.C;
        N_m1 = EqSs.N;
        Y_m1 = EqSs.Y;
        Dividend_m1 = EqSs.Dividend;
        ii_m1 = EqSs.ii - 1;
        wage_m1 = EqSs.wage;
        Pi_m1 = EqSs.Pi;
        tau_m1 = EqSs.tau;
    end
    
    if t<TransPeriods
        C_1 = x(iC,t+1);
        N_1 = x(iN,t+1);
        Y_1 = x(iY,t+1);
        Dividend_1 = x(iDividend,t+1);
        ii_1 = x(iii,t+1);
        wage_1 = x(iwage,t+1);
        Pi_1 = x(iPi,t+1);
        tau_1 = x(itau,t+1);
    else
        C_1 = EqSs.C;
        N_1 = EqSs.N;
        Y_1 = EqSs.Y;
        Dividend_1 = EqSs.Dividend;
        ii_1 = EqSs.ii - 1;
        wage_1 = EqSs.wage;
        Pi_1 = EqSs.Pi;
        tau_1 = EqSs.tau;
    end

    if isfield(Params,'peg_r')
        eta1 = eta1_t(t);
        eta2 = eta2_t(t);
        BetaShock = BetaShock_t(t);
        MShock = MShock_t(t);
        z = z_t(t);
        G = G_t(t);
        R0 = R0_t(t);
        
        equation_peg_r;
        F(:,t) = TEMP;
        
        % Evaluate Jacobian
        jac_peg_r;
        Coef(:,:,t) = TEMP;
        jacF_peg_r;
        Coef_F(:,:,t) = TEMP;
        jacB_peg_r;
        Coef_B(:,:,t) = TEMP;
    end
end

dX = solve_block(Coef_B,Coef,Coef_F,F);
end

