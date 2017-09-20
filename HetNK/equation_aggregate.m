function [F,J] = equation_aggregate(x,params,ss,wedge_c_t,wedge_n_t,betaShock_t,mShock_t,z_t,G_t,tax_unit,B)
x = reshape(x,8,[]);
F = zeros(size(x));

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
        C_m1 = ss.C;
        N_m1 = ss.N;
        Y_m1 = ss.Y;
        Dividend_m1 = ss.Dividend;
        ii_m1 = ss.ii - 1;
        wage_m1 = ss.wage;
        Pi_m1 = ss.Pi;
        tau_m1 = ss.tau;
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
        C_1 = ss.C;
        N_1 = ss.N;
        Y_1 = ss.Y;
        Dividend_1 = ss.Dividend;
        ii_1 = ss.ii - 1;
        wage_1 = ss.wage;
        Pi_1 = ss.Pi;
        tau_1 = ss.tau;
    end

    if isfield(params,'peg_r')
        eta1 = wedge_c_t(t);
        eta2 = wedge_n_t(t);
        BetaShock = betaShock_t(t);
        MShock = mShock_t(t);
        z = z_t(t);
        G = G_t(t);
        r0 = R0_t(t);
        
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

J = solve_block_sparse(Coef_B,Coef,Coef_F,F);
F = F(:);
end

