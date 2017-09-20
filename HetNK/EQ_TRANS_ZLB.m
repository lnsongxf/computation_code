function [EqTrans,EqTransAgg] = EQ_TRANS_ZLB(EqSs,Params)
v2struct(Params);

Metric = 1;
Iter = 0;

% Specify the sequence of unknowns
C_t = EqSs.C*ones(1,TransPeriods);
N_t = EqSs.N*ones(1,TransPeriods);
Y_t = EqSs.Y*ones(1,TransPeriods);
Dividend_t = EqSs.Dividend*ones(1,TransPeriods);
ii_t = (EqSs.II-1) * ones(1,TransPeriods);
ii_t(1:CommitPeriods) = 0;
% ii_t(:) = 0;
wage_t = EqSs.wage*ones(1,TransPeriods);
R_t = R0 * ones(1,TransPeriods);
% r_t = r0_t;
pRatio_t = EqSs.pRatio*ones(1,TransPeriods);
S_t = EqSs.S*ones(1,TransPeriods);
tau_t = EqSs.tau*ones(1,TransPeriods);
PA_t = EqSs.PA*ones(1,TransPeriods);
PB_t = EqSs.PB*ones(1,TransPeriods);
Pi_t = EqSs.Pi*ones(1,TransPeriods);

eta1_t = 1 ./ (Beta*R_t.* ([C_t(2:end),EqSs.C] ./ C_t).^(-Gamma));
eta2_t = C_t(1:TransPeriods).^(-Gamma).*wage_t ./ (N_t.^(1/Nu));

% Sequence of decision and distribution
EVPrime_t = repmat(EqSs.EVPrime, [1 1 TransPeriods+1]);
bPolicy_t = repmat(EqSs.bPolicy, [1 1 TransPeriods]);
nPolicy_t = repmat(EqSs.nPolicy, [1 1 TransPeriods]);
cPolicy_t = repmat(EqSs.cPolicy, [1 1 TransPeriods]);
Dist_t = repmat(EqSs.Dist, [1 1 TransPeriods+1]);

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
nVars = 11;
while (Metric>TolEqTrans)
    tic;
    % Given relative prices, solve decision and distribution
    for t=TransPeriods:-1:1
        Params.Beta = Beta*Params.BetaShock_t(t);
        [EVPrime_t(dimEV{:},t),bPolicy(dimEV{:},t),nPolicy(dimEV{:},t),cPolicy(dimEV{:},t)] = VFI(EVPrime_t(dimEV{:},t+1),Params,tau_t(t),Dividend_t(t),wage_t(t),R_t(t));
    end
    
    for t=1:1:TransPeriods
        Params.Beta = Beta*Params.BetaShock_t(t);
        [bPolicyDist,nPolicyDist,cPolicyDist] = EXPAND_TO_DIST(bPolicy(dimEV{:},t),nPolicy(dimEV{:},t),wage_t(t),R_t(t),tau_t(t),Dividend_t(t),Params);
        Dist_t(dimDist{:},t+1) = UPDATE_DIST(Dist_t(dimDist{:},t),bPolicyDist,Params);
        [N_t(t),~,~,~,C_t(t),~,~] = AGGREGATE_SS(bPolicyDist,nPolicyDist,cPolicyDist,Dist_t(dimDist{:},t),Params,z_t(t),G_t(t),R_t(t),wage_t(t));
    end
    
    % Calibrate a distorting parameter
    eta1_t_new = 1 ./ (Beta*Params.BetaShock_t.*R_t.* ([C_t(2:end),EqSs.C] ./ C_t).^(-Gamma));
    eta2_t_new = C_t(1:TransPeriods).^(-Gamma).*wage_t ./ (N_t.^(1/Nu));

    % Compute the block matrix
    Coef_B = zeros(nVars,nVars,TransPeriods);
    Coef = zeros(nVars,nVars,TransPeriods);
    Coef_F = zeros(nVars,nVars,TransPeriods);
    F = zeros(nVars,TransPeriods);
    
    x = [C_t;N_t;Y_t;Dividend_t;ii_t;wage_t;Pi_t;pRatio_t;S_t;PA_t;PB_t];
    
    MinorMetric = 1;
    while (MinorMetric>1e-8)
        for t=1:TransPeriods
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
            
            eta1 = eta1_t(t);
            eta2 = eta2_t(t);
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
    
    C_t_new = x(iC,:);
    N_t_new = x(iN,:);
    Y_t_new = x(iY,:);
    Dividend_t_new = x(iDividend,:);
    ii_t_new = x(iii,:);
    wage_t_new = x(iwage,:);
    Pi_t_new = x(iPi,:);
    pRatio_t_new = x(ipRatio,:);
    S_t_new = x(iS,:);
    PA_t_new = x(iPA,:);
    PB_t_new = x(iPB,:);
    
    R_t_new = (1+ii_t_new) ./ ([Pi_t_new(2:end) EqSs.Pi]);
    tau_t_new = (EqSs.B + G_t - EqSs.B./R_t_new) / EqSs.DistE(3);
    
    % Do a line search
    Metric = max(abs([
        tau_t(:) - tau_t_new(:);
        Dividend_t(:) - Dividend_t_new(:);
        wage_t(:) - wage_t_new(:);
        R_t(:) - R_t_new(:);
        eta1_t(:) - eta1_t_new(:);
        eta2_t(:) - eta2_t_new(:);
        ]));
    
    C_t = TransUpdateSpeed*C_t_new + (1-TransUpdateSpeed)*C_t;
    N_t = TransUpdateSpeed*N_t_new + (1-TransUpdateSpeed)*N_t;
    Y_t = TransUpdateSpeed*Y_t_new + (1-TransUpdateSpeed)*Y_t;
    Dividend_t = TransUpdateSpeed*Dividend_t_new + (1-TransUpdateSpeed)*Dividend_t;
    ii_t = TransUpdateSpeed*ii_t_new + (1-TransUpdateSpeed)*ii_t;
    tau_t = TransUpdateSpeed*tau_t_new + (1-TransUpdateSpeed)*tau_t;
    Pi_t = TransUpdateSpeed*Pi_t_new + (1-TransUpdateSpeed)*Pi_t;
    pRatio_t = TransUpdateSpeed*pRatio_t_new + (1-TransUpdateSpeed)*pRatio_t;
    S_t = TransUpdateSpeed*S_t_new + (1-TransUpdateSpeed)*S_t;
    PA_t = TransUpdateSpeed*PA_t_new + (1-TransUpdateSpeed)*PA_t;
    PB_t = TransUpdateSpeed*PB_t_new + (1-TransUpdateSpeed)*PB_t;
    Dividend_t = TransUpdateSpeed*Dividend_t_new + (1-TransUpdateSpeed)*Dividend_t;
    wage_t = TransUpdateSpeed*wage_t_new + (1-TransUpdateSpeed)*wage_t;
    R_t = TransUpdateSpeed*R_t_new + (1-TransUpdateSpeed)*R_t;
    
    if Metric<1e-2
        EtaUpdateSpeed = 1;
    end
    eta1_t = EtaUpdateSpeed*eta1_t_new + eta1_t.*(1-EtaUpdateSpeed);
    eta2_t = EtaUpdateSpeed*eta2_t_new + eta2_t.*(1-EtaUpdateSpeed);
    
    Iter = Iter+1;
    
    fprintf('Iter: %d. Metric: %g\n',Iter,Metric);
    toc;
end

EqTrans = v2struct(EVPrime_t,Dist_t,bPolicy_t,nPolicy_t,cPolicy_t,C_t,N_t,Y_t,Dividend_t,ii_t,wage_t,Pi_t,R_t,pRatio_t,S_t,tau_t,PA_t,PB_t);
EqTransAgg = v2struct(C_t,N_t,Y_t,Dividend_t,ii_t,wage_t,Pi_t,R_t,pRatio_t,S_t,tau_t,PA_t,PB_t);
end