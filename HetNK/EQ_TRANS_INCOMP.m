function [EqTrans,EqTransAgg] = EQ_TRANS_INCOMP(EqSs,Params)
v2struct(Params);

Metric = 1;
Iter = 0;

% Specify the sequence of unknowns
C_t = EqSs.C*ones(1,TransPeriods);
N_t = EqSs.N*ones(1,TransPeriods);
Y_t = EqSs.Y*ones(1,TransPeriods);
Dividend_t = EqSs.Dividend*ones(1,TransPeriods);
ii_t = EqSs.ii * ones(1,TransPeriods);
wage_t = EqSs.wage*ones(1,TransPeriods);
R_t = EqSs.R*ones(1,TransPeriods);
% R_t = 1.004*ones(1,TransPeriods);
R_t(19:21) = 1.00;
tau_t = EqSs.tau*ones(1,TransPeriods);
Pi_t = EqSs.Pi*ones(1,TransPeriods);

% Sequence of decision and distribution
EVPrime_t = repmat(EqSs.EVPrime, [1 1 TransPeriods+1]);
bPolicy_t = repmat(EqSs.bPolicy, [1 1 TransPeriods]);
nPolicy_t = repmat(EqSs.nPolicy, [1 1 TransPeriods]);
cPolicy_t = repmat(EqSs.cPolicy, [1 1 TransPeriods]);
Dist_t = repmat(EqSs.Dist, [1 1 TransPeriods+1]);

% Initial wedge
eta1_t = 1 ./ (Beta*BetaShock_t.*R_t.* ([C_t(2:end),EqSs.C] ./ C_t).^(-Gamma));
eta2_t = C_t(1:TransPeriods).^(-Gamma).*wage_t ./ (N_t.^(1/Nu));

iC = 1;
iN = 2;
iY = 3;
iDividend = 4;
iii = 5;
iwage = 6;
iPi = 7;
itau = 8;
nVars = 8;

tax_unit = EqSs.DistE(3);
B = EqSs.B;
x0 = [C_t;N_t;Y_t;Dividend_t;ii_t;wage_t;Pi_t;tau_t];
x = x0;

C_t_Het = C_t;
N_t_Het = N_t;

while (Metric>TolEqTrans)
    tic;
    % Given relative prices, solve decision and distribution
    for t=TransPeriods:-1:1
        [EVPrime_t(dimEV{:},t),bPolicy(dimEV{:},t),nPolicy(dimEV{:},t),cPolicy(dimEV{:},t)] = VFI(EVPrime_t(dimEV{:},t+1),Params,tau_t(t),Dividend_t(t),wage_t(t),R_t(t));
    end
    
    for t=1:1:TransPeriods
        [bPolicyDist,nPolicyDist,cPolicyDist] = EXPAND_TO_DIST(bPolicy(dimEV{:},t),nPolicy(dimEV{:},t),wage_t(t),R_t(t),tau_t(t),Dividend_t(t),Params);
        Dist_t(dimDist{:},t+1) = UPDATE_DIST(Dist_t(dimDist{:},t),bPolicyDist,Params);
        [N_t_Het(t),~,~,~,C_t_Het(t),~,~] = AGGREGATE_SS(bPolicyDist,nPolicyDist,cPolicyDist,Dist_t(dimDist{:},t),Params,z_t(t),G_t(t),R_t(t),wage_t(t));
    end
    
    % Calibrate a distorting parameter
    eta1_t_new = 1 ./ (Beta*BetaShock_t.*R_t.* ([C_t_Het(2:end),EqSs.C] ./ C_t_Het).^(-Gamma));
    eta2_t_new = C_t_Het(1:TransPeriods).^(-Gamma).*wage_t ./ (N_t_Het.^(1/Nu));
    
    % x = [C_t;N_t;Y_t;Dividend_t;ii_t;wage_t;Pi_t;tau_t];
    MinorMetric = 1;
    MinorIter = 0;
    step = 1;
    while (MinorMetric>1e-8)
        [dX,F] = get_dX(x,Params,EqSs,eta1_t,eta2_t,BetaShock_t,MShock_t,z_t,G_t,tax_unit,B);
        x = x+step*dX;
        
        MinorMetric = max(abs(F(:)));
        if mod(MinorIter,10) ==0
            fprintf('  MinorIter: %d, MinorMetric: %g\n',MinorIter,MinorMetric);
        end
        MinorIter = MinorIter+1;
    end
    fprintf('  MinorIter: %d, MinorMetric: %g\n',MinorIter,MinorMetric);

    C_t_new = x(iC,:);
    N_t_new = x(iN,:);
    Y_t_new = x(iY,:);
    Dividend_t_new = x(iDividend,:);
    ii_t_new = x(iii,:);
    wage_t_new = x(iwage,:);
    Pi_t_new = x(iPi,:);
    tau_t_new = x(itau,:);
    
    % Other variables
    R_t_new = (1+ii_t) ./ [Pi_t(2:end) EqSs.Pi];
    
    % Do a line search
    Metric = max(abs([
        tau_t(:) - tau_t_new(:);
        Dividend_t(:) - Dividend_t_new(:);
        wage_t(:) - wage_t_new(:);
        R_t(:) - R_t_new(:);
        eta1_t(:) - eta1_t_new(:);
        ]));
    
    C_t = TransUpdateSpeed*C_t_new + (1-TransUpdateSpeed)*C_t;
    N_t = TransUpdateSpeed*N_t_new + (1-TransUpdateSpeed)*N_t;
    Y_t = TransUpdateSpeed*Y_t_new + (1-TransUpdateSpeed)*Y_t;
    Dividend_t = TransUpdateSpeed*Dividend_t_new + (1-TransUpdateSpeed)*Dividend_t;
    ii_t = TransUpdateSpeed*ii_t_new + (1-TransUpdateSpeed)*ii_t;
    tau_t = TransUpdateSpeed*tau_t_new + (1-TransUpdateSpeed)*tau_t;
    Pi_t = TransUpdateSpeed*Pi_t_new + (1-TransUpdateSpeed)*Pi_t;
    Dividend_t = TransUpdateSpeed*Dividend_t_new + (1-TransUpdateSpeed)*Dividend_t;
    wage_t = TransUpdateSpeed*wage_t_new + (1-TransUpdateSpeed)*wage_t;
    R_t = TransUpdateSpeed*R_t_new + (1-TransUpdateSpeed)*R_t;
    
    eta1_t = EtaUpdateSpeed*eta1_t_new + eta1_t.*(1-EtaUpdateSpeed);
    eta2_t = EtaUpdateSpeed*eta2_t_new + eta2_t.*(1-EtaUpdateSpeed);
    
    Iter = Iter+1;
    
    fprintf('Iter: %d. Metric: %g\n',Iter,Metric);
    toc;
end

EqTrans = v2struct(EVPrime_t,Dist_t,bPolicy_t,nPolicy_t,cPolicy_t,C_t,N_t,Y_t,Dividend_t,ii_t,wage_t,Pi_t,R_t,tau_t);
EqTransAgg = v2struct(C_t,N_t,Y_t,Dividend_t,ii_t,wage_t,Pi_t,R_t,tau_t);
end