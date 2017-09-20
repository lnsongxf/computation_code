function [EqTrans,EqTransAgg] = EQ_TRANS(EqSs,Params)
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

while (Metric>TolEqTrans)
    tic;
    % Given relative prices, solve decision and distribution
    for t=TransPeriods:-1:1
        [EVPrime_t(dimEV{:},t),bPolicy(dimEV{:},t),nPolicy(dimEV{:},t),cPolicy(dimEV{:},t)] = VFI(EVPrime_t(dimEV{:},t+1),Params,tau_t(t),Dividend_t(t),wage_t(t),R_t(t));
    end
    
    for t=1:1:TransPeriods
        [bPolicyDist,nPolicyDist,cPolicyDist] = EXPAND_TO_DIST(bPolicy(dimEV{:},t),nPolicy(dimEV{:},t),wage_t(t),R_t(t),tau_t(t),Dividend_t(t),Params);
        Dist_t(dimDist{:},t+1) = UPDATE_DIST(Dist_t(dimDist{:},t),bPolicyDist,Params);
        [N_t(t),~,~,~,C_t(t),~,~] = AGGREGATE_SS(bPolicyDist,nPolicyDist,cPolicyDist,Dist_t(dimDist{:},t),Params,z_t(t),G_t(t),R_t(t),wage_t(t));
    end
    
    % Calibrate a distorting parameter
    eta1_t_new = 1 ./ (Beta*BetaShock_t.*R_t.* ([C_t(2:end),EqSs.C] ./ C_t).^(-Gamma));
    eta2_t_new = C_t(1:TransPeriods).^(-Gamma).*wage_t ./ (N_t.^(1/Nu));
    
    x = [C_t;N_t;Y_t;Dividend_t;ii_t;wage_t;Pi_t;tau_t];
    xx = fsolve(@(x) eq_nompeg_rotemberg(x,Params,EqSs,eta1_t,eta2_t,BetaShock_t,MShock_t,z_t,G_t,tax_unit,B), x(:), FsolveOptions);
    x = reshape(xx,nVars,TransPeriods);

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
    
    if Metric<1e-2
        EtaUpdateSpeed = 1;
    end
    eta1_t = EtaUpdateSpeed*eta1_t_new + eta1_t.*(1-EtaUpdateSpeed);
    eta2_t = EtaUpdateSpeed*eta2_t_new + eta2_t.*(1-EtaUpdateSpeed);
    
    Iter = Iter+1;
    
    fprintf('Iter: %d. Metric: %g\n',Iter,Metric);
    toc;
end

EqTrans = v2struct(EVPrime_t,Dist_t,bPolicy_t,nPolicy_t,cPolicy_t,C_t,N_t,Y_t,Dividend_t,ii_t,wage_t,Pi_t,R_t,tau_t);
EqTransAgg = v2struct(C_t,N_t,Y_t,Dividend_t,ii_t,wage_t,Pi_t,R_t,tau_t);
end