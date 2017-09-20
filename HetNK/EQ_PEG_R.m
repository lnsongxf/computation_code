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
eta1_t = 1 ./ (Beta*R_t.* ([C_t(2:end),EqSs.C] ./ C_t).^(-Gamma));
eta2_t = C_t(1:TransPeriods).^(-Gamma).*wage_t ./ (N_t.^(1/Nu));

nVars = 13;

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
    eta1_t = 1 ./ (Beta*R_t.* ([C_t(2:end),EqSs.C] ./ C_t).^(-Gamma));
    eta2_t = C_t(1:TransPeriods).^(-Gamma).*wage_t ./ (N_t.^(1/Nu));
    
    % Solve the system of equations for aggregates
    %{
    eq_resid = eq_aggregate(x,Params,EqSs.C,0,1,EqSs.PA,EqSs.PB,eta1_t,eta2_t,z_t,G_t,r0_t,EqSs.B,EqSs.DistE);
    max(abs(eq_resid))
    %}
    % x = fsolve(@(x) eq_aggregate(x,Params,EqSs.C,0,1,EqSs.PA,EqSs.PB,eta1_t,eta2_t,z_t,G_t,r0_t,EqSs.B,EqSs.DistE), x, options);
    
    % Compute the block matrix
    Coef_B = zeros(nVars,nVars,TransPeriods);
    Coef = zeros(nVars,nVars,TransPeriods);
    Coef_F = zeros(nVars,nVars,TransPeriods);
    F = zeros(nVars,TransPeriods);
    
    x = [C_t;N_t;Y_t;Dividend_t;ii_t;wage_t;Pi_t;R_t;pRatio_t;S_t;tau_t;PA_t;PB_t];
    MinorMetric = 1;
    while (MinorMetric>1e-8)
        for t=1:TransPeriods
            % Unpack things at current period
            % var =[C,N,Y,Dividend,ii,wage,Pi,r,pRatio,S,tau,PA,PB];
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

            eta1 = eta1_t(t);
            eta2 = eta2_t(t);
            

            Beta_C = Beta;
            z = z_t(t);
            G = G_t(t);
            B = EqSs.B;
            RFix = R0_t(t);
            DistE3 = EqSs.DistE(3);

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
        
        MinorMetric = max(abs(F(:)));
    end
    
    ptr = 0;
    
    %1 C
    C_t_new = x(ptr+1,:);
    ptr = ptr+1;
    
    %2 N
    N_t_new = x(ptr+1,:);
    ptr = ptr+1;
    
    %3 Y
    Y_t_new = x(ptr+1,:);
    ptr = ptr+1;
    
    %4 D
    Dividend_t_new = x(ptr+1,:);
    ptr = ptr+1;
    
    %5 i
    II_t_new = x(ptr+1,:);
    ptr = ptr+1;
    
    %6 wage
    wage_t_new = x(ptr+1,:);
    ptr = ptr+1;
    
    %7 Pi
    Pi_t_new = x(ptr+1,:);
    ptr = ptr+1;
    
    %8 r
    R_t_new = x(ptr+1,:);
    ptr = ptr+1;
    
    %9 pRatio
    pRatio_t_new = x(ptr+1,:);
    ptr = ptr+1;
    
    %10 S
    S_t_new = x(ptr+1,:);
    ptr = ptr+1;
    
    %11 tau
    tau_t_new = x(ptr+1,:);
    ptr = ptr+1;
    
    %12 PA
    PA_t_new = x(ptr+1,:);
    ptr = ptr+1;
    
    %13 PB
    PB_t_new = x(ptr+1,:);
    ptr = ptr+1;
    
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
    ii_t = TransUpdateSpeed*II_t_new + (1-TransUpdateSpeed)*ii_t;
    tau_t = TransUpdateSpeed*tau_t_new + (1-TransUpdateSpeed)*tau_t;
    Pi_t = TransUpdateSpeed*Pi_t_new + (1-TransUpdateSpeed)*Pi_t;
    pRatio_t = TransUpdateSpeed*pRatio_t_new + (1-TransUpdateSpeed)*pRatio_t;
    S_t = TransUpdateSpeed*S_t_new + (1-TransUpdateSpeed)*S_t;
    PA_t = TransUpdateSpeed*PA_t_new + (1-TransUpdateSpeed)*PA_t;
    PB_t = TransUpdateSpeed*PB_t_new + (1-TransUpdateSpeed)*PB_t;
    Dividend_t = TransUpdateSpeed*Dividend_t_new + (1-TransUpdateSpeed)*Dividend_t;
    wage_t = TransUpdateSpeed*wage_t_new + (1-TransUpdateSpeed)*wage_t;
    R_t = TransUpdateSpeed*R_t_new + (1-TransUpdateSpeed)*R_t;
    
    Iter = Iter+1;
    
    fprintf('Iter: %d. Metric: %g\n',Iter,Metric);
    toc;
end

EqTrans = v2struct(EVPrime_t,Dist_t,bPolicy_t,nPolicy_t,cPolicy_t,C_t,N_t,Y_t,Dividend_t,ii_t,wage_t,Pi_t,R_t,pRatio_t,S_t,tau_t,PA_t,PB_t);
EqTransAgg = v2struct(C_t,N_t,Y_t,Dividend_t,ii_t,wage_t,Pi_t,R_t,pRatio_t,S_t,tau_t,PA_t,PB_t);
end