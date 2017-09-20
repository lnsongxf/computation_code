function [EVPrime_new,bPolicy,nPolicy,cPolicy] = VFI(EVPrime,Params,tau,Dividend,wage,R)
Chi = Params.Chi;
Nu = Params.Nu;
Gamma = Params.Gamma;

BGrid = Params.BGrid;
BPts = Params.BPts;
BMesh = Params.BMesh;
TauEMesh = Params.TauEMesh;
EMesh = Params.EMesh;
BetaMesh = Params.BetaMesh;
BMax = Params.BMax;

BetaPts = Params.BetaPts;
EPts = Params.EPts;
EBetaTrans = Params.EBetaTrans;

sizeEVPrime = Params.sizeEVPrime;

MEX_SOLVE_N_BP0 = Params.MEX_SOLVE_N_BP0;
NumThreads = Params.NumThreads;

% Endogenous grid
Up = R*BetaMesh.*EVPrime;
% Compute labor suply
nPositive = (wage*EMesh.*Up/Chi).^Nu;
% Compute current b
BCurrent = Up.^(-1/Gamma) + BMesh/R - Dividend + tau*TauEMesh - wage*EMesh.*nPositive;

% Interpolate Policy
bPolicy = zeros(sizeEVPrime);
nPositivePolicy = zeros(sizeEVPrime);
for j=1:EPts*BetaPts
    bPolicy(j,:) = min(interp1(BCurrent(j,:),BGrid,BGrid,'pchip'),BMax);
    nPositivePolicy(j,:) = interp1(BCurrent(j,:),nPositive(j,:),BGrid,'pchip');
end

% Use golden x to compute labor when bp is zero

BudgetWithoutLabor = BMesh - tau*TauEMesh + Dividend;
nZero = zeros(sizeEVPrime);
%{
focForLabor = @(n) (BudgetWithoutLabor(:)+wage*EMesh(:).*n(:)).^(1-Gamma)/(1-Gamma) - Chi*n(:).^(1+1/Nu)/(1+1/Nu) ;
[nZero,focResid] = goldenx(focForLabor,zeros(numel(EVPrime),1),10*ones(numel(EVPrime),1),100,1e-12);
nZero = reshape(nZero,sizeEVPrime);
%}

% %{
NumProblems = numel(nZero);
MEX_TASK = MEX_SOLVE_N_BP0;
HetNKMex;
%}

% Use corresponding value
nPolicy = nPositivePolicy.*(bPolicy>0) + nZero.*(bPolicy<=0);
bPolicy(bPolicy<0) = 0;

% Envelope Theorem to compute EVPrime_new
cPolicy = BudgetWithoutLabor + wage*EMesh.*nPolicy - bPolicy/R;
VPrime_new = cPolicy.^(-Gamma);
EVPrime_new = EBetaTrans*VPrime_new;
end