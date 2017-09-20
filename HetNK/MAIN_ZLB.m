% Het NK
% Author: Wenlan Luo

Params.Incomplete = true;
Params = SETUP(Params);
Params = COMMON(Params);
load('EqSs');

% Zero other shock
Params.G_t = zeros(1,Params.TransPeriods);
Params.z_t = ones(1,Params.TransPeriods);

% A beta shock that lowers real interest rate
Params.CommitPeriods = 20;
Params.BetaShock_t = ones(1,Params.TransPeriods);
Params.BetaShock_t(1:33) = 1.001638;
Params.MShock_t = zeros(1,Params.TransPeriods);

[EqTrans,EqTransAgg] = EQ_TRANS_ZLB(EqSs,Params);

save('EqTransInComp_ZLB','EqTrans');
save('EqTransInComp_ZLB_Agg','EqTransAgg');