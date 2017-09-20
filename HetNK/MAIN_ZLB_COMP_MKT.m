% Het NK
% Author: Wenlan Luo

clear;

Params = SETUP;
Params = COMMON(Params);
Params = MODIFY_COMP_MKT(Params);

EqSsCompMkt = load('EqSsCompMkt');

% A beta shock that lowers real interest rate
Params.CommitPeriods = 20;
Params.BetaShock_t = ones(1,Params.TransPeriods);
Params.BetaShock_t(1:33) = 1.0014825;
Params.MShock_t = zeros(1,Params.TransPeriods);
Params.z_t = Params.z0*ones(1,Params.TransPeriods);
Params.G_t = zeros(1,Params.TransPeriods);

EqTrans = EQ_COMP_MKT_ZLB(Params,EqSsCompMkt.EqSs);

% save('EqSsCompMkt','EqSsCompMkt');
save('EqTransCompMkt_ZLB','EqTrans');