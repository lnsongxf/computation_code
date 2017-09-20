% Het NK
% Author: Wenlan Luo

clear;

Params = SETUP;
Params = COMMON(Params);
Params = MODIFY_COMP_MKT(Params);

EqSsCompMkt = load('EqSsCompMkt');

% A beta shock that lowers real interest rate
Params.PegPeriods = 0;
Params.BetaShock_t = ones(1,Params.TransPeriods);
% Params.BetaShock_t(1:33) = 1.0014825;
Params.MShock_t = zeros(1,Params.TransPeriods);
% Params.MShock_t(1:20) = -100;
Params.z_t = Params.z0*ones(1,Params.TransPeriods);
Params.G_t = zeros(1,Params.TransPeriods);
Params.G_t(28) = EqSsCompMkt.EqSs.Y*0.29;
% Params.z_t(10) = Params.z0*0.999;

% EqTrans = EQ_NOMPEG_COMP(Params,EqSsCompMkt.EqSs);
EqTrans = EQ_NOMPEG_COMP(Params,EqSsCompMkt.EqSs,EqTrans);

% save('EqSsCompMkt','EqSsCompMkt');
save('EqTransComp_Nompeg','EqTrans');