% Het NK
% Author: Wenlan Luo

Params.Incomplete = true;
Params.peg_r = true;
Params = SETUP(Params);
Params = COMMON(Params);
EqSsIncomp = load('EqSsIncomp.mat');

% Specify an exgeonous real interest rate
Params.R0_t = Params.R0 * ones(1,Params.TransPeriods);
Params.R0_t(20) = 1;

[EqTrans,EqTransAgg] = EQ_TRANS_INCOMP(EqSsIncomp.EqSs,Params);

save('EqTrans_IncompPegR','EqTrans');
save('EqTransAgg_IncompPegR','EqTransAgg');