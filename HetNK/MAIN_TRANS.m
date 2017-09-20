% Het NK
% Author: Wenlan Luo

Params = SETUP;
Params = COMMON(Params);
% Params.TransJacPattern = GET_JAC_PATTERN(Params);
load('EqSs');

% Specify an exgeonous real interest rate
Params.R0_t = Params.R0 * ones(1,Params.TransPeriods);
% Lower the interest rate at quarter 20
Params.R0_t(20) = 1;
Params.G_t = zeros(1,Params.TransPeriods);
Params.z_t = ones(1,Params.TransPeriods);

[EqTrans,EqTransAgg] = EQ_TRANS(EqSs,Params);

save('EqTrans','EqTrans');
save('EqTransAgg','EqTransAgg');