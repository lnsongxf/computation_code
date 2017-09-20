% Het NK
% Author: Wenlan Luo

Params.Complete = true;
Params = SETUP(Params);
Params = COMMON(Params);
Params = MODIFY_COMP_MKT(Params);

load('EqSs');

% Specify an exgeonous real interest rate
Params.R0_t = Params.R0 * ones(1,Params.TransPeriods);
% Lower the interest rate at quarter 20
Params.R0_t(20) = Params.R0 - 0.005;
Params.G_t = zeros(1,Params.TransPeriods);
Params.z_t = ones(1,Params.TransPeriods);

[EqSs, EqTrans] = EQ_COMP_MKT(Params,EqSs);

save('EqSsCompMkt','EqSs');
save('EqTransCompMkt','EqTrans');