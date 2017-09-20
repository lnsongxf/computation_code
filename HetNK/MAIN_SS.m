% Het NK
% Author: Wenlan Luo

Params = SETUP;
Params = COMMON(Params);
% load('EqSs');
% EqSs = EQ_SS(Params,EqSs);

EqSs = EQ_SS(Params);

save('EqSs','EqSs');