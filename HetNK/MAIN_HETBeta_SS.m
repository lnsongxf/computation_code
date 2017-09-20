% Het NK
% Author: Wenlan Luo

Params = SETUP_HETBeta;
Params = COMMON(Params);
% load('EqSs');
% EqSs = EQ_SS(Params,EqSs);

EqSs = EQ_SS(Params);

save('EqHetBetaSs','EqSs');