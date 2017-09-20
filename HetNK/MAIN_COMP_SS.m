% Het NK
% Complete Market Steady State
% Author: Wenlan Luo

Params.Complete = true;
Params = SETUP(Params);
EqSsIncomp = load('EqSsIncomp.mat');

EqSs = EQ_SS_COMP(Params,EqSsIncomp.EqSs);
save('EqSsComp','EqSs');

