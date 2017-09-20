% Het NK
% Incomplete Market Steady State
% Author: Wenlan Luo

Params.Incomplete = true;
Params = SETUP(Params);
Params = COMMON(Params);

EqSs = EQ_SS_INCOMP(Params);
%{
EqSs = EQ_SS_INCOMP(Params,EqSsIncomp.EqSs);
%}
save('EqSsIncomp','EqSs');
