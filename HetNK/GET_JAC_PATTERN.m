function TransJacPattern = GET_JAC_PATTERN(Params)
TransPeriods = Params.TransPeriods;

%% Transition path Jacobian pattern
TransJacPattern = sparse(13*TransPeriods,13*TransPeriods);

x = [1:13*TransPeriods]';

ptr = 0;
%1 C
C_t = x(ptr+1:ptr+TransPeriods);
ptr = ptr+TransPeriods;

%2 N
N_t = x(ptr+1:ptr+TransPeriods);
ptr = ptr+TransPeriods;

%3 Y
Y_t = x(ptr+1:ptr+TransPeriods);
ptr = ptr+TransPeriods;

%4 D
Dividend_t = x(ptr+1:ptr+TransPeriods);
ptr = ptr+TransPeriods;

%5 i
i_t = x(ptr+1:ptr+TransPeriods);
ptr = ptr+TransPeriods;

%6 wage
wage_t = x(ptr+1:ptr+TransPeriods);
ptr = ptr+TransPeriods;

%7 Pi
Pi_t = x(ptr+1:ptr+TransPeriods);
ptr = ptr+TransPeriods;

%8 r
r_t = x(ptr+1:ptr+TransPeriods);
ptr = ptr+TransPeriods;

%9 pRatio
pRatio_t = x(ptr+1:ptr+TransPeriods);
ptr = ptr+TransPeriods;

%10 S
S_t = x(ptr+1:ptr+TransPeriods);
ptr = ptr+TransPeriods;

%11 tau
tau_t = x(ptr+1:ptr+TransPeriods);
ptr = ptr+TransPeriods;

%12 PA
PA_t = x(ptr+1:ptr+TransPeriods);
ptr = ptr+TransPeriods;

%13 PB
PB_t = x(ptr+1:ptr+TransPeriods);
ptr = ptr+TransPeriods;


% Forward band
CurrentBand = speye(TransPeriods);
ForwardBand = spdiags(ones(TransPeriods,1),1,TransPeriods,TransPeriods);

ptr = 0;

%1
TransJacPattern(ptr+1:ptr+TransPeriods,C_t) = CurrentBand + ForwardBand;
TransJacPattern(ptr+1:ptr+TransPeriods,r_t) = CurrentBand;
ptr = ptr+TransPeriods;

%2
TransJacPattern(ptr+1:ptr+TransPeriods,[N_t]) = CurrentBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[C_t]) = CurrentBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[wage_t]) = CurrentBand;
ptr = ptr+TransPeriods;

%3
TransJacPattern(ptr+1:ptr+TransPeriods,[Y_t]) = CurrentBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[S_t]) = CurrentBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[N_t]) = CurrentBand;
ptr = ptr+TransPeriods;

%4
TransJacPattern(ptr+1:ptr+TransPeriods,[Y_t]) = CurrentBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[wage_t]) = CurrentBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[N_t]) = CurrentBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[Dividend_t]) = CurrentBand;
ptr = ptr+TransPeriods;

%5
TransJacPattern(ptr+1:ptr+TransPeriods,[Y_t]) = CurrentBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[C_t]) = CurrentBand;
ptr = ptr+TransPeriods;

%6
TransJacPattern(ptr+1:ptr+TransPeriods,[S_t]) = CurrentBand + ForwardBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[pRatio_t]) = ForwardBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[Pi_t]) = ForwardBand;
ptr = ptr+TransPeriods;

%7
TransJacPattern(ptr+1:ptr+TransPeriods,[Pi_t]) = CurrentBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[pRatio_t]) = CurrentBand;
ptr = ptr+TransPeriods;

%8
TransJacPattern(ptr+1:ptr+TransPeriods,[pRatio_t]) = CurrentBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[PA_t]) = CurrentBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[PB_t]) = CurrentBand;
ptr = ptr+TransPeriods;

%9
TransJacPattern(ptr+1:ptr+TransPeriods,[wage_t]) = CurrentBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[PA_t]) = CurrentBand + ForwardBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[Y_t]) = CurrentBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[Pi_t]) = ForwardBand;
ptr = ptr+TransPeriods;

%10
TransJacPattern(ptr+1:ptr+TransPeriods,[PB_t]) = CurrentBand + ForwardBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[Y_t]) = CurrentBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[Pi_t]) = ForwardBand;
ptr = ptr+TransPeriods;

%11
TransJacPattern(ptr+1:ptr+TransPeriods,[r_t]) = CurrentBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[Pi_t]) = ForwardBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[i_t]) = CurrentBand;
ptr = ptr+TransPeriods;

%12
TransJacPattern(ptr+1:ptr+TransPeriods,[r_t]) = CurrentBand;
ptr = ptr+TransPeriods;

%13
TransJacPattern(ptr+1:ptr+TransPeriods,[r_t]) = CurrentBand;
TransJacPattern(ptr+1:ptr+TransPeriods,[tau_t]) = CurrentBand;
ptr = ptr+TransPeriods;
end