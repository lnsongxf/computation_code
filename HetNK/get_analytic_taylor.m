function get_analytic_taylor
get_analytic(@get_symbol, '_taylor');
end

function [eq,jac,jacF,jacB] = get_symbol
warning('off');
syms C N Y Dividend ii wage Pi pRatio S PA PB
syms C_1 N_1 Y_1 Dividend_1 ii_1 wage_1 Pi_1 pRatio_1 S_1 PA_1 PB_1
syms C_m1 N_m1 Y_m1 Dividend_m1 ii_m1 wage_m1 Pi_m1 pRatio_m1 S_m1 PA_m1 PB_m1

vaR_m1 = [C_m1,N_m1,Y_m1,Dividend_m1,ii_m1,wage_m1,Pi_m1,pRatio_m1,S_m1,PA_m1,PB_m1];
var =[C,N,Y,Dividend,ii,wage,Pi,pRatio,S,PA,PB];
vaR_1 = [C_1,N_1,Y_1,Dividend_1,ii_1,wage_1,Pi_1,pRatio_1,S_1,PA_1,PB_1];

eq = sym('a',[11 1]);

ptr = 0;

eq(ptr+1) = 'Y*S - z*N';
ptr = ptr+1;

eq(ptr+1) = 'Y - C - G';
ptr = ptr+1;

eq(ptr+1) = '-ii + sin(R0-1 + Phi*(Pi-1) + MShock)';
ptr = ptr+1;

eq(ptr+1) = '-Pi + ( (1-Theta) / (1-Theta*pRatio^(1/(1-Mu))) )^(1-Mu)';
ptr = ptr+1;

eq(ptr+1) = '-pRatio*PB + PA';
ptr = ptr+1;

eq(ptr+1) = '-S + Theta*pRatio^(Mu/(1-Mu)) + Pi^(-Mu/(1-Mu))*(1-Theta)*S_m1';
ptr = ptr+1;

eq(ptr+1) = '-PB + Y + (1-Theta)*Beta*BetaShock*Pi_1^(-1/(1-Mu))*PB_1';
ptr = ptr+1;

eq(ptr+1) = '-PA + Mu*wage/z*Y + (1-Theta)*Beta*BetaShock*Pi_1^(-Mu/(1-Mu))*PA_1';
ptr = ptr+1;

eq(ptr+1) = 'Y - wage*N - Dividend';
ptr = ptr+1;

eq(ptr+1) = 'eta2*N^(1/Nu) - wage*C^(-Gamma)';
ptr = ptr+1;

eq(ptr+1) = 'eta1*Beta*BetaShock*(1+ii)/Pi_1*C_1^(-Gamma) - C^(-Gamma)';
ptr = ptr+1;



jac = jacobian(eq,var);
jacF = jacobian(eq,vaR_1);
jacB = jacobian(eq,vaR_m1);

% jac = sparse(jac);
% jacF = sparse(jacF);

warning('on');
end