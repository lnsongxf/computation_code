function get_analytic_peg_r
get_analytic(@get_symbol, '');
end

function [eq,jac,jacF,jacB] = get_symbol
warning('off');
syms C N Y Dividend II wage Pi R pRatio S tau PA PB
syms C_1 N_1 Y_1 Dividend_1 II_1 wage_1 Pi_1 R_1 pRatio_1 S_1 tau_1 PA_1 PB_1
syms C_m1 N_m1 Y_m1 Dividend_m1 II_m1 wage_m1 Pi_m1 R_m1 pRatio_m1 S_m1 tau_m1 PA_m1 PB_m1

vaR_m1 = [C_m1,N_m1,Y_m1,Dividend_m1,II_m1,wage_m1,Pi_m1,R_m1,pRatio_m1,S_m1,tau_m1,PA_m1,PB_m1];
var =[C,N,Y,Dividend,II,wage,Pi,R,pRatio,S,tau,PA,PB];
vaR_1 = [C_1,N_1,Y_1,Dividend_1,II_1,wage_1,Pi_1,R_1,pRatio_1,S_1,tau_1,PA_1,PB_1];

eq = sym('a',[13 1]);

ptr = 0;
%1
eq(ptr+1) = 'eta1*Beta_C*R*C_1^(-Gamma) - C^(-Gamma)';
ptr = ptr+1;

%2
eq(ptr+1) = 'eta2*N^(1/Nu) - wage*C^(-Gamma)';
ptr = ptr+1;

%3
eq(ptr+1) = 'Y*S - z*N';
ptr = ptr+1;

%4
eq(ptr+1) = 'Y - wage*N - Dividend';
ptr = ptr+1;

%5
eq(ptr+1) = 'Y - C - G';
ptr = ptr+1;

%6
eq(ptr+1) = 'Theta*pRatio^(Mu/(1-Mu)) + Pi^(-Mu/(1-Mu))*(1-Theta)*S_m1 - S';
ptr = ptr+1;

%7
eq(ptr+1) = 'Pi  - ( (1-Theta) / (1-Theta*pRatio^(1/(1-Mu))) )^(1-Mu)';
ptr = ptr+1;

%8
eq(ptr+1) = 'pRatio*PB - PA';
ptr = ptr+1;

%9
eq(ptr+1) = 'Mu*wage/z*Y + (1-Theta)*Beta*Pi_1^(-Mu/(1-Mu))*PA_1 - PA';
ptr = ptr+1;

%10
eq(ptr+1) = 'Y + (1-Theta)*Beta*Pi_1^(-1/(1-Mu))*PB_1 - PB';
ptr = ptr+1;

%11
eq(ptr+1) = 'R*Pi_1 - II';
ptr = ptr+1;

%12
eq(ptr+1) = 'R-RFix';
ptr = ptr+1;

%13
eq(ptr+1) = 'B + G - B/R - tau*DistE3';
ptr = ptr+1;


jac = jacobian(eq,var);
jacF = jacobian(eq,vaR_1);
jacB = jacobian(eq,vaR_m1);

% jac = sparse(jac);
% jacF = sparse(jacF);

warning('on');
end