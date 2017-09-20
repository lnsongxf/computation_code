function get_analytic_nom
get_analytic(@get_symbol_rotemberg, '_nom_rotemberg');
end

function [eq,jac,jacF,jacB] = get_symbol
warning('off');
syms C N Y Dividend ii wage Pi
syms C_1 N_1 Y_1 Dividend_1 ii_1 wage_1 Pi_1
syms C_m1 N_m1 Y_m1 Dividend_m1 ii_m1 wage_m1 Pi_m1

vaR_m1 = [C_m1,N_m1,Y_m1,Dividend_m1,ii_m1,wage_m1,Pi_m1];
var =[C,N,Y,Dividend,ii,wage,Pi];
vaR_1 = [C_1,N_1,Y_1,Dividend_1,ii_1,wage_1,Pi_1];

eq = sym('a',[length(var) 1]);

ptr = 0;

eq(ptr+1) = 'Y - z*N';
ptr = ptr+1;

eq(ptr+1) = 'Y - C - G';
ptr = ptr+1;

eq(ptr+1) = 'ii - ii0';
ptr = ptr+1;

eq(ptr+1) = 'Y - wage*N - Dividend';
ptr = ptr+1;

eq(ptr+1) = 'eta2*N^(1/Nu) - wage*C^(-Gamma)';
ptr = ptr+1;

eq(ptr+1) = 'eta1*Beta*BetaShock*(1+ii)/Pi_1*C_1^(-Gamma) - C^(-Gamma)';
ptr = ptr+1;

eq(ptr+1) = 'Pi_1/(1+ii) * Pi_1 * (Pi_1-1) * Y_1/Y + Eta/Theta*(wage/z - (Eta-1)/Eta) - Pi*(Pi-1)';
ptr = ptr+1;



jac = jacobian(eq,var);
jacF = jacobian(eq,vaR_1);
jacB = jacobian(eq,vaR_m1);

% jac = sparse(jac);
% jacF = sparse(jacF);

warning('on');
end