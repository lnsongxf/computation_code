function [equation,jac,jacF,jacB] = get_symbol_rotemberg(SPEC)
warning('off');
syms C N Y Dividend ii wage Pi tau
syms C_1 N_1 Y_1 Dividend_1 ii_1 wage_1 Pi_1 tau_1
syms C_m1 N_m1 Y_m1 Dividend_m1 ii_m1 wage_m1 Pi_m1 tau_m1
 
vaR_m1 = [C_m1,N_m1,Y_m1,Dividend_m1,ii_m1,wage_m1,Pi_m1,tau_m1];
var =[C,N,Y,Dividend,ii,wage,Pi,tau];
vaR_1 = [C_1,N_1,Y_1,Dividend_1,ii_1,wage_1,Pi_1,tau_1];

equation = sym('a',[length(var) 1]);

ptr = 0;

equation(ptr+1) = 'Y - z*N';
ptr = ptr+1;


equation(ptr+1) = 'Y - C - G - Theta/2*(Pi-1)^2*Y';
ptr = ptr+1;

switch SPEC
    case 'peg_r'
        equation(ptr+1) = '(1+ii)/Pi_1 - R0';
    otherwise
        error('No Spec found');
end
ptr = ptr+1;


equation(ptr+1) = 'Y - wage*N - Theta/2*(Pi-1)^2*Y - Dividend';
ptr = ptr+1;

equation(ptr+1) = 'eta2*N^(1/Nu) - wage*C^(-Gamma)';
ptr = ptr+1;

equation(ptr+1) = 'eta1*Beta*BetaShock*(1+ii)/Pi_1*C_1^(-Gamma) - C^(-Gamma)';
ptr = ptr+1;

equation(ptr+1) = 'Beta*BetaShock * Pi_1 * (Pi_1-1) * Y_1/Y + Eta/Theta*(wage/z - (Eta-1)/Eta) - Pi*(Pi-1)';
ptr = ptr+1;

equation(ptr+1) = 'B+G-B/( (1+ii)/Pi_1 ) - tau*tax_unit';
ptr = ptr+1;


jac = jacobian(equation,var);
jacF = jacobian(equation,vaR_1);
jacB = jacobian(equation,vaR_m1);

warning('on');
end