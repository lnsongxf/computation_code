function get_analytic(get_symbol,attach)
% Get Jacobian
[equation,jac,jacF,jacB] = get_symbol();

get_func_from(equation,attach);
get_func_from(jac,attach);
get_func_from(jacF,attach);
get_func_from(jacB,attach);
end

function func = get_func_from(func,attach)
func_str = func2str(matlabFunction(func));
for j=1:length(func_str)
    if func_str(j)==')'
        break;
    end
end
func_str(1:j) = [];
func = ['TEMP=' func_str ';'];
func = regexprep(func,'sin(','max(0,');
func = regexprep(func,'cos(','(0<');
fid = fopen(['generated_equation/' inputname(1),attach,'.m'], 'wt');
fprintf(fid, '%s\n', func);
fclose(fid);
clear([inputname(1),attach,'_func.m']);
end