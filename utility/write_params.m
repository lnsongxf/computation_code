function write_params(params)
fieldNames = fieldnames(params);
str = [];
LINEBREAK = ' ';
for i=1:numel(fieldNames)
    fieldName = fieldNames{i};
    line = [fieldName '=params.' fieldName ';'];
    str = [str line LINEBREAK];
end

fid = fopen(['read_params.m'], 'wt');
fprintf(fid, '%s\n', str);
fclose(fid);
clear(['read_params.m']);
end