
function writeToParamsPy(filename, fieldName, value)
% function writeToParamsPy(filename, fieldName, value)
% write a row in a params.py file. Will append to the end. 

if isstr(value)
    value = ['''' value '''']; % add a single quote on both sides
elseif value==true
    value = 'True';
elseif value==false;
    value = 'False';
else
    value = num2str(value);
end

fid = fopen(filename, 'a');

% here the newline comes before the field since params.py is not normally
% terminated with a newline
fprintf(fid, '\n%s = %s', fieldName, value);
fclose(fid);