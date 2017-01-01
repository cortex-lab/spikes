

function S = loadParamsPy(fn)
% Loads a phy-style "params.py" into a matlab struct. The params.py is
% python code but just a series of assignments, so most of it will run
% directly in matlab
%
% now based on text2struct from James Jun, 2016-12-31, modified slightly by
% N. Steinmetz. 

fid = fopen(fn, 'r');
mcFileMeta = textscan(fid, '%s%s', 'Delimiter', '=',  'ReturnOnError', false);
fclose(fid);
csName = mcFileMeta{1};
csValue = mcFileMeta{2};
S = struct();
for i=1:numel(csName)
    vcName1 = csName{i};
    if vcName1(1) == '~', vcName1(1) = []; end    
    try         
        if csValue{i}(1)==''''; % if first character of the value is a single quote
            % then it's a string and we can evaluate without adding a
            % single quote
            eval(sprintf('%s = %s;', vcName1, csValue{i}));
        else
            eval(sprintf('%s = ''%s'';', vcName1, csValue{i}));
        end
        eval(sprintf('num = str2double(%s);', vcName1));
        if ~isnan(num)
            eval(sprintf('%s = num;', vcName1));
        end
        eval(sprintf('S = setfield(S, ''%s'', %s);', vcName1, vcName1));
    catch
        fprintf('%s = %s error\n', csName{i}, csValue{i});
    end
end

% handle special case of true/false
fnames = fieldnames(S);
for f = 1:length(fnames)
    if strcmp(S.(fnames{f}), 'True') || strcmp(S.(fnames{f}), 'true')
        S.(fnames{f})=true;
    elseif strcmp(S.(fnames{f}), 'False') || strcmp(S.(fnames{f}), 'false')
        S.(fnames{f})=false;
    end
end


