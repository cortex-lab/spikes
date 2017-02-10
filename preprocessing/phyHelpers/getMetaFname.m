
function fn = getMetaFname(ksDir, rawDir, type)
% function fn = getMetaFname(ksDir, rawDir, type)
% 
% type is either 'LFP', or 'AP'
% expects that params.py is in ksDir, and 

paramsFn = fullfile(ksDir, 'params.py');
p = loadParamsPy(paramsFn);
[~, rawFn] = fileparts(p.dat_path);
if strcmp(rawFn(end-3:end), '_CAR')
    rawFn = rawFn(1:end-4);
end

switch type
    case 'LFP'
        metaFn = [rawFn(1:end-2) 'lf.meta']; % chop "ap"
    case 'AP'
        metaFn = [rawFn '.meta'];
    otherwise
        fprintf(1, 'supply type = ''AP'' or ''LFP''\n');
end

fn = fullfile(rawDir, metaFn);