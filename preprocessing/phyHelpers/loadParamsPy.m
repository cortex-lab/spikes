

function p = loadParamsPy(fn)
% Loads a phy-style "params.py" into a matlab struct. The params.py is
% python code but just a series of assignments, so most of it will run
% directly in matlab

% have to copy params.py into params.m because "run" won't try to run 
% something without .m extension
newFn = fullfile(fileparts(fn), 'params.m');
copyfile(fn, newFn); 

try
    % use evalc here so that the lack of semicolons doesn't produce
    % terminal output
    evalc('run(newFn)');
catch
    % have to put this in a try/catch because "False" will produce an
    % error. A future version of this code could replace True/False with
    % true/false first. 
end

clear fn newFn
q = whos();

for x = 1:length(q)    
    eval(sprintf('p.(q(x).name) = %s;', q(x).name));    
end


