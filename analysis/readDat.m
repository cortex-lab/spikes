
function dat = readDat(filename, numChans, varargin)
% takes the same amount of time to read through for just one channel, but
% you don't run into memory problems that way if the file is too big. 

requestedChan = [];
if ~isempty(varargin)
    requestedChan = varargin{1};
end

fid = fopen(filename);


if isempty(requestedChan)
    dat = fread(fid, [numChans, Inf], 'int16=>int16');
else

    if requestedChan>1
        % skip over the first samples of the other channels
        q = fread(fid, (requestedChan-1), 'int16=>int16'); 
    end
    dat = fread(fid, [1, Inf], 'int16=>int16', (numChans-1)*2); % skipping other channels
end

fclose(fid);