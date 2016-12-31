

function syncDat = extractSyncChannel(folder, numChans, syncChanIndex)
% extraChanIndices are 1-indexed

maxReadSize = 1e9;
dataFiles = dir(fullfile(folder,'*.lf.bin'));

for d = 1:length(dataFiles)    
    
    [~,fn] = fileparts(dataFiles(d).name);
    syncFname =  fullfile(folder, [fn '_sync.dat']);
    fidOut = fopen(syncFname, 'w');
    
    fprintf(1,' loading %s\n', dataFiles(d).name);
    
    fid = fopen(fullfile(folder, dataFiles(d).name), 'r');
    
    % skip over the first samples of the other channels
    q = fread(fid, (syncChanIndex-1), 'int16=>int16');
    
    nSamp = dataFiles(d).bytes/2/numChans;
    
    if nargout>0
        syncDat = zeros(1, nSamp, 'int16');
    end
    
    nBatch = floor(nSamp/maxReadSize);
    for b = 1:nBatch
        dat = fread(fid, [1, maxReadSize], 'int16=>int16', (numChans-1)*2); % skipping other channels
        fwrite(fidOut, dat, 'int16');
        if nargout>0
            syncDat((b-1)*maxReadSize+1:b*maxReadSize) = dat;
        end
    end
    
    % all the other samples
    dat = fread(fid, [1, Inf], 'int16=>int16', (numChans-1)*2); % skipping other channels
    fwrite(fidOut, dat, 'int16');
    if nargout>0
        syncDat(nBatch*maxReadSize+1:end) = dat;
    end
end

fclose(fid);
fclose(fidOut);   

disp(' done.')