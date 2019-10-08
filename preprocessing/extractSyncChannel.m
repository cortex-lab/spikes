

function syncDat = extractSyncChannel(folder, numChans, syncChanIndex)
% extraChanIndices are 1-indexed

dataFiles = dir(fullfile(folder,'*.lf.bin'));

for d = 1:length(dataFiles)    
    
    filename = fullfile(folder, dataFiles(d).name);
    syncDat = extractSyncChannelFromFile(filename, numChans, syncChanIndex);
    
end
 

disp(' done.')