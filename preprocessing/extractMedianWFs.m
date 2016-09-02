

function medWFs = extractMedianWFs(clu, st, Fs, datPath, dataType, dataSize, chanMap, gain)
% medWFs = extractMedianWFs(clu, st, cids, cgs, Fs, datPath, dataType, dataSize, chanMap, gain)
%
% - medWFs [nClusters, nChannels, nSamples] will be median waveforms of 
%    every cluster represented in clu
% -- note: if you desire only "good", try this:
% >> cluSub = clu(ismember(clu, cids(cgs==2))); % and same for spike times
%
% - clu [nSpikes,1] cluster identities
% - st [nSpikes,1] spike times (sec)
% - Fs [1,1] sampling frequency
% - datPath [string] filename
% - dataType [string] e.g. 'int16'
% - dataSize [2,1] e.g. [nChannelsInFile nSamples]
% - chanMap [vector] note: assume this is zero-indexed
% - gain [1,1] to convert to uV if desired, otherwise use 1

nWFsToLoad = 1000; % I find this to be enough, but 100 not enough. Could try other values.

% window is -0.5 to 1.25ms
wfWin = -round(0.5/1000*Fs):round(1.25/1000*Fs); nWFsamps = numel(wfWin);

nChInFile = dataSize(1);
nSamp = dataSize(2);
mmf = memmapfile(datPath, 'Format', {dataType, [nChInFile nSamp], 'x'});

cids = unique(clu);

nClu = length(cids);
nCh = length(chanMap);

medWFs = zeros(nClu, nCh, nWFsamps);

for q = 1:nClu
    fprintf(1, 'cluster %d (%d/%d)\n', cids(q), q, nClu);
    theseST = st(clu==cids(q));
    
    nWFsToLoad = min(nWFsToLoad, length(theseST));
    extractST = round(theseST(randperm(length(theseST), nWFsToLoad))*Fs);
    
    theseWF = zeros(nWFsToLoad, nCh, nWFsamps);
    for i=1:nWFsToLoad
        tempWF = mmf.Data.x(1:nChInFile,extractST(i)+wfWin(1):extractST(i)+wfWin(end));
        theseWF(i,:,:) = tempWF(chanMap+1,:);
    end
    
    
    medWF = squeeze(median(double(theseWF),1));
    medWFuV = medWF.*gain;
    
    medWFs(q,:,:) = medWFuV;
    
end