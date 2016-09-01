

function S = sparsePCs(pcFeat,pcFeatInd, spikeTemplates)
%function S = sparsePCs(pdFeat,pcFeatInd)
% make a sparse matrix version of the pc features
% for the moment, only the top PC

nPCchans = size(pcFeat,3);
nFeat = size(pcFeat,2);
nSpikes = size(pcFeat,1);
nTemplates = size(pcFeatInd,1);
nChannels = double(max(pcFeatInd(:))+1);


rowInds = repmat([1:nSpikes]', [nPCchans*nFeat,1]);
colIndsTemp = zeros(nSpikes*nPCchans,1);

for q = 1:nPCchans
    colIndsTemp((q-1)*nSpikes+1:q*nSpikes) = pcFeatInd(spikeTemplates+1,q)+1;
end

colInds = zeros(nSpikes*nPCchans*nFeat,1);
for thisFeat = 1:3
    colInds(((thisFeat-1)*nSpikes*nPCchans+1):thisFeat*nSpikes*nPCchans) = ...
        (colIndsTemp-1)*nFeat+thisFeat;
end
clear colIndsTemp

pcFeatRS = zeros(nSpikes*nPCchans*nFeat,1);
for thisFeat = 1:3
    pcFeatRS(((thisFeat-1)*nSpikes*nPCchans+1):thisFeat*nSpikes*nPCchans) = ...
        double(reshape(squeeze(pcFeat(:,thisFeat,:)), [nSpikes*nPCchans 1]));
end
S = sparse(rowInds, colInds, pcFeatRS, nSpikes, nChannels*nFeat);