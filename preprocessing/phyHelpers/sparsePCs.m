

function S = sparsePCs(pcFeat,pcFeatInd, spikeTemplates, varargin)
%function S = sparsePCs(pdFeat,pcFeatInd, spikeTemplates[, nPCsPerChan, nPCchans])
% make a sparse matrix version of the pc features
% - nPCsPerChan is the number of PCs per channel to keep, usually 3 (but use
% fewer for memory reasons)
% - nPCchans is number of channels of PC to keep for each spike (like maybe
% 12, or fewer for memory reasons

if ~isempty(varargin)
    nPCsPerChan = varargin{1};
    nPCchans = varargin{2};
else
    nPCsPerChan = size(pcFeat,2);
    nPCchans = size(pcFeat,3);
end

nPCchans = min(nPCchans, size(pcFeat,3));

if nPCchans<size(pcFeat,3)
    pcFeat = pcFeat(:,:,1:nPCchans);
    pcFeatInd = pcFeatInd(:,1:nPCchans);
end

nPCsPerChan = min(nPCsPerChan,size(pcFeat,2));

if nPCsPerChan<size(pcFeat,2)
    pcFeat = pcFeat(:,1:nPCsPerChan,:);
end

nSpikes = size(pcFeat,1);
nTemplates = size(pcFeatInd,1);
nChannels = double(max(pcFeatInd(:))+1);

rowInds = repmat([1:nSpikes]', [nPCchans*nPCsPerChan,1]);
colIndsTemp = zeros(nSpikes*nPCchans,1);

for q = 1:nPCchans
    colIndsTemp((q-1)*nSpikes+1:q*nSpikes) = pcFeatInd(spikeTemplates+1,q)+1;
end

colInds = zeros(nSpikes*nPCchans*nPCsPerChan,1);
for thisFeat = 1:nPCsPerChan
    colInds(((thisFeat-1)*nSpikes*nPCchans+1):thisFeat*nSpikes*nPCchans) = ...
        (colIndsTemp-1)*nPCsPerChan+thisFeat;
end
clear colIndsTemp

pcFeatRS = zeros(nSpikes*nPCchans*nPCsPerChan,1);
for thisFeat = 1:nPCsPerChan
    pcFeatRS(((thisFeat-1)*nSpikes*nPCchans+1):thisFeat*nSpikes*nPCchans) = ...
        double(reshape(squeeze(pcFeat(:,thisFeat,:)), [nSpikes*nPCchans 1]));
end
S = sparse(rowInds, colInds, pcFeatRS, nSpikes, nChannels*nPCsPerChan);